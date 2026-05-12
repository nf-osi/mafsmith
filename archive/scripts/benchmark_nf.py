#!/usr/bin/env python3
"""
Benchmark mafsmith vs vcf2maf.pl on NF Data Portal VCFs.

Downloads up to --n-vcfs Strelka somatic VCFs from the NF Data Portal Synapse
fileview (syn16858331), runs fastVEP annotation once per file, then times
vcf2maf.pl (--inhibit-vep) and mafsmith (--skip-annotation) on the annotated
output. Reports variants/sec, concordance, variant type distribution, and
total clock time. Results are saved to a JSON file for reproducibility.

Usage:
    python3 scripts/benchmark_nf.py                        # 100 VCFs, default output
    python3 scripts/benchmark_nf.py --n-vcfs 20            # quick test with 20 VCFs
    python3 scripts/benchmark_nf.py --results-dir /tmp/nf_bench --n-vcfs 50
    python3 scripts/benchmark_nf.py --report-only results/benchmark_nf_YYYY.json

Requirements:
    - synapseclient (pip install synapseclient)
    - conda env 'vcf2maf-env' with vcf2maf.pl
    - ~/.mafsmith/bin/fastvep
    - ~/.mafsmith/GRCh38/reference.fa + genes.gff3
    - mafsmith binary at target/release/mafsmith
"""

import argparse
import collections
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import synapseclient

# ── Paths ──────────────────────────────────────────────────────────────────────
REPO       = Path(__file__).parent.parent
MAFSMITH   = str(REPO / "target" / "release" / "mafsmith")
FASTVEP    = os.path.expanduser("~/.mafsmith/bin/fastvep")
REF_FASTA  = os.path.expanduser("~/.mafsmith/GRCh38/reference.fa")
GFF3       = os.path.expanduser("~/.mafsmith/GRCh38/genes.gff3")

# Saved vcf2maf.pl reference MAFs reuse directory (same as synapse_validate.py)
EXPECTED_DIR = REPO / "tests" / "fixtures" / "expected"

# NF Data Portal Synapse fileview entity ID
NF_FILEVIEW = "syn16858331"

# Variant types recognised as meaningful in output comparison
COMPARE_FIELDS = [
    "Hugo_Symbol", "Variant_Classification", "Variant_Type",
    "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
    "HGVSc", "HGVSp", "HGVSp_Short", "Transcript_ID", "Exon_Number",
    "t_depth", "t_ref_count", "t_alt_count",
    "n_depth", "n_ref_count", "n_alt_count",
]


# ── Helpers ────────────────────────────────────────────────────────────────────

def run(cmd, **kwargs):
    r = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if r.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(str(c) for c in cmd)}\n{r.stderr[-800:]}"
        )
    return r


def read_vcf_header(path):
    """Return (tumor_id, normal_id, has_chr_prefix, is_grch38)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    samples, is_grch38, has_chr = [], None, None
    with opener(path, "rt") as f:
        for line in f:
            if "GRCh38" in line or "hg38" in line:
                is_grch38 = True
            if ("GRCh37" in line or "hg19" in line or "b37" in line) and is_grch38 is None:
                is_grch38 = False
            if line.startswith("#CHROM"):
                cols = line.strip().split("\t")
                if len(cols) > 9:
                    samples = cols[9:]
            elif not line.startswith("#"):
                has_chr = line.startswith("chr")
                break
    tumor_id  = samples[0] if samples else None
    normal_id = samples[1] if len(samples) > 1 else None
    return tumor_id, normal_id, has_chr, is_grch38


def decompress_and_strip_chr(src, dst, has_chr):
    opener = gzip.open if str(src).endswith(".gz") else open
    with opener(src, "rt") as f, open(dst, "w") as out:
        for line in f:
            if has_chr and not line.startswith("#") and line.startswith("chr"):
                out.write(line[3:])
            else:
                out.write(line)


def parse_maf_for_comparison(path, skip_empty_chr=False):
    """Parse MAF into {(chr, pos): [row_dict, ...]}. Returns (dict, skipped_count)."""
    records: dict = {}
    skipped = 0
    with open(path) as f:
        header = None
        for line in f:
            if line.startswith("#version"):
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            chr_col = header.index("Chromosome") if "Chromosome" in header else 4
            pos_col = header.index("Start_Position") if "Start_Position" in header else 5
            if skip_empty_chr and not fields[chr_col]:
                skipped += 1
                continue
            key = (fields[chr_col], fields[pos_col])
            records.setdefault(key, []).append(dict(zip(header, fields)))
    return records, skipped


def normalize_hgvsp(s):
    for three, one in [
        ("Ala","A"),("Arg","R"),("Asn","N"),("Asp","D"),("Cys","C"),
        ("Gln","Q"),("Glu","E"),("Gly","G"),("His","H"),("Ile","I"),
        ("Leu","L"),("Lys","K"),("Met","M"),("Phe","F"),("Pro","P"),
        ("Ser","S"),("Thr","T"),("Trp","W"),("Tyr","Y"),("Val","V"),
        ("Ter","*"),("Sec","U"),("Pyl","O"),("Xaa","X"),
    ]:
        s = s.replace(three, one)
    return s


HGVS_COLS = {"HGVSp", "HGVSp_Short", "HGVSc"}
SV_ALTS   = {"<BND>", "<DEL>", "<DUP>", "<INV>", "<TRA>"}


def _field_diffs(r1, r2):
    diffs = 0
    for field in COMPARE_FIELDS:
        v1, v2 = r1.get(field, ""), r2.get(field, "")
        if field in HGVS_COLS:
            v1, v2 = normalize_hgvsp(v1), normalize_hgvsp(v2)
        if v1 != v2:
            diffs += 1
    return diffs


def compare_mafs(ms_path, vc_path):
    ms, _         = parse_maf_for_comparison(ms_path, skip_empty_chr=False)
    vc, sv_sec_vc = parse_maf_for_comparison(vc_path, skip_empty_chr=True)

    total_ms = sum(len(v) for v in ms.values())
    field_mismatches = 0
    variant_mismatches = 0  # rows with at least one field difference
    sv_secondary_ms = 0
    variant_types: dict = collections.Counter()
    matched_ms_keys: set = set()

    # Count variant types seen in mafsmith output
    for rows in ms.values():
        for r in rows:
            variant_types[r.get("Variant_Type", "?")] += 1
            variant_types[r.get("Variant_Classification", "?")] += 1

    for key in sorted(vc):
        vc_row = vc[key][0]
        if key not in ms:
            field_mismatches += 1
            variant_mismatches += 1
            continue
        matched_ms_keys.add(key)
        ms_row = min(ms[key], key=lambda r: _field_diffs(r, vc_row))
        diffs = _field_diffs(ms_row, vc_row)
        field_mismatches += diffs
        if diffs > 0:
            variant_mismatches += 1

    for key in set(ms) - matched_ms_keys:
        for row in ms[key]:
            if row.get("Tumor_Seq_Allele2", "") in SV_ALTS:
                sv_secondary_ms += 1
            else:
                field_mismatches += 1
                variant_mismatches += 1

    # Concordance metrics:
    # - variant_concordance: fraction of variants with zero field differences (most meaningful)
    # - field_concordance: fraction of compared field values that agree
    compared_variants = sum(len(v) for v in vc.values())
    variant_concordance = (compared_variants - variant_mismatches) / compared_variants * 100 if compared_variants else 0
    total_fields_compared = compared_variants * len(COMPARE_FIELDS)
    field_concordance = (total_fields_compared - field_mismatches) / total_fields_compared * 100 if total_fields_compared else 0

    return {
        "total_ms_rows": total_ms,
        "total_vc_rows": compared_variants,
        "sv_secondary_vc": sv_sec_vc,
        "sv_secondary_ms": sv_secondary_ms,
        "field_mismatches": field_mismatches,
        "variant_mismatches": variant_mismatches,
        "variant_concordance_pct": variant_concordance,
        "field_concordance_pct": field_concordance,
        "variant_types": dict(variant_types),
    }


# ── Core per-file processing ────────────────────────────────────────────────────

def process_one(syn_id, vcf_name, vcf_path, work_dir, force_vcf2maf=False):
    tumor_id, normal_id, has_chr, is_grch38 = read_vcf_header(vcf_path)

    if is_grch38 is False:
        return {"status": "skipped", "reason": "Not GRCh38"}
    if not tumor_id:
        return {"status": "skipped", "reason": "No sample columns"}

    # Decompress (fastVEP requires plain-text VCF)
    plain_vcf = os.path.join(work_dir, "input_plain.vcf")
    decompress_and_strip_chr(vcf_path, plain_vcf, has_chr or False)

    # ── fastVEP annotation (shared input for both tools) ──────────────────────
    fv_out = os.path.join(work_dir, "fastvep_annotated.vcf")
    t0 = time.perf_counter()
    run([FASTVEP, "annotate",
         "-i", plain_vcf, "-o", fv_out,
         "--fasta", REF_FASTA, "--gff3", GFF3,
         "--hgvs", "--output-format", "vcf"])
    fastvep_secs = time.perf_counter() - t0

    # Count annotated variants (data lines)
    n_variants = sum(1 for line in open(fv_out) if not line.startswith("#"))

    # ── vcf2maf.pl ────────────────────────────────────────────────────────────
    saved_ref = EXPECTED_DIR / f"{syn_id}_vcf2maf.maf"
    if not saved_ref.exists() or force_vcf2maf:
        EXPECTED_DIR.mkdir(parents=True, exist_ok=True)
        cmd_vc = ["conda", "run", "-n", "vcf2maf-env", "vcf2maf.pl",
                  "--input-vcf",  fv_out,
                  "--output-maf", str(saved_ref),
                  "--tumor-id",     tumor_id,
                  "--vcf-tumor-id", tumor_id,
                  "--ref-fasta",  REF_FASTA,
                  "--inhibit-vep"]
        if normal_id:
            cmd_vc += ["--normal-id", normal_id, "--vcf-normal-id", normal_id]
        t0 = time.perf_counter()
        run(cmd_vc)
        vc_secs = time.perf_counter() - t0
    else:
        vc_secs = None  # reused saved reference — timing not available

    # ── mafsmith ──────────────────────────────────────────────────────────────
    ms_maf = os.path.join(work_dir, "mafsmith.maf")
    cmd_ms = [MAFSMITH, "vcf2maf",
              "-i", fv_out, "-o", ms_maf,
              "--vcf-tumor-id", tumor_id,
              "--tumor-id",     tumor_id,
              "--skip-annotation"]
    if normal_id:
        cmd_ms += ["--vcf-normal-id", normal_id, "--normal-id", normal_id]
    t0 = time.perf_counter()
    run(cmd_ms)
    ms_secs = time.perf_counter() - t0

    # ── Compare ───────────────────────────────────────────────────────────────
    cmp = compare_mafs(ms_maf, str(saved_ref))

    return {
        "status":        "ok",
        "syn_id":        syn_id,
        "vcf_name":      vcf_name,
        "tumor_id":      tumor_id,
        "normal_id":     normal_id,
        "n_vcf_records": n_variants,
        "fastvep_secs":  fastvep_secs,
        "vc_secs":       vc_secs,
        "ms_secs":       ms_secs,
        **cmp,
    }


# ── Report generation ──────────────────────────────────────────────────────────

def print_report(results):
    ok      = [r for r in results if r.get("status") == "ok"]
    skipped = [r for r in results if r.get("status") == "skipped"]
    errors  = [r for r in results if r.get("status") == "error"]

    print("\n" + "=" * 70)
    print(f"NF PORTAL VCF BENCHMARK  —  {len(results)} files processed")
    print("=" * 70)
    print(f"OK: {len(ok)}  Skipped: {len(skipped)}  Errors: {len(errors)}")

    if not ok:
        return

    # Timing
    total_variants   = sum(r.get("total_ms_rows", 0) for r in ok)
    total_vcf_vars   = sum(r.get("n_vcf_records", 0) for r in ok)
    total_ms_secs    = sum(r.get("ms_secs", 0) for r in ok)
    total_fv_secs    = sum(r.get("fastvep_secs", 0) for r in ok)
    vc_timed         = [r for r in ok if r.get("vc_secs") is not None]
    total_vc_secs    = sum(r.get("vc_secs", 0) for r in vc_timed)

    print(f"\n{'─'*70}")
    print(f"TIMING  (conversion step only, VEP pre-run via fastVEP for both)")
    print(f"{'─'*70}")
    print(f"  VCF records processed:           {total_vcf_vars:>12,}")
    print(f"  MAF rows output (mafsmith):      {total_variants:>12,}")
    print(f"")
    print(f"  fastVEP total time:              {total_fv_secs:>10.1f}s")
    if total_fv_secs > 0:
        print(f"  fastVEP variants/sec:            {total_vcf_vars/total_fv_secs:>10,.0f}")
    print(f"")
    if vc_timed:
        print(f"  vcf2maf.pl total time:           {total_vc_secs:>10.1f}s  ({len(vc_timed)} files timed)")
        print(f"  vcf2maf.pl variants/sec:         {total_vcf_vars/total_vc_secs*len(ok)/len(vc_timed):>10,.0f}  (extrapolated)")
    else:
        print(f"  vcf2maf.pl: all runs used saved references (re-run with --force-vcf2maf to time)")
    print(f"")
    print(f"  mafsmith total time:             {total_ms_secs:>10.1f}s")
    if total_ms_secs > 0:
        print(f"  mafsmith variants/sec:           {total_variants/total_ms_secs:>10,.0f}")
    if vc_timed and total_ms_secs > 0:
        vc_extrap = total_vc_secs * len(ok) / len(vc_timed)
        speedup = vc_extrap / total_ms_secs
        print(f"")
        print(f"  Speedup (conversion):            {speedup:>10.1f}×  (mafsmith vs vcf2maf.pl)")
        fv_plus_ms = total_fv_secs + total_ms_secs
        fv_plus_vc = total_fv_secs + vc_extrap
        if fv_plus_vc > 0:
            print(f"  Speedup (full pipeline):         {fv_plus_vc/fv_plus_ms:>10.1f}×  (fastVEP+mafsmith vs fastVEP+vcf2maf.pl)")

    # Concordance
    total_field_mm   = sum(r.get("field_mismatches", 0) for r in ok)
    total_variant_mm = sum(r.get("variant_mismatches", 0) for r in ok)
    total_sv_sec_vc  = sum(r.get("sv_secondary_vc", 0) for r in ok)
    total_sv_sec_ms  = sum(r.get("sv_secondary_ms", 0) for r in ok)
    total_compared   = sum(r.get("total_vc_rows", 0) for r in ok)
    variant_conc_pct = (total_compared - total_variant_mm) / total_compared * 100 if total_compared else 0
    total_fields     = total_compared * len(COMPARE_FIELDS)
    field_conc_pct   = (total_fields - total_field_mm) / total_fields * 100 if total_fields else 0

    print(f"\n{'─'*70}")
    print(f"CONCORDANCE  (mafsmith vs vcf2maf.pl, {len(COMPARE_FIELDS)} fields per variant)")
    print(f"{'─'*70}")
    print(f"  Variants compared:               {total_compared:>12,}")
    print(f"  Variants perfectly concordant:   {total_compared-total_variant_mm:>12,}")
    print(f"  Variants with any difference:    {total_variant_mm:>12,}")
    print(f"  Variant concordance rate:        {variant_conc_pct:>11.4f}%  (rows with 0 field diffs)")
    print(f"  Field concordance rate:          {field_conc_pct:>11.4f}%  (individual field values)")
    print(f"  Total field mismatches:          {total_field_mm:>12,}")
    print(f"  SV secondary rows (vcf2maf bug): {total_sv_sec_vc:>12,}  (skipped: empty Chromosome)")
    print(f"  SV secondary rows (mafsmith):    {total_sv_sec_ms:>12,}  (skipped: correctly-keyed, no vcf2maf match)")

    if total_variant_mm > 0:
        top_mm = sorted(ok, key=lambda r: r.get("variant_mismatches", 0), reverse=True)[:5]
        print(f"\n  Top files with variant-level mismatches:")
        for r in top_mm:
            if r.get("variant_mismatches", 0) > 0:
                pct = r.get("variant_concordance_pct", 0)
                print(f"    [{r['syn_id']}] {r['vcf_name'][:50]}  {r['variant_mismatches']} variants ({pct:.2f}%)")

    # Variant type distribution
    all_vt: dict = collections.Counter()
    all_vc_dist: dict = collections.Counter()
    for r in ok:
        for k, v in r.get("variant_types", {}).items():
            if k in ("SNP","DNP","TNP","ONP","INS","DEL"):
                all_vt[k] += v
            else:
                all_vc_dist[k] += v

    print(f"\n{'─'*70}")
    print(f"VARIANT TYPE DISTRIBUTION  (mafsmith output)")
    print(f"{'─'*70}")
    for vt, cnt in sorted(all_vt.items(), key=lambda x: -x[1]):
        print(f"  {vt:<12} {cnt:>10,}  ({cnt/total_variants*100:.1f}%)")

    print(f"\n  Variant_Classification (top 15):")
    for vc, cnt in sorted(all_vc_dist.items(), key=lambda x: -x[1])[:15]:
        print(f"  {vc:<30} {cnt:>10,}  ({cnt/total_variants*100:.1f}%)")

    print(f"\n{'─'*70}")

    # Per-file table
    print(f"\nPER-FILE RESULTS  (all {len(ok)} successful runs):")
    print(f"{'Syn_ID':<14} {'Variants':>9} {'fastVEP':>8} {'mafsmith':>9} {'vcf2maf':>8} {'Conc%':>8}  File")
    print(f"{'-'*14} {'-'*9} {'-'*8} {'-'*9} {'-'*8} {'-'*8}  {'-'*40}")
    for r in ok:
        vc_str = f"{r['vc_secs']:.2f}s" if r.get("vc_secs") else "  cached"
        print(f"{r['syn_id']:<14} {r['n_vcf_records']:>9,} "
              f"{r['fastvep_secs']:>7.2f}s "
              f"{r['ms_secs']:>8.2f}s "
              f"{vc_str:>8}  "
              f"{r.get('variant_concordance_pct', 0):>7.4f}%  "
              f"{r['vcf_name'][:50]}")


# ── Main ────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--n-vcfs",      type=int, default=100,
                        help="Number of VCFs to benchmark (default: 100)")
    parser.add_argument("--results-dir", default=str(REPO / "results"),
                        help="Directory for JSON results output")
    parser.add_argument("--results-file",
                        help="Specific JSON output file (overrides --results-dir default)")
    parser.add_argument("--force-vcf2maf", action="store_true",
                        help="Re-run vcf2maf.pl even if saved reference exists")
    parser.add_argument("--report-only",
                        help="Load existing JSON results file and print report, skip processing")
    parser.add_argument("--skip-somatic", action="store_true",
                        help="Skip somatic VCFs (only process germline)")
    args = parser.parse_args()

    if args.report_only:
        with open(args.report_only) as f:
            data = json.load(f)
        print_report(data["results"])
        return

    os.makedirs(args.results_dir, exist_ok=True)
    from datetime import datetime
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_path = args.results_file or os.path.join(
        args.results_dir, f"benchmark_nf_{ts}.json"
    )

    syn = synapseclient.Synapse()
    syn.login(silent=True)

    print(f"Querying NF portal fileview ({NF_FILEVIEW}) for Strelka somatic VCFs...")
    # Prefer somatic SNV/indel VCFs over genome/variants VCFs (smaller, more relevant)
    query_parts = [
        f"SELECT id, name FROM {NF_FILEVIEW}",
        "WHERE (name LIKE '%.somatic_snvs.vcf.gz' OR name LIKE '%.somatic_indels.vcf.gz')",
        f"LIMIT {args.n_vcfs * 3}",  # fetch extra in case some are skipped
    ]
    results_df = syn.tableQuery(" ".join(query_parts)).asDataFrame()
    syn_ids = list(results_df["id"])[:args.n_vcfs * 2]
    vcf_names = dict(zip(results_df["id"], results_df["name"]))
    print(f"  Candidate VCFs: {len(syn_ids)}")

    results = []
    processed = 0
    start_wall = time.time()

    for syn_id in syn_ids:
        if processed >= args.n_vcfs:
            break
        vcf_name = vcf_names.get(syn_id, syn_id)
        print(f"\n[{syn_id}] {vcf_name}", flush=True)
        work_dir = tempfile.mkdtemp(dir="/tmp", prefix=f"mafsmith_bench_{syn_id}_")
        try:
            f = syn.get(syn_id, downloadFile=True, downloadLocation=work_dir)
            r = process_one(syn_id, vcf_name, f.path, work_dir,
                            force_vcf2maf=args.force_vcf2maf)
            results.append({"syn_id": syn_id, "vcf_name": vcf_name, **r})

            if r["status"] == "ok":
                processed += 1
                vmm = r.get("variant_mismatches", 0)
                vpct = r.get("variant_concordance_pct", 0)
                print(
                    f"    {r['n_vcf_records']:,} variants | "
                    f"fastVEP {r['fastvep_secs']:.1f}s | "
                    f"mafsmith {r['ms_secs']:.3f}s | "
                    f"{'✓' if vmm == 0 else '✗'} {vpct:.4f}% variant concordance",
                    flush=True,
                )
            else:
                print(f"    → {r['status'].upper()}: {r.get('reason','')}", flush=True)

        except Exception as e:
            print(f"    → ERROR: {e}", flush=True)
            results.append({"syn_id": syn_id, "vcf_name": vcf_name,
                             "status": "error", "reason": str(e)[:300]})
        finally:
            shutil.rmtree(work_dir, ignore_errors=True)
            shutil.rmtree(os.path.expanduser("~/.synapseCache"), ignore_errors=True)

    elapsed = time.time() - start_wall
    print(f"\nTotal wall time: {elapsed:.0f}s ({elapsed/60:.1f} min)")

    # Save JSON results
    payload = {
        "timestamp": ts,
        "n_requested": args.n_vcfs,
        "n_processed": processed,
        "wall_secs": elapsed,
        "results": results,
    }
    with open(results_path, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"Results saved to: {results_path}")

    print_report(results)
    sys.exit(0 if all(r.get("variant_mismatches", 0) == 0 for r in results if r.get("status") == "ok") else 1)


if __name__ == "__main__":
    main()
