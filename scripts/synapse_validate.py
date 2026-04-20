#!/usr/bin/env python3
"""
Validate mafsmith against vcf2maf.pl on Synapse VCFs.

Usage:
    # Compare the 4 previously-failing IDs (default):
    python3 scripts/synapse_validate.py

    # Compare specific IDs:
    python3 scripts/synapse_validate.py syn21296193 syn21296175

    # Run vcf2maf.pl even if a saved reference already exists:
    python3 scripts/synapse_validate.py --force-vcf2maf

vcf2maf.pl outputs are saved to tests/fixtures/expected/{syn_id}_vcf2maf.maf
so subsequent runs skip vcf2maf.pl entirely (which is very slow).
"""

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import synapseclient

REPO = Path(__file__).parent.parent
EXPECTED_DIR = REPO / "tests" / "fixtures" / "expected"

MAFSMITH  = str(REPO / "target" / "release" / "mafsmith")
FASTVEP   = os.path.expanduser("~/.mafsmith/bin/fastvep")
REF_FASTA = os.path.expanduser("~/.mafsmith/GRCh38/reference.fa")
GFF3      = os.path.expanduser("~/.mafsmith/GRCh38/genes.gff3")

DEFAULT_IDS = [
    "syn21296193",
    "syn21296175",
    "syn21296184",
    "syn31625234",
]

COMPARE_FIELDS = [
    "Hugo_Symbol",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "HGVSc",
    "HGVSp",
    "HGVSp_Short",
    "Transcript_ID",
    "Exon_Number",
    "t_depth",
    "t_ref_count",
    "t_alt_count",
    "n_depth",
    "n_ref_count",
    "n_alt_count",
]


def run(cmd, **kwargs):
    r = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if r.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(str(c) for c in cmd)}\n{r.stderr[-800:]}"
        )
    return r


def read_vcf_header(path):
    """Return (sample_names, has_chr_prefix, is_grch38)."""
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
    return samples, has_chr, is_grch38


def strip_chr(src, dst):
    opener = gzip.open if str(src).endswith(".gz") else open
    with opener(src, "rt") as f, open(dst, "w") as out:
        for line in f:
            if line.startswith("#"):
                out.write(line)
            elif line.startswith("chr"):
                out.write(line[3:])
            else:
                out.write(line)


def is_sv_only(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    has_svtype, has_snv = False, False
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) > 7:
                if "SVTYPE" in cols[7]:
                    has_svtype = True
                if not cols[4].startswith("<"):
                    has_snv = True
            if has_snv:
                break
    return has_svtype and not has_snv


def parse_maf(path, skip_empty_chr=False):
    """Parse MAF into {(chr, pos): [row_dict, ...]}.

    skip_empty_chr: drop rows where Chromosome is blank (vcf2maf.pl secondary-row bug).
    Returns (records_dict, skipped_count).
    """
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
    three_to_one = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "Ter": "*", "Sec": "U", "Pyl": "O", "Xaa": "X",
    }
    for three, one in three_to_one.items():
        s = s.replace(three, one)
    return s


HGVS_COLS = {"HGVSp", "HGVSp_Short", "HGVSc"}


def _field_diffs(row1, row2):
    """Count mismatched fields between two rows."""
    diffs = 0
    for field in COMPARE_FIELDS:
        v1 = row1.get(field, "")
        v2 = row2.get(field, "")
        if field in HGVS_COLS:
            v1, v2 = normalize_hgvsp(v1), normalize_hgvsp(v2)
        if v1 != v2:
            diffs += 1
    return diffs


def _best_ms_row(ms_rows, vc_row):
    """Pick the mafsmith row (from a list of candidates at the same key) that best
    matches the vcf2maf row.  Used to resolve collisions where mafsmith outputs both a
    primary row and a secondary SV row at the same (chr, pos) key."""
    return min(ms_rows, key=lambda r: _field_diffs(r, vc_row))


def compare_mafs(ms_path, vc_path):
    # Skip vcf2maf.pl rows with empty Chromosome — vcf2maf.pl has a known bug where
    # it emits secondary SV breakpoint rows with Chromosome="" and Start_Position="".
    # mafsmith correctly fills in the actual partner chr/pos for these rows.
    ms, _ = parse_maf(ms_path, skip_empty_chr=False)
    vc, sv_secondary_vc = parse_maf(vc_path, skip_empty_chr=True)

    # Count distinct primary rows in mafsmith (keys that have at least one row)
    total_ms_rows = sum(len(v) for v in ms.values())

    if not ms and not vc:
        return {"variants": 0, "mismatches": 0, "sv_secondary_vc": 0,
                "sv_secondary_ms": 0, "details": []}

    mismatches, details = 0, []
    sv_secondary_ms = 0
    matched_ms_keys: set = set()

    # Iterate over vcf2maf keys first to find matches in mafsmith.
    for key in sorted(vc):
        vc_rows = vc[key]
        vc_row = vc_rows[0]  # vcf2maf should have at most one row per key after filtering
        if key not in ms:
            details.append(f"MISSING in mafsmith: {key}")
            mismatches += 1
            continue
        matched_ms_keys.add(key)
        # For keys with multiple mafsmith rows (primary + secondary collision),
        # pick the best-matching row to compare against vcf2maf.
        ms_row = _best_ms_row(ms[key], vc_row)
        for field in COMPARE_FIELDS:
            v1 = ms_row.get(field, "")
            v2 = vc_row.get(field, "")
            if field in HGVS_COLS:
                v1, v2 = normalize_hgvsp(v1), normalize_hgvsp(v2)
            if v1 != v2:
                details.append(f"{key[0]}:{key[1]} {field}: ms={v1!r} vc={v2!r}")
                mismatches += 1

    # mafsmith rows not in vcf2maf: secondary SV rows (correctly keyed, no vcf2maf match)
    for key in sorted(set(ms) - matched_ms_keys):
        for row in ms[key]:
            tsa2 = row.get("Tumor_Seq_Allele2", "")
            sv_alts = {"<BND>", "<DEL>", "<DUP>", "<INV>", "<TRA>"}
            if tsa2 in sv_alts:
                sv_secondary_ms += 1
            else:
                details.append(f"MISSING in vcf2maf:  {key}")
                mismatches += 1

    return {
        "variants": total_ms_rows,
        "mismatches": mismatches,
        "sv_secondary_vc": sv_secondary_vc,
        "sv_secondary_ms": sv_secondary_ms,
        "details": details,
    }


def process_one(syn_id, vcf_path, work_dir, force_vcf2maf=False, max_variants=None):
    samples, has_chr, is_grch38 = read_vcf_header(vcf_path)
    if is_grch38 is False:
        return {"status": "skipped", "reason": "Not GRCh38"}
    if not samples:
        return {"status": "skipped", "reason": "No sample columns"}
    if is_sv_only(vcf_path):
        return {"status": "skipped", "reason": "SV-only VCF"}

    tumor_id  = samples[0]
    normal_id = samples[1] if len(samples) > 1 else None

    # fastVEP requires plain-text VCF; decompress gzip/bgzf first.
    # Also strip chr prefix if needed (GFF3/FASTA use bare chr names).
    # Optionally subsample to max_variants variant lines.
    plain_vcf = os.path.join(work_dir, "input_plain.vcf")
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt") as fin, open(plain_vcf, "w") as fout:
        variant_count = 0
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
            else:
                if max_variants is not None and variant_count >= max_variants:
                    break
                if has_chr and line.startswith("chr"):
                    line = line[3:]
                fout.write(line)
                variant_count += 1
    if max_variants is not None:
        print(f"    Subsampled to first {variant_count:,} variants", flush=True)
    vcf_for_fv = plain_vcf

    # Annotate with fastVEP
    fv_out = os.path.join(work_dir, "fastvep_annotated.vcf")
    run([FASTVEP, "annotate",
         "-i", vcf_for_fv, "-o", fv_out,
         "--fasta", REF_FASTA, "--gff3", GFF3,
         "--hgvs", "--output-format", "vcf"])

    # vcf2maf.pl reference — generate once and save; reuse on subsequent runs
    saved_ref = EXPECTED_DIR / f"{syn_id}_vcf2maf.maf"
    if not saved_ref.exists() or force_vcf2maf:
        print(f"    Running vcf2maf.pl → saving to {saved_ref.name} ...", flush=True)
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
        run(cmd_vc)
        print(f"    vcf2maf.pl done: {saved_ref}", flush=True)
    else:
        print(f"    Using saved vcf2maf.pl reference: {saved_ref.name}", flush=True)

    # Run mafsmith
    ms_maf = os.path.join(work_dir, "mafsmith.maf")
    cmd_ms = [MAFSMITH, "vcf2maf",
              "-i", fv_out, "-o", ms_maf,
              "--vcf-tumor-id", tumor_id,
              "--tumor-id",     tumor_id,
              "--skip-annotation", "--strict"]
    if normal_id:
        cmd_ms += ["--vcf-normal-id", normal_id, "--normal-id", normal_id]
    run(cmd_ms)

    result = compare_mafs(ms_maf, str(saved_ref))
    result["status"]    = "ok"
    result["tumor_id"]  = tumor_id
    result["normal_id"] = normal_id
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("syn_ids", nargs="*", default=DEFAULT_IDS,
                        help="Synapse IDs to validate (default: 4 previously-failing IDs)")
    parser.add_argument("--force-vcf2maf", action="store_true",
                        help="Re-run vcf2maf.pl even if a saved reference already exists")
    parser.add_argument("--max-variants", type=int, default=None,
                        help="Subsample to the first N variant lines (useful for large VCFs)")
    args = parser.parse_args()

    syn = synapseclient.Synapse()
    syn.login(silent=True)

    results = []
    for syn_id in args.syn_ids:
        print(f"\n[{syn_id}]", flush=True)
        work_dir = tempfile.mkdtemp(dir="/tmp", prefix=f"mafsmith_val_{syn_id}_")
        try:
            f = syn.get(syn_id, downloadFile=True, downloadLocation=work_dir)
            print(f"    {Path(f.path).name}", flush=True)
            r = process_one(syn_id, f.path, work_dir, force_vcf2maf=args.force_vcf2maf,
                            max_variants=args.max_variants)
            results.append({"syn_id": syn_id, **r})

            if r["status"] == "ok":
                mm = r["mismatches"]
                sv_vc = r.get("sv_secondary_vc", 0)
                sv_ms = r.get("sv_secondary_ms", 0)
                print(f"    → {'✓ OK' if mm == 0 else '✗ FAIL'}: "
                      f"{r['variants']} variants, {mm} mismatches "
                      f"(SV secondary rows: vcf2maf={sv_vc} skipped, mafsmith={sv_ms} skipped)",
                      flush=True)
                for d in r["details"][:10]:
                    print(f"      {d}", flush=True)
                if len(r["details"]) > 10:
                    print(f"      ... and {len(r['details']) - 10} more", flush=True)
            else:
                print(f"    → SKIPPED: {r.get('reason', '')}", flush=True)

        except Exception as e:
            print(f"    → ERROR: {e}", flush=True)
            results.append({"syn_id": syn_id, "status": "error", "reason": str(e)[:300]})
        finally:
            shutil.rmtree(work_dir, ignore_errors=True)
            # Purge entire Synapse cache to reclaim disk space between files.
            shutil.rmtree(os.path.expanduser("~/.synapseCache"), ignore_errors=True)

    print("\n" + "=" * 60)
    ok      = [r for r in results if r.get("status") == "ok"]
    skipped = [r for r in results if r.get("status") == "skipped"]
    errors  = [r for r in results if r.get("status") == "error"]
    print(f"OK: {len(ok)}  Skipped: {len(skipped)}  Errors: {len(errors)}")

    total_vars = sum(r.get("variants", 0) for r in ok)
    total_mm   = sum(r.get("mismatches", 0) for r in ok)
    print(f"Variants compared: {total_vars:,}")
    print(f"Mismatches:        {total_mm}")
    if total_vars:
        rate = (total_vars - total_mm) / total_vars * 100
        print(f"Match rate:        {rate:.4f}%")

    if total_mm > 0:
        print("\nAll mismatches:")
        for r in ok:
            for d in r.get("details", []):
                print(f"  [{r['syn_id']}] {d}")

    sys.exit(0 if total_mm == 0 and not errors else 1)


if __name__ == "__main__":
    main()
