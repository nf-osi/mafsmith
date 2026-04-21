#!/usr/bin/env python3
"""
vcf_diff.py — run mafsmith+VEP and vcf2maf.pl+VEP on the same VCF and
produce a structured diff focused purely on conversion logic differences.

Both pipelines use the same Ensembl VEP installation and cache, so
annotation differences are eliminated and only conversion bugs remain.

Input: a Synapse entity ID (auto-downloaded) or a local VCF / VCF.gz path.

Usage:
    python3 scripts/vcf_diff.py syn12345678
    python3 scripts/vcf_diff.py /path/to/file.vcf.gz
    python3 scripts/vcf_diff.py syn12345678 --max-variants 5000
    python3 scripts/vcf_diff.py syn12345678 --tumor-id TUMOR --normal-id NORMAL
    python3 scripts/vcf_diff.py syn12345678 --output-dir /tmp/my_diff
    python3 scripts/vcf_diff.py syn12345678 --no-vcf2maf  # mafsmith only

Outputs (saved to --output-dir, default /tmp/vcf_diff_<id>/):
    mafsmith.maf    mafsmith + VEP MAF
    vcf2maf.maf     vcf2maf.pl + VEP MAF
    diff.txt        structured per-variant diff for Claude
"""

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
MAFSMITH  = REPO / "target" / "release" / "mafsmith"
REF_FASTA = Path.home() / ".mafsmith" / "hg38" / "reference.fa"
VEP_CACHE = Path.home() / ".vep"

# Only these fields are compared. They reflect VCF→MAF conversion logic and
# should be identical regardless of which transcript/gene model was chosen.
CONVERSION_FIELDS = {
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "Start_Position",
    "End_Position",
    "t_depth",
    "t_ref_count",
    "t_alt_count",
    "n_depth",
    "n_ref_count",
    "n_alt_count",
    "Matched_Norm_Sample_Barcode",
    "Tumor_Sample_Barcode",
}

SV_ALTS = {"<BND>", "<DEL>", "<DUP>", "<INV>", "<TRA>"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd, env_extra=None, label=None):
    env = os.environ.copy()
    if env_extra:
        env.update(env_extra)
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        tag = label or " ".join(str(c) for c in cmd[:2])
        raise RuntimeError(
            f"{tag} failed (exit {r.returncode}):\n"
            f"  stderr: {r.stderr[-1500:]}"
        )
    return r


def read_vcf_header(path):
    """Return (sample_names: list[str], has_chr_prefix: bool)."""
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    samples, has_chr = [], None
    with opener(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("#CHROM"):
                cols = line.strip().split("\t")
                samples = cols[9:] if len(cols) > 9 else []
            elif not line.startswith("#"):
                has_chr = line.startswith("chr")
                break
    return samples, bool(has_chr)


def decompress_vcf(src, dst, max_variants=None):
    """Write a plain-text VCF to dst, optionally capping at max_variants."""
    opener = gzip.open if str(src).endswith((".gz", ".bgz")) else open
    n = 0
    with opener(src, "rt", errors="replace") as fin, open(dst, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
            else:
                if max_variants is not None and n >= max_variants:
                    break
                fout.write(line)
                n += 1
    return n


def find_tool(name):
    which = shutil.which(name)
    if which:
        return Path(which)
    for candidate in [
        Path.home() / "miniconda3/envs/vcf2maf-env/bin" / name,
        Path.home() / "anaconda3/envs/vcf2maf-env/bin" / name,
    ]:
        if candidate.exists():
            return candidate
    return None


def read_maf(path):
    rows, header = [], None
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            rows.append(dict(zip(header, fields)))
    return rows


def variant_key(row):
    return (
        row.get("Chromosome", ""),
        row.get("Start_Position", ""),
        row.get("Reference_Allele", ""),
        row.get("Tumor_Seq_Allele2", ""),
    )


# ---------------------------------------------------------------------------
# Comparison logic
# ---------------------------------------------------------------------------

def compare_mafs(ms_rows, vc_rows):
    """Return a structured diff dict comparing only CONVERSION_FIELDS."""
    ms_idx = defaultdict(list)
    for r in ms_rows:
        if r.get("Chromosome"):
            ms_idx[variant_key(r)].append(r)

    vc_idx = {}
    vc_empty_chr = 0
    for r in vc_rows:
        if not r.get("Chromosome"):
            vc_empty_chr += 1
            continue
        vc_idx[variant_key(r)] = r

    only_in_vc, only_in_ms, diffs_by_variant = [], [], []
    matched_keys = set()

    for k, vc_row in sorted(vc_idx.items()):
        if k not in ms_idx:
            only_in_vc.append((k, vc_row))
            continue
        matched_keys.add(k)

        ms_candidates = ms_idx[k]
        if len(ms_candidates) == 1:
            ms_row = ms_candidates[0]
        else:
            # Pick the candidate with fewest conversion-field differences.
            ms_row = min(ms_candidates, key=lambda r: sum(
                1 for c in CONVERSION_FIELDS if r.get(c, "") != vc_row.get(c, "")
            ))

        diffs = []
        for col in sorted(CONVERSION_FIELDS):
            ms_val = ms_row.get(col, "")
            vc_val = vc_row.get(col, "")
            if ms_val != vc_val:
                diffs.append({"field": col, "ms_val": ms_val, "vc_val": vc_val})

        if diffs:
            diffs_by_variant.append({"key": k, "diffs": diffs})

    # mafsmith-only rows (not SV secondary)
    for k in sorted(set(ms_idx) - matched_keys):
        for row in ms_idx[k]:
            if row.get("Tumor_Seq_Allele2", "") not in SV_ALTS:
                only_in_ms.append((k, row))

    return {
        "ms_total": len(ms_rows),
        "vc_total": len(vc_rows),
        "vc_empty_chr": vc_empty_chr,
        "matched": len(matched_keys),
        "only_in_vc": only_in_vc,
        "only_in_ms": only_in_ms,
        "diffs_by_variant": diffs_by_variant,
    }


# ---------------------------------------------------------------------------
# Report writer
# ---------------------------------------------------------------------------

def write_diff_report(diff, path, input_label, tumor_id, normal_id, run_desc):
    SEP  = "=" * 72
    SEP2 = "-" * 72

    lines = [
        SEP,
        " VCF Diff Report",
        SEP,
        f" Input:    {input_label}",
        f" Run:      {run_desc}",
        f" Tumor:    {tumor_id}" + (f"   Normal: {normal_id}" if normal_id else ""),
        f" Compared: {', '.join(sorted(CONVERSION_FIELDS))}",
        SEP, "",
    ]

    # Row counts
    lines += [
        "ROW COUNTS",
        f"  vcf2maf.pl rows:       {diff['vc_total']:,}",
        f"  mafsmith rows:         {diff['ms_total']:,}",
        f"  Matched (by key):      {diff['matched']:,}",
    ]
    if diff["vc_empty_chr"]:
        lines.append(f"  vcf2maf SV empty-chr:  {diff['vc_empty_chr']:,}  (known vcf2maf.pl bug — skipped)")
    lines.append("")

    # Missing row warnings
    if diff["only_in_vc"]:
        lines.append(f"  WARNING: {len(diff['only_in_vc'])} variant(s) in vcf2maf MISSING from mafsmith:")
        for k, _ in diff["only_in_vc"][:20]:
            lines.append(f"     {k[0]}:{k[1]}  {k[2]}>{k[3]}")
        if len(diff["only_in_vc"]) > 20:
            lines.append(f"     ... and {len(diff['only_in_vc']) - 20} more")
        lines.append("")

    if diff["only_in_ms"]:
        lines.append(f"  WARNING: {len(diff['only_in_ms'])} variant(s) in mafsmith MISSING from vcf2maf:")
        for k, _ in diff["only_in_ms"][:20]:
            lines.append(f"     {k[0]}:{k[1]}  {k[2]}>{k[3]}")
        if len(diff["only_in_ms"]) > 20:
            lines.append(f"     ... and {len(diff['only_in_ms']) - 20} more")
        lines.append("")

    # Field mismatch summary
    field_counts: dict[str, int] = defaultdict(int)
    for v in diff["diffs_by_variant"]:
        for d in v["diffs"]:
            field_counts[d["field"]] += 1

    matched = diff["matched"]
    pct = lambda n: f"{n / matched * 100:.1f}%" if matched else "—"

    lines += [
        "FIELD MISMATCH SUMMARY",
        f"  {'Field':<38} {'Count':>8}  {'%':>7}",
        "  " + SEP2,
    ]
    for field, count in sorted(field_counts.items(), key=lambda x: -x[1]):
        lines.append(f"  {field:<38} {count:>8}  {pct(count):>7}")
    if not field_counts:
        lines.append("  (no differences — tools agree on all conversion fields)")
    lines.append("")

    # Per-variant details
    n_diffs = len(diff["diffs_by_variant"])
    lines += [
        SEP,
        f" CONVERSION DIFFERENCES ({n_diffs} variant(s))",
        " Fields that depend on VCF parsing, not which transcript was chosen.",
        SEP, "",
    ]
    if diff["diffs_by_variant"]:
        for i, v in enumerate(diff["diffs_by_variant"], 1):
            k = v["key"]
            lines.append(f"[{i}/{n_diffs}]  {k[0]}:{k[1]}  {k[2]}>{k[3]}")
            for d in sorted(v["diffs"], key=lambda x: x["field"]):
                lines += [
                    f"  {d['field']}:",
                    f"    vcf2maf:   {d['vc_val']!r}",
                    f"    mafsmith:  {d['ms_val']!r}",
                ]
            lines.append("")
    else:
        lines += ["  (none — perfect agreement)", ""]

    report = "\n".join(lines) + "\n"
    with open(path, "w") as f:
        f.write(report)
    return report


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("input",
                   help="Synapse ID (syn...) or path to VCF / VCF.gz")
    p.add_argument("--tumor-id",     default=None,
                   help="Tumor sample name for MAF (auto-detected from header)")
    p.add_argument("--normal-id",    default=None,
                   help="Normal sample name for MAF (auto-detected from header)")
    p.add_argument("--vcf-tumor-id", default=None,
                   help="Tumor column name in VCF (if different from --tumor-id)")
    p.add_argument("--vcf-normal-id", default=None,
                   help="Normal column name in VCF (if different from --normal-id)")
    p.add_argument("--max-variants", type=int, default=None,
                   help="Subsample to first N variant lines (for quick tests)")
    p.add_argument("--output-dir",   type=Path, default=None,
                   help="Where to write MAFs and diff.txt (default: /tmp/vcf_diff_<id>/)")
    p.add_argument("--mafsmith",     type=Path, default=MAFSMITH)
    p.add_argument("--ref-fasta",    type=Path, default=REF_FASTA)
    p.add_argument("--vep-cache",    type=Path, default=VEP_CACHE)
    p.add_argument("--vep-forks",    type=int,  default=16,
                   help="VEP forks for both mafsmith and vcf2maf.pl [default: 16]")
    p.add_argument("--strict",        action="store_true",
                   help="Pass --strict to mafsmith (match vcf2maf.pl AD-count behavior)")
    p.add_argument("--no-vcf2maf",   action="store_true",
                   help="Skip vcf2maf.pl — run mafsmith only (no diff produced)")
    p.add_argument("--keep-work",    action="store_true",
                   help="Keep temp working directory after finishing")
    args = p.parse_args()

    # Validate tools / paths
    for attr, label in [("mafsmith", "mafsmith binary"),
                         ("ref_fasta", "reference FASTA")]:
        path = getattr(args, attr)
        if not path.exists():
            sys.exit(f"ERROR: {label} not found at {path}")

    vcf2maf = None if args.no_vcf2maf else find_tool("vcf2maf.pl")
    if not args.no_vcf2maf and not vcf2maf:
        sys.exit("ERROR: vcf2maf.pl not found. Install via conda or pass --no-vcf2maf.")
    vep = None if args.no_vcf2maf else find_tool("vep")

    # Resolve input
    input_label = args.input
    work_dir = Path(tempfile.mkdtemp(prefix="vcf_diff_work_"))

    try:
        if args.input.startswith("syn"):
            print(f"Downloading {args.input} from Synapse...", flush=True)
            import synapseclient
            syn = synapseclient.Synapse()
            syn.login(silent=True)
            ent = syn.get(args.input, downloadFile=True,
                          downloadLocation=str(work_dir))
            vcf_path = Path(ent.path)
            input_label = f"{args.input} ({vcf_path.name})"
            print(f"  → {vcf_path.name}", flush=True)
        else:
            vcf_path = Path(args.input).resolve()
            if not vcf_path.exists():
                sys.exit(f"ERROR: file not found: {vcf_path}")

        # Auto-detect sample IDs
        samples, has_chr = read_vcf_header(vcf_path)
        print(f"VCF samples: {samples}  chr-prefix: {has_chr}", flush=True)

        # Auto-detect tumor/normal order: if first sample looks like normal (name contains
        # "normal" or "germline") and second looks like tumor, swap them.
        if not args.vcf_tumor_id and len(samples) >= 2:
            s0, s1 = samples[0].lower(), samples[1].lower()
            is_normal_first = any(w in s0 for w in ("normal", "germline")) and any(w in s1 for w in ("tumor", "somatic"))
            detected_tumor, detected_normal = (samples[1], samples[0]) if is_normal_first else (samples[0], samples[1])
        else:
            detected_tumor  = samples[0] if samples else "TUMOR"
            detected_normal = samples[1] if len(samples) > 1 else ""
        vcf_tumor_id  = args.vcf_tumor_id  or detected_tumor
        vcf_normal_id = args.vcf_normal_id or detected_normal
        tumor_id  = args.tumor_id  or vcf_tumor_id
        normal_id = args.normal_id or vcf_normal_id

        print(f"Tumor:  vcf={vcf_tumor_id!r}  maf={tumor_id!r}", flush=True)
        if normal_id:
            print(f"Normal: vcf={vcf_normal_id!r}  maf={normal_id!r}", flush=True)

        # Output directory
        if args.output_dir:
            out_dir = args.output_dir
        else:
            tag = args.input.replace("/", "_").lstrip(".")
            out_dir = Path("/tmp") / f"vcf_diff_{tag}"
        out_dir.mkdir(parents=True, exist_ok=True)

        ms_maf    = out_dir / "mafsmith.maf"
        vc_maf    = out_dir / "vcf2maf.maf"
        diff_path = out_dir / "diff.txt"

        # Decompressed VCF for vcf2maf.pl (kept in its own subdir so VEP temp
        # files don't collide with anything else in work_dir)
        vc_work = work_dir / "vcf2maf_work"
        vc_work.mkdir()
        plain_vcf = vc_work / "input.vcf"

        print(f"\nPreparing plain VCF...", flush=True)
        n_variants = decompress_vcf(vcf_path, plain_vcf, args.max_variants)
        if args.max_variants:
            print(f"  Subsampled to first {n_variants:,} variants", flush=True)
        else:
            print(f"  {n_variants:,} variants", flush=True)

        # --- mafsmith + VEP ---
        print(f"\nRunning mafsmith + VEP (--vep-forks {args.vep_forks})...", flush=True)
        ms_cmd = [
            str(args.mafsmith), "vcf2maf",
            "-i", str(plain_vcf),
            "-o", str(ms_maf),
            "--vcf-tumor-id", vcf_tumor_id,
            "--tumor-id",     tumor_id,
            "--genome", "grch38",
            "--ref-fasta", str(args.ref_fasta),
            "--annotator", "vep",
            "--vep-data",  str(args.vep_cache),
            "--vep-forks", str(args.vep_forks),
        ]
        if normal_id:
            ms_cmd += ["--vcf-normal-id", vcf_normal_id, "--normal-id", normal_id]
        if vep:
            ms_cmd += ["--vep-path", str(vep)]
        if args.strict:
            ms_cmd += ["--strict"]
        run(ms_cmd, label="mafsmith")
        print(f"  → {ms_maf}", flush=True)

        if not args.no_vcf2maf and vcf2maf:
            # --- vcf2maf.pl + VEP ---
            print(f"\nRunning vcf2maf.pl + VEP 115 (--vep-forks {args.vep_forks})...",
                  flush=True)

            colocated_perl = vcf2maf.parent / "perl"
            perl = str(colocated_perl) if colocated_perl.exists() else "perl"

            vc_cmd = [
                perl, str(vcf2maf),
                "--input-vcf",     str(plain_vcf),
                "--output-maf",    str(vc_maf),
                "--tumor-id",      tumor_id,
                "--vcf-tumor-id",  vcf_tumor_id,
                "--ref-fasta",     str(args.ref_fasta),
                "--ncbi-build",    "GRCh38",
                "--cache-version", "115",
                "--vep-forks",     str(args.vep_forks),
                "--vep-data",      str(args.vep_cache),
            ]
            if normal_id:
                vc_cmd += ["--normal-id", normal_id, "--vcf-normal-id", vcf_normal_id]
            if vep:
                vc_cmd += ["--vep-path", str(vep.parent)]
            # Point samtools/tabix at the conda env so PATH order doesn't matter.
            for tool in ("samtools", "tabix"):
                tp = vcf2maf.parent / tool
                if tp.exists():
                    vc_cmd += [f"--{tool}-exec", str(tp)]
            run(vc_cmd, label="vcf2maf.pl")
            print(f"  → {vc_maf}", flush=True)

            # --- Compare ---
            print("\nComparing MAFs...", flush=True)
            ms_rows = read_maf(str(ms_maf))
            vc_rows = read_maf(str(vc_maf))
            diff = compare_mafs(ms_rows, vc_rows)

            run_desc = (
                f"mafsmith{'(--strict)' if args.strict else ''}+VEP vs vcf2maf.pl+VEP (--vep-forks {args.vep_forks})"
                + (f"; max-variants {args.max_variants:,}" if args.max_variants else "")
            )
            write_diff_report(diff, diff_path, input_label,
                              tumor_id, normal_id, run_desc)

            n_diffs = len(diff["diffs_by_variant"])

            print(f"\n{'='*60}")
            print(f"Matched rows:          {diff['matched']:,}")
            print(f"Conversion diffs:      {n_diffs:,}  ← likely fixable bugs")
            print(f"Missing in mafsmith:   {len(diff['only_in_vc']):,}")
            print(f"Missing in vcf2maf:    {len(diff['only_in_ms']):,}")
            print(f"\nOutputs written to {out_dir}/")
            print(f"  mafsmith.maf")
            print(f"  vcf2maf.maf")
            print(f"  diff.txt   ← feed this to Claude for edge-case analysis")
        else:
            print(f"\nOutputs written to {out_dir}/")
            print(f"  mafsmith.maf")

    finally:
        if args.keep_work:
            print(f"\nWork dir kept: {work_dir}", flush=True)
        else:
            shutil.rmtree(work_dir, ignore_errors=True)
        if args.input.startswith("syn"):
            shutil.rmtree(os.path.expanduser("~/.synapseCache"), ignore_errors=True)


if __name__ == "__main__":
    main()
