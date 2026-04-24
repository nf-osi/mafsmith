#!/usr/bin/env python3
"""
Compare mafsmith maf2maf vs maf2maf.pl on a MAF file.

mafsmith uses fastVEP for annotation; maf2maf.pl uses Ensembl VEP 112.
Comparison is on conversion fields (same set as vcf_diff.py), excluding
Variant_Classification (annotation-policy-dependent).

Usage:
  python3 maf2maf_diff.py input.maf \\
      --genome grch38 \\
      --output-dir /path/to/out \\
      --work-dir /tmp/maf2maf_work \\
      --mafsmith /path/to/mafsmith \\
      --vep-forks 4
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

CONDA_BIN = Path.home() / "miniconda3/envs/vcf2maf-env/bin"
VEP_SHARE_CANDIDATES = sorted(Path(CONDA_BIN).parent.glob("share/ensembl-vep-*"))
VEP_SHARE = VEP_SHARE_CANDIDATES[-1] if VEP_SHARE_CANDIDATES else None
VEP_BIN = CONDA_BIN / "vep"
VEP_DATA = Path.home() / ".vep"

REF_FASTA_BY_GENOME = {
    "grch38": Path.home() / ".mafsmith/GRCh38/reference.fa",
    "grch37": Path.home() / ".mafsmith/GRCh37/reference.fa",
}
NCBI_BUILD_BY_GENOME = {
    "grch38": "GRCh38",
    "grch37": "GRCh37",
}

# Conversion fields (exclude annotation-policy fields)
CONVERSION_FIELDS = {
    "Chromosome", "Start_Position", "End_Position", "Strand",
    "Variant_Type", "Reference_Allele",
    "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
    "t_depth", "t_ref_count", "t_alt_count",
    "n_depth", "n_ref_count", "n_alt_count",
}


def perl_env():
    env = os.environ.copy()
    env["PATH"] = f"{CONDA_BIN}:{env.get('PATH', '')}"
    if VEP_SHARE:
        env["PERL5LIB"] = f"{VEP_SHARE}:{VEP_SHARE}/modules:{env.get('PERL5LIB', '')}"
    return env


def parse_maf(path: Path):
    rows = []
    col_idx = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            if col_idx is None:
                col_idx = {c: i for i, c in enumerate(line.split("\t"))}
                continue
            cols = line.split("\t")
            row = {}
            for k, i in col_idx.items():
                row[k] = cols[i] if i < len(cols) else ""
            rows.append(row)
    return rows, col_idx


def maf_key(row):
    return (
        row.get("Chromosome", ""),
        row.get("Start_Position", ""),
        row.get("Reference_Allele", ""),
        row.get("Tumor_Seq_Allele2", ""),
        row.get("Tumor_Sample_Barcode", ""),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="Input MAF file")
    p.add_argument("--genome", default="grch38", choices=["grch38", "grch37"])
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--work-dir", type=Path, default=Path("/tmp/maf2maf_diff_work"))
    p.add_argument("--mafsmith", type=Path, required=True)
    p.add_argument("--ref-fasta", type=Path)
    p.add_argument("--vep-forks", type=int, default=4)
    args = p.parse_args()

    ref_fasta = args.ref_fasta or REF_FASTA_BY_GENOME.get(args.genome)
    ncbi_build = NCBI_BUILD_BY_GENOME.get(args.genome, "GRCh38")

    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    maf_path = Path(args.input)

    # Run mafsmith maf2maf (fastVEP)
    ms_maf = args.work_dir / "ms_maf2maf.maf"
    cmd = [str(args.mafsmith), "maf2maf",
           "-i", str(maf_path),
           "-o", str(ms_maf),
           "--genome", args.genome]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: mafsmith maf2maf failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    # Run maf2maf.pl (VEP 112)
    pl_maf = args.work_dir / "pl_maf2maf.maf"
    maf2maf_pl = CONDA_BIN / "maf2maf.pl"
    cmd = ["perl", str(maf2maf_pl),
           "--input-maf", str(maf_path),
           "--output-maf", str(pl_maf),
           "--ncbi-build", ncbi_build,
           "--cache-version", "112",
           "--ref-fasta", str(ref_fasta),
           "--vep-path", str(CONDA_BIN),
           "--vep-data", str(VEP_DATA),
           "--vep-forks", str(args.vep_forks)]
    result = subprocess.run(cmd, capture_output=True, text=True, env=perl_env())
    if result.returncode != 0:
        print(f"ERROR: maf2maf.pl failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    # Compare
    ms_rows, ms_cols = parse_maf(ms_maf)
    pl_rows, pl_cols = parse_maf(pl_maf)

    ms_by_key = {maf_key(r): r for r in ms_rows}
    pl_by_key = {maf_key(r): r for r in pl_rows}

    ms_only = sorted(set(ms_by_key) - set(pl_by_key))
    pl_only = sorted(set(pl_by_key) - set(ms_by_key))
    common = set(ms_by_key) & set(pl_by_key)

    # Field-level diffs on CONVERSION_FIELDS
    field_diffs = []
    cmp_fields = CONVERSION_FIELDS & (set(ms_cols or {}) | set(pl_cols or {}))
    for key in sorted(common):
        ms_r = ms_by_key[key]
        pl_r = pl_by_key[key]
        for field in sorted(cmp_fields):
            mv = ms_r.get(field, "")
            pv = pl_r.get(field, "")
            if mv != pv:
                field_diffs.append((key, field, mv, pv))

    n_diffs = len(ms_only) + len(pl_only) + len(field_diffs)

    diff_path = args.output_dir / "diff.txt"
    with open(diff_path, "w") as f:
        f.write(f"maf2maf COMPARISON\n")
        f.write(f"  mafsmith records:   {len(ms_rows)}\n")
        f.write(f"  maf2maf.pl records: {len(pl_rows)}\n")
        f.write(f"  In common:          {len(common)}\n")
        f.write(f"  mafsmith only:      {len(ms_only)}\n")
        f.write(f"  maf2maf.pl only:    {len(pl_only)}\n")
        f.write(f"  Field differences:  {len(field_diffs)}\n")
        f.write(f"\nCONVERSION DIFFERENCES ({n_diffs})\n")
        if field_diffs:
            f.write("\nField differences (first 30):\n")
            for key, fld, mv, pv in field_diffs[:30]:
                f.write(f"  {key[0]}:{key[1]} {key[2]}>{key[3]} | {fld}: ms={mv!r} pl={pv!r}\n")

    print(f"maf2maf: {len(ms_rows)} mafsmith, {len(pl_rows)} maf2maf.pl, "
          f"{len(common)} common, CONVERSION DIFFERENCES ({n_diffs})")
    return 0 if n_diffs == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
