#!/usr/bin/env python3
"""
Compare mafsmith vcf2vcf vs vcf2vcf.pl on a VCF file.

Compares normalized variant keys (CHROM, POS, REF, ALT) after filtering.
Both tools filter non-PASS variants; mafsmith also drops ref-only records.

Usage:
  python3 vcf2vcf_diff.py input.vcf \\
      --genome grch38 \\
      --output-dir /path/to/out \\
      --work-dir /tmp/vcf2vcf_work \\
      --mafsmith /path/to/mafsmith
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

CONDA_BIN = Path.home() / "miniconda3/envs/vcf2maf-env/bin"
BCFTOOLS = CONDA_BIN / "bcftools"
VEP_SHARE = sorted(Path(CONDA_BIN).parent.glob("share/ensembl-vep-*"))[-1] if list(Path(CONDA_BIN).parent.glob("share/ensembl-vep-*")) else None
REF_FASTA_BY_GENOME = {
    "grch38": Path.home() / ".mafsmith/GRCh38/reference.fa",
    "grch37": Path.home() / ".mafsmith/GRCh37/reference.fa",
    "grcm39": Path.home() / ".mafsmith/grcm39/reference.fa",
}


def perl_env():
    env = os.environ.copy()
    env["PATH"] = f"{CONDA_BIN}:{env.get('PATH', '')}"
    if VEP_SHARE:
        env["PERL5LIB"] = f"{VEP_SHARE}:{VEP_SHARE}/modules:{env.get('PERL5LIB', '')}"
    return env


def parse_vcf_keys(path: Path) -> list[tuple]:
    keys = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                continue
            keys.append((cols[0], int(cols[1]), cols[3].upper(), cols[4].upper()))
    return keys


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="Input VCF (or VCF.gz) file")
    p.add_argument("--genome", default="grch38", choices=["grch38", "grch37", "grcm39"])
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--work-dir", type=Path, default=Path("/tmp/vcf2vcf_diff_work"))
    p.add_argument("--mafsmith", type=Path, required=True)
    p.add_argument("--ref-fasta", type=Path)
    p.add_argument("--vcf-tumor-id", default="")
    p.add_argument("--vcf-normal-id", default="")
    p.add_argument("--max-variants", type=int, default=2000, help="Subsample input to this many variants")
    args = p.parse_args()

    ref_fasta = args.ref_fasta or REF_FASTA_BY_GENOME.get(args.genome)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Possibly subsample the input VCF
    input_vcf = Path(args.input)
    if args.max_variants:
        sampled = args.work_dir / "input_sampled.vcf"
        # Use bcftools to sample, or just head the file
        cmd = [str(BCFTOOLS), "view", "-h", str(input_vcf)]
        header = subprocess.run(cmd, capture_output=True, text=True).stdout
        cmd2 = [str(BCFTOOLS), "view", "-H", str(input_vcf)]
        body = subprocess.run(cmd2, capture_output=True, text=True).stdout
        lines = [l for l in body.splitlines() if l][:args.max_variants]
        with open(sampled, "w") as f:
            f.write(header)
            f.write("\n".join(lines) + "\n")
        input_vcf = sampled

    # Run mafsmith vcf2vcf
    ms_vcf = args.work_dir / "ms_vcf2vcf.vcf"
    cmd = [str(args.mafsmith), "vcf2vcf",
           "-i", str(input_vcf),
           "-o", str(ms_vcf),
           "--genome", args.genome]
    if args.vcf_tumor_id:
        cmd += ["--vcf-tumor-id", args.vcf_tumor_id]
    if args.vcf_normal_id:
        cmd += ["--vcf-normal-id", args.vcf_normal_id]
    if ref_fasta:
        cmd += ["--ref-fasta", str(ref_fasta)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: mafsmith vcf2vcf failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    # Run vcf2vcf.pl
    pl_vcf = args.work_dir / "pl_vcf2vcf.vcf"
    vcf2vcf_pl = CONDA_BIN / "vcf2vcf.pl"
    cmd = ["perl", str(vcf2vcf_pl),
           "--input-vcf", str(input_vcf),
           "--output-vcf", str(pl_vcf)]
    if args.vcf_tumor_id:
        cmd += ["--vcf-tumor-id", args.vcf_tumor_id]
    if args.vcf_normal_id:
        cmd += ["--vcf-normal-id", args.vcf_normal_id]
    if ref_fasta:
        cmd += ["--ref-fasta", str(ref_fasta)]
    result = subprocess.run(cmd, capture_output=True, text=True, env=perl_env())
    if result.returncode != 0:
        print(f"ERROR: vcf2vcf.pl failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    # Compare variant keys
    ms_keys = set(parse_vcf_keys(ms_vcf))
    pl_keys = set(parse_vcf_keys(pl_vcf))

    ms_only = sorted(ms_keys - pl_keys)
    pl_only = sorted(pl_keys - ms_keys)
    common = ms_keys & pl_keys
    n_diffs = len(ms_only) + len(pl_only)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    diff_path = args.output_dir / "diff.txt"
    with open(diff_path, "w") as f:
        f.write(f"vcf2vcf COMPARISON\n")
        f.write(f"  mafsmith records:   {len(ms_keys)}\n")
        f.write(f"  vcf2vcf.pl records: {len(pl_keys)}\n")
        f.write(f"  In common:          {len(common)}\n")
        f.write(f"  mafsmith only:      {len(ms_only)}\n")
        f.write(f"  vcf2vcf.pl only:    {len(pl_only)}\n")
        f.write(f"\nCONVERSION DIFFERENCES ({n_diffs})\n")
        if ms_only:
            f.write("\nVariants in mafsmith only (first 20):\n")
            for k in ms_only[:20]:
                f.write(f"  {k[0]}:{k[1]} {k[2]}>{k[3]}\n")
        if pl_only:
            f.write("\nVariants in vcf2vcf.pl only (first 20):\n")
            for k in pl_only[:20]:
                f.write(f"  {k[0]}:{k[1]} {k[2]}>{k[3]}\n")

    print(f"vcf2vcf: {len(ms_keys)} mafsmith, {len(pl_keys)} vcf2vcf.pl, "
          f"{len(common)} common, CONVERSION DIFFERENCES ({n_diffs})")
    return 0 if n_diffs == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
