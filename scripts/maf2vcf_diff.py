#!/usr/bin/env python3
"""
Compare mafsmith maf2vcf vs maf2vcf.pl on a MAF file.

Usage:
  python3 maf2vcf_diff.py input.maf \\
      --genome grch38 \\
      --output-dir /path/to/out \\
      --work-dir /tmp/maf2vcf_work \\
      --mafsmith /path/to/mafsmith \\
      --ref-fasta /path/to/reference.fa

If --input-vcf is provided, the MAF is generated from the VCF first using
mafsmith vcf2maf --skip-annotation, and the generated MAF is used as input.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

CONDA_BIN = Path.home() / "miniconda3/envs/vcf2maf-env/bin"
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


def parse_vcf(path: Path) -> list[dict]:
    records = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            records.append({
                "chrom": cols[0],
                "pos": int(cols[1]),
                "ref": cols[3],
                "alt": cols[4],
                "format": cols[8],
                "samples": cols[9:],
            })
    return records


def vcf_key(rec):
    return (rec["chrom"], rec["pos"], rec["ref"], rec["alt"])


def run_mafsmith_maf2vcf(mafsmith, maf_path, out_vcf, genome, ref_fasta=None):
    cmd = [str(mafsmith), "maf2vcf",
           "-i", str(maf_path),
           "-o", str(out_vcf),
           "--genome", genome]
    if ref_fasta:
        cmd += ["--ref-fasta", str(ref_fasta)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"mafsmith maf2vcf failed:\n{result.stderr}")


def run_maf2vcf_pl(maf_path, out_vcf, out_dir, ref_fasta):
    maf2vcf_pl = CONDA_BIN / "maf2vcf.pl"
    cmd = ["perl", str(maf2vcf_pl),
           "--input-maf", str(maf_path),
           "--output-dir", str(out_dir),
           "--output-vcf", str(out_vcf),
           "--ref-fasta", str(ref_fasta)]
    result = subprocess.run(cmd, capture_output=True, text=True, env=perl_env())
    if result.returncode != 0:
        raise RuntimeError(f"maf2vcf.pl failed:\n{result.stderr}")
    return result.stderr  # may contain WARNING about skipped variants


def count_skipped(stderr_text: str) -> int:
    m = re.search(r"skipped.*?(\d+)", stderr_text, re.IGNORECASE)
    return int(m.group(1)) if m else 0


def compare_vcfs(ms_vcf: Path, pl_vcf: Path, out_dir: Path) -> dict:
    ms_recs = parse_vcf(ms_vcf)
    pl_recs = parse_vcf(pl_vcf)

    ms_by_key = {vcf_key(r): r for r in ms_recs}
    pl_by_key = {vcf_key(r): r for r in pl_recs}

    ms_only = sorted(set(ms_by_key) - set(pl_by_key))
    pl_only = sorted(set(pl_by_key) - set(ms_by_key))
    common = set(ms_by_key) & set(pl_by_key)

    # Check GT differences among common records
    gt_diffs = []
    for key in sorted(common):
        ms_r = ms_by_key[key]
        pl_r = pl_by_key[key]
        ms_gt = ms_r["samples"][0].split(":")[0] if ms_r["samples"] else "."
        pl_gt = pl_r["samples"][0].split(":")[0] if pl_r["samples"] else "."
        if ms_gt != pl_gt:
            gt_diffs.append((key, ms_gt, pl_gt))

    # Check AD:DP differences among common records (optional)
    depth_diffs = []
    for key in sorted(common):
        ms_r = ms_by_key[key]
        pl_r = pl_by_key[key]
        ms_vals = {k: v for k, v in zip(ms_r["format"].split(":"), (ms_r["samples"][0].split(":") if ms_r["samples"] else []))}
        pl_vals = {k: v for k, v in zip(pl_r["format"].split(":"), (pl_r["samples"][0].split(":") if pl_r["samples"] else []))}
        for field in ("AD", "DP"):
            mv = ms_vals.get(field, ".")
            pv = pl_vals.get(field, ".")
            if mv != pv and mv not in (".", ".,.") and pv not in (".", ".,."):
                depth_diffs.append((key, field, mv, pv))
                break

    # Write diff summary
    out_dir.mkdir(parents=True, exist_ok=True)
    diff_path = out_dir / "diff.txt"
    with open(diff_path, "w") as f:
        n_ms = len(ms_recs)
        n_pl = len(pl_recs)
        n_common = len(common)
        n_gt_diffs = len(gt_diffs)
        n_depth_diffs = len(depth_diffs)
        n_diffs = len(ms_only) + len(pl_only) + n_gt_diffs

        f.write(f"maf2vcf COMPARISON\n")
        f.write(f"  mafsmith records:   {n_ms}\n")
        f.write(f"  maf2vcf.pl records: {n_pl}\n")
        f.write(f"  In common:          {n_common}\n")
        f.write(f"  mafsmith only:      {len(ms_only)}\n")
        f.write(f"  maf2vcf.pl only:    {len(pl_only)}\n")
        f.write(f"  GT differences:     {n_gt_diffs}\n")
        f.write(f"  AD/DP differences:  {n_depth_diffs}\n")
        f.write(f"\nCONVERSION DIFFERENCES ({n_diffs})\n")
        if ms_only:
            f.write("\nVariants in mafsmith only (first 20):\n")
            for k in ms_only[:20]:
                f.write(f"  {k[0]}:{k[1]} {k[2]}>{k[3]}\n")
        if pl_only:
            f.write("\nVariants in maf2vcf.pl only (first 20):\n")
            for k in pl_only[:20]:
                f.write(f"  {k[0]}:{k[1]} {k[2]}>{k[3]}\n")
        if gt_diffs:
            f.write("\nGT differences (first 20):\n")
            for k, mg, pg in gt_diffs[:20]:
                f.write(f"  {k[0]}:{k[1]} {k[2]}>{k[3]}: mafsmith={mg} pl={pg}\n")

    return {
        "ms_records": n_ms,
        "pl_records": n_pl,
        "common": n_common,
        "ms_only": len(ms_only),
        "pl_only": len(pl_only),
        "gt_diffs": n_gt_diffs,
        "depth_diffs": n_depth_diffs,
        "total_diffs": n_diffs,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="Input MAF file (or Synapse ID to download)")
    p.add_argument("--input-vcf", help="If set, generate MAF from this VCF first (mafsmith vcf2maf --skip-annotation)")
    p.add_argument("--genome", default="grch38", choices=["grch38", "grch37", "grcm39"])
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--work-dir", type=Path, default=Path("/tmp/maf2vcf_diff_work"))
    p.add_argument("--mafsmith", type=Path, required=True)
    p.add_argument("--ref-fasta", type=Path)
    p.add_argument("--max-variants", type=int, default=2000)
    p.add_argument("--vcf-tumor-id")
    p.add_argument("--vcf-normal-id")
    p.add_argument("--tumor-id")
    p.add_argument("--normal-id")
    args = p.parse_args()

    ref_fasta = args.ref_fasta or REF_FASTA_BY_GENOME.get(args.genome)
    if not ref_fasta or not Path(ref_fasta).exists():
        print(f"ERROR: reference FASTA not found at {ref_fasta}", file=sys.stderr)
        sys.exit(1)

    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    maf_path = Path(args.input)

    # If an input VCF is provided, generate MAF first
    if args.input_vcf:
        maf_path = args.work_dir / "input_noanno.maf"
        cmd = [str(args.mafsmith), "vcf2maf",
               "-i", args.input_vcf,
               "-o", str(maf_path),
               "--genome", args.genome,
               "--skip-annotation"]
        if args.vcf_tumor_id:
            cmd += ["--vcf-tumor-id", args.vcf_tumor_id]
        if args.vcf_normal_id:
            cmd += ["--vcf-normal-id", args.vcf_normal_id]
        if args.tumor_id:
            cmd += ["--tumor-id", args.tumor_id]
        if args.normal_id:
            cmd += ["--normal-id", args.normal_id]
        subprocess.run(cmd, check=True, capture_output=True)

        # Optionally trim to max_variants
        if args.max_variants:
            trimmed = args.work_dir / "input_trimmed.maf"
            with open(maf_path) as fin, open(trimmed, "w") as fout:
                header_written = False
                count = 0
                for line in fin:
                    if line.startswith("#") or not header_written:
                        fout.write(line)
                        if not line.startswith("#"):
                            header_written = True
                        continue
                    if count >= args.max_variants:
                        break
                    fout.write(line)
                    count += 1
            maf_path = trimmed

    # Run mafsmith maf2vcf
    ms_vcf = args.work_dir / "ms_maf2vcf.vcf"
    run_mafsmith_maf2vcf(args.mafsmith, maf_path, ms_vcf, args.genome, ref_fasta)

    # Run maf2vcf.pl
    pl_out_dir = args.work_dir / "pl_maf2vcf_dir"
    pl_out_dir.mkdir(exist_ok=True)
    pl_vcf = args.work_dir / "pl_maf2vcf.vcf"
    pl_stderr = run_maf2vcf_pl(maf_path, pl_vcf, pl_out_dir, ref_fasta)

    # Compare
    stats = compare_vcfs(ms_vcf, pl_vcf, args.output_dir)

    # Print summary
    n_diffs = stats["total_diffs"]
    print(f"maf2vcf: {stats['ms_records']} mafsmith, {stats['pl_records']} maf2vcf.pl, "
          f"{stats['common']} common, CONVERSION DIFFERENCES ({n_diffs})")
    if pl_stderr and "WARNING" in pl_stderr:
        print(f"  [NOTE] maf2vcf.pl: {pl_stderr.strip().splitlines()[0][:120]}")

    return 0 if n_diffs == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
