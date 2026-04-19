#!/usr/bin/env python3
"""
Throughput benchmark: mafsmith vs vcf2maf.pl

Measures wall-clock variants/second for both tools on the same pre-annotated VCF,
isolating conversion speed from annotation. Also benchmarks the full annotation
pipeline including fastVEP.

Usage:
    python3 scripts/benchmark.py INPUT.vcf [OPTIONS]

Options:
    --tumor-id ID        VCF tumor sample column name
    --normal-id ID       VCF normal sample column name
    --ref-fasta FILE     Reference FASTA (required for vcf2maf.pl; used for SV lookup)
    --gff3 FILE          GFF3 annotation for fastVEP
    --fastvep PATH       Path to fastVEP binary
    --mafsmith PATH      Path to mafsmith binary
    --iterations N       Timing iterations per tool (default: 3)
    --conda-env NAME     Conda environment containing vcf2maf.pl (default: vcf2maf-env)
    --no-vcf2maf         Skip vcf2maf.pl comparison (only time mafsmith)
"""

import argparse
import gzip
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MAFSMITH = REPO_ROOT / "target" / "release" / "mafsmith"
DEFAULT_FASTVEP = Path.home() / ".mafsmith" / "bin" / "fastvep"
DEFAULT_REF_FASTA = Path.home() / ".mafsmith" / "GRCh38" / "reference.fa"
DEFAULT_GFF3 = Path.home() / ".mafsmith" / "GRCh38" / "genes.gff3"


def is_gzip(path: Path) -> bool:
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def open_vcf(path: Path):
    return gzip.open(path, "rt") if is_gzip(path) else open(path, "rt", errors="replace")


def detect_samples(vcf_path: Path) -> tuple[str, str]:
    with open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                cols = line.strip().split("\t")
                samples = cols[9:] if len(cols) > 9 else []
                return (samples[0] if samples else "TUMOR",
                        samples[1] if len(samples) > 1 else "NORMAL")
    return "TUMOR", "NORMAL"


def strip_chr_prefix(src: Path, dst: Path) -> int:
    """Write plain-text VCF with chr stripped from chromosome column. Returns variant count."""
    count = 0
    with open_vcf(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
            else:
                if line.startswith("chr"):
                    line = line[3:]
                fout.write(line)
                count += 1
    return count


def count_variants(maf_path: Path) -> int:
    n = 0
    with open(maf_path) as f:
        for line in f:
            if not line.startswith("#") and not line.startswith("Hugo_Symbol"):
                n += 1
    return n


def time_cmd(cmd: list, label: str, iterations: int, output_template: str) -> dict:
    """Run cmd `iterations` times, substituting {i} in output_template. Returns timing stats."""
    elapsed = []
    for i in range(iterations):
        out_path = output_template.format(i=i)
        full_cmd = [str(x).replace("{out}", out_path) for x in cmd]
        t0 = time.monotonic()
        result = subprocess.run(full_cmd, capture_output=True)
        t1 = time.monotonic()
        dur = t1 - t0
        elapsed.append(dur)
        status = "✓" if result.returncode == 0 else f"✗ (exit {result.returncode})"
        print(f"    {label} run {i+1}/{iterations}: {dur:.3f}s {status}")
        if result.returncode != 0 and result.stderr:
            print(f"      stderr: {result.stderr.decode()[-200:]}", file=sys.stderr)
    mean = sum(elapsed) / len(elapsed)
    return {"mean": mean, "times": elapsed, "success": all(t > 0 for t in elapsed)}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input_vcf", type=Path)
    p.add_argument("--tumor-id", default="")
    p.add_argument("--normal-id", default="")
    p.add_argument("--ref-fasta", type=Path, default=DEFAULT_REF_FASTA)
    p.add_argument("--gff3", type=Path, default=DEFAULT_GFF3)
    p.add_argument("--fastvep", type=Path, default=DEFAULT_FASTVEP)
    p.add_argument("--mafsmith", type=Path, default=DEFAULT_MAFSMITH)
    p.add_argument("--iterations", type=int, default=3)
    p.add_argument("--conda-env", default="vcf2maf-env")
    p.add_argument("--no-vcf2maf", action="store_true")
    args = p.parse_args()

    if not args.mafsmith.exists():
        sys.exit(f"ERROR: mafsmith binary not found at {args.mafsmith}\n  Run: cargo build --release")

    if not args.tumor_id or not args.normal_id:
        t, n = detect_samples(args.input_vcf)
        args.tumor_id = args.tumor_id or t
        args.normal_id = args.normal_id or n

    print("=" * 67)
    print(" mafsmith vs vcf2maf.pl throughput benchmark")
    print("=" * 67)
    print(f" Input:        {args.input_vcf}")
    print(f" Tumor ID:     {args.tumor_id}")
    print(f" Normal ID:    {args.normal_id}")
    print(f" Iterations:   {args.iterations}")
    print(f" mafsmith:     {args.mafsmith}")
    print("=" * 67)

    tmpdir = Path(tempfile.mkdtemp(prefix="mafsmith_bench_"))
    try:
        # ── Prepare VCF ─────────────────────────────────────────────────────────
        print("\n[1/4] Preparing VCF...")
        plain_vcf = tmpdir / "plain.vcf"
        variant_count = strip_chr_prefix(args.input_vcf, plain_vcf)
        print(f"    Variants in input: {variant_count:,}")

        # ── Annotate with fastVEP ────────────────────────────────────────────────
        annotated_vcf = tmpdir / "annotated.vcf"
        annot_time = None
        if args.fastvep.exists() and args.gff3.exists() and args.ref_fasta.exists():
            print("\n[2/4] Running fastVEP annotation (once, shared)...")
            t0 = time.monotonic()
            result = subprocess.run(
                [str(args.fastvep), "annotate",
                 "-i", str(plain_vcf),
                 "-o", str(annotated_vcf),
                 "--fasta", str(args.ref_fasta),
                 "--gff3", str(args.gff3),
                 "--hgvs",
                 "--output-format", "vcf"],
                capture_output=True)
            annot_time = time.monotonic() - t0
            if result.returncode != 0:
                print(f"    WARNING: fastVEP exited {result.returncode} — using unannotated VCF")
                annotated_vcf = plain_vcf
            else:
                ann_count = sum(1 for l in open(annotated_vcf) if not l.startswith("#"))
                print(f"    Annotation complete: {ann_count:,} variants in {annot_time:.2f}s "
                      f"({ann_count / annot_time:,.0f} variants/s)")
        else:
            print("\n[2/4] fastVEP not found — using unannotated VCF (all variants → Targeted_Region)")
            annotated_vcf = plain_vcf

        # ── Benchmark mafsmith ───────────────────────────────────────────────────
        print(f"\n[3/4] Benchmarking mafsmith (--skip-annotation, {args.iterations} runs)...")
        ms_cmd = [str(args.mafsmith), "vcf2maf",
                  "--input-vcf", str(annotated_vcf),
                  "--output-maf", "{out}",
                  "--genome", "grch38",
                  "--vcf-tumor-id", args.tumor_id, "--tumor-id", args.tumor_id,
                  "--vcf-normal-id", args.normal_id, "--normal-id", args.normal_id,
                  "--skip-annotation"]
        if args.ref_fasta.exists():
            ms_cmd += ["--ref-fasta", str(args.ref_fasta)]

        ms_stats = time_cmd(ms_cmd, "mafsmith", args.iterations,
                            str(tmpdir / "ms_{i}.maf"))
        ms_vars = count_variants(tmpdir / f"ms_{args.iterations - 1}.maf")

        # ── Benchmark vcf2maf.pl ─────────────────────────────────────────────────
        vc_stats = None
        vc_vars = None
        if not args.no_vcf2maf and args.ref_fasta.exists():
            print(f"\n[4/4] Benchmarking vcf2maf.pl (--inhibit-vep, {args.iterations} runs)...")
            vc_cmd = ["conda", "run", "-n", args.conda_env, "vcf2maf.pl",
                      "--input-vcf", str(annotated_vcf),
                      "--output-maf", "{out}",
                      "--inhibit-vep",
                      "--tumor-id", args.tumor_id, "--vcf-tumor-id", args.tumor_id,
                      "--normal-id", args.normal_id, "--vcf-normal-id", args.normal_id,
                      "--ref-fasta", str(args.ref_fasta)]
            vc_stats = time_cmd(vc_cmd, "vcf2maf.pl", args.iterations,
                                str(tmpdir / "vc_{i}.maf"))
            vc_last = tmpdir / f"vc_{args.iterations - 1}.maf"
            if vc_last.exists():
                vc_vars = count_variants(vc_last)
        elif args.no_vcf2maf:
            print("\n[4/4] vcf2maf.pl comparison skipped (--no-vcf2maf)")
        else:
            print("\n[4/4] vcf2maf.pl skipped (no --ref-fasta)")

        # ── Summary ──────────────────────────────────────────────────────────────
        print()
        print("=" * 67)
        print(" RESULTS")
        print("=" * 67)
        print(f" Input variants:     {variant_count:,}")
        if annot_time:
            print(f" fastVEP annotation: {annot_time:.2f}s  ({variant_count / annot_time:,.0f} variants/s)")
        print()
        print(f" {'Tool':<20} {'Mean time':>10}  {'Variants/s':>12}  {'MAF rows':>10}")
        print(f" {'-'*20} {'-'*10}  {'-'*12}  {'-'*10}")

        ms_mean = ms_stats["mean"]
        ms_vps = ms_vars / ms_mean if ms_mean > 0 else 0
        print(f" {'mafsmith':<20} {ms_mean:>9.3f}s  {ms_vps:>12,.0f}  {ms_vars:>10,}")

        if vc_stats:
            vc_mean = vc_stats["mean"]
            vc_vars = vc_vars or 0
            vc_vps = vc_vars / vc_mean if vc_mean > 0 else 0
            print(f" {'vcf2maf.pl':<20} {vc_mean:>9.3f}s  {vc_vps:>12,.0f}  {vc_vars:>10,}")
            if ms_mean > 0 and vc_mean > 0:
                speedup = vc_mean / ms_mean
                print()
                print(f" Speedup: mafsmith is {speedup:.1f}x faster than vcf2maf.pl (conversion only)")
                if annot_time:
                    ms_total = ms_mean + annot_time
                    vc_total = vc_mean + annot_time
                    full_speedup = vc_total / ms_total
                    print(f"          Full pipeline (incl. fastVEP): {full_speedup:.1f}x faster")

        print()
        print(" Individual run times:")
        print(f"   mafsmith:   {', '.join(f'{t:.3f}s' for t in ms_stats['times'])}")
        if vc_stats:
            print(f"   vcf2maf.pl: {', '.join(f'{t:.3f}s' for t in vc_stats['times'])}")

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
