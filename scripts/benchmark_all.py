#!/usr/bin/env python3
"""
Full benchmark suite: mafsmith vs vcf2maf.pl across all benchmark datasets.

Covers:
  - GIAB germline (HG001–HG007): single-sample GRCh38 benchmarks
  - GIAB HG008 somatic: paired T/N (MuTect2, Strelka2 SNV, Strelka2 INDEL)
  - SEQC2 HCC1395: paired T/N (MuTect2, Strelka)

For each dataset, mafsmith is timed at both 1 core (RAYON_NUM_THREADS=1) and
the default (all cores), with --iterations runs each. vcf2maf.pl is run once.

Usage:
    python3 scripts/benchmark_all.py [OPTIONS]

Options:
    --bench-vcfs DIR     Root directory containing benchmark VCFs
                         (default: ~/bench_vcfs)
    --ref-fasta FILE     chr-prefixed GRCh38 FASTA for T/N and germline datasets
                         (default: ~/.mafsmith/hg38/reference.fa)
    --mafsmith PATH      mafsmith binary (default: target/release/mafsmith)
    --vcf2maf PATH       vcf2maf.pl binary (default: auto-detected from PATH)
    --iterations N       Timing iterations per tool/thread-count (default: 3)
    --no-vcf2maf         Skip vcf2maf.pl comparison
    --datasets LIST      Comma-separated subset to run: germline,hg008,seqc2
                         (default: all available)
    --output FILE        Write markdown results table to this file
                         (default: results/benchmark_<mode>_YYYYMMDD.md)
    --annotated          Run full annotation pipeline (fastVEP + mafsmith vs
                         VEP + vcf2maf.pl). Requires VEP cache at ~/.vep/
                         and chr-prefixed GFF3 at ~/.mafsmith/hg38/genes.gff3.
    --gff3 FILE          Chr-prefixed GFF3 for mafsmith annotated mode
                         (default: ~/.mafsmith/hg38/genes.gff3)
    --vep-path PATH      Dir containing VEP binary (auto-detected if absent)
    --vep-cache DIR      VEP cache directory (default: ~/.vep)
"""

import argparse
import gzip
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_BENCH_VCFS = Path.home() / "bench_vcfs"
DEFAULT_REF_FASTA   = Path.home() / ".mafsmith" / "hg38" / "reference.fa"
DEFAULT_MAFSMITH    = REPO_ROOT / "target" / "release" / "mafsmith"
# Chr-prefixed GFF3 for annotated mode (created from ~/.mafsmith/GRCh38/genes.gff3.gz)
DEFAULT_GFF3_CHR    = Path.home() / ".mafsmith" / "hg38" / "genes.gff3"
DEFAULT_VEP_CACHE   = Path.home() / ".vep"


@dataclass
class Dataset:
    group: str          # "germline", "hg008", "seqc2"
    label: str          # display name for table
    vcf: Path
    vcf_tumor_id: str
    tumor_id: str
    vcf_normal_id: str  # empty string for single-sample
    normal_id: str      # empty string for single-sample


def build_datasets(bench_vcfs: Path) -> list[Dataset]:
    g = bench_vcfs / "GIAB_germline"
    h = bench_vcfs / "GIAB_HG008_somatic"
    s = bench_vcfs / "SEQC2_HCC1395"

    germline_samples = [
        ("HG001 (NA12878)", "HG001"),
        ("HG002 (NA24385)", "HG002"),
        ("HG003 (NA24149)", "HG003"),
        ("HG004 (NA24143)", "HG004"),
        ("HG005 (NA24631)", "HG005"),
        ("HG006 (NA24694)", "HG006"),
        ("HG007 (NA24695)", "HG007"),
    ]
    datasets = []
    for label, sid in germline_samples:
        vcf = g / f"{sid}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
        datasets.append(Dataset("germline", label, vcf, sid, sid, "", ""))

    datasets += [
        Dataset("hg008", "HG008 MuTect2",
                h / "HG008-T--HG008-N.mutect2.vcf.gz",
                "HG008-T", "HG008-T", "HG008-N", "HG008-N"),
        Dataset("hg008", "HG008 Strelka2 SNV",
                h / "HG008-T--HG008-N.snv.strelka2.vcf.gz",
                "HG008-T", "HG008-T", "HG008-N", "HG008-N"),
        Dataset("hg008", "HG008 Strelka2 INDEL",
                h / "HG008-T--HG008-N.indel.strelka2.vcf.gz",
                "HG008-T", "HG008-T", "HG008-N", "HG008-N"),
        Dataset("seqc2", "SEQC2 MuTect2",
                s / "WGS_FD_1.bwa.muTect2.vcf.gz",
                "TUMOR", "HCC1395", "NORMAL", "HCC1395BL"),
        Dataset("seqc2", "SEQC2 Strelka",
                s / "WGS_FD_1.bwa.strelka.vcf.gz",
                "TUMOR", "HCC1395", "NORMAL", "HCC1395BL"),
    ]
    return datasets


def count_variants(vcf: Path) -> int:
    opener = gzip.open if vcf.suffix == ".gz" else open
    n = 0
    with opener(vcf, "rt", errors="replace") as f:
        for line in f:
            if not line.startswith("#"):
                n += 1
    return n


def decompress_vcf(src: Path, dst: Path) -> None:
    with gzip.open(src, "rt", errors="replace") as fin, open(dst, "w") as fout:
        for line in fin:
            fout.write(line)


def mafsmith_cmd(ms: Path, ds: Dataset, ref: Path, out: Path, threads: int,
                 annotated: bool = False, gff3: Optional[Path] = None) -> list:
    cmd = [str(ms), "vcf2maf",
           "-i", str(ds.vcf),
           "-o", str(out),
           "--vcf-tumor-id", ds.vcf_tumor_id,
           "--tumor-id", ds.tumor_id,
           "--genome", "grch38",
           "--ref-fasta", str(ref)]
    if not annotated:
        cmd += ["--skip-annotation"]
    elif gff3:
        cmd += ["--gff3", str(gff3)]
    if ds.vcf_normal_id:
        cmd += ["--vcf-normal-id", ds.vcf_normal_id, "--normal-id", ds.normal_id]
    return cmd


def vcf2maf_cmd(vcf2maf: Path, ds: Dataset, ref: Path, plain_vcf: Path, out: Path,
                annotated: bool = False, vep_path: Optional[Path] = None,
                vep_cache: Optional[Path] = None) -> list:
    # Prefer the perl co-located with vcf2maf.pl to avoid shebang resolving to system perl.
    colocated_perl = vcf2maf.parent / "perl"
    perl = str(colocated_perl) if colocated_perl.exists() else "perl"
    cmd = [perl, str(vcf2maf),
           "--input-vcf", str(plain_vcf),
           "--output-maf", str(out),
           "--tumor-id", ds.tumor_id,
           "--vcf-tumor-id", ds.vcf_tumor_id,
           "--ref-fasta", str(ref)]
    if not annotated:
        cmd += ["--inhibit-vep"]
    else:
        if vep_path:
            cmd += ["--vep-path", str(vep_path)]
        if vep_cache:
            cmd += ["--vep-data", str(vep_cache)]
        # Use offline cache (no internet required) and assembly GRCh38
        cmd += ["--ncbi-build", "GRCh38", "--cache-version", "115"]
    if ds.vcf_normal_id:
        cmd += ["--normal-id", ds.normal_id, "--vcf-normal-id", ds.vcf_normal_id]
    # Point samtools/tabix to the conda env bin so PATH ordering doesn't matter.
    # vcf2maf.pl accepts --samtools-exec and --tabix-exec (not --bgzip-exec).
    env_bin = vcf2maf.parent
    for tool in ("samtools", "tabix"):
        tool_path = env_bin / tool
        if tool_path.exists():
            cmd += [f"--{tool}-exec", str(tool_path)]
    return cmd


def time_cmd(cmd: list, env: Optional[dict] = None) -> float:
    full_env = os.environ.copy()
    if env:
        full_env.update(env)
    t0 = time.monotonic()
    r = subprocess.run(cmd, capture_output=True, env=full_env)
    t1 = time.monotonic()
    if r.returncode != 0:
        print(f"    WARNING: exit {r.returncode}: {r.stderr.decode()[-200:]}", file=sys.stderr)
    return t1 - t0


@dataclass
class Result:
    dataset: Dataset
    n_variants: int
    ms1_times: list = field(default_factory=list)   # 1-core runs
    ms16_times: list = field(default_factory=list)  # all-core runs
    vc_times: list = field(default_factory=list)    # vcf2maf.pl runs


def run_benchmark(ds: Dataset, ms: Path, vcf2maf: Optional[Path], ref: Path,
                  iterations: int, tmpdir: Path,
                  annotated: bool = False, gff3: Optional[Path] = None,
                  vep_path: Optional[Path] = None,
                  vep_cache: Optional[Path] = None) -> Result:
    out_ms  = tmpdir / "ms_out.maf"
    out_vc  = tmpdir / "vc_out.maf"
    plain   = tmpdir / "plain.vcf"

    print(f"  Counting variants...", end=" ", flush=True)
    n = count_variants(ds.vcf)
    print(f"{n:,}")

    res = Result(dataset=ds, n_variants=n)

    # mafsmith — 1 core
    print(f"  mafsmith 1-core  ({iterations} runs):", end="", flush=True)
    for i in range(iterations):
        t = time_cmd(mafsmith_cmd(ms, ds, ref, out_ms, 1, annotated=annotated, gff3=gff3),
                     env={"RAYON_NUM_THREADS": "1"})
        res.ms1_times.append(t)
        print(f"  {t:.3f}s", end="", flush=True)
    print(f"  → mean {statistics.mean(res.ms1_times):.3f}s")

    # mafsmith — all cores
    ncores = os.cpu_count() or 1
    print(f"  mafsmith {ncores}-core ({iterations} runs):", end="", flush=True)
    for i in range(iterations):
        t = time_cmd(mafsmith_cmd(ms, ds, ref, out_ms, ncores, annotated=annotated, gff3=gff3))
        res.ms16_times.append(t)
        print(f"  {t:.3f}s", end="", flush=True)
    print(f"  → mean {statistics.mean(res.ms16_times):.3f}s")

    # vcf2maf.pl — multiple runs (needs decompressed input)
    if vcf2maf:
        print(f"  vcf2maf.pl       ({iterations} runs):", end="", flush=True)
        decompress_vcf(ds.vcf, plain)
        for i in range(iterations):
            t = time_cmd(vcf2maf_cmd(vcf2maf, ds, ref, plain, out_vc,
                                     annotated=annotated, vep_path=vep_path,
                                     vep_cache=vep_cache))
            res.vc_times.append(t)
            print(f"  {t:.3f}s", end="", flush=True)
        plain.unlink(missing_ok=True)
        print(f"  → mean {statistics.mean(res.vc_times):.3f}s")

    for f in [out_ms, out_vc]:
        f.unlink(missing_ok=True)

    return res


def find_vep(override: Optional[Path]) -> Optional[Path]:
    if override:
        return override if override.exists() else None
    for candidate in [
        shutil.which("vep"),
        Path.home() / "miniconda3/envs/vcf2maf-env/bin/vep",
        Path.home() / "anaconda3/envs/vcf2maf-env/bin/vep",
    ]:
        if candidate and Path(candidate).exists():
            return Path(candidate)
    return None


def find_vcf2maf(override: Optional[Path]) -> Optional[Path]:
    if override:
        return override if override.exists() else None
    # Check PATH then common conda locations
    for candidate in [
        shutil.which("vcf2maf.pl"),
        Path.home() / "miniconda3/envs/vcf2maf-env/bin/vcf2maf.pl",
        Path.home() / "anaconda3/envs/vcf2maf-env/bin/vcf2maf.pl",
    ]:
        if candidate and Path(candidate).exists():
            return Path(candidate)
    return None


def format_table(results: list[Result], ncores: int) -> str:
    lines = []
    lines.append(f"| Dataset | Variants | ms 1-core (s) | ms {ncores}-core (s) | vcf2maf.pl (s) | Speedup (1-core) | Speedup ({ncores}-core) |")
    lines.append(f"|---------|----------|---------------|-----------------|----------------|-----------------|--------------|")

    for r in results:
        ms1  = statistics.mean(r.ms1_times)
        ms16 = statistics.mean(r.ms16_times)
        vc   = statistics.mean(r.vc_times) if r.vc_times else None

        sp1  = f"{vc/ms1:.1f}×"  if vc else "—"
        sp16 = f"{vc/ms16:.1f}×" if vc else "—"
        vc_s = f"{vc:.3f}"        if vc else "—"

        lines.append(
            f"| {r.dataset.label} | {r.n_variants:,} | {ms1:.3f} | {ms16:.3f}"
            f" | {vc_s} | {sp1} | {sp16} |"
        )

    # Mean row for speedups where vcf2maf ran
    have_vc = [r for r in results if r.vc_times]
    if have_vc:
        sp1_vals  = [statistics.mean(r.vc_times) / statistics.mean(r.ms1_times)  for r in have_vc]
        sp16_vals = [statistics.mean(r.vc_times) / statistics.mean(r.ms16_times) for r in have_vc]
        lines.append(
            f"| **Mean** | | | | | **{statistics.mean(sp1_vals):.1f}×** | **{statistics.mean(sp16_vals):.1f}×** |"
        )

    return "\n".join(lines)


def write_report(results: list[Result], ncores: int, path: Path,
                 ref: Path, mafsmith: Path, vcf2maf: Optional[Path]) -> None:
    today = date.today().isoformat()
    sections = []

    sections.append(f"# mafsmith benchmark results\n")
    sections.append(
        f"**Date:** {today}  \n"
        f"**Instance:** {os.uname().nodename}  \n"
        f"**mafsmith:** `{mafsmith}`  \n"
        f"**vcf2maf.pl:** `{vcf2maf or 'not run'}`  \n"
        f"**Reference:** `{ref}`  \n"
    )

    groups = {}
    for r in results:
        groups.setdefault(r.dataset.group, []).append(r)

    group_titles = {"germline": "## Germline (GIAB HG001–HG007)",
                    "hg008":    "## Somatic T/N — GIAB HG008",
                    "seqc2":    "## Somatic T/N — SEQC2 HCC1395"}

    for g, title in group_titles.items():
        if g not in groups:
            continue
        sections.append(f"\n{title}\n")
        sections.append(format_table(groups[g], ncores))

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(sections) + "\n")
    print(f"\nResults written to {path}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bench-vcfs",  type=Path, default=DEFAULT_BENCH_VCFS)
    p.add_argument("--ref-fasta",   type=Path, default=DEFAULT_REF_FASTA)
    p.add_argument("--mafsmith",    type=Path, default=DEFAULT_MAFSMITH)
    p.add_argument("--vcf2maf",     type=Path, default=None)
    p.add_argument("--iterations",  type=int,  default=3)
    p.add_argument("--no-vcf2maf",  action="store_true")
    p.add_argument("--datasets",    default="germline,hg008,seqc2",
                   help="Comma-separated subset: germline,hg008,seqc2")
    p.add_argument("--output",      type=Path, default=None)
    p.add_argument("--annotated",   action="store_true",
                   help="Run with annotation enabled (fastVEP for mafsmith, VEP for vcf2maf.pl)")
    p.add_argument("--gff3",        type=Path, default=DEFAULT_GFF3_CHR,
                   help="Chr-prefixed GFF3 for mafsmith annotation mode "
                        f"(default: {DEFAULT_GFF3_CHR})")
    p.add_argument("--vep-path",    type=Path, default=None,
                   help="Path to VEP binary dir (for annotated mode; auto-detected if absent)")
    p.add_argument("--vep-cache",   type=Path, default=DEFAULT_VEP_CACHE,
                   help=f"Path to VEP cache dir (default: {DEFAULT_VEP_CACHE})")
    args = p.parse_args()

    ms = args.mafsmith
    if not ms.exists():
        sys.exit(f"ERROR: mafsmith not found at {ms}\n  Run: cargo build --release")

    ref = args.ref_fasta
    if not ref.exists():
        sys.exit(f"ERROR: reference FASTA not found at {ref}")

    vcf2maf = None if args.no_vcf2maf else find_vcf2maf(args.vcf2maf)
    if not args.no_vcf2maf and not vcf2maf:
        print("WARNING: vcf2maf.pl not found — skipping vcf2maf comparison", file=sys.stderr)

    vep_bin_dir: Optional[Path] = None
    if args.annotated:
        vep_exec = find_vep(args.vep_path)
        if vep_exec:
            vep_bin_dir = vep_exec.parent
        else:
            print("WARNING: vep not found — vcf2maf.pl will try PATH", file=sys.stderr)
        if not args.gff3.exists():
            sys.exit(
                f"ERROR: GFF3 not found at {args.gff3}.\n"
                "  For annotated mode with chr-prefixed VCFs, create it with:\n"
                "  zcat ~/.mafsmith/GRCh38/genes.gff3.gz | "
                "sed 's/^##sequence-region   /##sequence-region   chr/; /^#/!s/^/chr/' "
                "> ~/.mafsmith/hg38/genes.gff3"
            )

    requested = set(args.datasets.split(","))
    all_datasets = build_datasets(args.bench_vcfs)
    datasets = [d for d in all_datasets
                if d.group in requested and d.vcf.exists()]

    if not datasets:
        sys.exit("ERROR: no datasets found — check --bench-vcfs path")

    missing = [d for d in all_datasets if d.group in requested and not d.vcf.exists()]
    if missing:
        print(f"WARNING: {len(missing)} dataset(s) not found, skipping:", file=sys.stderr)
        for d in missing:
            print(f"  {d.vcf}", file=sys.stderr)

    ncores = os.cpu_count() or 1
    mode_tag = "annotated" if args.annotated else "conversion-only"
    out_path = args.output or (REPO_ROOT / "results" / f"benchmark_{mode_tag}_{date.today().isoformat()}.md")

    print("=" * 65)
    print(f" mafsmith benchmark suite")
    print("=" * 65)
    print(f" Mode:        {mode_tag}")
    print(f" Datasets:    {', '.join(sorted(requested))}")
    print(f" Iterations:  {args.iterations}")
    print(f" Cores:       1 and {ncores}")
    print(f" mafsmith:    {ms}")
    print(f" vcf2maf.pl:  {vcf2maf or 'skipped'}")
    print(f" Reference:   {ref}")
    if args.annotated:
        print(f" GFF3:        {args.gff3}")
        print(f" VEP cache:   {args.vep_cache}")
    print(f" Output:      {out_path}")
    print("=" * 65)

    results = []
    tmpdir = Path(tempfile.mkdtemp(prefix="mafsmith_bench_"))
    try:
        for ds in datasets:
            print(f"\n[{ds.group}] {ds.label}")
            r = run_benchmark(
                ds, ms, vcf2maf, ref, args.iterations, tmpdir,
                annotated=args.annotated,
                gff3=args.gff3 if args.annotated else None,
                vep_path=vep_bin_dir if args.annotated else None,
                vep_cache=args.vep_cache if args.annotated else None,
            )
            results.append(r)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    print("\n" + "=" * 65)
    print(" SUMMARY")
    print("=" * 65)
    for r in results:
        ms1  = statistics.mean(r.ms1_times)
        ms16 = statistics.mean(r.ms16_times)
        vc_mean = statistics.mean(r.vc_times) if r.vc_times else None
        sp1  = f"{vc_mean/ms1:.1f}×"  if vc_mean else "—"
        sp16 = f"{vc_mean/ms16:.1f}×" if vc_mean else "—"
        print(f"  {r.dataset.label:<30}  1-core {ms1:.3f}s  {ncores}-core {ms16:.3f}s  speedup {sp1} / {sp16}")

    write_report(results, ncores, out_path, ref, ms, vcf2maf)


if __name__ == "__main__":
    main()
