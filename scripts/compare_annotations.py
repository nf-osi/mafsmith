#!/usr/bin/env python3
"""
Compare CSQ annotations produced by fastVEP vs VEP on the same input VCF.

Usage:
    python3 scripts/compare_annotations.py \
        --vcf bench_vcfs/GIAB_HG008_somatic/HG008-T--HG008-N.mutect2.vcf.gz \
        --ref-fasta ~/.mafsmith/hg38/reference.fa \
        --gff3 ~/.mafsmith/hg38/genes.gff3 \
        --vep-cache ~/.vep \
        [--max-variants N]

Output: per-field mismatch summary + example mismatches.
Does NOT modify any source files.
"""

import argparse
import gzip
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

MAFSMITH_DATA = Path.home() / ".mafsmith"
DEFAULT_REF_FASTA = MAFSMITH_DATA / "hg38" / "reference.fa"
DEFAULT_GFF3 = MAFSMITH_DATA / "hg38" / "genes.gff3"
DEFAULT_VEP_CACHE = Path.home() / ".vep"
FASTVEP_BIN = MAFSMITH_DATA / "bin" / "fastvep"

# CSQ fields to compare (subset most relevant to mafsmith conversion)
COMPARE_FIELDS = [
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "CANONICAL",
    "STRAND",
]


def find_fastvep() -> Path:
    if FASTVEP_BIN.exists():
        return FASTVEP_BIN
    which = shutil.which("fastvep")
    if which:
        return Path(which)
    sys.exit("ERROR: fastvep not found. Run `mafsmith fetch` or add to PATH.")


def find_vep() -> Path:
    for candidate in [
        Path.home() / "miniconda3/envs/vcf2maf-env/bin/vep",
        Path.home() / "anaconda3/envs/vcf2maf-env/bin/vep",
        shutil.which("vep"),
    ]:
        if candidate and Path(candidate).exists():
            return Path(candidate)
    sys.exit("ERROR: vep not found. Install via conda or add to PATH.")


def decompress_vcf(src: Path, dst: Path, max_variants: int) -> None:
    opener = gzip.open if src.suffix == ".gz" else open
    n = 0
    with opener(src, "rt", errors="replace") as fin, open(dst, "w") as fout:
        for line in fin:
            fout.write(line)
            if not line.startswith("#"):
                n += 1
                if n >= max_variants:
                    break


def run_fastvep(fastvep: Path, vcf: Path, ref: Path, gff3: Path, out: Path) -> None:
    cmd = [str(fastvep), "annotate",
           "-i", str(vcf),
           "-o", str(out),
           "--fasta", str(ref),
           "--gff3", str(gff3),
           "--hgvs",
           "--output-format", "vcf"]
    r = subprocess.run(cmd, capture_output=True)
    if r.returncode != 0:
        print(f"fastVEP stderr:\n{r.stderr.decode()[-500:]}", file=sys.stderr)
        sys.exit(f"ERROR: fastVEP failed with exit {r.returncode}")


def run_vep(vep: Path, vcf: Path, ref: Path, cache: Path, out: Path) -> None:
    # VEP needs its co-located perl
    perl = vep.parent / "perl"
    perl_bin = str(perl) if perl.exists() else "perl"
    cmd = [perl_bin, str(vep),
           "--input_file", str(vcf),
           "--output_file", str(out),
           "--format", "vcf",
           "--vcf",
           "--offline",
           "--cache",
           "--dir_cache", str(cache),
           "--assembly", "GRCh38",
           "--hgvs",
           "--symbol",
           "--canonical",
           "--biotype",
           "--numbers",
           "--fork", "4",
           "--fasta", str(ref),
           "--force_overwrite",
           "--no_stats"]
    r = subprocess.run(cmd, capture_output=True)
    if r.returncode != 0:
        print(f"VEP stderr:\n{r.stderr.decode()[-500:]}", file=sys.stderr)
        sys.exit(f"ERROR: VEP failed with exit {r.returncode}")


def parse_vcf_csq(vcf_path: Path) -> dict:
    """Parse annotated VCF, return {(chrom, pos, ref, alt): csq_string}."""
    csq_fmt = []
    variants = {}
    with open(vcf_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("##INFO=<ID=CSQ"):
                m = re.search(r"Format: ([^\"]+)", line)
                if m:
                    csq_fmt = m.group(1).rstrip('"').split("|")
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            info = cols[7]
            csq = ""
            for part in info.split(";"):
                if part.startswith("CSQ="):
                    csq = part[4:]
                    break
            key = (chrom, pos, ref, alt)
            variants[key] = (csq, csq_fmt)
    return variants


def parse_transcripts(csq_str: str, csq_fmt: list) -> list:
    """Split multi-transcript CSQ into list of field dicts."""
    transcripts = []
    for entry in csq_str.split(","):
        fields = entry.split("|")
        d = {csq_fmt[i]: fields[i] if i < len(fields) else "" for i in range(len(csq_fmt))}
        transcripts.append(d)
    return transcripts


def pick_canonical(transcripts: list) -> dict:
    """Pick canonical transcript (CANONICAL=YES), else first."""
    for t in transcripts:
        if t.get("CANONICAL") == "YES":
            return t
    return transcripts[0] if transcripts else {}


def compare(fv_vcf: Path, vep_vcf: Path) -> None:
    fv_data  = parse_vcf_csq(fv_vcf)
    vep_data = parse_vcf_csq(vep_vcf)

    common = set(fv_data) & set(vep_data)
    only_fv = set(fv_data) - set(vep_data)
    only_vep = set(vep_data) - set(fv_data)

    print(f"\n=== Annotation comparison: fastVEP vs VEP ===")
    print(f"Variants in fastVEP output:  {len(fv_data):,}")
    print(f"Variants in VEP output:      {len(vep_data):,}")
    print(f"Variants in both:            {len(common):,}")
    if only_fv:
        print(f"Only in fastVEP:             {len(only_fv):,}")
    if only_vep:
        print(f"Only in VEP:                 {len(only_vep):,}")

    # Per-field mismatch counts (on canonical transcript)
    field_mismatches: Counter = Counter()
    examples: dict = defaultdict(list)

    for key in sorted(common):
        fv_csq, fv_fmt = fv_data[key]
        vep_csq, vep_fmt = vep_data[key]

        if not fv_csq or not vep_csq:
            continue

        fv_t  = pick_canonical(parse_transcripts(fv_csq, fv_fmt))
        vep_t = pick_canonical(parse_transcripts(vep_csq, vep_fmt))

        for field in COMPARE_FIELDS:
            fv_val  = fv_t.get(field, "")
            vep_val = vep_t.get(field, "")
            if fv_val != vep_val:
                field_mismatches[field] += 1
                if len(examples[field]) < 3:
                    examples[field].append((key, fv_val, vep_val))

    print(f"\n--- Per-field mismatch summary (canonical transcript) ---")
    print(f"{'Field':<25} {'Mismatches':>12} {'%':>8}")
    print("-" * 48)
    for field in COMPARE_FIELDS:
        n = field_mismatches.get(field, 0)
        pct = 100 * n / len(common) if common else 0.0
        marker = "  ← differs" if n > 0 else ""
        print(f"  {field:<23} {n:>12,} {pct:>7.1f}%{marker}")

    if examples:
        print(f"\n--- Example mismatches (up to 3 per field) ---")
        for field, exs in examples.items():
            print(f"\n  {field}:")
            for (chrom, pos, ref, alt), fv_val, vep_val in exs:
                print(f"    {chrom}:{pos} {ref}>{alt}")
                print(f"      fastVEP: {fv_val!r}")
                print(f"      VEP:     {vep_val!r}")
    else:
        print("\n✓ No mismatches found in compared fields.")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf",         type=Path, required=True)
    p.add_argument("--ref-fasta",   type=Path, default=DEFAULT_REF_FASTA)
    p.add_argument("--gff3",        type=Path, default=DEFAULT_GFF3)
    p.add_argument("--vep-cache",   type=Path, default=DEFAULT_VEP_CACHE)
    p.add_argument("--vep-path",    type=Path, default=None,
                   help="Path to VEP binary (auto-detected if absent)")
    p.add_argument("--max-variants", type=int, default=1000,
                   help="Compare only the first N variants [default: 1000]")
    args = p.parse_args()

    fastvep = find_fastvep()
    vep = args.vep_path or find_vep()

    if not args.ref_fasta.exists():
        sys.exit(f"ERROR: reference FASTA not found: {args.ref_fasta}")
    if not args.gff3.exists():
        sys.exit(f"ERROR: GFF3 not found: {args.gff3}")
    if not args.vep_cache.exists():
        sys.exit(f"ERROR: VEP cache not found: {args.vep_cache}")

    with tempfile.TemporaryDirectory(prefix="annot_cmp_") as td:
        td = Path(td)
        plain_vcf = td / "input.vcf"
        fv_out    = td / "fastvep.vcf"
        vep_out   = td / "vep.vcf"

        print(f"Extracting first {args.max_variants} variants from {args.vcf.name}...")
        decompress_vcf(args.vcf, plain_vcf, args.max_variants)
        n = sum(1 for l in open(plain_vcf) if not l.startswith("#"))
        print(f"  {n} variants extracted")

        print(f"Running fastVEP...")
        run_fastvep(fastvep, plain_vcf, args.ref_fasta, args.gff3, fv_out)
        print(f"  Done → {fv_out.name}")

        print(f"Running VEP 115...")
        run_vep(vep, plain_vcf, args.ref_fasta, args.vep_cache, vep_out)
        print(f"  Done → {vep_out.name}")

        compare(fv_out, vep_out)


if __name__ == "__main__":
    main()
