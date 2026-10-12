#!/usr/bin/env python3
"""Build the generated inputs of the vcf2maf.pl parity suite (tests/fixtures/parity/).

The parity fixtures sit on a tiny synthetic genome so the reference FASTA can be
committed alongside them. This script writes two things:

  annotated.vcf   CSQ-annotated records from tests/fixtures/integration_annotated.vcf,
                  moved onto the synthetic chromosomes so the classification/HGVS
                  columns are covered too. SV records get CHR2 added when it is
                  missing (vcf2maf.pl reads only CHR2, never the BND ALT, so without
                  it its secondary breakpoint row has an empty Chromosome), and the
                  BND that carries its partner only in ALT is left out.
  reference.fa    One sequence per contig used by any parity VCF, random bases
  (+ .fai)        (fixed seed) overwritten with every record's REF allele, so the
                  flanking_bps column and REF checks are consistent.

The hand-written caller fixtures (*.vcf other than annotated.vcf) are inputs to
this script, not outputs. Re-run after editing any of them, then regenerate the
expected MAFs with scripts/generate_parity_expected.sh.
"""
import random
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PARITY = REPO / "tests/fixtures/parity"
SOURCE = REPO / "tests/fixtures/integration_annotated.vcf"
CONTIG_LEN = 4000

# Source chromosome -> synthetic chromosome for annotated.vcf.
# The source mixes "chr1" and "1" (its CHR2 values), so names are canonicalized first.
CHROM_MAP = {"chr21": "chr3", "chr1": "chr4", "chr2": "chr5"}


def canon(chrom):
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def build_annotated():
    lines = SOURCE.read_text().splitlines()
    header = [l for l in lines if l.startswith("##") and not l.startswith("##contig")]
    header = [l for l in header if not l.startswith("##reference")]
    colnames = next(l for l in lines if l.startswith("#CHROM")).split("\t")
    records = [l.split("\t") for l in lines if not l.startswith("#")]

    # Every coordinate a record refers to (POS, END, BND partner), per chromosome,
    # gets a new position 100 bp apart in source order.
    coords = {}

    def note(chrom, pos):
        coords.setdefault(canon(chrom), set()).add(int(pos))

    kept = []
    for r in records:
        info = dict(kv.split("=", 1) for kv in r[7].split(";") if "=" in kv)
        bnd = re.search(r"[\]\[]([^:\]\[]+):(\d+)[\]\[]", r[4])
        if info.get("SVTYPE") and "END" not in info:
            continue  # partner only in the ALT notation — see module docstring
        note(r[0], r[1])
        if "END" in info:
            note(info.get("CHR2", bnd.group(1) if bnd else r[0]), info["END"])
        kept.append(r)

    remap = {
        (chrom, p): (CHROM_MAP[chrom], 100 * (i + 1))
        for chrom, ps in coords.items()
        for i, p in enumerate(sorted(ps))
    }

    out = []
    for r in kept:
        r = list(r)
        chrom, pos = remap[(canon(r[0]), int(r[1]))]
        info_kv = r[7].split(";")
        info = dict(kv.split("=", 1) for kv in info_kv if "=" in kv)
        if "END" in info:
            # A BND ALT is rewritten to point at the remapped CHR2:END.
            bnd = re.search(r"([\]\[])([^:\]\[]+):(\d+)([\]\[])", r[4])
            src2 = info.get("CHR2", bnd.group(2) if bnd else r[0])
            chr2, end = remap[(canon(src2), int(info["END"]))]
            new_kv = []
            for kv in info_kv:
                k = kv.split("=", 1)[0]
                if k == "END":
                    kv = f"END={end}"
                elif k == "CHR2":
                    kv = f"CHR2={chr2}"
                new_kv.append(kv)
            if "CHR2" not in info:
                new_kv.insert(new_kv.index(f"END={end}"), f"CHR2={chr2}")
            r[7] = ";".join(new_kv)
            if bnd:
                r[4] = r[4].replace(bnd.group(0), f"{bnd.group(1)}{chr2}:{end}{bnd.group(4)}")
        r[0], r[1] = chrom, str(pos)
        out.append(r)

    contigs = sorted({c for c, _ in remap.values()})
    text = ["##fileformat=VCFv4.2"] if header[0] != "##fileformat=VCFv4.2" else []
    text += header[:1] + [f"##contig=<ID={c},length={CONTIG_LEN}>" for c in contigs] + header[1:]
    text.append("\t".join(colnames))
    text += ["\t".join(r) for r in out]
    (PARITY / "annotated.vcf").write_text("\n".join(text) + "\n")


def build_reference():
    rng = random.Random(20261008)
    seqs = {}
    for vcf in sorted(PARITY.glob("*.vcf")):
        for line in vcf.read_text().splitlines():
            if line.startswith("##contig"):
                cid = re.search(r"ID=([^,>]+)", line).group(1)
                seqs.setdefault(cid, [rng.choice("ACGT") for _ in range(CONTIG_LEN)])
    owner = {}
    for vcf in sorted(PARITY.glob("*.vcf")):
        for line in vcf.read_text().splitlines():
            if line.startswith("#"):
                continue
            chrom, pos, _, ref = line.split("\t")[:4]
            for i, base in enumerate(ref.upper()):
                key = (chrom, int(pos) + i)
                prev = owner.get(key)
                if prev and prev[0] != base:
                    sys.exit(f"REF conflict at {chrom}:{key[1]}: {prev[1]} vs {vcf.name}")
                owner[key] = (base, vcf.name)
                seqs[chrom][key[1] - 1] = base

    fa = PARITY / "reference.fa"
    fai = []
    offset = 0
    with fa.open("w") as f:
        for chrom in sorted(seqs):
            seq = "".join(seqs[chrom])
            hdr = f">{chrom}\n"
            f.write(hdr)
            offset += len(hdr)
            fai.append(f"{chrom}\t{len(seq)}\t{offset}\t60\t61")
            for i in range(0, len(seq), 60):
                chunk = seq[i : i + 60] + "\n"
                f.write(chunk)
                offset += len(chunk)
    (PARITY / "reference.fa.fai").write_text("\n".join(fai) + "\n")


if __name__ == "__main__":
    build_annotated()
    build_reference()
    print(f"Wrote {PARITY / 'annotated.vcf'} and {PARITY / 'reference.fa'}(.fai)")
