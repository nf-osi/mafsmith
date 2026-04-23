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
VEP_CACHE = Path.home() / ".vep"

REF_FASTA_BY_GENOME = {
    "grch38": Path.home() / ".vep" / "homo_sapiens" / "112_GRCh38" / "Homo_sapiens.GRCh38.dna.toplevel.chr.fa.gz",
    "grch37": Path.home() / ".vep" / "homo_sapiens" / "112_GRCh37" / "Homo_sapiens.GRCh37.dna.toplevel.fa.gz",
    "grcm39": Path.home() / ".mafsmith" / "grcm39" / "reference.fa",
}
NCBI_BUILD_BY_GENOME = {
    "grch38": "GRCh38",
    "grch37": "GRCh37",
    "grcm39": "GRCm39",
}

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


def decompress_vcf(src, dst, max_variants=None, add_chr_prefix=False, seed=42,
                   skip_ref_blocks=False):
    """Write a plain-text VCF to dst, optionally reservoir-sampling max_variants.

    When max_variants is set, uses Knuth reservoir sampling (Algorithm R) so
    the selected variants are spread across the whole genome rather than
    biased towards early chromosomes. The sample is written in original
    genomic order (sorted by chrom/pos) so VEP/vcf2maf receive a valid VCF.

    add_chr_prefix: prepend 'chr' to chromosome names in data lines and
    ##contig headers (for VCFs using Ensembl-style names without prefix).

    skip_ref_blocks: drop GVCF reference-block records (ALT is '.' or a
    symbolic allele like '<NON_REF>'/'<*>') before passing to tools.
    Use for all-sites genome.vcf files where reference blocks vastly
    outnumber actual variants.
    """
    import random
    opener = gzip.open if str(src).endswith((".gz", ".bgz")) else open

    def _is_ref_block(line):
        parts = line.split("\t", 5)
        if len(parts) < 5:
            return False
        alt = parts[4]
        return alt == "." or alt.startswith("<")

    def _fix(line):
        if add_chr_prefix:
            line = "chr" + line
        return line

    def _chrom_key(line):
        parts = line.split("\t")
        chrom = parts[0].lstrip("chr") if parts else ""
        try:
            pos = int(parts[1]) if len(parts) > 1 else 0
        except (ValueError, IndexError):
            pos = 0
        _special = {"X": 23, "Y": 24, "MT": 25, "M": 25}
        try:
            return (0, int(chrom), pos)
        except ValueError:
            return (0, _special.get(chrom, 99), pos) if chrom in _special else (1, chrom, pos)

    header_lines = []
    rng = random.Random(seed)

    if max_variants is None:
        n = 0
        with opener(src, "rt", errors="replace") as fin, open(dst, "w") as fout:
            for line in fin:
                if line.startswith("#"):
                    if add_chr_prefix and line.startswith("##contig=<ID="):
                        line = line.replace("##contig=<ID=", "##contig=<ID=chr", 1)
                    fout.write(line)
                else:
                    if skip_ref_blocks and _is_ref_block(line):
                        continue
                    fout.write(_fix(line))
                    n += 1
        return n

    # Reservoir sampling — single pass, O(k) memory
    reservoir = []
    n_total = 0
    with opener(src, "rt", errors="replace") as fin:
        for line in fin:
            if line.startswith("#"):
                if add_chr_prefix and line.startswith("##contig=<ID="):
                    line = line.replace("##contig=<ID=", "##contig=<ID=chr", 1)
                header_lines.append(line)
            else:
                if skip_ref_blocks and _is_ref_block(line):
                    continue
                data = _fix(line)
                n_total += 1
                if len(reservoir) < max_variants:
                    reservoir.append(data)
                else:
                    j = rng.randint(0, n_total - 1)
                    if j < max_variants:
                        reservoir[j] = data

    # Sort by genomic position so VEP/vcf2maf receive a valid ordered VCF
    reservoir.sort(key=_chrom_key)

    with open(dst, "w") as fout:
        for line in header_lines:
            fout.write(line)
        for line in reservoir:
            fout.write(line)

    return len(reservoir)


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

def _iter_maf_sparse(path, need_fields):
    """Yield sparse row dicts (only need_fields + variant-key columns) from a MAF file."""
    need = set(need_fields) | {"Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"}
    col_idx = {}
    header_parsed = False
    with open(str(path)) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if not header_parsed:
                col_idx = {c: i for i, c in enumerate(cols) if c in need}
                header_parsed = True
                continue
            yield {c: cols[i] if i < len(cols) else "" for c, i in col_idx.items()}


def _emit_diff(ms_row, vc_row, cmp_fields, diffs_by_variant):
    diffs = [
        {"field": col, "ms_val": ms_row.get(col, ""), "vc_val": vc_row.get(col, "")}
        for col in sorted(cmp_fields)
        if ms_row.get(col, "") != vc_row.get(col, "")
    ]
    if diffs:
        diffs_by_variant.append({"key": variant_key(ms_row), "diffs": diffs})


def compare_mafs_streaming(ms_path, vc_path, cmp_fields=None):
    """Memory-efficient MAF comparison: streams both files simultaneously.

    Assumes both MAFs are in the same genomic order (always true when both
    tools processed the same sorted input VCF). Falls back to a pending-row
    dict only when order diverges (SV secondary rows, etc.) — normally empty.

    Returns the same structure as compare_mafs().
    """
    from itertools import zip_longest

    if cmp_fields is None:
        cmp_fields = CONVERSION_FIELDS

    ms_iter = _iter_maf_sparse(ms_path, cmp_fields)
    vc_iter = _iter_maf_sparse(vc_path, cmp_fields)

    ms_total = vc_total = vc_empty_chr = matched_count = 0
    only_in_vc = []
    only_in_ms_rows = []
    diffs_by_variant = []

    # Pending rows that didn't match in expected position (normally empty).
    ms_pending = defaultdict(list)  # key -> [row, ...]
    vc_pending = {}                 # key -> row

    _sentinel = object()

    for ms_raw, vc_raw in zip_longest(ms_iter, vc_iter, fillvalue=_sentinel):
        ms_row = None if ms_raw is _sentinel else ms_raw
        vc_row = None if vc_raw is _sentinel else vc_raw

        if ms_row is not None:
            ms_total += 1
        if vc_row is not None:
            vc_total += 1
            if not vc_row.get("Chromosome"):
                vc_empty_chr += 1
                vc_row = None

        ms_key = variant_key(ms_row) if ms_row else None
        vc_key = variant_key(vc_row) if vc_row else None

        if ms_key and vc_key and ms_key == vc_key:
            matched_count += 1
            _emit_diff(ms_row, vc_row, cmp_fields, diffs_by_variant)
            continue

        if ms_row:
            if ms_key in vc_pending:
                matched_count += 1
                _emit_diff(ms_row, vc_pending.pop(ms_key), cmp_fields, diffs_by_variant)
            else:
                ms_pending[ms_key].append(ms_row)

        if vc_row:
            if vc_key in ms_pending:
                ms_candidates = ms_pending[vc_key]
                ms_match = ms_candidates.pop(0)
                if not ms_candidates:
                    del ms_pending[vc_key]
                matched_count += 1
                _emit_diff(ms_match, vc_row, cmp_fields, diffs_by_variant)
            else:
                vc_pending[vc_key] = vc_row

    for k, rows in sorted(ms_pending.items()):
        for row in rows:
            if row.get("Tumor_Seq_Allele2", "") not in SV_ALTS:
                only_in_ms_rows.append((k, row))

    for k, row in sorted(vc_pending.items()):
        only_in_vc.append((k, row))

    return {
        "ms_total": ms_total,
        "vc_total": vc_total,
        "vc_empty_chr": vc_empty_chr,
        "matched": matched_count,
        "only_in_vc": only_in_vc,
        "only_in_ms": only_in_ms_rows,
        "diffs_by_variant": diffs_by_variant,
    }


def compare_mafs(ms_rows, vc_rows, cmp_fields=None):
    """Return a structured diff dict comparing only cmp_fields (default: CONVERSION_FIELDS)."""
    if cmp_fields is None:
        cmp_fields = CONVERSION_FIELDS

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
                1 for c in cmp_fields if r.get(c, "") != vc_row.get(c, "")
            ))

        diffs = []
        for col in sorted(cmp_fields):
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

def write_diff_report(diff, path, input_label, tumor_id, normal_id, run_desc, cmp_fields=None):
    if cmp_fields is None:
        cmp_fields = CONVERSION_FIELDS
    SEP  = "=" * 72
    SEP2 = "-" * 72

    lines = [
        SEP,
        " VCF Diff Report",
        SEP,
        f" Input:    {input_label}",
        f" Run:      {run_desc}",
        f" Tumor:    {tumor_id}" + (f"   Normal: {normal_id}" if normal_id else ""),
        f" Compared: {', '.join(sorted(cmp_fields))}",
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
    p.add_argument("--genome",        default="grch38",
                   choices=["grch38", "grch37", "grcm39"],
                   help="Reference genome assembly [default: grch38]")
    p.add_argument("--mafsmith",     type=Path, default=MAFSMITH)
    p.add_argument("--ref-fasta",    type=Path, default=None,
                   help="Reference FASTA (auto-detected from --genome if omitted)")
    p.add_argument("--vep-cache",    type=Path, default=VEP_CACHE)
    p.add_argument("--vep-forks",    type=int,  default=16,
                   help="VEP forks for both mafsmith and vcf2maf.pl [default: 16]")
    p.add_argument("--strict",        action="store_true",
                   help="Pass --strict to mafsmith (match vcf2maf.pl AD-count behavior)")
    p.add_argument("--inhibit-vep",  action="store_true",
                   help="Skip VEP annotation (mafsmith --skip-annotation / vcf2maf.pl --inhibit-vep); "
                        "Variant_Classification is excluded from the diff")
    p.add_argument("--no-vcf2maf",   action="store_true",
                   help="Skip vcf2maf.pl — run mafsmith only (no diff produced)")
    p.add_argument("--no-keep-mafs", action="store_true",
                   help="Delete MAF files after writing diff.txt (saves disk space in batch runs)")
    p.add_argument("--skip-ref-blocks", action="store_true",
                   help="Drop GVCF reference-block records (ALT='.' or symbolic <...>) before "
                        "running tools. Use for all-sites genome.vcf files to avoid processing "
                        "the ~95%% of records that are non-variant reference blocks.")
    p.add_argument("--keep-work",    action="store_true",
                   help="Keep temp working directory after finishing")
    p.add_argument("--work-dir",     type=Path, default=None,
                   help="Parent directory for temp working files (default: system /tmp). "
                        "Use a real disk path to avoid RAM-backed tmpfs pressure on "
                        "instances where /tmp is a tmpfs.")
    args = p.parse_args()

    # Resolve ref_fasta from genome if not given explicitly
    if args.ref_fasta is None:
        args.ref_fasta = REF_FASTA_BY_GENOME.get(args.genome)

    # Validate tools / paths
    if not args.mafsmith.exists():
        sys.exit(f"ERROR: mafsmith binary not found at {args.mafsmith}")
    if args.ref_fasta and not args.ref_fasta.exists():
        print(f"WARNING: ref FASTA not found at {args.ref_fasta} — running without --fasta", flush=True)
        args.ref_fasta = None

    vcf2maf = None if args.no_vcf2maf else find_tool("vcf2maf.pl")
    if not args.no_vcf2maf and not vcf2maf:
        sys.exit("ERROR: vcf2maf.pl not found. Install via conda or pass --no-vcf2maf.")
    vep = None if (args.no_vcf2maf or args.inhibit_vep) else find_tool("vep")

    # Resolve input
    input_label = args.input
    if args.work_dir:
        args.work_dir.mkdir(parents=True, exist_ok=True)
    work_dir = Path(tempfile.mkdtemp(prefix="vcf_diff_work_", dir=args.work_dir))

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

        # GRCh38/GRCm39 reference FASTAs use chr-prefix; add it if the VCF lacks it.
        add_chr = (not has_chr) and args.genome in ("grch38", "grcm39")
        if add_chr:
            print("  Note: adding chr prefix (VCF uses Ensembl-style names)", flush=True)

        print(f"\nPreparing plain VCF...", flush=True)
        n_variants = decompress_vcf(vcf_path, plain_vcf, args.max_variants,
                                    add_chr_prefix=add_chr,
                                    skip_ref_blocks=args.skip_ref_blocks)
        if args.max_variants:
            print(f"  Reservoir-sampled {n_variants:,} variants (genome-wide)", flush=True)
        elif args.skip_ref_blocks:
            print(f"  {n_variants:,} variants (reference blocks excluded)", flush=True)
        else:
            print(f"  {n_variants:,} variants", flush=True)

        ncbi_build = NCBI_BUILD_BY_GENOME.get(args.genome, "GRCh38")

        # Fields excluded in --inhibit-vep mode because they depend on annotation.
        cmp_fields = CONVERSION_FIELDS - {"Variant_Classification"} if args.inhibit_vep else CONVERSION_FIELDS

        # --- mafsmith ---
        if args.inhibit_vep:
            print(f"\nRunning mafsmith --skip-annotation...", flush=True)
            ms_cmd = [
                str(args.mafsmith), "vcf2maf",
                "-i", str(plain_vcf),
                "-o", str(ms_maf),
                "--vcf-tumor-id", vcf_tumor_id,
                "--tumor-id",     tumor_id,
                "--genome", args.genome,
                "--skip-annotation",
            ]
        else:
            print(f"\nRunning mafsmith + VEP (--vep-forks {args.vep_forks})...", flush=True)
            ms_cmd = [
                str(args.mafsmith), "vcf2maf",
                "-i", str(plain_vcf),
                "-o", str(ms_maf),
                "--vcf-tumor-id", vcf_tumor_id,
                "--tumor-id",     tumor_id,
                "--genome", args.genome,
                "--annotator", "vep",
                "--vep-data",  str(args.vep_cache),
                "--vep-forks", str(args.vep_forks),
            ]
        if args.ref_fasta and not args.inhibit_vep:
            ms_cmd += ["--ref-fasta", str(args.ref_fasta)]
        if normal_id:
            ms_cmd += ["--vcf-normal-id", vcf_normal_id, "--normal-id", normal_id]
        if not args.inhibit_vep and vep:
            ms_cmd += ["--vep-path", str(vep)]
        if args.strict:
            ms_cmd += ["--strict"]
        run(ms_cmd, label="mafsmith")
        print(f"  → {ms_maf}", flush=True)

        if not args.no_vcf2maf and vcf2maf:
            # --- vcf2maf.pl ---
            colocated_perl = vcf2maf.parent / "perl"
            perl = str(colocated_perl) if colocated_perl.exists() else "perl"

            if args.inhibit_vep:
                print(f"\nRunning vcf2maf.pl --inhibit-vep...", flush=True)
                vc_cmd = [
                    perl, str(vcf2maf),
                    "--input-vcf",    str(plain_vcf),
                    "--output-maf",   str(vc_maf),
                    "--tumor-id",     tumor_id,
                    "--vcf-tumor-id", vcf_tumor_id,
                    "--ncbi-build",   ncbi_build,
                    "--inhibit-vep",
                ]
            else:
                print(f"\nRunning vcf2maf.pl + VEP 115 (--vep-forks {args.vep_forks})...",
                      flush=True)
                vc_cmd = [
                    perl, str(vcf2maf),
                    "--input-vcf",     str(plain_vcf),
                    "--output-maf",    str(vc_maf),
                    "--tumor-id",      tumor_id,
                    "--vcf-tumor-id",  vcf_tumor_id,
                    "--ncbi-build",    ncbi_build,
                    "--cache-version", "115",
                    "--vep-forks",     str(args.vep_forks),
                    "--vep-data",      str(args.vep_cache),
                ]
            if args.ref_fasta:
                vc_cmd += ["--ref-fasta", str(args.ref_fasta)]
            if normal_id:
                vc_cmd += ["--normal-id", normal_id, "--vcf-normal-id", vcf_normal_id]
            if not args.inhibit_vep and vep:
                vc_cmd += ["--vep-path", str(vep.parent)]
            # Point samtools/tabix at the conda env so PATH order doesn't matter.
            for tool in ("samtools", "tabix"):
                tp = vcf2maf.parent / tool
                if tp.exists():
                    vc_cmd += [f"--{tool}-exec", str(tp)]
            run(vc_cmd, label="vcf2maf.pl")
            print(f"  → {vc_maf}", flush=True)

            # --- Compare (streaming: O(1) memory, no full MAF load) ---
            print("\nComparing MAFs...", flush=True)
            diff = compare_mafs_streaming(str(ms_maf), str(vc_maf), cmp_fields=cmp_fields)

            mode_tag = "--inhibit-vep" if args.inhibit_vep else f"+VEP ({ncbi_build}, --vep-forks {args.vep_forks})"
            strict_tag = "(--strict)" if args.strict else ""
            run_desc = (
                f"mafsmith{strict_tag} vs vcf2maf.pl {mode_tag}"
                + (f"; max-variants {args.max_variants:,}" if args.max_variants else "")
            )
            write_diff_report(diff, diff_path, input_label,
                              tumor_id, normal_id, run_desc, cmp_fields=cmp_fields)

            if args.no_keep_mafs:
                ms_maf.unlink(missing_ok=True)
                vc_maf.unlink(missing_ok=True)

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
