#!/usr/bin/env python3
"""
Compare two MAF files column-by-column.

Usage:
    python scripts/compare_maf.py reference.maf mafsmith_output.maf [options]

Matches rows by (Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2).
Reports mismatches for a configurable set of columns.
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

# Columns that must match exactly — derived purely from VCF input, not VEP database.
EXACT_COLUMNS = [
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "Transcript_ID",
    "Exon_Number",
    "t_depth",
    "t_ref_count",
    "t_alt_count",
    "n_depth",
    "n_ref_count",
    "n_alt_count",
    "Matched_Norm_Sample_Barcode",
    "Tumor_Sample_Barcode",
]

# Columns we compare but accept minor formatting differences (e.g. "p.Ter" vs "p.*").
LENIENT_COLUMNS = [
    "HGVSc",
    "HGVSp",
    "HGVSp_Short",
]

# Columns we skip entirely (version-dependent annotation data).
SKIP_COLUMNS = {
    "dbSNP_RS",           # rsID assignment changes between VEP releases
    "DOMAINS",            # randomly-ordered comma list (vcf2maf excludes this in its own tests)
    "SIFT",               # database version dependent
    "PolyPhen",           # database version dependent
    "all_effects",        # contains full transcript list — complex comparison
    "Entrez_Gene_Id",     # not always populated
    "gnomADe_AF",
    "gnomADe_AFR_AF",
    "gnomADe_AMR_AF",
    "gnomADe_ASJ_AF",
    "gnomADe_EAS_AF",
    "gnomADe_FIN_AF",
    "gnomADe_NFE_AF",
    "gnomADe_OTH_AF",
    "gnomADe_SAS_AF",
    "AF", "AFR_AF", "AMR_AF", "ASN_AF", "EAS_AF", "EUR_AF", "SAS_AF",
}


def read_maf(path: Path) -> tuple[list[str], list[dict]]:
    """Read a MAF file, returning (headers, rows). Skips # comment lines."""
    headers = []
    rows = []
    with open(path, newline="") as fh:
        for line in fh:
            stripped = line.rstrip("\n")
            if stripped.startswith("#"):
                continue
            if not headers:
                headers = stripped.split("\t")
                continue
            vals = stripped.split("\t")
            row = dict(zip(headers, vals))
            rows.append(row)
    return headers, rows


def variant_key(row: dict) -> tuple:
    return (
        row.get("Chromosome", ""),
        row.get("Start_Position", ""),
        row.get("Reference_Allele", ""),
        row.get("Tumor_Seq_Allele2", ""),
    )


def normalize_hgvsp(s: str) -> str:
    """Normalize HGVSp: p.Ter → p.*, Ala → A, etc., for lenient comparison.

    Keep in sync with shorten_hgvsp / THREE_TO_ONE in src/annotation/csq.rs.
    Note: Rust has Xle→J which is absent here (no practical impact for comparison).
    """
    three_to_one = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "Ter": "*", "Sec": "U", "Pyl": "O", "Xaa": "X",
    }
    result = s
    for three, one in three_to_one.items():
        result = result.replace(three, one)
    return result


def compare(ref_path: Path, got_path: Path, exact_cols, lenient_cols, skip_cols) -> dict:
    _, ref_rows = read_maf(ref_path)
    _, got_rows = read_maf(got_path)

    # Index reference rows by variant key
    ref_by_key: dict[tuple, dict] = {}
    for row in ref_rows:
        k = variant_key(row)
        if k in ref_by_key:
            print(f"WARNING: duplicate key in reference MAF: {k}", file=sys.stderr)
        ref_by_key[k] = row

    exact_mismatches = defaultdict(list)
    lenient_mismatches = defaultdict(list)
    missing_in_ref = []
    missing_in_got = []

    for got_row in got_rows:
        k = variant_key(got_row)
        if k not in ref_by_key:
            missing_in_ref.append(k)
            continue
        ref_row = ref_by_key[k]

        for col in exact_cols:
            if col in skip_cols:
                continue
            got_val = got_row.get(col, "")
            ref_val = ref_row.get(col, "")
            if got_val != ref_val:
                exact_mismatches[col].append({
                    "key": k,
                    "got": got_val,
                    "ref": ref_val,
                })

        for col in lenient_cols:
            if col in skip_cols:
                continue
            got_val = normalize_hgvsp(got_row.get(col, ""))
            ref_val = normalize_hgvsp(ref_row.get(col, ""))
            if got_val != ref_val:
                lenient_mismatches[col].append({
                    "key": k,
                    "got": got_row.get(col, ""),
                    "ref": ref_row.get(col, ""),
                })

    got_keys = {variant_key(r) for r in got_rows}
    for ref_row in ref_rows:
        k = variant_key(ref_row)
        if k not in got_keys:
            missing_in_got.append(k)

    return {
        "ref_count": len(ref_rows),
        "got_count": len(got_rows),
        "exact_mismatches": dict(exact_mismatches),
        "lenient_mismatches": dict(lenient_mismatches),
        "missing_in_ref": missing_in_ref,
        "missing_in_got": missing_in_got,
    }


def print_report(result: dict):
    print(f"\n{'='*60}")
    print(f"Reference rows:  {result['ref_count']}")
    print(f"mafsmith rows:   {result['got_count']}")

    if result["missing_in_ref"]:
        print(f"\n⚠  {len(result['missing_in_ref'])} variants in mafsmith output NOT in reference:")
        for k in result["missing_in_ref"][:10]:
            print(f"   {k}")

    if result["missing_in_got"]:
        print(f"\n⚠  {len(result['missing_in_got'])} variants in reference NOT in mafsmith output:")
        for k in result["missing_in_got"][:10]:
            print(f"   {k}")

    exact = result["exact_mismatches"]
    if not exact:
        print("\n✓  All exact columns match!")
    else:
        total = sum(len(v) for v in exact.values())
        print(f"\n✗  {total} exact-column mismatches across {len(exact)} columns:")
        for col, mismatches in sorted(exact.items()):
            print(f"\n  [{col}] ({len(mismatches)} mismatches)")
            for m in mismatches[:5]:
                print(f"    {m['key']}")
                print(f"      ref: {m['ref']!r}")
                print(f"      got: {m['got']!r}")
            if len(mismatches) > 5:
                print(f"    ... and {len(mismatches)-5} more")

    lenient = result["lenient_mismatches"]
    if lenient:
        total = sum(len(v) for v in lenient.values())
        print(f"\n~  {total} lenient-column differences (may be formatting only):")
        for col, mismatches in sorted(lenient.items()):
            print(f"\n  [{col}] ({len(mismatches)} differences)")
            for m in mismatches[:3]:
                print(f"    {m['key']}")
                print(f"      ref: {m['ref']!r}")
                print(f"      got: {m['got']!r}")

    total_exact = sum(len(v) for v in exact.values())
    print(f"\n{'='*60}")
    return total_exact == 0 and not result["missing_in_ref"] and not result["missing_in_got"]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reference", type=Path, help="Reference MAF (from vcf2maf)")
    parser.add_argument("mafsmith", type=Path, help="mafsmith output MAF")
    parser.add_argument("--extra-exact", nargs="*", default=[], help="Additional columns to compare exactly")
    parser.add_argument("--skip", nargs="*", default=[], help="Additional columns to skip")
    args = parser.parse_args()

    skip = SKIP_COLUMNS | set(args.skip)
    exact = EXACT_COLUMNS + args.extra_exact

    result = compare(args.reference, args.mafsmith, exact, LENIENT_COLUMNS, skip)
    ok = print_report(result)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
