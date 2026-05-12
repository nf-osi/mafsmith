#!/usr/bin/env bash
# Compare mafsmith vs vcf2maf on the real-world Strelka VCF (syn29349669).
#
# Prerequisites:
#   - Rust installed + `cargo build --release` done
#   - `mafsmith fetch --genome grch38` done (installs fastVEP + downloads GRCh38 data)
#   - vcf2maf + VEP available (or the vcf2maf Docker image)
#   - tests/fixtures/real_world_subset.vcf present (extracted from syn29349669)
#
# Usage:
#   bash scripts/run_real_world_comparison.sh
#
# To use the full VCF instead of the 500-variant subset:
#   REAL_WORLD_VCF=tests/fixtures/real_world.vcf bash scripts/run_real_world_comparison.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
FIXTURES="$REPO_ROOT/tests/fixtures"
EXPECTED="$FIXTURES/expected"

MAFSMITH="${MAFSMITH:-$REPO_ROOT/target/release/mafsmith}"
if [[ ! -f "$MAFSMITH" ]]; then
    MAFSMITH="$REPO_ROOT/target/debug/mafsmith"
fi

REAL_WORLD_VCF="${REAL_WORLD_VCF:-$FIXTURES/real_world_subset.vcf}"
REFERENCE_MAF="$EXPECTED/real_world_vcf2maf.maf"
TMPDIR_OUT="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_OUT"' EXIT

echo "==> Input VCF: $REAL_WORLD_VCF"
echo "    $(grep -v '^#' "$REAL_WORLD_VCF" | wc -l) variants"

# ── Step 1: generate vcf2maf reference output ─────────────────────────────────
if [[ ! -f "$REFERENCE_MAF" ]]; then
    echo ""
    echo "==> Generating vcf2maf reference MAF..."
    echo "    (Set VCF2MAF, VEP_PATH, VEP_CACHE, REF_FASTA env vars as needed)"

    VEP_PATH="${VEP_PATH:-/usr/local/bin}"
    VEP_CACHE="${VEP_CACHE:-$HOME/.vep}"
    REF_FASTA="${REF_FASTA:-$HOME/.mafsmith/GRCh38/reference.fa}"
    VCF2MAF="${VCF2MAF:-vcf2maf.pl}"

    mkdir -p "$EXPECTED"
    perl "$VCF2MAF" \
        --vep-path "$VEP_PATH" \
        --vep-data "$VEP_CACHE" \
        --vep-overwrite \
        --ncbi-build GRCh38 \
        --input-vcf "$REAL_WORLD_VCF" \
        --output-maf "$REFERENCE_MAF" \
        --ref-fasta "$REF_FASTA" \
        --vcf-tumor-id "cNF00.10a"

    echo "==> Reference MAF: $REFERENCE_MAF"
fi

# ── Step 2: run mafsmith ──────────────────────────────────────────────────────
MAFSMITH_OUT="$TMPDIR_OUT/mafsmith.maf"
echo ""
echo "==> Running mafsmith vcf2maf..."

"$MAFSMITH" vcf2maf \
    --input-vcf "$REAL_WORLD_VCF" \
    --output-maf "$MAFSMITH_OUT" \
    --genome grch38 \
    --vcf-tumor-id "cNF00.10a" \
    --tumor-id "cNF00.10a"

echo "    Wrote $(grep -v '^#' "$MAFSMITH_OUT" | wc -l) rows"

# ── Step 3: compare ───────────────────────────────────────────────────────────
if [[ -f "$REFERENCE_MAF" ]]; then
    echo ""
    echo "==> Comparing outputs..."
    python3 "$SCRIPT_DIR/compare_maf.py" "$REFERENCE_MAF" "$MAFSMITH_OUT"
else
    echo ""
    echo "-- No vcf2maf reference MAF available for comparison."
    echo "   Run with vcf2maf installed, or set REFERENCE_MAF=/path/to/ref.maf"
    echo "   mafsmith output: $MAFSMITH_OUT"
    cp "$MAFSMITH_OUT" "$EXPECTED/real_world_mafsmith.maf"
    echo "   Saved to: $EXPECTED/real_world_mafsmith.maf"
fi
