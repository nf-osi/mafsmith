#!/usr/bin/env bash
# Run mafsmith vcf2maf on test fixtures and compare against reference MAF output.
#
# Usage:
#   bash scripts/run_comparison.sh [--real-world /path/to/input.vcf]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
FIXTURES="$REPO_ROOT/tests/fixtures"
EXPECTED="$FIXTURES/expected"

MAFSMITH="${MAFSMITH:-$REPO_ROOT/target/release/mafsmith}"
if [[ ! -f "$MAFSMITH" ]]; then
    MAFSMITH="$REPO_ROOT/target/debug/mafsmith"
fi

if [[ ! -f "$MAFSMITH" ]]; then
    echo "ERROR: mafsmith binary not found. Run 'cargo build --release' first."
    exit 1
fi

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

echo "==> mafsmith binary: $MAFSMITH"

# ── Test 1: pre-annotated VCF (no external tools needed) ─────────────────────
echo ""
echo "==> Test 1: pre-annotated VCF fixture (fast, no fastVEP needed)"

"$MAFSMITH" vcf2maf \
    --input-vcf "$FIXTURES/test_b38_annotated.vcf" \
    --output-maf "$TMPDIR/test1.maf" \
    --genome grch38 \
    --tumor-id TUMOR \
    --normal-id NORMAL \
    --skip-annotation

python3 "$SCRIPT_DIR/compare_maf.py" \
    "$EXPECTED/test_b38_key_cols.tsv" \
    "$TMPDIR/test1.maf" && echo "✓ Test 1 PASSED" || echo "✗ Test 1 FAILED"

# ── Test 2: full annotation pipeline vs vcf2maf reference ────────────────────
REFERENCE_MAF="$EXPECTED/test_b38_vcf2maf.maf"
if [[ -f "$REFERENCE_MAF" ]]; then
    echo ""
    echo "==> Test 2: full annotation pipeline vs vcf2maf reference"

    "$MAFSMITH" vcf2maf \
        --input-vcf "$FIXTURES/test_b38.vcf" \
        --output-maf "$TMPDIR/test2.maf" \
        --genome grch38 \
        --tumor-id TUMOR \
        --normal-id NORMAL

    python3 "$SCRIPT_DIR/compare_maf.py" \
        "$REFERENCE_MAF" \
        "$TMPDIR/test2.maf" && echo "✓ Test 2 PASSED" || echo "✗ Test 2 FAILED"
else
    echo ""
    echo "-- Test 2 SKIPPED: $REFERENCE_MAF not found."
    echo "   Run scripts/generate_expected.sh to create the reference MAF."
fi

# ── Test 3: real-world VCF ────────────────────────────────────────────────────
REAL_WORLD_VCF="${REAL_WORLD_VCF:-$FIXTURES/real_world.vcf}"
REAL_WORLD_REF="$EXPECTED/real_world_vcf2maf.maf"

if [[ -f "$REAL_WORLD_VCF" ]]; then
    echo ""
    echo "==> Test 3: real-world VCF ($REAL_WORLD_VCF)"

    "$MAFSMITH" vcf2maf \
        --input-vcf "$REAL_WORLD_VCF" \
        --output-maf "$TMPDIR/test3.maf" \
        --genome grch38

    if [[ -f "$REAL_WORLD_REF" ]]; then
        python3 "$SCRIPT_DIR/compare_maf.py" \
            "$REAL_WORLD_REF" \
            "$TMPDIR/test3.maf" && echo "✓ Test 3 PASSED" || echo "✗ Test 3 FAILED"
    else
        echo "   No reference MAF for real-world VCF — mafsmith output at $TMPDIR/test3.maf"
        echo "   Row count: $(grep -v '^#' "$TMPDIR/test3.maf" | wc -l)"
    fi
else
    echo ""
    echo "-- Test 3 SKIPPED: set REAL_WORLD_VCF=/path/to/input.vcf to enable."
fi

echo ""
echo "Done."
