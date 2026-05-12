#!/usr/bin/env bash
# Timing benchmark: mafsmith vs vcf2maf.pl
#
# Measures wall-clock throughput (variants/second) for both tools on the same
# pre-annotated VCF, isolating conversion speed from annotation differences.
#
# Usage:
#   bash scripts/benchmark.sh [OPTIONS] INPUT.vcf
#
# Options:
#   --tumor-id ID       VCF tumor sample ID (default: first sample)
#   --normal-id ID      VCF normal sample ID (default: second sample)
#   --ref-fasta FILE    Reference FASTA (required for vcf2maf.pl; optional for mafsmith)
#   --gff3 FILE         GFF3 annotation (required for fastVEP annotation)
#   --fastvep PATH      Path to fastVEP binary
#   --iterations N      Repeat each timing N times (default: 3)
#   --full              Also benchmark full annotation pipeline (not just conversion)
#   -h, --help          Show this message

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

MAFSMITH="${MAFSMITH:-$REPO_ROOT/target/release/mafsmith}"
FASTVEP="${FASTVEP:-$HOME/.mafsmith/bin/fastvep}"
REF_FASTA="${REF_FASTA:-$HOME/.mafsmith/GRCh38/reference.fa}"
GFF3="${GFF3:-$HOME/.mafsmith/GRCh38/genes.gff3}"
VCF2MAF="${VCF2MAF:-vcf2maf.pl}"
CONDA_ENV="${CONDA_ENV:-vcf2maf-env}"

TUMOR_ID=""
NORMAL_ID=""
ITERATIONS=3
FULL_PIPELINE=0
INPUT_VCF=""

usage() {
    grep '^#' "$0" | head -20 | sed 's/^# \{0,2\}//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tumor-id)   TUMOR_ID="$2";  shift 2 ;;
        --normal-id)  NORMAL_ID="$2"; shift 2 ;;
        --ref-fasta)  REF_FASTA="$2"; shift 2 ;;
        --gff3)       GFF3="$2";      shift 2 ;;
        --fastvep)    FASTVEP="$2";   shift 2 ;;
        --iterations) ITERATIONS="$2"; shift 2 ;;
        --full)       FULL_PIPELINE=1; shift ;;
        -h|--help)    usage ;;
        -*)           echo "Unknown option: $1"; exit 1 ;;
        *)            INPUT_VCF="$1"; shift ;;
    esac
done

if [[ -z "$INPUT_VCF" ]]; then
    echo "Usage: benchmark.sh [OPTIONS] INPUT.vcf"
    exit 1
fi

if [[ ! -f "$MAFSMITH" ]]; then
    echo "ERROR: mafsmith binary not found at $MAFSMITH"
    echo "  Run: cargo build --release"
    exit 1
fi

if ! command -v conda >/dev/null 2>&1 && ! command -v "$VCF2MAF" >/dev/null 2>&1; then
    echo "ERROR: vcf2maf.pl not found (conda env '$CONDA_ENV' or in PATH)"
    exit 1
fi

TMPDIR_BM="$(mktemp -d -t mafsmith_bench_XXXX)"
trap 'rm -rf "$TMPDIR_BM"' EXIT

# ── Portable high-resolution timestamp (works on macOS and Linux) ─────────────
now_s() { python3 -c "import time; print(f'{time.monotonic():.6f}')"; }

# ── Detect sample IDs ──────────────────────────────────────────────────────────
if [[ -z "$TUMOR_ID" || -z "$NORMAL_ID" ]]; then
    CHROM_LINE=$(zcat "$INPUT_VCF" 2>/dev/null || cat "$INPUT_VCF")
    CHROM_LINE=$(echo "$CHROM_LINE" | grep '^#CHROM' | head -1)
    samples=($(echo "$CHROM_LINE" | cut -f10- | tr '\t' ' '))
    TUMOR_ID="${TUMOR_ID:-${samples[0]:-TUMOR}}"
    NORMAL_ID="${NORMAL_ID:-${samples[1]:-NORMAL}}"
fi

echo "==================================================================="
echo " mafsmith vs vcf2maf.pl throughput benchmark"
echo "==================================================================="
echo " Input:      $INPUT_VCF"
echo " Tumor ID:   $TUMOR_ID"
echo " Normal ID:  $NORMAL_ID"
echo " Iterations: $ITERATIONS"
echo "==================================================================="
echo ""

# ── Step 1: Annotate with fastVEP (shared, run once) ─────────────────────────
ANNOTATED="$TMPDIR_BM/annotated.vcf"
PLAIN_VCF="$TMPDIR_BM/plain.vcf"

echo "[1/4] Preparing VCF (strip chr prefix if needed)..."
# Strip chr prefix to match bare FASTA chromosome names
if zcat "$INPUT_VCF" 2>/dev/null | grep -v '^#' | head -1 | grep -q '^chr'; then
    (zcat "$INPUT_VCF" 2>/dev/null || cat "$INPUT_VCF") | \
        awk 'BEGIN{OFS="\t"} /^#/{print; next} {$1=substr($1,4); print}' \
        > "$PLAIN_VCF"
else
    (zcat "$INPUT_VCF" 2>/dev/null || cat "$INPUT_VCF") > "$PLAIN_VCF"
fi

VARIANT_COUNT=$(grep -v '^#' "$PLAIN_VCF" | wc -l | tr -d ' ')
echo "    Variants in input VCF: $VARIANT_COUNT"

echo "[2/4] Running fastVEP annotation (once, shared)..."
if [[ -f "$FASTVEP" ]] && [[ -f "$GFF3" ]] && [[ -f "$REF_FASTA" ]]; then
    T_ANNOT_START=$(now_s)
    "$FASTVEP" annotate \
        -i "$PLAIN_VCF" \
        -o "$ANNOTATED" \
        --fasta "$REF_FASTA" \
        --gff3 "$GFF3" \
        --hgvs \
        --output-format vcf 2>/dev/null
    T_ANNOT_END=$(now_s)
    ANNOT_SECS=$(python3 -c "print(f'{$T_ANNOT_END - $T_ANNOT_START:.2f}')")
    ANNOT_VARS=$(grep -v '^#' "$ANNOTATED" | wc -l | tr -d ' ')
    echo "    Annotation complete: $ANNOT_VARS variants in ${ANNOT_SECS}s"
else
    echo "    WARNING: fastVEP/GFF3/FASTA not found — using unannotated VCF"
    echo "             (results will show Targeted_Region for all variants)"
    ANNOTATED="$PLAIN_VCF"
    ANNOT_SECS="N/A"
fi

# ── Helper: time a command N times, report mean ───────────────────────────────
time_cmd() {
    local label="$1"
    local n="$2"
    shift 2
    local total=0
    local elapsed
    local all_times=()
    for i in $(seq 1 "$n"); do
        local t0 t1
        t0=$(now_s)
        "$@" 2>/dev/null || { echo "    WARNING: command failed (exit $?)"; }
        t1=$(now_s)
        elapsed=$(python3 -c "print(f'{$t1 - $t0:.3f}')")
        all_times+=("${elapsed}s")
        total=$(python3 -c "print(f'{$total + $elapsed:.6f}')")
    done
    MEAN_TIME=$(python3 -c "print(f'{$total / $n:.3f}')")
    echo "    $label: mean ${MEAN_TIME}s over $n runs (${all_times[*]})"
}

count_variants() {
    grep -v '^#' "$1" | wc -l | tr -d ' '
}

# ── Step 3: Time mafsmith (conversion only, --skip-annotation) ────────────────
echo ""
echo "[3/4] Benchmarking mafsmith (--skip-annotation)..."
MS_OUT="$TMPDIR_BM/mafsmith.maf"
MS_ARGS=("$MAFSMITH" vcf2maf
    --input-vcf "$ANNOTATED"
    --output-maf "$MS_OUT"
    --genome grch38
    --vcf-tumor-id "$TUMOR_ID"
    --tumor-id "$TUMOR_ID"
    --skip-annotation)
[[ -n "$NORMAL_ID" ]] && MS_ARGS+=(--vcf-normal-id "$NORMAL_ID" --normal-id "$NORMAL_ID")
[[ -f "$REF_FASTA" ]] && MS_ARGS+=(--ref-fasta "$REF_FASTA")

time_cmd "mafsmith" "$ITERATIONS" "${MS_ARGS[@]}"
MS_TIME="$MEAN_TIME"
MS_VARS=$(count_variants "$MS_OUT")

# ── Step 4: Time vcf2maf.pl (--inhibit-vep) ────────────────────────────────────
echo ""
echo "[4/4] Benchmarking vcf2maf.pl (--inhibit-vep)..."
VC_OUT="$TMPDIR_BM/vcf2maf.maf"
if [[ -f "$REF_FASTA" ]]; then
    VC_ARGS=(conda run -n "$CONDA_ENV" vcf2maf.pl
        --input-vcf "$ANNOTATED"
        --output-maf "$VC_OUT"
        --inhibit-vep
        --tumor-id "$TUMOR_ID"
        --vcf-tumor-id "$TUMOR_ID"
        --ref-fasta "$REF_FASTA")
    [[ -n "$NORMAL_ID" ]] && VC_ARGS+=(--normal-id "$NORMAL_ID" --vcf-normal-id "$NORMAL_ID")

    time_cmd "vcf2maf.pl" "$ITERATIONS" "${VC_ARGS[@]}"
    VC_TIME="$MEAN_TIME"
    VC_VARS=$(count_variants "$VC_OUT")
else
    echo "    SKIPPED (no --ref-fasta)"
    VC_TIME="N/A"
    VC_VARS="N/A"
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "==================================================================="
echo " RESULTS"
echo "==================================================================="
echo " Input variants:   $VARIANT_COUNT"
echo " Annotated VCF:    $([ "$ANNOTATED" != "$PLAIN_VCF" ] && echo "$ANNOT_SECS s" || echo "N/A (no fastVEP)")"
echo ""
printf " %-20s %10s %12s\n" "Tool" "Time (s)" "Variants/s"
printf " %-20s %10s %12s\n" "----" "--------" "----------"

if [[ "$MS_TIME" != "N/A" ]]; then
    MS_VPS=$(echo "scale=0; $MS_VARS / $MS_TIME" | bc 2>/dev/null || echo "N/A")
    printf " %-20s %10s %12s\n" "mafsmith" "$MS_TIME" "$MS_VPS"
fi
if [[ "$VC_TIME" != "N/A" ]]; then
    VC_VPS=$(echo "scale=0; $VC_VARS / $VC_TIME" | bc 2>/dev/null || echo "N/A")
    printf " %-20s %10s %12s\n" "vcf2maf.pl" "$VC_TIME" "$VC_VPS"
fi

if [[ "$MS_TIME" != "N/A" && "$VC_TIME" != "N/A" ]]; then
    SPEEDUP=$(echo "scale=1; $VC_TIME / $MS_TIME" | bc 2>/dev/null || echo "N/A")
    echo ""
    echo " Speedup: mafsmith is ${SPEEDUP}x faster than vcf2maf.pl"
fi

echo ""
echo " MAF rows written:"
echo "   mafsmith:   $MS_VARS"
echo "   vcf2maf.pl: $VC_VARS"

if [[ "$FULL_PIPELINE" -eq 1 ]] && [[ -f "$FASTVEP" ]]; then
    echo ""
    echo "==================================================================="
    echo " FULL PIPELINE (including fastVEP annotation)"
    echo "==================================================================="
    echo " fastVEP annotation: ${ANNOT_SECS}s"

    if [[ "$MS_TIME" != "N/A" ]]; then
        MS_TOTAL=$(echo "scale=3; $MS_TIME + $ANNOT_SECS" | bc)
        echo " mafsmith total:     ${MS_TOTAL}s"
    fi
    if [[ "$VC_TIME" != "N/A" ]]; then
        VC_TOTAL=$(echo "scale=3; $VC_TIME + $ANNOT_SECS" | bc)
        echo " vcf2maf.pl total:   ${VC_TOTAL}s"
        if [[ "$MS_TIME" != "N/A" ]]; then
            FULL_SPEEDUP=$(echo "scale=1; $VC_TOTAL / $MS_TOTAL" | bc 2>/dev/null || echo "N/A")
            echo " Full pipeline speedup: ${FULL_SPEEDUP}x"
        fi
    fi
fi

echo ""
echo "Done."
