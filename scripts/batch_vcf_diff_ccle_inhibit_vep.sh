#!/bin/bash
# Batch vcf_diff (--inhibit-vep) for all 802 DepMap CCLE WGS VCFs
# Single-sample GATK Mutect2 hg38; sample ID = cell line name from VCF column.
#
# Usage: batch_vcf_diff_ccle_inhibit_vep.sh [PARALLEL [MAX_VARIANTS]]
#   PARALLEL:     number of concurrent VCFs (default: 4; use 1 for safe sequential run)
#   MAX_VARIANTS: subsample each VCF to N variants (default: none = full VCF)
#                 WARNING: full VCFs have 4-5M variants; keep PARALLEL=1 without this.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VCF_DIR="/home/ssm-user/bench_vcfs/CCLE_WGS/vcfs"
OUT_BASE="/home/ssm-user/bench_vcfs/CCLE_WGS/diffs_inhibit_vep"
WORK_DIR="/home/ssm-user/bench_vcfs/CCLE_WGS/work"
MAFSMITH="${SCRIPT_DIR}/../target/release/mafsmith"
PARALLEL="${1:-4}"
MAX_VARIANTS="${2:-}"

run_one() {
    vcf="$1"
    cds_id="$(basename "$vcf" .vcf.gz)"
    outdir="${OUT_BASE}/${cds_id}"

    if [ -f "${outdir}/diff.txt" ]; then
        echo "[SKIP] ${cds_id} already done"
        return 0
    fi

    mkdir -p "$outdir"

    mv_arg=()
    [ -n "$MAX_VARIANTS" ] && mv_arg=(--max-variants "$MAX_VARIANTS")

    python3 "${SCRIPT_DIR}/vcf_diff.py" "$vcf" \
        --genome grch38 \
        --output-dir "$outdir" \
        --work-dir "${WORK_DIR}" \
        "${mv_arg[@]}" \
        --strict \
        --inhibit-vep \
        --no-keep-mafs \
        --mafsmith "$MAFSMITH" \
        > "${outdir}/run.log" 2>&1

    diffs=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "${outdir}/diff.txt" 2>/dev/null || echo "?")
    echo "[DONE] ${cds_id}: ${diffs} conversion diffs"
}

export OUT_BASE WORK_DIR SCRIPT_DIR MAFSMITH MAX_VARIANTS
export -f run_one

mkdir -p "$OUT_BASE" "$WORK_DIR"
echo "Starting CCLE --inhibit-vep batch: $(ls "${VCF_DIR}"/*.vcf.gz | wc -l) VCFs, ${PARALLEL} parallel${MAX_VARIANTS:+, max-variants ${MAX_VARIANTS}}"
ls "${VCF_DIR}"/*.vcf.gz | xargs -P "${PARALLEL}" -I{} bash -c 'run_one "$@"' _ {}
echo "CCLE --inhibit-vep batch complete."

echo ""
echo "=== SUMMARY ==="
total=0
nonzero=0
for f in "${OUT_BASE}"/*/diff.txt; do
    n=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "$f" 2>/dev/null || echo 0)
    total=$((total + n))
    [ "$n" -gt 0 ] && nonzero=$((nonzero + 1))
done
samples=$(ls -d "${OUT_BASE}"/*/  2>/dev/null | wc -l)
echo "Samples processed: ${samples}"
echo "Samples with diffs: ${nonzero}"
echo "Total conversion diffs: ${total}"
