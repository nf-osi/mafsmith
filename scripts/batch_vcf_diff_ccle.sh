#!/bin/bash
# Batch vcf_diff for all 802 CCLE WGS VCFs — 8 parallel, --vep-forks 2
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VCF_DIR="/home/ssm-user/bench_vcfs/CCLE_WGS"
OUT_BASE="/home/ssm-user/bench_vcfs/CCLE_WGS/diffs"
MAFSMITH="${SCRIPT_DIR}/../target/release/mafsmith"

run_one() {
    vcf="$1"
    cds_id="$(basename "$vcf" .vcf.gz)"
    outdir="${OUT_BASE}/${cds_id}"

    if [ -f "${outdir}/diff.txt" ]; then
        echo "[SKIP] ${cds_id} already done"
        return 0
    fi

    sample="$(zcat "$vcf" | awk '/^#CHROM/{print $NF; exit}')"
    mkdir -p "$outdir"

    python3 "${SCRIPT_DIR}/vcf_diff.py" "$vcf" \
        --vcf-tumor-id "$sample" --tumor-id "$sample" \
        --output-dir "$outdir" \
        --max-variants 3000 \
        --strict \
        --vep-forks 2 \
        --mafsmith "$MAFSMITH" \
        > "${outdir}/run.log" 2>&1

    diffs=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "${outdir}/diff.txt" 2>/dev/null || echo "?")
    echo "[DONE] ${cds_id} (${sample}): ${diffs} conversion diffs"
}

export OUT_BASE SCRIPT_DIR MAFSMITH
export -f run_one

mkdir -p "$OUT_BASE"
echo "Starting CCLE batch: $(ls ${VCF_DIR}/*.vcf.gz | wc -l) VCFs, 8 parallel"
ls "${VCF_DIR}"/*.vcf.gz | xargs -P 8 -I{} bash -c 'run_one "$@"' _ {}
echo "CCLE batch complete."
