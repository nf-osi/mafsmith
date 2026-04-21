#!/bin/bash
# Batch vcf_diff for all 1,902 PCAWG consensus SNV/MNV VCFs — 8 parallel, --vep-forks 2
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VCF_DIR="/home/ssm-user/bench_vcfs/PCAWG_consensus/vcfs/snv_mnv"
OUT_BASE="/home/ssm-user/bench_vcfs/PCAWG_consensus/diffs"
MAFSMITH="${SCRIPT_DIR}/../target/release/mafsmith"

run_one() {
    vcf="$1"
    uuid="$(basename "$vcf" .consensus.20160830.somatic.snv_mnv.vcf.gz)"
    outdir="${OUT_BASE}/${uuid}"

    if [ -f "${outdir}/diff.txt" ]; then
        echo "[SKIP] ${uuid} already done"
        return 0
    fi

    mkdir -p "$outdir"

    python3 "${SCRIPT_DIR}/vcf_diff.py" "$vcf" \
        --genome grch37 \
        --tumor-id "$uuid" \
        --output-dir "$outdir" \
        --strict \
        --vep-forks 2 \
        --mafsmith "$MAFSMITH" \
        > "${outdir}/run.log" 2>&1

    diffs=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "${outdir}/diff.txt" 2>/dev/null || echo "?")
    echo "[DONE] ${uuid}: ${diffs} conversion diffs"
}

export OUT_BASE SCRIPT_DIR MAFSMITH
export -f run_one

mkdir -p "$OUT_BASE"
echo "Starting PCAWG batch: $(ls ${VCF_DIR}/*.vcf.gz | wc -l) VCFs, 8 parallel"
ls "${VCF_DIR}"/*.vcf.gz | xargs -P 8 -I{} bash -c 'run_one "$@"' _ {}
echo "PCAWG batch complete."
