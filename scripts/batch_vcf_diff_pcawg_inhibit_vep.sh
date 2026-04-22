#!/bin/bash
# Batch vcf_diff (--inhibit-vep) for all 1,902 PCAWG consensus SNV/MNV VCFs — 16 parallel
# Compares mafsmith --skip-annotation vs vcf2maf.pl --inhibit-vep (no VEP annotation).
# Variant_Classification is excluded; depth/allele/coordinate fields are compared.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VCF_DIR="/home/ssm-user/bench_vcfs/PCAWG_consensus/vcfs/snv_mnv"
OUT_BASE="/home/ssm-user/bench_vcfs/PCAWG_consensus/diffs_inhibit_vep"
MAFSMITH="${SCRIPT_DIR}/../target/release/mafsmith"
PARALLEL="${1:-4}"

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
        --inhibit-vep \
        --mafsmith "$MAFSMITH" \
        > "${outdir}/run.log" 2>&1

    diffs=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "${outdir}/diff.txt" 2>/dev/null || echo "?")
    echo "[DONE] ${uuid}: ${diffs} conversion diffs"
}

export OUT_BASE SCRIPT_DIR MAFSMITH
export -f run_one

mkdir -p "$OUT_BASE"
echo "Starting PCAWG --inhibit-vep batch: $(ls "${VCF_DIR}"/*.vcf.gz | wc -l) VCFs, ${PARALLEL} parallel"
ls "${VCF_DIR}"/*.vcf.gz | xargs -P "${PARALLEL}" -I{} bash -c 'run_one "$@"' _ {}
echo "PCAWG --inhibit-vep batch complete."

# Print summary: total diffs across all samples
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
