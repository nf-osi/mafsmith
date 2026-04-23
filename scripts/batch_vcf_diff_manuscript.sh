#!/bin/bash
# vcf_diff (--inhibit-vep, full VCF, no sampling) for non-Synapse manuscript datasets:
#   GIAB germline HG001-HG007, GIAB HG008 somatic, SEQC2 HCC1395
# Sequential (PARALLEL=1) by default — GIAB germline files are ~4M variants each.
# Usage: batch_vcf_diff_manuscript.sh [PARALLEL]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BENCH="/home/ssm-user/bench_vcfs"
OUT_BASE="${BENCH}/manuscript_diffs"
WORK_DIR="${BENCH}/manuscript_work"
MAFSMITH="${SCRIPT_DIR}/../target/release/mafsmith"
PARALLEL="${1:-1}"

mkdir -p "$OUT_BASE" "$WORK_DIR"

run_one() {
    local label="$1" vcf="$2" outdir="$3" extra_args="${4:-}"
    if [ -f "${outdir}/diff.txt" ]; then
        echo "[SKIP] ${label} already done"
        return 0
    fi
    mkdir -p "$outdir"
    # shellcheck disable=SC2086
    python3 "${SCRIPT_DIR}/vcf_diff.py" "$vcf" \
        --genome grch38 \
        --output-dir "$outdir" \
        --work-dir "${WORK_DIR}" \
        --strict \
        --inhibit-vep \
        --no-keep-mafs \
        --mafsmith "$MAFSMITH" \
        $extra_args \
        > "${outdir}/run.log" 2>&1
    diffs=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "${outdir}/diff.txt" 2>/dev/null || echo "?")
    echo "[DONE] ${label}: ${diffs} conversion diffs"
}
export -f run_one
export SCRIPT_DIR WORK_DIR MAFSMITH

# ── GIAB germline HG001-HG007 (single-sample, GRCh38) ────────────────────────
G="${BENCH}/GIAB_germline"
for sid in HG001 HG002 HG003 HG004 HG005 HG006 HG007; do
    vcf="${G}/${sid}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    printf '%s\t%s\t%s\t%s\n' \
        "GIAB germline ${sid}" "$vcf" "${OUT_BASE}/GIAB_germline_${sid}" \
        "--vcf-tumor-id ${sid} --tumor-id ${sid}"
done |
xargs -P "$PARALLEL" -I{} bash -c '
    IFS=$'"'"'\t'"'"' read -r label vcf outdir extra <<< "$1"
    run_one "$label" "$vcf" "$outdir" "$extra"
' _ {}

# ── GIAB HG008 somatic (paired T/N, normal listed first in VCF) ──────────────
H="${BENCH}/GIAB_HG008_somatic"
TN_HG008="--vcf-tumor-id HG008-T --tumor-id HG008-T --vcf-normal-id HG008-N --normal-id HG008-N"
for spec in \
    "GIAB_HG008_mutect2:${H}/HG008-T--HG008-N.mutect2.vcf.gz" \
    "GIAB_HG008_strelka2_snv:${H}/HG008-T--HG008-N.snv.strelka2.vcf.gz" \
    "GIAB_HG008_strelka2_indel:${H}/HG008-T--HG008-N.indel.strelka2.vcf.gz"
do
    label="${spec%%:*}"
    vcf="${spec#*:}"
    run_one "$label" "$vcf" "${OUT_BASE}/${label}" "$TN_HG008" &
    [ "$PARALLEL" -eq 1 ] && wait
done
wait

# ── SEQC2 HCC1395 (paired T/N, samples named TUMOR/NORMAL) ───────────────────
S="${BENCH}/SEQC2_HCC1395"
TN_SEQC2="--vcf-tumor-id TUMOR --tumor-id HCC1395 --vcf-normal-id NORMAL --normal-id HCC1395BL"
for spec in \
    "SEQC2_mutect2:${S}/WGS_FD_1.bwa.muTect2.vcf.gz" \
    "SEQC2_strelka:${S}/WGS_FD_1.bwa.strelka.vcf.gz" \
    "SEQC2_somaticSniper:${S}/WGS_FD_1.bwa.somaticSniper.vcf.gz"
do
    label="${spec%%:*}"
    vcf="${spec#*:}"
    run_one "$label" "$vcf" "${OUT_BASE}/${label}" "$TN_SEQC2" &
    [ "$PARALLEL" -eq 1 ] && wait
done
wait

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "=== SUMMARY ==="
total=0; nonzero=0
for f in "${OUT_BASE}"/*/diff.txt; do
    [ -f "$f" ] || continue
    n=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "$f" 2>/dev/null || echo 0)
    label=$(basename "$(dirname "$f")")
    total=$((total + n))
    [ "$n" -gt 0 ] && nonzero=$((nonzero + 1))
    echo "  ${label}: ${n} diffs"
done
echo "Datasets with diffs: ${nonzero} / $(ls -d "${OUT_BASE}"/*/ 2>/dev/null | wc -l)"
echo "Total conversion diffs: ${total}"
