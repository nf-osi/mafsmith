#!/bin/bash
# Batch validation of maf2vcf, vcf2vcf, and maf2maf subcommands
# against their respective Perl reference tools.
#
# Runs SEQUENTIALLY to avoid OOM.
# Results written to $OUT_BASE/{maf2vcf,vcf2vcf,maf2maf}/<label>/diff.txt
#
# Prerequisites:
#   - mafsmith binary built at ../target/release/mafsmith
#   - VEP 112 caches at ~/.vep/homo_sapiens/112_GRCh38/ and ~/...112_GRCh37/
#   - Reference FASTAs at ~/.mafsmith/GRCh38/reference.fa and ~/.mafsmith/GRCh37/reference.fa
#   - Bench VCFs at /home/ssm-user/bench_vcfs/
#
# Usage: bash batch_subcommand_diff.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_BASE="/home/ssm-user/bench_vcfs/subcommand_diffs"
WORK_DIR="/home/ssm-user/bench_vcfs/subcommand_work"
MAFSMITH="${SCRIPT_DIR}/../target/release/mafsmith"
PYTHON="/home/ssm-user/miniconda3/bin/python3"
BENCH="/home/ssm-user/bench_vcfs"

MAX_VARIANTS=2000
VEP_FORKS=4

mkdir -p "$OUT_BASE" "$WORK_DIR"

# ── maf2vcf: generate MAF from VCF, compare back-conversion ──────────────────

run_maf2vcf_grch38() {
    local label="$1" vcf="$2"
    local outdir="${OUT_BASE}/maf2vcf/${label}"
    shift 2
    local extra_args=("$@")

    if [ -f "${outdir}/diff.txt" ]; then
        echo "[SKIP] maf2vcf/${label} already done"
        return 0
    fi
    mkdir -p "$outdir"

    "$PYTHON" "${SCRIPT_DIR}/maf2vcf_diff.py" dummy_placeholder \
        --input-vcf "$vcf" \
        --genome grch38 \
        --output-dir "$outdir" \
        --work-dir "${WORK_DIR}/maf2vcf_${label}" \
        --mafsmith "$MAFSMITH" \
        --max-variants "$MAX_VARIANTS" \
        "${extra_args[@]}" \
        > "${outdir}/run.log" 2>&1 && rc=0 || rc=$?

    diffs=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "${outdir}/diff.txt" 2>/dev/null || echo "?")
    if [ $rc -ne 0 ]; then
        echo "[FAIL] maf2vcf/${label}: exit ${rc} — see ${outdir}/run.log"
    else
        echo "[DONE] maf2vcf/${label}: ${diffs} diffs"
    fi
}

run_maf2vcf_grch37() {
    local label="$1" vcf="$2"
    local outdir="${OUT_BASE}/maf2vcf/${label}"
    shift 2
    local extra_args=("$@")

    if [ -f "${outdir}/diff.txt" ]; then
        echo "[SKIP] maf2vcf/${label} already done"
        return 0
    fi
    mkdir -p "$outdir"

    "$PYTHON" "${SCRIPT_DIR}/maf2vcf_diff.py" dummy_placeholder \
        --input-vcf "$vcf" \
        --genome grch37 \
        --output-dir "$outdir" \
        --work-dir "${WORK_DIR}/maf2vcf_${label}" \
        --mafsmith "$MAFSMITH" \
        --max-variants "$MAX_VARIANTS" \
        "${extra_args[@]}" \
        > "${outdir}/run.log" 2>&1 && rc=0 || rc=$?

    diffs=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "${outdir}/diff.txt" 2>/dev/null || echo "?")
    if [ $rc -ne 0 ]; then
        echo "[FAIL] maf2vcf/${label}: exit ${rc} — see ${outdir}/run.log"
    else
        echo "[DONE] maf2vcf/${label}: ${diffs} diffs"
    fi
}

# ── vcf2vcf: compare VCF normalization ───────────────────────────────────────

run_vcf2vcf_grch38() {
    local label="$1" vcf="$2"
    local outdir="${OUT_BASE}/vcf2vcf/${label}"
    shift 2
    local extra_args=("$@")

    if [ -f "${outdir}/diff.txt" ]; then
        echo "[SKIP] vcf2vcf/${label} already done"
        return 0
    fi
    mkdir -p "$outdir"

    "$PYTHON" "${SCRIPT_DIR}/vcf2vcf_diff.py" "$vcf" \
        --genome grch38 \
        --output-dir "$outdir" \
        --work-dir "${WORK_DIR}/vcf2vcf_${label}" \
        --mafsmith "$MAFSMITH" \
        --max-variants "$MAX_VARIANTS" \
        "${extra_args[@]}" \
        > "${outdir}/run.log" 2>&1 && rc=0 || rc=$?

    diffs=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "${outdir}/diff.txt" 2>/dev/null || echo "?")
    if [ $rc -ne 0 ]; then
        echo "[FAIL] vcf2vcf/${label}: exit ${rc} — see ${outdir}/run.log"
    else
        echo "[DONE] vcf2vcf/${label}: ${diffs} diffs"
    fi
}

run_vcf2vcf_grch37() {
    local label="$1" vcf="$2"
    local outdir="${OUT_BASE}/vcf2vcf/${label}"
    shift 2
    local extra_args=("$@")

    if [ -f "${outdir}/diff.txt" ]; then
        echo "[SKIP] vcf2vcf/${label} already done"
        return 0
    fi
    mkdir -p "$outdir"

    "$PYTHON" "${SCRIPT_DIR}/vcf2vcf_diff.py" "$vcf" \
        --genome grch37 \
        --output-dir "$outdir" \
        --work-dir "${WORK_DIR}/vcf2vcf_${label}" \
        --mafsmith "$MAFSMITH" \
        --max-variants "$MAX_VARIANTS" \
        "${extra_args[@]}" \
        > "${outdir}/run.log" 2>&1 && rc=0 || rc=$?

    diffs=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "${outdir}/diff.txt" 2>/dev/null || echo "?")
    if [ $rc -ne 0 ]; then
        echo "[FAIL] vcf2vcf/${label}: exit ${rc} — see ${outdir}/run.log"
    else
        echo "[DONE] vcf2vcf/${label}: ${diffs} diffs"
    fi
}

# ── maf2maf: compare round-trip re-annotation ────────────────────────────────

run_maf2maf_grch38() {
    local label="$1" vcf="$2"
    local outdir="${OUT_BASE}/maf2maf/${label}"
    shift 2
    local extra_args=("$@")

    if [ -f "${outdir}/diff.txt" ]; then
        echo "[SKIP] maf2maf/${label} already done"
        return 0
    fi
    mkdir -p "$outdir"

    # Generate input MAF first
    local work="${WORK_DIR}/maf2maf_${label}"
    mkdir -p "$work"
    local maf="${work}/input.maf"

    "$MAFSMITH" vcf2maf -i "$vcf" -o "$maf" \
        --genome grch38 --skip-annotation \
        "${extra_args[@]}" \
        >> "${outdir}/run.log" 2>&1 || true

    # Trim to max variants
    local maf_trimmed="${work}/input_trimmed.maf"
    { head -1 "$maf"; tail -n +2 "$maf" | head -"$MAX_VARIANTS"; } > "$maf_trimmed" || true

    "$PYTHON" "${SCRIPT_DIR}/maf2maf_diff.py" "$maf_trimmed" \
        --genome grch38 \
        --output-dir "$outdir" \
        --work-dir "$work" \
        --mafsmith "$MAFSMITH" \
        --vep-forks "$VEP_FORKS" \
        >> "${outdir}/run.log" 2>&1 && rc=0 || rc=$?

    diffs=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "${outdir}/diff.txt" 2>/dev/null || echo "?")
    if [ $rc -ne 0 ]; then
        echo "[FAIL] maf2maf/${label}: exit ${rc} — see ${outdir}/run.log"
    else
        echo "[DONE] maf2maf/${label}: ${diffs} diffs"
    fi
}

run_maf2maf_grch37() {
    local label="$1" vcf="$2"
    local outdir="${OUT_BASE}/maf2maf/${label}"
    shift 2
    local extra_args=("$@")

    if [ -f "${outdir}/diff.txt" ]; then
        echo "[SKIP] maf2maf/${label} already done"
        return 0
    fi
    mkdir -p "$outdir"

    local work="${WORK_DIR}/maf2maf_${label}"
    mkdir -p "$work"
    local maf="${work}/input.maf"

    "$MAFSMITH" vcf2maf -i "$vcf" -o "$maf" \
        --genome grch37 --skip-annotation \
        "${extra_args[@]}" \
        >> "${outdir}/run.log" 2>&1 || true

    local maf_trimmed="${work}/input_trimmed.maf"
    { head -1 "$maf"; tail -n +2 "$maf" | head -"$MAX_VARIANTS"; } > "$maf_trimmed" || true

    "$PYTHON" "${SCRIPT_DIR}/maf2maf_diff.py" "$maf_trimmed" \
        --genome grch37 \
        --output-dir "$outdir" \
        --work-dir "$work" \
        --mafsmith "$MAFSMITH" \
        --vep-forks "$VEP_FORKS" \
        >> "${outdir}/run.log" 2>&1 && rc=0 || rc=$?

    diffs=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "${outdir}/diff.txt" 2>/dev/null || echo "?")
    if [ $rc -ne 0 ]; then
        echo "[FAIL] maf2maf/${label}: exit ${rc} — see ${outdir}/run.log"
    else
        echo "[DONE] maf2maf/${label}: ${diffs} diffs"
    fi
}

echo "=== Subcommand validation (maf2vcf, vcf2vcf, maf2maf) ==="
echo "Max variants: $MAX_VARIANTS per dataset"
echo "VEP forks: $VEP_FORKS"
echo "Sequential execution"
echo "Started: $(date)"
echo ""

# ── maf2vcf validation ────────────────────────────────────────────────────────
echo "--- maf2vcf validation ---"

S="${BENCH}/SEQC2_HCC1395"
run_maf2vcf_grch38 "SEQC2_mutect2" \
    "${S}/WGS_FD_1.bwa.muTect2.vcf.gz" \
    --vcf-tumor-id TUMOR --tumor-id HCC1395 \
    --vcf-normal-id NORMAL --normal-id HCC1395BL

run_maf2vcf_grch38 "SEQC2_strelka" \
    "${S}/WGS_FD_1.bwa.strelka.vcf.gz" \
    --vcf-tumor-id TUMOR --tumor-id HCC1395 \
    --vcf-normal-id NORMAL --normal-id HCC1395BL

run_maf2vcf_grch38 "SEQC2_somaticSniper" \
    "${S}/WGS_FD_1.bwa.somaticSniper.vcf.gz" \
    --vcf-tumor-id TUMOR --tumor-id HCC1395 \
    --vcf-normal-id NORMAL --normal-id HCC1395BL

H="${BENCH}/GIAB_HG008_somatic"
run_maf2vcf_grch38 "GIAB_HG008_mutect2" \
    "${H}/HG008-T--HG008-N.mutect2.vcf.gz" \
    --vcf-tumor-id HG008-T --tumor-id HG008-T \
    --vcf-normal-id HG008-N --normal-id HG008-N

G="${BENCH}/GIAB_germline"
run_maf2vcf_grch38 "GIAB_germline_HG001" \
    "${G}/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    --vcf-tumor-id HG001 --tumor-id HG001

VCF_DIR="${BENCH}/PCAWG_consensus/vcfs/snv_mnv"
run_maf2vcf_grch37 "PCAWG_0009b464" \
    "${VCF_DIR}/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20160830.somatic.snv_mnv.vcf.gz" \
    --tumor-id "0009b464-b376-4fbc-8a56-da538269a02f"

# ── vcf2vcf validation ────────────────────────────────────────────────────────
echo ""
echo "--- vcf2vcf validation ---"

run_vcf2vcf_grch38 "SEQC2_mutect2" \
    "${S}/WGS_FD_1.bwa.muTect2.vcf.gz" \
    --vcf-tumor-id TUMOR --vcf-normal-id NORMAL

run_vcf2vcf_grch38 "SEQC2_strelka" \
    "${S}/WGS_FD_1.bwa.strelka.vcf.gz" \
    --vcf-tumor-id TUMOR --vcf-normal-id NORMAL

run_vcf2vcf_grch38 "SEQC2_somaticSniper" \
    "${S}/WGS_FD_1.bwa.somaticSniper.vcf.gz" \
    --vcf-tumor-id TUMOR --vcf-normal-id NORMAL

run_vcf2vcf_grch38 "GIAB_HG008_mutect2" \
    "${H}/HG008-T--HG008-N.mutect2.vcf.gz" \
    --vcf-tumor-id HG008-T --vcf-normal-id HG008-N

run_vcf2vcf_grch38 "GIAB_germline_HG001" \
    "${G}/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    --vcf-tumor-id HG001

run_vcf2vcf_grch37 "PCAWG_0009b464" \
    "${VCF_DIR}/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20160830.somatic.snv_mnv.vcf.gz"

# ── maf2maf validation ────────────────────────────────────────────────────────
echo ""
echo "--- maf2maf validation (fastVEP vs VEP 112) ---"

run_maf2maf_grch38 "SEQC2_mutect2" \
    "${S}/WGS_FD_1.bwa.muTect2.vcf.gz" \
    --vcf-tumor-id TUMOR --tumor-id HCC1395 \
    --vcf-normal-id NORMAL --normal-id HCC1395BL

run_maf2maf_grch38 "GIAB_HG008_mutect2" \
    "${H}/HG008-T--HG008-N.mutect2.vcf.gz" \
    --vcf-tumor-id HG008-T --tumor-id HG008-T \
    --vcf-normal-id HG008-N --normal-id HG008-N

run_maf2maf_grch37 "PCAWG_0009b464" \
    "${VCF_DIR}/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20160830.somatic.snv_mnv.vcf.gz" \
    --tumor-id "0009b464-b376-4fbc-8a56-da538269a02f"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "=== SUMMARY ==="
echo "Completed: $(date)"

for mode in maf2vcf vcf2vcf maf2maf; do
    total=0; nonzero=0; failed=0; ndatasets=0
    for f in "${OUT_BASE}/${mode}"/*/diff.txt; do
        [ -f "$f" ] || continue
        ndatasets=$((ndatasets + 1))
        n=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "$f" 2>/dev/null || echo 0)
        label=$(basename "$(dirname "$f")")
        total=$((total + n))
        if [ "$n" -gt 0 ]; then
            nonzero=$((nonzero + 1))
            echo "  [DIFF] ${mode}/${label}: ${n} diffs"
        else
            echo "  [OK]   ${mode}/${label}: 0 diffs"
        fi
    done
    for d in "${OUT_BASE}/${mode}"/*/; do
        [ -d "$d" ] || continue
        [ -f "${d}/diff.txt" ] || { failed=$((failed + 1)); echo "  [FAIL] ${mode}/$(basename $d)"; }
    done
    echo "  ${mode}: ${ndatasets} datasets, ${nonzero} with diffs, ${failed} failed, ${total} total diffs"
    echo ""
done
