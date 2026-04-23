#!/bin/bash
# VEP-mode Table 1 validation: run mafsmith+VEP and vcf2maf.pl+VEP on all Table 1 datasets.
# Both tools use the same Ensembl VEP 112 cache so annotation differences are eliminated.
# Runs SEQUENTIALLY to avoid OOM.
#
# Prerequisites:
#   - VEP 112 caches extracted at ~/.vep/homo_sapiens/112_GRCh38/ and ~/...112_GRCh37/
#     (download + extract: see download_and_extract_vep_caches.sh)
#   - VEP 112 binary: ~/miniconda3/envs/vcf2maf-env/bin/vep
#   - synapseclient (for Synapse dataset downloads)
#   - Local VCFs in /home/ssm-user/bench_vcfs/{GIAB_germline,GIAB_HG008_somatic,
#                                                SEQC2_HCC1395,PCAWG_consensus,CCLE_WGS}
#
# Usage: bash batch_vcf_diff_vep.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_BASE="/home/ssm-user/bench_vcfs/vep_diffs"
WORK_DIR="/home/ssm-user/bench_vcfs/vep_work"
MAFSMITH="${SCRIPT_DIR}/../target/release/mafsmith"
PYTHON="/home/ssm-user/miniconda3/bin/python3"
BENCH="/home/ssm-user/bench_vcfs"

MAX_VARIANTS=2000
VEP_FORKS=4

mkdir -p "$OUT_BASE" "$WORK_DIR"

run_grch38() {
    local label="$1" input="$2"
    local outdir="${OUT_BASE}/${label}"
    shift 2
    local extra_args=("$@")

    if [ -f "${outdir}/diff.txt" ]; then
        echo "[SKIP] ${label} already done"
        return 0
    fi
    mkdir -p "$outdir"

    "$PYTHON" "${SCRIPT_DIR}/vcf_diff.py" "$input" \
        --genome grch38 \
        --output-dir "$outdir" \
        --work-dir "${WORK_DIR}" \
        --strict \
        --max-variants "$MAX_VARIANTS" \
        --vep-forks "$VEP_FORKS" \
        --no-keep-mafs \
        --mafsmith "$MAFSMITH" \
        "${extra_args[@]}" \
        > "${outdir}/run.log" 2>&1
    rc=$?

    diffs=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "${outdir}/diff.txt" 2>/dev/null || echo "?")
    if [ $rc -ne 0 ]; then
        echo "[FAIL] ${label}: exit ${rc} — see ${outdir}/run.log"
    else
        echo "[DONE] ${label}: ${diffs} conversion diffs"
    fi
}

run_grch37() {
    local label="$1" input="$2"
    local outdir="${OUT_BASE}/${label}"
    shift 2
    local extra_args=("$@")

    if [ -f "${outdir}/diff.txt" ]; then
        echo "[SKIP] ${label} already done"
        return 0
    fi
    mkdir -p "$outdir"

    "$PYTHON" "${SCRIPT_DIR}/vcf_diff.py" "$input" \
        --genome grch37 \
        --output-dir "$outdir" \
        --work-dir "${WORK_DIR}" \
        --strict \
        --max-variants "$MAX_VARIANTS" \
        --vep-forks "$VEP_FORKS" \
        --no-keep-mafs \
        --mafsmith "$MAFSMITH" \
        "${extra_args[@]}" \
        > "${outdir}/run.log" 2>&1
    rc=$?

    diffs=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "${outdir}/diff.txt" 2>/dev/null || echo "?")
    if [ $rc -ne 0 ]; then
        echo "[FAIL] ${label}: exit ${rc} — see ${outdir}/run.log"
    else
        echo "[DONE] ${label}: ${diffs} conversion diffs"
    fi
}

echo "=== VEP-mode Table 1 validation ==="
echo "Max variants: $MAX_VARIANTS per dataset"
echo "VEP forks: $VEP_FORKS"
echo "Sequential execution (no parallel)"
echo "Started: $(date)"
echo ""

# ── Synapse datasets (GRCh38) ─────────────────────────────────────────────────
echo "--- Synapse datasets (GRCh38) ---"

run_grch38 "deepvariant_VA01"     "syn31624545"
run_grch38 "mutect2_single_VA01"  "syn31624525"
run_grch38 "freebayes_VA01"       "syn31624535"
run_grch38 "strelka2_germline_variants_VA05" "syn31624939"
run_grch38 "strelka2_germline_genome_VA01"   "syn31624637" "--skip-ref-blocks"
run_grch38 "strelka2_somatic_indels"         "syn68172710"
run_grch38 "varscan2_somatic"                "syn6840402"
run_grch38 "vardict_paired_tn"               "syn6039268"

# ── GIAB HG008 somatic (GRCh38, paired T/N) ──────────────────────────────────
echo ""
echo "--- GIAB HG008 somatic (GRCh38) ---"

H="${BENCH}/GIAB_HG008_somatic"
run_grch38 "GIAB_HG008_mutect2" \
    "${H}/HG008-T--HG008-N.mutect2.vcf.gz" \
    --vcf-tumor-id HG008-T --tumor-id HG008-T \
    --vcf-normal-id HG008-N --normal-id HG008-N

run_grch38 "GIAB_HG008_strelka2_snv" \
    "${H}/HG008-T--HG008-N.snv.strelka2.vcf.gz" \
    --vcf-tumor-id HG008-T --tumor-id HG008-T \
    --vcf-normal-id HG008-N --normal-id HG008-N

run_grch38 "GIAB_HG008_strelka2_indel" \
    "${H}/HG008-T--HG008-N.indel.strelka2.vcf.gz" \
    --vcf-tumor-id HG008-T --tumor-id HG008-T \
    --vcf-normal-id HG008-N --normal-id HG008-N

# ── SEQC2 HCC1395 (GRCh38, paired T/N, samples TUMOR/NORMAL) ─────────────────
echo ""
echo "--- SEQC2 HCC1395 somatic (GRCh38) ---"

S="${BENCH}/SEQC2_HCC1395"
run_grch38 "SEQC2_mutect2" \
    "${S}/WGS_FD_1.bwa.muTect2.vcf.gz" \
    --vcf-tumor-id TUMOR --tumor-id HCC1395 \
    --vcf-normal-id NORMAL --normal-id HCC1395BL

run_grch38 "SEQC2_strelka" \
    "${S}/WGS_FD_1.bwa.strelka.vcf.gz" \
    --vcf-tumor-id TUMOR --tumor-id HCC1395 \
    --vcf-normal-id NORMAL --normal-id HCC1395BL

run_grch38 "SEQC2_somaticSniper" \
    "${S}/WGS_FD_1.bwa.somaticSniper.vcf.gz" \
    --vcf-tumor-id TUMOR --tumor-id HCC1395 \
    --vcf-normal-id NORMAL --normal-id HCC1395BL

# ── GIAB germline benchmarks (GRCh38, 3 representative samples) ──────────────
echo ""
echo "--- GIAB germline benchmarks (GRCh38, 3 samples) ---"

G="${BENCH}/GIAB_germline"
for sid in HG001 HG002 HG003; do
    vcf="${G}/${sid}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    run_grch38 "GIAB_germline_${sid}" "$vcf" \
        --vcf-tumor-id "$sid" --tumor-id "$sid"
done

# ── DepMap CCLE WGS (GRCh38, 3 representative samples) ───────────────────────
echo ""
echo "--- DepMap CCLE WGS (GRCh38, 3 samples) ---"

for cid in CDS-NMp6b2 CDS-7nLbDx CDS-rBDWH4; do
    vcf="${BENCH}/CCLE_WGS/vcfs/${cid}.vcf.gz"
    if [ -f "$vcf" ]; then
        run_grch38 "CCLE_${cid}" "$vcf"
    else
        echo "[SKIP] CCLE ${cid}: VCF not found at $vcf"
    fi
done

# ── ICGC PCAWG consensus (GRCh37, 3 representative samples) ──────────────────
echo ""
echo "--- ICGC PCAWG consensus (GRCh37, 3 samples) ---"

VCF_DIR="${BENCH}/PCAWG_consensus/vcfs/snv_mnv"
for uuid in \
    "0009b464-b376-4fbc-8a56-da538269a02f" \
    "003819bc-c415-4e76-887c-931d60ed39e7" \
    "0040b1b6-b07a-4b6e-90ef-133523eaf412"
do
    vcf="${VCF_DIR}/${uuid}.consensus.20160830.somatic.snv_mnv.vcf.gz"
    if [ -f "$vcf" ]; then
        run_grch37 "PCAWG_${uuid:0:8}" "$vcf" --tumor-id "$uuid"
    else
        echo "[SKIP] PCAWG ${uuid}: VCF not found"
    fi
done

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "=== SUMMARY ==="
echo "Completed: $(date)"
total=0; nonzero=0; failed=0
for f in "${OUT_BASE}"/*/diff.txt; do
    [ -f "$f" ] || continue
    n=$(grep -oP "CONVERSION DIFFERENCES \(\K[0-9]+" "$f" 2>/dev/null || echo 0)
    label=$(basename "$(dirname "$f")")
    total=$((total + n))
    if [ "$n" -gt 0 ]; then
        nonzero=$((nonzero + 1))
        echo "  [DIFF] ${label}: ${n} diffs"
    else
        echo "  [OK]   ${label}: 0 diffs"
    fi
done
for d in "${OUT_BASE}"/*/; do
    [ -f "${d}/diff.txt" ] || { failed=$((failed + 1)); echo "  [FAIL] $(basename $d)"; }
done
total_datasets=$(ls -d "${OUT_BASE}"/*/ 2>/dev/null | wc -l)
echo ""
echo "Datasets processed: $((total_datasets - failed)) / ${total_datasets}"
echo "Datasets with diffs: ${nonzero}"
echo "Failed runs: ${failed}"
echo "Total conversion diffs: ${total}"
