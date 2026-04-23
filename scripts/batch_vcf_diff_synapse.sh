#!/bin/bash
# vcf_diff (--inhibit-vep, full VCF, no sampling) for Synapse validation datasets
# from the mafsmith manuscript (NF-OSI project syn16858331).
# Excludes: syn64156972 (removed from manuscript), syn4988483 (not a VCF).
# Usage: batch_vcf_diff_synapse.sh [PARALLEL]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_BASE="/home/ssm-user/bench_vcfs/synapse_diffs"
WORK_DIR="/home/ssm-user/bench_vcfs/synapse_work"
MAFSMITH="${SCRIPT_DIR}/../target/release/mafsmith"
PYTHON="${PYTHON:-/home/ssm-user/miniconda3/bin/python3}"
PARALLEL="${1:-1}"

mkdir -p "$OUT_BASE" "$WORK_DIR"

run_one() {
    local syn_id="$1" label="$2" extra_args="${3:-}"
    local outdir="${OUT_BASE}/${label}"

    if [ -f "${outdir}/diff.txt" ]; then
        echo "[SKIP] ${label} already done"
        return 0
    fi
    mkdir -p "$outdir"

    # shellcheck disable=SC2086
    "$PYTHON" "${SCRIPT_DIR}/vcf_diff.py" "$syn_id" \
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
export SCRIPT_DIR OUT_BASE WORK_DIR MAFSMITH PYTHON

# Each entry: "syn_id|label|extra_args"
# T/N IDs inferred from filenames; Strelka2 somatic always uses TUMOR/NORMAL columns.
datasets=(
    # DeepVariant single-sample gVCF (VA01)
    "syn31624545|deepvariant_VA01|"

    # GATK MuTect2 single-sample (VA01)
    "syn31624525|mutect2_VA01|"

    # FreeBayes single-sample (VA01)
    "syn31624535|freebayes_VA01|"

    # Strelka2 germline variants.vcf (VA05)
    "syn31624939|strelka2_germline_variants_VA05|"

    # Strelka2 germline genome.vcf (VA01) — 267M variants (2 GB); sample 20k to avoid OOM
    "syn31624637|strelka2_germline_genome_VA01|--sample 20000"

    # Strelka2 somatic indels — paired T/N, Strelka2 uses TUMOR/NORMAL columns
    "syn68172710|strelka2_somatic_indels|"

    # SV callers (Manta/DELLY) — vcf2maf.pl crashes on symbolic ALT / BND notation;
    # SV datasets are validated separately via synapse_validate.py, not vcf_diff.py.

    # VarScan2 somatic — paired T/N
    "syn6840402|varscan2_somatic|"

    # VarDict paired T/N — normal listed first (patient_11_normal_* / patient_11_tumor_*)
    "syn6039268|vardict_paired_tn|"
)

printf '%s\n' "${datasets[@]}" |
xargs -P "$PARALLEL" -I{} bash -c '
    IFS="|" read -r syn_id label extra <<< "$1"
    run_one "$syn_id" "$label" "$extra"
' _ {}

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
