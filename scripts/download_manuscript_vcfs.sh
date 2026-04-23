#!/bin/bash
# Download non-Synapse VCFs used in the mafsmith manuscript (validation + benchmark tables).
# Targets: GIAB germline HG001-HG007, GIAB HG008 somatic, SEQC2 HCC1395.
# Output layout matches what benchmark_all.py and batch_vcf_diff_manuscript.sh expect.
set -euo pipefail

BENCH="${1:-/home/ssm-user/bench_vcfs}"
GIAB_G="${BENCH}/GIAB_germline"
GIAB_H="${BENCH}/GIAB_HG008_somatic"
SEQC2="${BENCH}/SEQC2_HCC1395"

mkdir -p "$GIAB_G" "$GIAB_H" "$SEQC2"

GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release"
SEQC2_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/analysis/SNVs/vcfs/WGS"
HG008_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/NYGC-somatic-pipeline_20240412/GRCh38-GIABv3"

fetch() {
    local url="$1" dst="$2"
    local fname
    fname="$(basename "$dst")"
    if [ -f "$dst" ]; then
        echo "[SKIP] $fname (already present)"
        return 0
    fi
    echo "[DL]   $fname"
    wget -q --show-progress -O "${dst}.tmp" "$url" && mv "${dst}.tmp" "$dst"
    echo "[OK]   $fname ($(du -sh "$dst" | cut -f1))"
}

# GIAB germline HG001-HG007 (GRCh38, NIST v4.2.1)
fetch "${GIAB_BASE}/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
      "${GIAB_G}/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
fetch "${GIAB_BASE}/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
      "${GIAB_G}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
fetch "${GIAB_BASE}/AshkenazimTrio/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
      "${GIAB_G}/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
fetch "${GIAB_BASE}/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
      "${GIAB_G}/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
fetch "${GIAB_BASE}/ChineseTrio/HG005_NA24631_son/latest/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
      "${GIAB_G}/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
fetch "${GIAB_BASE}/ChineseTrio/HG006_NA24694_father/latest/GRCh38/HG006_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
      "${GIAB_G}/HG006_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
fetch "${GIAB_BASE}/ChineseTrio/HG007_NA24695_mother/latest/GRCh38/HG007_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
      "${GIAB_G}/HG007_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

# GIAB HG008 somatic (GRCh38, NYGC pipeline, paired T/N: normal listed first in VCF)
fetch "${HG008_BASE}/HG008-T--HG008-N.mutect2.vcf.gz"     "${GIAB_H}/HG008-T--HG008-N.mutect2.vcf.gz"
fetch "${HG008_BASE}/HG008-T--HG008-N.snv.strelka2.vcf.gz" "${GIAB_H}/HG008-T--HG008-N.snv.strelka2.vcf.gz"
fetch "${HG008_BASE}/HG008-T--HG008-N.indel.strelka2.vcf.gz" "${GIAB_H}/HG008-T--HG008-N.indel.strelka2.vcf.gz"

# SEQC2 HCC1395 (GRCh38, paired T/N: samples named TUMOR/NORMAL)
fetch "${SEQC2_BASE}/WGS_FD_1.bwa.muTect2.vcf.gz"       "${SEQC2}/WGS_FD_1.bwa.muTect2.vcf.gz"
fetch "${SEQC2_BASE}/WGS_FD_1.bwa.strelka.vcf.gz"        "${SEQC2}/WGS_FD_1.bwa.strelka.vcf.gz"
fetch "${SEQC2_BASE}/WGS_FD_1.bwa.somaticSniper.vcf.gz"  "${SEQC2}/WGS_FD_1.bwa.somaticSniper.vcf.gz"

echo ""
echo "Downloads complete. Files:"
ls -lh "$GIAB_G"/*.vcf.gz "$GIAB_H"/*.vcf.gz "$SEQC2"/*.vcf.gz
