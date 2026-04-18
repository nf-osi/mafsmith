#!/usr/bin/env bash
# Generate reference MAF output using vcf2maf for comparison against mafsmith.
#
# Prerequisites:
#   - vcf2maf installed (perl vcf2maf.pl)
#   - VEP installed and cache available
#   - samtools in PATH
#   - GRCh38 chr21 reference FASTA at $REF_FASTA
#   - VEP cache at $VEP_CACHE
#
# Usage:
#   VEP_PATH=/path/to/vep \
#   VEP_CACHE=/path/to/vep_cache \
#   REF_FASTA=/path/to/chr21.fa \
#   bash scripts/generate_expected.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
FIXTURES="$REPO_ROOT/tests/fixtures"
EXPECTED="$FIXTURES/expected"

VEP_PATH="${VEP_PATH:-/usr/local/bin}"
VEP_CACHE="${VEP_CACHE:-$HOME/.vep}"
REF_FASTA="${REF_FASTA:-$FIXTURES/Homo_sapiens.GRCh38.dna.chromosome.21.fa}"
VCF2MAF="${VCF2MAF:-vcf2maf.pl}"

mkdir -p "$EXPECTED"

echo "==> Generating reference MAF for test_b38.vcf..."
perl "$VCF2MAF" \
    --vep-path "$VEP_PATH" \
    --vep-data "$VEP_CACHE" \
    --vep-overwrite \
    --ncbi-build GRCh38 \
    --input-vcf "$FIXTURES/test_b38.vcf" \
    --output-maf "$EXPECTED/test_b38_vcf2maf.maf" \
    --ref-fasta "$REF_FASTA" \
    --tumor-id TUMOR \
    --normal-id NORMAL \
    --vcf-tumor-id TUMOR \
    --vcf-normal-id NORMAL

echo "==> Reference MAF written to $EXPECTED/test_b38_vcf2maf.maf"

# If a real-world VCF is available, generate its reference MAF too.
if [[ -f "$FIXTURES/real_world.vcf" ]]; then
    echo "==> Generating reference MAF for real_world.vcf..."
    perl "$VCF2MAF" \
        --vep-path "$VEP_PATH" \
        --vep-data "$VEP_CACHE" \
        --vep-overwrite \
        --ncbi-build GRCh38 \
        --input-vcf "$FIXTURES/real_world.vcf" \
        --output-maf "$EXPECTED/real_world_vcf2maf.maf" \
        --ref-fasta "$REF_FASTA"
    echo "==> Done: $EXPECTED/real_world_vcf2maf.maf"
fi
