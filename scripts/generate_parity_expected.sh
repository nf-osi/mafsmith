#!/usr/bin/env bash
# Regenerate the vcf2maf.pl reference MAFs for the parity suite
# (tests/fixtures/parity/expected/), which tests/parity_test.rs compares
# mafsmith against.
#
# Runs vcf2maf.pl --inhibit-vep on every fixture in tests/fixtures/parity/manifest.tsv
# with the same sample IDs and --retain-* options the test passes to mafsmith.
# No VEP or network access is needed: annotated.vcf already carries CSQ, and the
# other fixtures are converted unannotated.
#
# Requirements: vcf2maf.pl, samtools and tabix. By default all three are taken from the
# vcf2maf-env conda environment; override with VCF2MAF=/path/to/vcf2maf.pl and
# SAMTOOLS=/path/to/samtools, TABIX=/path/to/tabix. The committed references were made with
# vcf2maf 1.6.22 (bioconda).
#
# Usage: bash scripts/generate_parity_expected.sh
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PARITY="$REPO/tests/fixtures/parity"
ENV_BIN="${CONDA_PREFIX_VCF2MAF:-$HOME/miniconda3/envs/vcf2maf-env/bin}"
VCF2MAF="${VCF2MAF:-$ENV_BIN/vcf2maf.pl}"
SAMTOOLS="${SAMTOOLS:-$ENV_BIN/samtools}"
TABIX="${TABIX:-$ENV_BIN/tabix}"
PERL="${PERL:-$ENV_BIN/perl}"
[[ -x "$PERL" ]] || PERL=perl

mkdir -p "$PARITY/expected"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

while IFS=$'\t' read -r name tumor normal info fmt ann; do
    [[ -z "$name" || "$name" == \#* ]] && continue
    args=(--inhibit-vep --input-vcf "$PARITY/$name.vcf" --output-maf "$PARITY/expected/$name.vcf2maf.maf"
          --ref-fasta "$PARITY/reference.fa" --ncbi-build GRCh38 --samtools-exec "$SAMTOOLS" --tabix-exec "$TABIX"
          --tmp-dir "$TMP" --tumor-id "$tumor" --vcf-tumor-id "$tumor")
    [[ "$normal" != "-" ]] && args+=(--normal-id "$normal" --vcf-normal-id "$normal")
    [[ "$info" != "-" ]] && args+=(--retain-info "$info")
    [[ "$fmt" != "-" ]] && args+=(--retain-fmt "$fmt")
    [[ "$ann" != "-" ]] && args+=(--retain-ann "$ann")
    echo "==> $name"
    "$PERL" "$VCF2MAF" "${args[@]}" 2> "$TMP/$name.log" || { cat "$TMP/$name.log"; exit 1; }
    grep -v "^STATUS" "$TMP/$name.log" | sed 's/^/    /' || true
done < "$PARITY/manifest.tsv"
