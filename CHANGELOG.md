# Changelog

All notable changes to mafsmith will be documented here.

## [0.1.0] — 2026-05-14

Initial release.

### Features

- `vcf2maf`: VCF to MAF conversion with fastVEP or Ensembl VEP annotation. Validated to 0 conversion-field differences against `vcf2maf.pl` in `--strict` mode across 21 caller types and 87.2 million variants.
- `maf2vcf`: MAF to VCF reconstruction, validated against `maf2vcf.pl` across 6 datasets.
- `vcf2vcf`: VCF normalization and sample column selection, validated against `vcf2vcf.pl`.
- `maf2maf`: MAF re-annotation via fastVEP, validated against `maf2maf.pl`.
- `fetch`: one-command download of reference FASTA and GFF3 data for GRCh38, GRCh37, and GRCm39.
- `--strict` mode for bit-for-bit output compatibility with `vcf2maf.pl`.
- Rayon-based parallel processing; 47.6× single-core and 79.4× 16-core speedup over `vcf2maf.pl` on GIAB germline benchmarks.
- Support for structural variant ALTs (BND, DEL, DUP, INV) including secondary breakpoint rows.
- Caller-specific depth extraction for GATK, Strelka2, VarScan2, VarDict, SomaticSniper, DeepVariant, FreeBayes, Ion Torrent, and DKFZ SNV callers.
