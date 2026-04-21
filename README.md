# mafsmith

Fast, self-contained VCF↔MAF converter with embedded variant annotation.

mafsmith is a Rust reimplementation of [vcf2maf](https://github.com/mskcc/vcf2maf), using [fastVEP](https://github.com/Ensembl/ensembl-vep) for annotation instead of the Perl VEP stack. It produces MAF files that match `vcf2maf.pl --inhibit-vep` output field-for-field.

> **Status:** mafsmith is research software provided as-is. The `vcf2maf` subcommand has been validated against real-world VCFs from multiple callers (see below), but users should verify outputs independently for their use case. The `maf2vcf`, `maf2maf`, and `vcf2vcf` subcommands are experimental and largely untested — use with caution.

## Acknowledgements

mafsmith is a rust-based adaptation of the design, field conventions, and edge-case handling from [vcf2maf](https://github.com/mskcc/vcf2maf) by [@ckandoth](https://github.com/ckandoth) (ORCID: [0000-0002-1345-3573](https://orcid.org/0000-0002-1345-3573)).

## Performance

Benchmarked on 7 GIAB NIST v4.2.1 GRCh38 samples (HG001–HG007, ~3.9–4.0M variants each),
conversion step only (`--skip-annotation` / `--inhibit-vep`), on an AWS c6a.4xlarge
(AMD EPYC 7R13, 16 vCPU, 30 GiB RAM):

| Tool | Mean time | Variants/s |
|------|-----------|-----------|
| `mafsmith --skip-annotation` | 7.050 ± 0.371 s | ~559,000 |
| `vcf2maf.pl --inhibit-vep` | 559.073 ± 17.016 s | ~7,000 |

**79.4× faster** (range 74.3–84.1× across 7 samples) for the conversion step.
Full end-to-end pipeline speedup (including fastVEP annotation) to be reported separately.

See [`results/conversion_benchmark_giab_grch38.md`](results/conversion_benchmark_giab_grch38.md)
for full per-sample data, cost savings, and carbon estimates.

## Installation

### Prerequisites

- Rust toolchain (`cargo`)
- Python ≥ 3.8 (for `mafsmith fetch`)

### Build

```bash
cargo build --release
# binary at target/release/mafsmith
```

### Download reference data

```bash
# GRCh38 (default)
mafsmith fetch

# GRCh37 or mouse (GRCm39)
mafsmith fetch --genome grch37
mafsmith fetch --genome grcm39

# Multiple genomes
mafsmith fetch --genome grch38,grch37

# Use existing files instead of downloading
mafsmith fetch --gff3 /path/to/genes.gff3 --ref-fasta /path/to/ref.fa
```

Reference data and the fastVEP binary are stored in `~/.mafsmith/` by default.

## Usage

### vcf2maf — Convert VCF to MAF

```bash
# Single-sample VCF
mafsmith vcf2maf \
  -i tumor.vcf.gz \
  -o output.maf \
  --tumor-id SAMPLE_ID

# Paired tumor/normal VCF
mafsmith vcf2maf \
  -i paired.vcf.gz \
  -o output.maf \
  --vcf-tumor-id TUMOR \
  --tumor-id TUMOR \
  --vcf-normal-id NORMAL \
  --normal-id NORMAL

# Skip annotation (VCF already has CSQ field from a prior fastVEP run)
mafsmith vcf2maf \
  -i annotated.vcf \
  -o output.maf \
  --skip-annotation

# GRCh37
mafsmith vcf2maf -i input.vcf -o output.maf --genome grch37
```

#### Key options

| Flag | Description |
|------|-------------|
| `--genome` | Reference assembly: `grch38` (default), `grch37`, `grcm39` |
| `--tumor-id` | MAF `Tumor_Sample_Barcode` (defaults to VCF sample column name) |
| `--normal-id` | MAF `Matched_Norm_Sample_Barcode` |
| `--vcf-tumor-id` | Sample column name in VCF for tumor |
| `--vcf-normal-id` | Sample column name in VCF for normal |
| `--custom-enst` | File of preferred Ensembl transcript IDs (one per line) |
| `--retain-ann` | Comma-separated CSQ field names to pass through to MAF |
| `--skip-annotation` | Use existing CSQ annotations in VCF (skips fastVEP) |
| `--strict` | Match `vcf2maf.pl` exactly for truncated AD arrays (outputs `.` for depth fields instead of partial counts) |
| `--min-hom-vaf` | VAF threshold for inferring homozygous-alt genotype (default: 0.7) |

### maf2vcf, maf2maf, vcf2vcf — Experimental subcommands

> These subcommands are experimental and have not been validated. Use at your own risk.

```bash
# Convert MAF back to VCF
mafsmith maf2vcf -i input.maf -o output.vcf

# Reannotate an existing MAF (round-trips through fastVEP)
mafsmith maf2maf -i input.maf -o reannotated.maf

# Normalize a VCF
mafsmith vcf2vcf -i input.vcf -o normalized.vcf
```

## Compatibility with vcf2maf.pl

mafsmith targets field-for-field agreement with `vcf2maf.pl --inhibit-vep` (same fastVEP-annotated input). Validated to 0 conversion-field mismatches in `--strict` mode across the following caller types and datasets:

| Caller | VCF type | Source |
|--------|----------|--------|
| DeepVariant 1.2.0 | Single-sample gVCF (GT=`0/0`/`./.'`), with and without `VAF` field | syn31624545; syn4988483 (VA02, VA06) |
| GATK MuTect2 | Single-sample GRCh38 | syn64156972; syn31624525 |
| GATK MuTect2 (paired T/N) | Paired tumor/normal (`GT:AD:AF:DP:F1R2:F2R1:SB` FORMAT) | GIAB HG008; SEQC2 HCC1395 |
| FreeBayes | Single-sample | syn31624535 |
| Strelka2 germline | Single-sample (`variants.vcf` and `genome.vcf` formats) | syn31624939; syn31624637 |
| Strelka2 somatic SNVs | Paired T/N (per-base `AU/CU/GU/TU` depth fields) | GIAB HG008; SEQC2 HCC1395 |
| Strelka2 somatic indels | Paired T/N (`TAR`/`TIR` depth fields) | syn68172710; GIAB HG008 |
| SV callers (Manta/DELLY) | SV-only (BND, DEL, DUP, INV symbolic ALTs) | syn21296193 |
| VarScan2 somatic | Paired T/N (`RD`+`AD` FORMAT) | syn6840402 |
| VarDict | Paired T/N (`RD` strand-bias field coexists with `AD`) | syn6039268 |
| SomaticSniper | Paired T/N (`DP4`+`BCOUNT` FORMAT, no `AD`) | SEQC2 HCC1395 |
| GIAB germline benchmarks (HG001–HG007) | Multi-caller consensus (`ADALL` field), GRCh38 | NIST v4.2.1 |
| COSMIC v103 | Annotation database VCF (no sample columns), GRCh38 | COSMIC GenomeScreensMutant; NonCodingVariants |

### Known intentional differences

- **`--strict` mode off (default)**: when a caller emits a truncated `AD` array (fewer values than `REF + all ALTs`), mafsmith still extracts available depth counts. Use `--strict` to output `.` for those fields instead, matching `vcf2maf.pl` exactly.
- **SV secondary rows**: mafsmith emits secondary breakpoint rows with the actual partner chromosome/position. `vcf2maf.pl` emits them with an empty `Chromosome` field (a known bug).
- **Multi-allelic tie-breaking**: when two ALTs have equal depth, tie-breaking may differ from `vcf2maf.pl` for a small number of variants (~4 per 50k-variant file).
- **Transcript selection at gene boundaries**: for variants near 5′/3′ UTR–flank and Intron/RNA boundaries, mafsmith and vcf2maf.pl may select different canonical transcripts, affecting `Variant_Classification` for ~2–5 variants per dataset. This reflects different gene-model versions rather than a conversion bug.

## Supported genomes

| Assembly | Species |
|----------|---------|
| GRCh38 | Human (hg38) |
| GRCh37 | Human (hg19/b37) |
| GRCm39 | Mouse |

## License

Apache 2.0
