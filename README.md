# mafsmith

Fast, self-contained VCF↔MAF converter with embedded variant annotation.

mafsmith is a Rust reimplementation of [vcf2maf](https://github.com/mskcc/vcf2maf), using [fastVEP](https://github.com/Ensembl/ensembl-vep) for annotation instead of the Perl VEP stack. It produces MAF files that match `vcf2maf.pl --inhibit-vep` output field-for-field.

> **Status:** mafsmith is research software provided as-is. The `vcf2maf` subcommand has been validated against real-world VCFs from multiple callers (see below), but users should verify outputs independently for their use case. The `maf2vcf`, `maf2maf`, and `vcf2vcf` subcommands are experimental and largely untested — use with caution.

## Acknowledgements

mafsmith builds directly on the design, field conventions, and edge-case handling documented in [vcf2maf](https://github.com/mskcc/vcf2maf). Without the years of work by [@ckandoth](https://github.com/ckandoth) (ORCID: [0000-0002-1345-3573](https://orcid.org/0000-0002-1345-3573)) in understanding and encoding the full complexity of the VCF and MAF specifications, this package would not exist.

## Performance

Benchmarked on a 6,292-variant SV VCF (post-annotation, conversion step only):

| Tool | Mean time | Variants/s |
|------|-----------|-----------|
| `mafsmith --skip-annotation` | ~0.06 s | ~104,000 |
| `vcf2maf.pl --inhibit-vep` | ~13.9 s | ~132 |

**~229× faster** for conversion. Full pipeline speedup including annotation: ~4.8×.

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

mafsmith targets field-for-field agreement with `vcf2maf.pl --inhibit-vep` (same fastVEP-annotated input). Validated to 0 mismatches across 20,000 variants each for the following caller types:

| Caller | VCF type | Synapse example |
|--------|----------|----------------|
| DRAGEN RefCall | Single-sample, GT=`0/0`/`./.'` | syn31624545 |
| MuTect2 | Single-sample GRCh38 | syn64156972 |
| FreeBayes | Single-sample | syn31624535 |
| Strelka2 | Paired tumor/normal | syn31624939 |
| Strelka2 somatic indels | Paired tumor/normal | syn68172710 |
| SV callers (Manta/DELLY) | SV-only | syn21296193 |

### Known intentional differences

- **`--strict` mode off (default)**: when a caller emits a truncated `AD` array (fewer values than `REF + all ALTs`), mafsmith still extracts available depth counts. Use `--strict` to output `.` for those fields instead, matching `vcf2maf.pl` exactly.
- **SV secondary rows**: mafsmith emits secondary breakpoint rows with the actual partner chromosome/position. `vcf2maf.pl` emits them with an empty `Chromosome` field (a known bug).
- **Multi-allelic tie-breaking**: when two ALTs have equal depth, tie-breaking may differ from `vcf2maf.pl` for a small number of variants (~4 per 50k-variant file).

## Supported genomes

| Assembly | Species |
|----------|---------|
| GRCh38 | Human (hg38) |
| GRCh37 | Human (hg19/b37) |
| GRCm39 | Mouse |

## License

Apache 2.0
