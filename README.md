# mafsmith

Fast, self-contained VCF↔MAF converter with embedded variant annotation.

mafsmith is a Rust reimplementation of [vcf2maf](https://github.com/mskcc/vcf2maf), using [fastVEP](https://github.com/Ensembl/ensembl-vep) for annotation. All four conversion subcommands (`vcf2maf`, `maf2vcf`, `vcf2vcf`, `maf2maf`) have been validated to 0 conversion-field differences against their reference Perl counterparts across 15+ caller types and 6 representative datasets.

## Acknowledgements

mafsmith is a Rust-based adaptation of the design, field conventions, and edge-case handling from [vcf2maf](https://github.com/mskcc/vcf2maf) by [@ckandoth](https://github.com/ckandoth) (ORCID: [0000-0002-1345-3573](https://orcid.org/0000-0002-1345-3573)).

## Performance

Benchmarked on 7 GIAB NIST v4.2.1 GRCh38 samples (HG001–HG007, ~3.9–4.0M variants each),
conversion step only (`--skip-annotation` / `--inhibit-vep`), on an AWS c6a.4xlarge
(AMD EPYC 7R13, 16 vCPU, 30 GiB RAM):

| Tool | Mean time | Variants/s | Speedup |
|------|-----------|-----------|---------|
| `mafsmith` (16-core) | 7.050 ± 0.371 s | ~559,000 | **79.4×** |
| `mafsmith` (1-core) | 11.755 ± 0.462 s | ~335,000 | **47.6×** |
| `vcf2maf.pl --inhibit-vep` | 559.073 ± 17.016 s | ~7,000 | — |

Full end-to-end pipeline (mafsmith + fastVEP vs. vcf2maf.pl + VEP): **25–83×** faster depending on VEP fork count.

See [`results/conversion_benchmark_giab_grch38.md`](results/conversion_benchmark_giab_grch38.md) for full per-sample data, cost savings, and carbon estimates.

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

# GRCh37
mafsmith fetch --genome grch37

# Mouse (GRCm39)
mafsmith fetch --genome grcm39

# Multiple genomes
mafsmith fetch --genome grch38,grch37

# Use existing files instead of downloading
mafsmith fetch --gff3 /path/to/genes.gff3 --ref-fasta /path/to/ref.fa
```

Reference data and the fastVEP binary are stored in `~/.mafsmith/` by default.

---

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

# Skip annotation (VCF already has CSQ fields from a prior fastVEP run)
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
| `--strict` | Match `vcf2maf.pl` exactly: output `.` for depth fields when AD arrays are truncated |
| `--min-hom-vaf` | VAF threshold for inferring homozygous-alt genotype (default: 0.7) |

---

### maf2vcf — Convert MAF to VCF

Reconstructs a standards-compliant VCF from a MAF file. Recovers allele representations, multi-allelic sites, genotype strings (`GT`), and allele depth fields (`GT:AD:DP` when depth columns are present) from MAF columns. Anchor bases for indels are looked up from the reference FASTA.

```bash
mafsmith maf2vcf \
  -i input.maf \
  -o output.vcf \
  --genome grch38
```

| Flag | Description |
|------|-------------|
| `--genome` | Reference assembly (default: `grch38`) |
| `--ref-fasta` | Override reference FASTA path |
| `--per-tn-vcfs` | Write one VCF per tumor/normal pair to a directory |

---

### vcf2vcf — Normalize a VCF

Passes a VCF through with FORMAT field standardization and optional sample column selection. Non-PASS variants and multi-allelic ALTs are preserved. Only ref-only records (ALT=`.`) are dropped.

```bash
mafsmith vcf2vcf \
  -i input.vcf.gz \
  -o normalized.vcf \
  --vcf-tumor-id TUMOR \
  --vcf-normal-id NORMAL
```

| Flag | Description |
|------|-------------|
| `--vcf-tumor-id` | Tumor sample column to select |
| `--vcf-normal-id` | Normal sample column to select |
| `--genome` | Reference assembly (default: `grch38`) |

---

### maf2maf — Reannotate a MAF

Re-annotates an existing MAF by internally converting it to VCF, running fastVEP, and converting back to MAF. Equivalent to `maf2maf.pl` with VEP replaced by fastVEP.

```bash
mafsmith maf2maf \
  -i input.maf \
  -o reannotated.maf \
  --genome grch38
```

| Flag | Description |
|------|-------------|
| `--genome` | Reference assembly (default: `grch38`) |
| `--custom-enst` | File of preferred Ensembl transcript IDs |
| `--fastvep-path` | Override fastVEP binary path |

---

## Validation

### vcf2maf

Validated to 0 conversion-field mismatches in `--strict` mode against `vcf2maf.pl --inhibit-vep` across the following caller types:

| Caller | VCF type | Source |
|--------|----------|--------|
| DeepVariant 1.2.0 | Single-sample gVCF (`GT=0/0`/`./.`), with and without `VAF` field | syn31624545; syn4988483 |
| GATK MuTect2 | Single-sample GRCh38 | syn64156972; syn31624525 |
| GATK MuTect2 (paired T/N) | Paired tumor/normal (`GT:AD:AF:DP:F1R2:F2R1:SB`) | GIAB HG008; SEQC2 HCC1395 |
| FreeBayes | Single-sample | syn31624535 |
| Strelka2 germline | `variants.vcf` and `genome.vcf` formats | syn31624939; syn31624637 |
| Strelka2 somatic SNVs | Paired T/N (`AU`/`CU`/`GU`/`TU` depth fields) | GIAB HG008; SEQC2 HCC1395 |
| Strelka2 somatic indels | Paired T/N (`TAR`/`TIR` depth fields) | syn68172710; GIAB HG008 |
| SV callers (Manta/DELLY) | SV-only (BND, DEL, DUP, INV symbolic ALTs) | syn21296193 |
| VarScan2 somatic | Paired T/N (`RD`+`AD` FORMAT) | syn6840402 |
| VarDict | Paired T/N (`RD` strand-bias field alongside `AD`) | syn6039268 |
| SomaticSniper | Paired T/N (`DP4`+`BCOUNT` FORMAT, no `AD`) | SEQC2 HCC1395 |
| DKFZ SNV caller | Paired T/N (`GT:DP:DP4`, GRCh37, no chr prefix) | ICGC PCAWG |
| GIAB germline benchmarks (HG001–HG007) | Multi-caller consensus (`ADALL` field) | NIST v4.2.1 |
| ICGC PCAWG consensus SNV/MNV | Consensus VCF (depth in INFO, no FORMAT), GRCh37 | 1,902 samples |
| DepMap CCLE WGS | GATK MuTect2 single-sample, hg38 | 802 samples |
| COSMIC v103 | Annotation database VCF (no sample columns) | COSMIC |

When run with the same Ensembl VEP 112 annotation cache, mafsmith produces 0 conversion differences versus `vcf2maf.pl` across 23 representative datasets (GRCh38 and GRCh37).

### maf2vcf, vcf2vcf, maf2maf

Validated to 0 conversion-field differences against `maf2vcf.pl`, `vcf2vcf.pl`, and `maf2maf.pl` across 6 datasets (2,000 variants each):

| Dataset | Genome | `maf2vcf` | `vcf2vcf` | `maf2maf` |
|---------|--------|-----------|-----------|-----------|
| SEQC2 HCC1395, GATK MuTect2 | GRCh38 | 0 diffs | 0 diffs | 0 diffs |
| SEQC2 HCC1395, Strelka2 somatic | GRCh38 | 0 diffs | 0 diffs | — |
| SEQC2 HCC1395, SomaticSniper | GRCh38 | 0 diffs | 0 diffs | — |
| GIAB HG008, GATK MuTect2 | GRCh38 | 0 diffs | 0 diffs | 0 diffs |
| GIAB HG001 germline benchmark | GRCh38 | 0 diffs | 0 diffs | — |
| PCAWG consensus (0009b464) | GRCh37 | 0 diffs | 0 diffs | 0 diffs |

---

## Known intentional differences from vcf2maf.pl

- **Default vs. `--strict` mode**: when a caller emits a truncated `AD` array (fewer values than `REF + all ALTs`), mafsmith extracts available depth counts by default. Use `--strict` to output `.` for those fields instead, exactly matching `vcf2maf.pl`.
- **SV secondary rows**: mafsmith emits secondary breakpoint rows with the actual partner chromosome and position. `vcf2maf.pl` leaves these fields blank (a known bug in the reference implementation).
- **Multi-allelic tie-breaking**: when two ALTs have identical depth, tie-breaking may differ from `vcf2maf.pl` for a small number of variants (~4 per 50k-variant file).
- **Transcript selection at gene boundaries**: for variants near 5′/3′ UTR–flank and Intron/RNA boundaries, mafsmith and `vcf2maf.pl` may select different canonical transcripts, affecting `Variant_Classification` for ~2–5 variants per dataset. This reflects different gene-model versions rather than a conversion bug.

---

## Supported genomes

| Assembly | Species |
|----------|---------|
| GRCh38 | Human (hg38) |
| GRCh37 | Human (hg19/b37) |
| GRCm39 | Mouse |

## License

Apache 2.0
