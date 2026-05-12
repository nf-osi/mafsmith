# scripts/

Validation, benchmarking, and data-prep scripts. See `results/README.md` for
the mapping from manuscript tables to specific entry points here.

## Validation (Table 1 / Table 2)

| Script | Purpose |
|--------|---------|
| `vcf_diff.py` | Core field-by-field comparison of two MAF files. Used by all batch validators below. Handles SV secondary rows, multi-allelic AD layouts, caller-specific allele-depth conventions. |
| `synapse_validate.py` | Wraps `vcf_diff.py` for Synapse-hosted VCFs (handles auth, downloads, expected-MAF generation). |
| `batch_vcf_diff_synapse.sh` | Runs the NF Data Portal validation cohort (Table 1 rows: DeepVariant, FreeBayes, Strelka2 germline, VarScan2, VarDict, Strelka2 somatic indels). |
| `batch_vcf_diff_ccle.sh` | DepMap CCLE WGS validation (Table 1: GATK MuTect2 single-sample). |
| `batch_vcf_diff_ccle_inhibit_vep.sh` | Same as above but passes `--inhibit-vep` to vcf2maf.pl, isolating conversion-logic differences from VEP version differences. |
| `batch_vcf_diff_pcawg.sh` | ICGC PCAWG cell-line validation (Table 1: DKFZ SNV/MNV; GRCh37). |
| `batch_vcf_diff_pcawg_inhibit_vep.sh` | PCAWG run with `--inhibit-vep`. |
| `batch_vcf_diff_manuscript.sh` | Consolidates per-caller runs from a single fixed set of manuscript datasets. |
| `batch_vcf_diff_vep.sh` | Cross-engine validation (mafsmith + VEP 112 vs. vcf2maf.pl + VEP 112, same indexed cache). Produces the 0-difference result reported in the Results section. |
| `batch_subcommand_diff.sh` | Drives Table 2: validates `maf2vcf` / `vcf2vcf` / `maf2maf` against their Perl references. Calls `maf2vcf_diff.py`, `vcf2vcf_diff.py`, `maf2maf_diff.py`. |
| `maf2vcf_diff.py`, `vcf2vcf_diff.py`, `maf2maf_diff.py` | Subcommand-specific diff utilities used by `batch_subcommand_diff.sh`. |

## Benchmarking (Tables 3 – 5)

| Script | Purpose |
|--------|---------|
| `benchmark_all.py` | Full benchmark harness. Drives both conversion-only (Tables 3, 4) and annotated-pipeline (Table 5) modes. Selects datasets via `--datasets germline,hg008,seqc2` and modes via `--annotated`. |
| `compare_annotations.py` | Compares fastVEP vs. Ensembl VEP CSQ output on the same input VCF; produces the transcript-level / consequence-level concordance numbers cited in the Results section. |
| `../benches/vcf2maf.rs` | Criterion microbenchmarks (`cargo bench`), unrelated to manuscript tables. |

## Data prep / fixtures

| Script | Purpose |
|--------|---------|
| `download_manuscript_vcfs.sh` | Pulls the VCFs listed in Table 1 / Data Availability from their source archives (Synapse, GIAB FTP, SEQC2 FTP, COSMIC, ICGC, DepMap). |
| `generate_expected.sh` | Regenerates the `tests/fixtures/expected/` MAFs from vcf2maf.pl. Used when fixture VCFs change. |

## Generated outputs

| Path | Purpose |
|------|---------|
| `../results/compute_savings.py` | Lives in `results/` (not here) because it ingests the timing tables from `results/` and emits Tables 6 – 9 there. Run via `python results/compute_savings.py`. |

Older or one-off scripts now live in `archive/scripts/`.
