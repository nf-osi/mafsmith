# results/

Validation and benchmark results referenced by the manuscript. Each manuscript
table maps to a result file in this directory and (where applicable) a script
in `scripts/` or `benches/` that produced it.

## Table → result file → code

| Manuscript table | What it reports | Result file | Generating code |
|------------------|-----------------|-------------|-----------------|
| **Table 1** | mafsmith vs vcf2maf.pl conversion validation across 15 caller types, 87.2 M variants. | Per-caller diff summaries in this directory plus `pcawg_inhibit_vep_validation.md` (PCAWG subset). | `scripts/batch_vcf_diff_synapse.sh` (NF Data Portal); `scripts/batch_vcf_diff_ccle.sh` and `batch_vcf_diff_ccle_inhibit_vep.sh` (DepMap CCLE); `scripts/batch_vcf_diff_pcawg.sh` and `batch_vcf_diff_pcawg_inhibit_vep.sh` (PCAWG); `scripts/batch_vcf_diff_manuscript.sh` (consolidates per-caller runs). Each batch script calls `scripts/vcf_diff.py` for the field-by-field comparison. |
| **Table 2** | `maf2vcf` / `vcf2vcf` / `maf2maf` subcommand validation against the Perl references. | `pcawg_inhibit_vep_validation.md` (subcommand-by-subcommand counts inline). | `scripts/batch_subcommand_diff.sh` (entry point) which calls `scripts/maf2vcf_diff.py`, `scripts/vcf2vcf_diff.py`, and `scripts/maf2maf_diff.py`. |
| **Table 3** | Conversion-only timing on GIAB HG001–HG007 (7 samples × 3 iterations each). | `conversion_benchmark_giab_grch38.md`. | `scripts/benchmark_all.py --datasets germline --iterations 3`. |
| **Table 4** | Conversion-only timing on five somatic T/N datasets (GIAB HG008 + SEQC2 HCC1395). | `somatic_tn_benchmark.md` (conversion-only section). | `scripts/benchmark_all.py --datasets hg008,seqc2 --iterations 3`. |
| **Table 5** | Full annotated-pipeline timing (mafsmith + fastVEP vs. vcf2maf.pl + VEP 115 `--vep-forks 16`) on the same five datasets. | `somatic_tn_benchmark.md` (annotated-pipeline section). | `scripts/benchmark_all.py --datasets hg008,seqc2 --annotated --iterations 1`. |
| **Tables 6, 7** | Per-sample and at-scale cost / energy / CO₂ savings for the conversion step. | Derived from `conversion_benchmark_giab_grch38.md` and AWS / EPA eGRID inputs hard-coded in the script. | `compute_savings.py` (`python results/compute_savings.py`). |
| **Tables 8, 9** | Per-run and at-scale cost / energy / CO₂ savings for the full annotated pipeline. | Derived from `somatic_tn_benchmark.md` plus the same AWS / EPA eGRID inputs. | `compute_savings.py`. |

Annotation concordance numbers in the Results section (fastVEP vs. VEP 115,
~26% transcript-level discordance / 87.6% consequence-level agreement) come
from `scripts/compare_annotations.py`.

## File-by-file

| File | Purpose |
|------|---------|
| `conversion_benchmark_giab_grch38.md` | Full per-sample timing for the GIAB germline benchmarks (Table 3). |
| `somatic_tn_benchmark.md` | Per-sample timing for the somatic T/N benchmarks; covers both the conversion-only (Table 4) and annotated (Table 5) modes. |
| `pcawg_inhibit_vep_validation.md` | PCAWG-specific validation subset plus the field-by-field counts that back Table 2. |
| `compute_savings.py` | Single script that ingests the timing tables above and emits Tables 6–9. |
| `README.md` | This file. |

Older intermediate drafts of the same benchmarks live in `archive/results/`.
