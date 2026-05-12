# archive/

Superseded scripts, older intermediate result files, and untracked scratch
artifacts that are not part of the current mafsmith package, validation
pipeline, or manuscript. Kept under version control so historical commits
remain reproducible.

## What's here

### `scripts/`

| File | Why archived |
|------|--------------|
| `benchmark.py` | Older single-file Python timing harness; superseded by `scripts/benchmark_all.py`. |
| `benchmark.sh` | Bash equivalent of the above; flagged in CLAUDE.md as less reliable on macOS. |
| `compare_maf.py` | Lightweight MAF column-diff helper; superseded by `scripts/vcf_diff.py` (per-caller-aware, handles SV secondary rows, etc.). |
| `run_comparison.sh` | Early entry-point that ran mafsmith against fixtures and compared MAF output. Replaced by the integration tests in `tests/integration_test.rs`. |
| `run_real_world_comparison.sh` | Older Strelka-specific comparison script; replaced by the `batch_vcf_diff_*.sh` family. |
| `benchmark_nf.py` | One-off NF Data Portal benchmark used for an early pass. The current NF validation lives in `scripts/batch_vcf_diff_synapse.sh`. |

### `results/`

| File | Why archived |
|------|--------------|
| `benchmark_2026-04-20.md` | Early draft of the somatic T/N timing table. Superseded by `results/somatic_tn_benchmark.md` (same date, polished form). |
| `benchmark_annotated_2026-04-20.md` | Early draft of the annotated-pipeline timing table. Superseded by the manuscript Table 5 / `results/somatic_tn_benchmark.md`. |
| `benchmark_annotated_2026-04-21.md` | Second draft of the same. |
| `benchmark_nf_20260418_234745.json` | Raw JSON output from `benchmark_nf.py`. |

### `tests-fixtures/`

| File | Why archived |
|------|--------------|
| `test_b38_vep113.vcf_summary.html` | VEP's auto-generated run summary; not consumed by any test. |
| `test_b38_vep113.vcf_warnings.txt` | VEP's auto-generated warnings file; not consumed by any test. |
| `integration_annotated.split.vcf` | Scratch split of `integration_annotated.vcf` from a debugging session; never referenced. |

### Top level

| File | Why archived |
|------|--------------|
| `profile.json.gz` | Output from a one-off `cargo flamegraph` / profiler run. |
