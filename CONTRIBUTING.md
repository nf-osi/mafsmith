# Contributing to mafsmith

## Reporting a caller discrepancy

If you find a variant where mafsmith produces different output from vcf2maf.pl:

1. Confirm the difference is reproducible with `--strict` mode. Some differences in default mode are intentional (see "Known intentional differences" in the README).
2. Isolate the minimal VCF record that reproduces it — a single-record VCF is ideal.
3. Open a GitHub issue with:
   - The VCF record (or a synthetic equivalent that reproduces the bug)
   - The mafsmith output field(s) that differ
   - The vcf2maf.pl output for the same record
   - The caller and VCF FORMAT fields involved

## Adding support for a new caller or edge case

1. Add a representative VCF record to `tests/fixtures/integration_annotated.vcf` that exercises the new behavior.
2. Fix the conversion logic and confirm `tests/fixtures/expected/` matches the expected MAF output.
3. Add a row to the scenario table in `CLAUDE.md` describing the record, expected classification, and what bug or case it covers.
4. Run the integration test suite to confirm no regressions:

```bash
cargo test
```

5. If you have access to a real-world VCF from this caller, run the field-by-field diff script to validate at scale:

```bash
python scripts/vcf_diff.py \
  --mafsmith target/release/mafsmith \
  --vcf2maf vcf2maf.pl \
  --vcf /path/to/caller.vcf.gz \
  --genome grch38
```

## Pull request checklist

- [ ] Integration fixture updated if a new scenario was added
- [ ] `CLAUDE.md` scenario table updated
- [ ] `CHANGELOG.md` entry added under `[Unreleased]`
- [ ] `cargo test` passes
- [ ] `cargo clippy` passes
- [ ] `cargo fmt --check` passes

## Running tests locally

```bash
# Unit and integration tests (fast; uses pre-annotated fixture)
cargo test

# Full integration test with the annotated fixture
target/release/mafsmith vcf2maf \
  -i tests/fixtures/integration_annotated.vcf \
  -o /tmp/out.maf \
  --vcf-tumor-id test_tumor --tumor-id test_tumor \
  --vcf-normal-id test_normal --normal-id test_normal \
  --genome grch38 --skip-annotation
```
