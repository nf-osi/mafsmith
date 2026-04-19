# mafsmith — development notes for Claude

## Integration test fixture

`tests/fixtures/integration_annotated.vcf` is a pre-annotated (post-fastVEP) VCF designed to
cover every conversion scenario that has required a bug fix. Run it with `--skip-annotation`
to exercise conversion logic in isolation:

```bash
target/release/mafsmith vcf2maf \
  -i tests/fixtures/integration_annotated.vcf \
  -o /tmp/out.maf \
  --vcf-tumor-id test_tumor --tumor-id test_tumor \
  --vcf-normal-id test_normal --normal-id test_normal \
  --genome grch38 --skip-annotation
```

### Scenarios covered and why each was added

| Record | What it tests | Bug that prompted it |
|--------|--------------|----------------------|
| 21:34792292 A>C | SNV missense, multi-transcript CSQ → picks MANE/canonical | baseline |
| 21:34880555 A>ACCTCTT | Insertion, Splice_Site classification | baseline |
| 21:34880650 T>TG | Frameshift insertion | baseline |
| 21:41479182 GT>G | Frameshift deletion | baseline |
| 21:41479219 C>T | Synonymous SNV | baseline |
| 21:44237111 C>T | Stop gained | baseline |
| 21:34880300 G>A,C (GT=0/0) | Multi-allelic GVCF site: picks highest-depth alt (C, 25 reads) over A (3 reads) | mafsmith was defaulting to first alt instead of highest-depth alt in GVCF hom-ref calls |
| 1:1000000 `<DEL>` CHR2+END | SV DEL: emits secondary breakpoint row at END position | mafsmith only emitted 1 MAF row per SV; vcf2maf.pl emits 2 |
| 1:2000000 `<DUP:TANDEM>` | ALT normalized to `<DUP>` (SVTYPE=DUP) | `<DUP:TANDEM>` was passed through verbatim; vcf2maf.pl rewrites to `<SVTYPE>` |
| 1:3000000 `G]2:5000000]` BND (CHR2=2) | BND secondary row on chr2; HGVSc greedy-stripped to `5000000]` | HGVSc was stripped at first colon; vcf2maf.pl uses last colon (greedy `s/^.*://`) |
| 1:4000000 `T]2:6000000]` BND (no CHR2) | BND secondary row with empty chromosome (Manta-style, no CHR2 in INFO) | same missing-secondary-row bug |
| 1:5000000 `<INV>` CHR2+END | INV secondary row | same missing-secondary-row bug |
| 1:6000000 `<DEL>` multi-consequence | `feature_truncation&splice_acceptor_variant` → `Splice_Site` (not `Targeted_Region`) | consequences were evaluated in VEP output order; `feature_truncation` appeared first and short-circuited; fix: sort by severity before classifying |

### Known remaining differences vs vcf2maf.pl (not fixed, not fixable)

- **BND with no CHR2 and no END**: vcf2maf.pl emits a secondary row at `('', '')`. mafsmith cannot generate this row without a position. One per typical somaticSV file.
- **Multi-allelic GVCF tie**: when two ALTs have identical depth, tool-specific tie-breaking differs. Affects ~4 variants in large GVCF files.
- **MISSING in vcf2maf**: occasionally mafsmith emits a secondary SV row that vcf2maf.pl does not (e.g. when CHR2 is present but vcf2maf.pl's perl regex fails). These are correct mafsmith rows.

## Performance

Benchmarked on a 305K-variant annotated VCF (578 MB):

| Tool | Mean time | Variants/s |
|------|-----------|-----------|
| mafsmith (`--skip-annotation`) | ~6.2s | ~49,000 |
| vcf2maf.pl (`--inhibit-vep`) | ~183s | ~1,665 |

**~30× faster** for the conversion step. Full pipeline speedup (including fastVEP) is lower but still substantial since fastVEP itself is fast.

## Benchmarking scripts

- `scripts/benchmark.py` — Python timing harness; recommended
- `scripts/benchmark.sh` — bash alternative (less reliable on macOS)
- `benches/vcf2maf.rs` — Criterion microbenchmarks (`cargo bench`)
