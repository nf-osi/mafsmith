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

| Record | Variant_Classification | Variant_Type | What it tests / bug it covers |
|--------|----------------------|--------------|-------------------------------|
| 21:34792292 A>C | Missense_Mutation | SNP | Multi-transcript CSQ → canonical MANE transcript selected |
| 21:34880555 A>ACCTCTT | Splice_Site | INS | Insertion classified correctly |
| 21:34880650 T>TG | Frame_Shift_Ins | INS | Frameshift insertion |
| 21:41479182 GT>G | Frame_Shift_Del | DEL | Frameshift deletion |
| 21:41479219 C>T | Silent | SNP | Synonymous SNV |
| 21:44237111 C>T | Nonsense_Mutation | SNP | Stop gained |
| 21:34880300 G>A,C GT=0/0 | Missense_Mutation | SNP | Multi-allelic GVCF: picks highest-depth alt (C, 25 reads) not A (3 reads). Bug: mafsmith defaulted to ALT[0] instead of highest-depth alt for hom-ref GT calls |
| 1:1000000 `<DEL>` CHR2+END | Intron | INS | SV DEL emits secondary breakpoint row. Bug: only 1 row was emitted |
| 1:2000000 `<DUP:TANDEM>` | Intron | INS | ALT normalized to `<DUP>`. Bug: `<DUP:TANDEM>` was passed through verbatim |
| 1:3000000 `G]2:5000000]` BND CHR2=2 | Intron | INS | BND secondary row on chr2; HGVSc=`5000000]` (greedy last-colon strip). Bug: used first colon |
| 1:4000000 `T]2:6000000]` BND no CHR2/END in INFO | Intron | INS | BND secondary row parsed from ALT notation. Bug: only INFO was checked; fix: parse `chr:pos` from `]chr:pos]` or `[chr:pos[` in ALT when CHR2 absent |
| 1:5000000 `<INV>` CHR2+END | Intron | INS | INV secondary breakpoint row |
| 1:6000000 `<DEL>` multi-consequence | Splice_Site | INS | `feature_truncation&splice_acceptor_variant` → Splice_Site. Bug: first consequence short-circuited; fix: sort by severity |
| 21:45000000 GGCT>G | In_Frame_Del | DEL | 3 bp in-frame deletion |
| 21:45100000 G>GCAG | In_Frame_Ins | INS | 3 bp in-frame insertion |
| 21:45200000 T>C | Nonstop_Mutation | SNP | stop_lost |
| 21:45300000 A>G | Translation_Start_Site | SNP | start_lost |
| 21:45400000 A>G | Splice_Region | SNP | splice_region_variant (also tests multi-consequence sort) |
| 21:45500000 C>T | 3'UTR | SNP | 3_prime_UTR_variant |
| 21:45600000 G>A | 5'UTR | SNP | 5_prime_UTR_variant |
| 21:45700000 T>C | 5'Flank | SNP | upstream_gene_variant |
| 21:45800000 A>G | 3'Flank | SNP | downstream_gene_variant |
| 21:45900000 G>A | Intron | SNP | intron_variant (SNV, not SV) |
| 21:46000000 C>T | IGR | SNP | intergenic_variant |
| 21:46100000 G>A | RNA | SNP | non_coding_transcript_exon_variant (lncRNA biotype) |
| 21:46200000 AC>GT | Missense_Mutation | DNP | Dinucleotide polymorphism Variant_Type |

### Known remaining differences vs vcf2maf.pl

- **Multi-allelic GVCF tie**: when two ALTs have identical depth, tool-specific tie-breaking may differ. Affects ~4 variants in large GVCF files.
- **MISSING in vcf2maf**: occasionally mafsmith emits a secondary SV row that vcf2maf.pl does not (e.g. when CHR2 is present but vcf2maf.pl's perl regex fails). These are correct mafsmith rows.

## Performance

Benchmarked on a 6,292-variant SV VCF (post-fastVEP annotation):

| Tool | Mean time | Variants/s |
|------|-----------|-----------|
| mafsmith (`--skip-annotation`) | ~0.06s | ~104,000 |
| vcf2maf.pl (`--inhibit-vep`) | ~13.9s | ~132 |

**~229× faster** for the conversion step. Full pipeline speedup (including fastVEP): ~4.8×.

Previously benchmarked on a 305K-variant annotated VCF (578 MB): mafsmith ~6.2s / ~49,000 variants/s (pre-optimization). On the `optimize/perf` branch this improves to ~5× due to rayon parallel processing and allocation reduction.

## Benchmarking scripts

- `scripts/benchmark.py` — Python timing harness; recommended
- `scripts/benchmark.sh` — bash alternative (less reliable on macOS)
- `benches/vcf2maf.rs` — Criterion microbenchmarks (`cargo bench`)
