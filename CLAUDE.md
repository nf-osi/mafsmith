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
| 21:46300000 A>G GT=0/1 AD=2,18 | Missense_Mutation | SNP | Het call with VAF=0.9 ≥ 0.7: Tumor_Seq_Allele1 must stay REF. Bug: VAF override fired on explicit GT=0/1; fix: suppress override when GT has both ref (0) and alt (>0) indices |

### Edge case bugs found and fixed via real-VCF comparison

| Bug | Root cause | Fix | Verified on |
|-----|-----------|-----|-------------|
| Ion Torrent multi-allelic: `t_alt_count` wrong for GT=0/2 | `AO` is alt-only (no REF entry at index 0). For VCF allele index N, correct AO element is `AO[N-1]`. Code always took `AO[0]`. | `depth.rs` case 6: use `ao.split(',').nth(alt_vcf_idx.saturating_sub(1))` | syn20443868 (147→0 diffs) |
| `protein_altering_variant` mis-classified as Frame_Shift for in-frame indels with "-" allele | `inframe` computed as `"-".len() - alt.len() = 1 - 9 = 8; 8%3≠0` instead of `0 - 9 = 9; 9%3==0`. Normalized insertions use "-" placeholder (len=1) not empty string (len=0). | `consequence.rs`: treat `"-"` as len=0 in inframe check | syn5614682 (1 diff fixed) |
| Truncated AD alt-index OOB: `t_ref_count`/`t_alt_count` serialized as `''` instead of `'.'` | For multi-allelic records where AD has fewer values than alts (e.g. AD="0,2" with 3 alleles and effective_alt_vcf_idx=2), alt_count=None was serialized as `''`. vcf2maf.pl outputs `'.'`. | `vcf2maf.rs`: added `tumor_alt_oor`/`normal_alt_oor` flags when AD has real values but alt_vcf_idx is OOB; required `ad_len >= 2` to avoid false-positive on VarScan single-value AD | syn5555584, syn5553158 |
| VarScan somatic: all depth fields outputting `'.'`, all TSA1 defaulting to REF | Two bugs: (1) `ad_oor` check fired on VarScan's single-value AD (ad_len=1 <= alt_vcf_idx=1), masking VarScan depth extraction. (2) Depth-based hom_alt TSA1 fallback used `ad_vals.get(alt_vcf_idx)` which fails for VarScan's AD[0] = alt_count. | (1) `vcf2maf.rs` `ad_oor`: add `ad_len >= 2` guard. (2) Replace AD-index hom_alt logic with `extract_depth()` call that handles all callers | syn6840402 (933→0 diffs) |
| Normal sample wrong alt depth for multi-allelic records | `extract_depth` for normal sample always used `alt_vcf_idx=1` (first ALT) regardless of which alt was selected. For GT=1/2 records the normal sample's AD index for the effective alt was never used. | `vcf2maf.rs` line 368: pass `effective_alt` and `effective_alt_vcf_idx` to normal `extract_depth` (same as tumor) | syn5553155 (n_alt_count diffs fixed) |
| SomaticSniper: t_alt_count/n_alt_count off by 1–many | `DP4` alt_fwd+alt_rev counts ALL non-ref reads at a site, not just the specific ALT allele. For T>C records with background G reads, DP4 over-counts. vcf2maf.pl uses `BCOUNT` (per-base A,C,G,T counts) for per-allele ref/alt counts. Both fields are always present in SomaticSniper VCFs. | `depth.rs` case 5: when BCOUNT is present alongside DP4, use BCOUNT[ref_base]/BCOUNT[alt_base] for ref/alt counts; use DP4 sum for total depth only | WGS_FD_1.bwa.somaticSniper.vcf.gz (270→0 diffs) |

### Known remaining differences vs vcf2maf.pl

- **5'UTR vs 5'Flank** (chr1:3857566, chr1:3857574): Transcript selection picks different canonical isoform. Affects ~2 variants per dataset. Under investigation.
- **3'UTR vs 3'Flank** (e.g. chr1:9729848, chr1:944194, chr1:53087567): mafsmith picks one gene's canonical transcript (3'UTR), vcf2maf picks an adjacent gene's transcript (3'Flank). Root cause is same: different gene/transcript selection at gene-boundary regions. Affects ~1–3 variants per dataset.
- **IGR vs 3'Flank / 5'Flank** (e.g. chr7:150676690): vcf2maf selects an Ensembl regulatory feature annotation (ENSR-prefixed ID) over a nearby transcript, yielding `IGR`. mafsmith selects the transcript and correctly classifies as `3'Flank` or `5'Flank`. Root cause: vcf2maf.pl does not filter regulatory-feature CSQ entries from transcript selection. Affects ~1 variant per genome-wide 3,000-variant sample.
- **Intron vs RNA** (chr1:~809967-827588 cluster, and occasionally elsewhere): When the only protein_coding transcript at a site is non-canonical, `select_transcript_light` falls through to pick an lncRNA canonical transcript. vcf2maf.pl reports `Intron` (protein_coding biotype), mafsmith reports `RNA` (lncRNA biotype). The chr1 ~818 kb cluster appears consistently across every dataset tested (CCLE WGS, SEQC2, GIAB HG008 etc.) — typically 5–6 affected variants per 3,000-variant window. Under investigation.
- **Multi-allelic truncated AD (non-strict mode)**: For records with N ALTs where AD has fewer than N+1 values, vcf2maf.pl outputs `'.'` for ref/alt counts even when the selected ALT's AD index is valid. mafsmith (non-strict) extracts the valid counts. In `--strict` mode mafsmith matches vcf2maf.pl. Affects ~3–5% of GATK paired T/N variants.
- **Multi-allelic GVCF tie**: when two ALTs have identical depth, tool-specific tie-breaking may differ. Affects ~4 variants in large GVCF files.
- **SV secondary rows**: mafsmith correctly emits secondary SV breakpoint rows with actual partner chromosome/position. vcf2maf.pl emits them with empty Chromosome (a vcf2maf.pl bug). The `synapse_validate.py` comparison script handles this with best-match selection for key collisions and skips empty-chromosome vcf2maf rows.
- **Unrecognized symbolic ALTs**: records with `<INS>`, `<CNV>` or other non-BND/TRA/DEL/DUP/INV symbolic ALTs are dropped (matching vcf2maf.pl behavior).

### Validated callers (20k variants each, 0 mismatches in --strict mode)

| Synapse ID | Filename pattern | Caller | VCF type |
|------------|-----------------|--------|----------|
| syn31624545 | VA01.vcf.gz | DeepVariant 1.2.0 | single-sample gVCF, GT=`./.' or 0/0 |
| syn64156972 | Hg38_2_024.mt2.Full.vcf | MuTect2 | single-sample GRCh38 |
| syn31624535 | FreeBayes_VA01.vcf.gz | FreeBayes | single-sample |
| syn31624939 | Strelka_VA05_variants.vcf.gz | Strelka2 germline | single-sample variants.vcf (no reference blocks) |
| syn68172710 | patient*.strelka.somatic_indels.vcf.gz | Strelka2 somatic | paired T/N, indel-only |
| syn21296193 | somaticSV.vcf.gz | SV caller (Manta/DELLY) | SV-only |
| syn6840402 | patient_5_*.VarScan_somatic.snp.snpeff.filtered.vcf.gz | VarScan2 somatic | paired T/N, SNP-only |
| syn6039268 | patient_11_normal_*_patient_11_tumor_*.snpeff.vcf.gz | VarDict paired T/N | paired T/N; `RD` FORMAT field (strand bias) coexists with standard `AD` |
| syn31624525 | Mutect2_filtered_VA01.vcf.gz | GATK Mutect2 4.x somatic | single-sample, `GT:AD:AF:DP:F1R2:F2R1:SB` FORMAT |
| syn31624637 | Strelka_VA01_genome.vcf.gz | Strelka2 germline | single-sample genome.vcf, `GT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL` FORMAT |
| HG008-T--HG008-N.mutect2.vcf.gz | GIAB HG008 somatic | GATK Mutect2 4.2.5 paired T/N | normal listed first, has `FAD` (fragment allele depths) FORMAT field |
| HG008-T--HG008-N.snv.strelka2.vcf.gz | GIAB HG008 somatic | Strelka2 2.9.3 somatic SNV | per-base allele counts (AU/CU/GU/TU), no AD field |
| HG008-T--HG008-N.indel.strelka2.vcf.gz | GIAB HG008 somatic | Strelka2 2.9.3 somatic indel | TIR/TAR allele count format, no AD field |
| WGS_FD_1.bwa.muTect2.vcf.gz | SEQC2 HCC1395 | GATK Mutect2 paired T/N | samples named TUMOR/NORMAL |
| WGS_FD_1.bwa.strelka.vcf.gz | SEQC2 HCC1395 | Strelka2 somatic | normal listed first, samples named NORMAL/TUMOR |
| VA02.vcf.gz, VA06.vcf.gz | SMMART DeepVariant | DeepVariant 1.2.0 | single-sample gVCF with `VAF` FORMAT field (absent in VA01); GT=`./.' or 0/0 |
| HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | GIAB germline benchmark | Multi-caller consensus (ADALL field) | uses `ADALL` (all-dataset allele depths) in addition to `AD` |
| HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | GIAB germline benchmark | Multi-caller consensus | same format as HG001 |
| WGS_FD_1.bwa.somaticSniper.vcf.gz | SEQC2 HCC1395 | SomaticSniper paired T/N | `DP4`+`BCOUNT` FORMAT (no `AD`); samples named TUMOR/NORMAL |
| HCC1143.dkfz-snvCalling.somatic.snv_mnv.vcf.gz | ICGC PCAWG (DKFZ cell line) | DKFZ SNV caller | GRCh37, `GT:DP:DP4` FORMAT (no `AD`, no `BCOUNT`); samples named CONTROL/TUMOR; no chr prefix |
| CDS-NMp6b2.vcf.gz, CDS-7nLbDx.vcf.gz, CDS-rBDWH4.vcf.gz (3 of 802) | DepMap CCLE WGS (s3://depmap-omics-ccle) | GATK Mutect2 single-sample | hg38, bcftools-normalized; sample column = cell line name; 802 VCFs total (~163 GB); 0 conversion diffs across 3 samples tested |

## Performance

All benchmarks on AWS c6a.4xlarge (AMD EPYC 7R13, 16 vCPU, 30 GiB RAM). Full per-sample data in `results/`.

### Conversion-only (`--skip-annotation` / `--inhibit-vep`)

**Somatic T/N datasets** (5 datasets: GIAB HG008 Mutect2/Strelka2, SEQC2 Mutect2/Strelka; 3 iterations each):

| Mode | Mean speedup | Throughput |
|------|-------------|-----------|
| mafsmith 1-core vs vcf2maf.pl | 39.1× | 255,807 vs 6,568 variants/s |
| mafsmith 16-core vs vcf2maf.pl | 69.5× | 454,661 vs 6,568 variants/s |

**GIAB germline benchmark** (HG001–HG007, 27.5M total variants; 3 iterations each):

| Mode | Mean speedup | Throughput |
|------|-------------|-----------|
| mafsmith 1-core vs vcf2maf.pl | 79.4× | 558,914 vs 7,036 variants/s |

### Annotated pipeline (mafsmith+fastVEP vs vcf2maf.pl+VEP 115)

Same 5 somatic T/N datasets; fastVEP 16-core Rayon vs VEP `--vep-forks N`:

| VEP forks | Speedup (1-core) | Speedup (16-core) |
|-----------|-----------------|-------------------|
| 16 (symmetric) | 25.3× | 65.5× |
| 4 (typical deployment) | 31.5× | 83.4× |

Detailed tables in `results/somatic_tn_benchmark.md` and `results/conversion_benchmark_giab_grch38.md`.

## Benchmarking scripts

- `scripts/benchmark_all.py` — Full benchmark suite (conversion-only and annotated modes)
  - `--datasets germline,hg008,seqc2 --iterations 3` — conversion-only, 3 runs each
  - `--datasets hg008,seqc2 --annotated --iterations 1` — annotated (fastVEP/VEP), 1 run each
  - Requires `~/.mafsmith/hg38/genes.gff3` (chr-prefixed) for annotated mode; create with:
    `zcat ~/.mafsmith/GRCh38/genes.gff3.gz | sed 's/^##sequence-region   /##sequence-region   chr/; /^#/!s/^/chr/' > ~/.mafsmith/hg38/genes.gff3`
  - Requires VEP 115 GRCh38 indexed cache at `~/.vep/` for vcf2maf.pl annotated mode:
    1. Download: `curl -# https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz -o ~/.vep/homo_sapiens_vep_115_GRCh38.tar.gz`
    2. Extract: `tar -xzf ~/.vep/homo_sapiens_vep_115_GRCh38.tar.gz -C ~/.vep/`
  - Requires Perl module fix (if not done): `conda install -n vcf2maf-env -c conda-forge perl-compress-raw-zlib`
- `scripts/compare_annotations.py` — Compare fastVEP vs VEP CSQ output on the same VCF
- `scripts/benchmark.py` — Older single-file Python timing harness
- `scripts/benchmark.sh` — bash alternative (less reliable on macOS)
- `benches/vcf2maf.rs` — Criterion microbenchmarks (`cargo bench`)
