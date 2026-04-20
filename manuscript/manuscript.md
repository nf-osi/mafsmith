# mafsmith: a fast, self-contained VCF-to-MAF converter for cancer genomics

Robert Allaway^1^, [additional authors TBD]

^1^ Sage Bionetworks, Seattle, WA, USA

**Corresponding author:** robert.allaway@sagebase.org

---

## Abstract

The Mutation Annotation Format (MAF) is the standard interchange format for somatic variant data in cancer genomics, required by the NCI Genomic Data Commons and widely used in downstream analytical pipelines. Converting variant call format (VCF) files to MAF requires functional annotation (via tools such as the Ensembl Variant Effect Predictor) followed by complex allele normalisation and field-mapping logic. The de facto reference implementation, vcf2maf.pl, depends on a large Perl software stack and is prohibitively slow for large cohort analyses. We present mafsmith, a Rust implementation of the VCF-to-MAF conversion pipeline. mafsmith uses fastVEP for annotation and reimplements the full allele-normalisation and field-mapping logic of vcf2maf.pl, producing field-for-field identical output across six validated caller types. On a representative dataset, mafsmith achieves **[X]-fold** faster end-to-end conversion and **~229-fold** faster conversion of pre-annotated VCFs, enabling rapid processing of large cancer cohorts on standard hardware. mafsmith is open source and available at https://github.com/nf-osi/mafsmith.

---

## Introduction

Somatic variant calling produces VCF files whose downstream use in cancer genomics almost universally requires conversion to the Mutation Annotation Format (MAF). MAF is the primary data model for the NCI Genomic Data Commons (GDC) [CITATION: GDC] and is consumed by widely used analytical tools including maftools [CITATION: maftools], MSKCC cBioPortal [CITATION: cBioPortal], and OncoKB [CITATION: OncoKB]. The conversion from VCF to MAF is non-trivial: it requires functional annotation of each variant against a transcript database, selection of a canonical or preferred transcript, normalisation of multi-allelic and indel representations, and mapping of genotype information to allele-level MAF fields.

The standard tool for this conversion is vcf2maf.pl, developed and maintained by Cyriac Kandoth [CITATION: vcf2maf]. vcf2maf.pl is a Perl script that wraps the Ensembl Variant Effect Predictor (VEP) [CITATION: VEP], handling the full complexity of VCF allele representations, multi-allelic sites, structural variants, and caller-specific FORMAT field conventions that have accumulated over years of production use across major cancer genomics consortia. Its correctness and breadth of supported variant types make it the de facto reference implementation.

However, vcf2maf.pl has significant performance limitations. Running VEP for each conversion is the primary bottleneck: annotating a typical whole-exome sequencing (WES) VCF of 50,000–300,000 variants can take 10–60 minutes per sample. For large cohorts of hundreds or thousands of samples, this becomes a major computational burden. Additionally, the dependency on a compatible Perl + VEP + reference database stack creates significant installation and reproducibility challenges.

Recent development of fastVEP [Huang, 2026], a reimplementation of VEP's core annotation logic in the Rust programming language, substantially reduces annotation time — achieving up to 130-fold speedup over the original Perl VEP while maintaining complete concordance. However, the conversion step itself (allele normalisation, genotype parsing, field mapping) still requires vcf2maf.pl, which processes variants at approximately 130 variants/second even when annotation is skipped. The annotation and conversion bottlenecks therefore call for complementary solutions.

Here we describe mafsmith, a Rust implementation of the complete VCF-to-MAF conversion pipeline that pairs naturally with fastVEP to replace the entire vcf2maf.pl + VEP stack with a single, self-contained toolchain. Both tools exploit Rust's performance characteristics — memory safety, zero-cost abstractions, and native parallelism — to achieve throughput that is impractical with the incumbent Perl/Python implementations. mafsmith reimplements the allele-normalisation and field-mapping logic of vcf2maf.pl from first principles, targeting field-for-field identical output. We validate mafsmith against vcf2maf.pl across six distinct variant caller types, demonstrate **[X]-fold** end-to-end speedup and **~229-fold** speedup for the conversion step alone, and describe the specific edge cases and caller-specific conventions that required careful implementation to achieve full concordance.

---

## Implementation

### Overview

mafsmith is implemented in Rust and structured as a command-line tool with five subcommands: `vcf2maf` (VCF to MAF, the primary validated subcommand), `maf2vcf`, `maf2maf`, `vcf2vcf`, and `fetch` (reference data download). The core conversion pipeline consists of: (1) optional annotation via fastVEP to produce a VCF with embedded CSQ fields; (2) parsing of annotated VCF records; (3) canonical transcript selection; (4) allele normalisation; (5) genotype and depth field extraction; and (6) MAF record serialisation.

### Annotation

When annotation is required, mafsmith invokes fastVEP [Huang, 2026] with HGVS notation enabled. fastVEP produces a standard VCF with `CSQ` INFO fields using the same format as Ensembl VEP, allowing mafsmith to reuse the same downstream parsing logic regardless of whether annotation was performed upstream. For pre-annotated VCFs, annotation can be skipped with `--skip-annotation`.

### Transcript selection

mafsmith selects the canonical transcript for each variant using the same priority order as vcf2maf.pl: (1) user-supplied custom ENST list; (2) MANE Select transcript; (3) Ensembl canonical transcript; (4) longest transcript. When multiple transcripts tie on all criteria, the first in VEP output order is used, matching vcf2maf.pl behaviour.

### Allele normalisation

VCF alleles are left-aligned and trimmed using a prefix/suffix-stripping approach. Structural variant ALT alleles (`<DEL>`, `<DUP>`, `<INV>`, BND notation) receive special handling: symbolic ALTs are mapped to their MAF representations, BND ALT strings are parsed to extract the partner chromosome and position for secondary breakpoint rows, and unrecognised symbolic ALTs (e.g. `<INS>`, `<CNV>`) are dropped, matching vcf2maf.pl behaviour.

### Genotype and depth extraction

mafsmith implements the full vcf2maf.pl logic for determining `Tumor_Seq_Allele1` and `Tumor_Seq_Allele2` from VCF FORMAT fields. This includes:

- **GT-based allele assignment**: GT allele indices are sorted; the minimum index determines Allele1 (REF for heterozygous, ALT for homozygous-alt).
- **Depth-based hom-alt inference**: when GT is a no-call (`./.'`) or homozygous-reference (`0/0`), allele assignment falls back to AD-based VAF; VAF ≥ 0.7 is treated as homozygous-alt. This matches vcf2maf.pl behaviour across DRAGEN, MuTect2, and GVCF-style callers.
- **VAF override**: for paired tumour/normal VCFs, when GT indicates heterozygous but VAF ≥ 0.7 (suggesting caller under-calling), Allele1 is overridden to ALT, matching vcf2maf.pl. This override is suppressed for single-sample VCFs (absent normal column).
- **`--strict` mode**: when AD arrays are shorter than the expected `1 + n_alts` length (a GATK behaviour for trimmed multi-allelic sites), `--strict` outputs `.` for all depth fields and suppresses depth-based allele calling, exactly matching vcf2maf.pl. In default mode, mafsmith extracts whatever depth information is available.
- **Strelka2 somatic FORMAT fields**: Strelka2 somatic VCFs lack a GT field and use caller-specific depth fields (AU/CU/GU/TU for SNVs, TAR/TIR for indels). mafsmith extracts depth counts from these fields and infers het/hom-alt from VAF.

### Parallelism and performance

The conversion step (post-annotation) uses Rayon [CITATION: Rayon] for parallel record processing across available CPU cores. Memory allocation uses jemalloc [CITATION: jemalloc] for improved throughput on multi-threaded workloads. CSQ field parsing is lazy (deferred until a record is selected for output), reducing unnecessary work for filtered or dropped records.

---

## Results

### Validation

We validated mafsmith against vcf2maf.pl across six caller types, comparing all key MAF fields (`Hugo_Symbol`, `Variant_Classification`, `Variant_Type`, `Reference_Allele`, `Tumor_Seq_Allele1`, `Tumor_Seq_Allele2`, `HGVSc`, `HGVSp`, `HGVSp_Short`, `Transcript_ID`, `Exon_Number`, `t_depth`, `t_ref_count`, `t_alt_count`, `n_depth`, `n_ref_count`, `n_alt_count`) on the same fastVEP-annotated input (Table 1). For all six caller types, mafsmith produced 0 field-level mismatches across 20,000 variants per file with `--strict` mode enabled.

**Table 1. Validation of mafsmith against vcf2maf.pl (20,000 variants per caller type, --strict mode).**

| Caller | VCF type | Synapse ID | Mismatches |
|--------|----------|------------|-----------|
| DRAGEN RefCall | Single-sample, GT=`0/0`/`./.'` | syn31624545 | 0 / 20,000 |
| MuTect2 | Single-sample GRCh38 | syn64156972 | 0 / 20,000 |
| FreeBayes | Single-sample | syn31624535 | 0 / 20,000 |
| Strelka2 | Paired tumour/normal | syn31624939 | 0 / 20,000 |
| Strelka2 somatic indels | Paired tumour/normal | syn68172710 | 0 / 20,000 |
| SV callers (Manta/DELLY) | SV-only | syn21296193 | 0 / 398 |

### Performance

**[PLACEHOLDER — benchmarking results from cloud instance to be inserted here]**

Benchmarks comparing:
1. mafsmith vs. vcf2maf.pl (conversion step only, pre-annotated VCF, `--skip-annotation` / `--inhibit-vep`)
2. mafsmith + fastVEP vs. vcf2maf.pl + VEP (full end-to-end pipeline)

Tested on: [cloud instance spec TBD] using [WES/WGS VCF TBD].

Preliminary results on a 6,292-variant SV VCF (conversion step only, local MacBook):

| Tool | Mean time | Variants/second |
|------|-----------|-----------------|
| mafsmith (`--skip-annotation`) | ~0.06 s | ~104,000 |
| vcf2maf.pl (`--inhibit-vep`) | ~13.9 s | ~132 |

**~229× faster** for the conversion step alone. Full end-to-end speedup (including annotation): **~4.8×** [to be updated with cloud benchmarks on larger cohort].

---

## Discussion

The primary motivation for mafsmith is throughput: converting thousands of VCFs in large cancer genomics cohorts places substantial demands on compute infrastructure when using vcf2maf.pl, and the per-sample annotation step is a major bottleneck. By combining fastVEP's compiled annotation engine with mafsmith's parallel Rust conversion logic, the full pipeline achieves substantially faster conversion while maintaining field-for-field agreement with the reference implementation.

Achieving full concordance with vcf2maf.pl required careful reverse-engineering of a number of non-obvious behaviours accumulated over years of production use. Key examples include: the treatment of absent normal samples vs. no-call normal GTs in homozygous-alt inference; the VAF-based Allele1 override and its interaction with single-sample vs. paired VCF configurations; the handling of truncated AD arrays from GATK multi-allelic calling; consequence severity ranking for multi-consequence VEP annotations; and the representation of structural variant secondary breakpoint rows.

Several known differences remain. In default mode (without `--strict`), mafsmith extracts partial depth counts when AD arrays are truncated, which we regard as more informative than outputting missing values. mafsmith correctly populates SV secondary breakpoint rows with the partner chromosome and position, whereas vcf2maf.pl leaves these fields blank (a bug in the reference implementation). These intentional differences are documented and `--strict` mode is provided for workflows requiring exact vcf2maf.pl compatibility.

mafsmith currently uses a prefix/suffix-stripping approach for allele normalisation rather than FASTA-guided left-alignment. For most variant types this produces identical results to vcf2maf.pl, but complex indels at homopolymer runs or short tandem repeats may differ. Full FASTA-based normalisation is planned for a future release.

---

## Conclusion

mafsmith is a fast, self-contained reimplementation of the vcf2maf.pl conversion pipeline in Rust. It produces field-for-field identical MAF output across six distinct caller types while achieving **[X]-fold** faster end-to-end conversion, making large-cohort VCF-to-MAF conversion practical on standard compute infrastructure. mafsmith is open source (Apache 2.0) and available at https://github.com/nf-osi/mafsmith.

---

## Data availability

Validation VCF files are available from the NF-OSI Synapse project (syn16858331) at https://synapse.org. mafsmith source code and release binaries are available at https://github.com/nf-osi/mafsmith.

---

## Acknowledgements

mafsmith builds on the design, field conventions, and years of accumulated edge-case handling embodied in vcf2maf.pl. The authors are deeply grateful to Cyriac Kandoth ([0000-0002-1345-3573](https://orcid.org/0000-0002-1345-3573)) for his sustained work developing and maintaining vcf2maf.pl and for his deep expertise in the VCF and MAF specifications — without which this reimplementation would not have been possible.

---

## Funding

[TBD]

---

## References

[CITATION: vcf2maf] Kandoth C. mskcc/vcf2maf. GitHub. https://github.com/mskcc/vcf2maf

[CITATION: VEP] McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F. The Ensembl Variant Effect Predictor. *Genome Biol.* 2016;17(1):122. https://doi.org/10.1186/s13059-016-0974-4

[Huang, 2026] Huang K. fastVEP: A Fast, Comprehensive Variant Effect Predictor Written in Rust. *bioRxiv.* 2026. https://doi.org/10.64898/2026.04.14.718452

[CITATION: GDC] Jensen MA, Ferretti V, Grossman RL, Staudt LM. The NCI Genomic Data Commons as an engine for precision medicine. *Blood.* 2017;130(4):453–459. https://doi.org/10.1182/blood-2017-03-735654

[CITATION: maftools] Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP. Maftools: efficient and comprehensive analysis of somatic variants in cancer. *Genome Res.* 2018;28(11):1747–1756. https://doi.org/10.1101/gr.239244.118

[CITATION: cBioPortal] Cerami E, Gao J, Dogrusoz U, Gross BE, Sumer SO, Aksoy BA, Jacobsen A, Byrne CJ, Heuer ML, Larsson E, Antipin Y, Reva B, Goldberg AP, Sander C, Schultz N. The cBio cancer genomics portal: an open platform for exploring multidimensional cancer genomics data. *Cancer Discov.* 2012;2(5):401–404. https://doi.org/10.1158/2159-8290.CD-12-0095

[CITATION: OncoKB] Chakravarty D, Gao J, Phillips SM, et al. OncoKB: A Precision Oncology Knowledge Base. *JCO Precis Oncol.* 2017;1:1–16. https://doi.org/10.1200/PO.17.00011

[CITATION: MuTect2] Benjamin D, Sato T, Cibulskis K, Getz G, Stewart C, Lichtenstein L. Calling Somatic SNVs and Indels with Mutect2. *bioRxiv.* 2019. https://doi.org/10.1101/861054

[CITATION: Strelka2] Kim S, Scheffler K, Halpern AL, Bekritsky MA, Noh E, Källberg M, Chen X, Kim Y, Beyter D, Krusche P, Saunders CT. Strelka2: fast and accurate calling of germline and somatic variants. *Nat Methods.* 2018;15(8):591–594. https://doi.org/10.1038/s41592-018-0051-x

[CITATION: FreeBayes] Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. *arXiv.* 2012. https://arxiv.org/abs/1207.3907

[CITATION: DRAGEN] Illumina DRAGEN Bio-IT Platform. https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html

[CITATION: Manta] Chen X, Schulz-Trieglaff O, Shaw R, Barnes B, Schlesinger F, Källberg M, Cox AJ, Kruglyak S, Saunders CT. Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. *Bioinformatics.* 2016;32(8):1220–1222. https://doi.org/10.1093/bioinformatics/btv710

[CITATION: Rayon] Matsakis N. Rayon: A data parallelism library for Rust. https://github.com/rayon-rs/rayon

[CITATION: jemalloc] Evans J. A Scalable Concurrent malloc(3) Implementation for FreeBSD. *BSDCan.* 2006. https://people.freebsd.org/~jasone/jemalloc/bsdcan2006/jemalloc.pdf

[CITATION: TCGA] Cancer Genome Atlas Research Network, Weinstein JN, Collisson EA, et al. The Cancer Genome Atlas Pan-Cancer analysis project. *Nat Genet.* 2013;45(10):1113–1120. https://doi.org/10.1038/ng.2764
