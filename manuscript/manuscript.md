# mafsmith: a fast, self-contained VCF-to-MAF converter for cancer genomics

Robert Allaway^1^, [additional authors TBD]

^1^ Sage Bionetworks, Seattle, WA, USA

**Corresponding author:** robert.allaway@sagebase.org

---

## Abstract

The Mutation Annotation Format (MAF) is the standard interchange format for somatic variant data in cancer genomics, required by the NCI Genomic Data Commons and widely used in downstream analytical pipelines. Converting variant call format (VCF) files to MAF requires functional annotation (via tools such as the Ensembl Variant Effect Predictor) followed by complex allele normalisation and field-mapping logic. The de facto reference implementation, vcf2maf.pl, depends on a large Perl software stack and is prohibitively slow for large cohort analyses. We present mafsmith, a Rust implementation of the VCF-to-MAF conversion pipeline. mafsmith uses fastVEP for annotation and reimplements the full allele-normalisation and field-mapping logic of vcf2maf.pl, producing field-for-field identical output across thirteen validated caller types and formats spanning germline, somatic, and structural variant VCFs. Benchmarked on seven GIAB reference samples totalling 27.5 million variants, mafsmith achieves **79.4-fold** faster conversion of pre-annotated VCFs (range 74.3–84.1×), enabling rapid processing of large cancer cohorts on standard hardware. mafsmith is open source and available at https://github.com/nf-osi/mafsmith.

---

## Introduction

Somatic variant calling produces VCF files whose downstream use in cancer genomics almost universally requires conversion to the Mutation Annotation Format (MAF). MAF is the primary data model for the NCI Genomic Data Commons (GDC) [CITATION: GDC] and is consumed by widely used analytical tools including maftools [CITATION: maftools], MSKCC cBioPortal [CITATION: cBioPortal], and OncoKB [CITATION: OncoKB]. The conversion from VCF to MAF is non-trivial: it requires functional annotation of each variant against a transcript database, selection of a canonical or preferred transcript, normalisation of multi-allelic and indel representations, and mapping of genotype information to allele-level MAF fields.

The standard tool for this conversion is vcf2maf.pl, developed and maintained by Cyriac Kandoth [CITATION: vcf2maf]. vcf2maf.pl is a Perl script that wraps the Ensembl Variant Effect Predictor (VEP) [CITATION: VEP], handling the full complexity of VCF allele representations, multi-allelic sites, structural variants, and caller-specific FORMAT field conventions that have accumulated over years of production use across major cancer genomics consortia. Its correctness and breadth of supported variant types make it the de facto reference implementation.

However, vcf2maf.pl has significant performance limitations. Running VEP for each conversion is the primary bottleneck: annotating a typical whole-exome sequencing (WES) VCF of 50,000–300,000 variants can take 10–60 minutes per sample. For large cohorts of hundreds or thousands of samples, this becomes a major computational burden. Additionally, the dependency on a compatible Perl + VEP + reference database stack creates significant installation and reproducibility challenges.

Recent development of fastVEP [Huang, 2026], a reimplementation of VEP's core annotation logic in the Rust programming language, substantially reduces annotation time — achieving up to 130-fold speedup over the original Perl VEP while maintaining complete concordance. However, the conversion step itself (allele normalisation, genotype parsing, field mapping) still requires vcf2maf.pl, which processes variants at approximately 130 variants/second even when annotation is skipped. The annotation and conversion bottlenecks therefore call for complementary solutions.

Here we describe mafsmith, a Rust implementation of the complete VCF-to-MAF conversion pipeline that pairs naturally with fastVEP to replace the entire vcf2maf.pl + VEP stack with a single, self-contained toolchain. Both tools exploit Rust's performance characteristics — memory safety, zero-cost abstractions, and native parallelism — to achieve throughput that is impractical with the incumbent Perl/Python implementations. mafsmith reimplements the allele-normalisation and field-mapping logic of vcf2maf.pl from first principles, targeting field-for-field identical output. We validate mafsmith against vcf2maf.pl across thirteen distinct variant caller types and formats spanning germline, somatic, and structural variant VCFs, demonstrate **79.4-fold** speedup for the conversion step alone (range 74.3–84.1× across seven GIAB reference samples), and describe the specific edge cases and caller-specific conventions that required careful implementation to achieve full concordance.

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

We validated mafsmith against vcf2maf.pl across thirteen caller types and VCF formats, comparing all key MAF fields (`Hugo_Symbol`, `Variant_Classification`, `Variant_Type`, `Reference_Allele`, `Tumor_Seq_Allele1`, `Tumor_Seq_Allele2`, `HGVSc`, `HGVSp`, `HGVSp_Short`, `Transcript_ID`, `Exon_Number`, `t_depth`, `t_ref_count`, `t_alt_count`, `n_depth`, `n_ref_count`, `n_alt_count`) on the same fastVEP-annotated input (Table 1). For all caller types, mafsmith produced 0 conversion-field mismatches in `--strict` mode across a total of 57.6 million variants. Remaining differences between tools are restricted to `Variant_Classification` for ~2–5 variants per dataset at gene-boundary regions where tools select different canonical transcripts, reflecting different Ensembl gene model versions rather than conversion logic.

**Table 1. Validation of mafsmith against vcf2maf.pl (full variant sets per dataset; --strict mode; 0 conversion-field mismatches in all cases; 57.6 million variants compared in total).**

| Caller | VCF type / FORMAT fields | Source | Variants compared |
|--------|--------------------------|--------|-------------------|
| DeepVariant 1.2.0 | Single-sample gVCF (GT=`0/0`/`./.'`) | syn31624545 | 2,936,426 |
| GATK MuTect2 | Single-sample GRCh38 | syn31624525 | 751,548 |
| GATK MuTect2 (paired T/N) | `GT:AD:AF:DP:F1R2:F2R1:SB` FORMAT | GIAB HG008; SEQC2 HCC1395 | 277,645; 271,945 |
| FreeBayes | Single-sample | syn31624535 | 2,858,303 |
| Strelka2 germline | Single-sample `variants.vcf` and `genome.vcf` | syn31624939; syn31624637 | 5,841,306; 3,486,411 |
| Strelka2 somatic SNVs | Paired T/N; per-base `AU/CU/GU/TU` depth fields | GIAB HG008; SEQC2 HCC1395 | 1,562,847; 2,191,720 |
| Strelka2 somatic indels | Paired T/N; `TAR`/`TIR` depth fields | syn68172710; GIAB HG008 | 23,375; 293,719 |
| SV callers (Manta/DELLY) | SV-only; BND, DEL, DUP, INV symbolic ALTs | syn21296193 | 398 (all SVs) |
| VarScan2 somatic | Paired T/N; `RD`+`AD` FORMAT | syn6840402 | 59,618 |
| VarDict (paired T/N) | `RD` strand-bias field coexists with standard `AD` | syn6039268 | 9,303,064 |
| SomaticSniper | Paired T/N; `DP4`+`BCOUNT` FORMAT (no `AD`) | SEQC2 HCC1395 | 164,704 |
| GIAB germline benchmarks | Multi-caller consensus; `ADALL` field; GRCh38 | HG001–HG007 (NIST v4.2.1) | 3,839,315–4,048,342 each (27,529,001 total) |
| COSMIC v103 | Annotation database VCF (no sample columns); GRCh38 | GenomeScreensMutant; NonCodingVariants | 3,000 each |

### Performance

We benchmarked mafsmith against vcf2maf.pl on the conversion step in isolation, passing the same raw VCF to each tool (`mafsmith --skip-annotation` and `vcf2maf.pl --inhibit-vep`) without any prior annotation. Benchmarks were run on an AWS c6a.4xlarge instance (AMD EPYC 7R13, 16 vCPU, 30 GiB RAM) using all seven GIAB NIST v4.2.1 GRCh38 benchmark VCFs (HG001–HG007; 3.84–4.05 million variants per file; 27.5 million variants total). To distinguish algorithmic speedup from parallelism, mafsmith was timed at both 1 core (`RAYON_NUM_THREADS=1`, 3 runs per sample) and 16 cores (default Rayon, 1 run per sample); vcf2maf.pl is single-threaded.

**Table 2. Conversion-only benchmark: mafsmith vs vcf2maf.pl on GIAB NIST v4.2.1 GRCh38.**

| Sample | Variants | mafsmith 1-core (s) | mafsmith 16-core (s) | vcf2maf.pl (s) | Speedup (1-core) | Speedup (16-core) |
|--------|----------|---------------------|----------------------|----------------|------------------|-------------------|
| HG001 (NA12878) | 3,893,341 | 11.579 | 7.398 | 549.988 | 47.5× | 74.3× |
| HG002 (NA24385) | 4,048,342 | 12.314 | 7.399 | 579.662 | 47.1× | 78.3× |
| HG003 (NA24149) | 4,000,097 | 12.084 | 7.194 | 573.716 | 47.5× | 79.7× |
| HG004 (NA24143) | 4,031,346 | 11.994 | 7.054 | 573.687 | 47.8× | 81.3× |
| HG005 (NA24631) | 3,856,856 | 11.986 | 7.210 | 558.302 | 46.6× | 77.4× |
| HG006 (NA24694) | 3,839,315 | 11.071 | 6.669 | 537.916 | 48.6× | 80.7× |
| HG007 (NA24695) | 3,859,704 | 11.256 | 6.424 | 540.243 | 48.0× | 84.1× |
| **Mean ± SD** | | **11.755 ± 0.462 s** | **7.050 ± 0.371 s** | **559.073 ± 17.016 s** | **47.6× ± 0.6×** | **79.4× ± 3.1×** |

Even on a single core, mafsmith achieved a mean throughput of 334,802 variants/s — a **47.6-fold speedup** over vcf2maf.pl (range 46.6–48.6×), confirming that the performance advantage is primarily algorithmic rather than a product of parallelism. With all 16 cores, throughput increased to 558,914 variants/s (**79.4-fold speedup**, range 74.3–84.1×), a 1.67× parallel scaling factor. Both speedups were highly consistent across samples (CV 1.3% and 3.9%, respectively).

The practical impact of this speedup scales directly with cohort size. On the same instance type (c6a.4xlarge, $0.612/hr on-demand), the reduction in instance time translates to a cost saving of approximately $0.094 per sample and a reduction of 2.3 g CO₂e per sample (location-based, EPA eGRID 2022 Virginia grid mix; see `results/conversion_benchmark_giab_grch38.md` for full methodology). For a cohort of 10,000 samples this represents approximately $938 in compute cost and 23 kg CO₂e avoided. The full end-to-end pipeline speedup (mafsmith + fastVEP vs. vcf2maf.pl + VEP 115) across the same five somatic T/N datasets is reported in Table 4. In a symmetric 16-core comparison (fastVEP 16 Rayon threads vs. VEP `--vep-forks 16`), the full pipeline achieved a mean single-core speedup of **25.3×** (range 17.8–32.4×) and a 16-core speedup of **65.5×** (range 39.6–101.7×). With VEP at its default fork count of 4, speedups increase to 31.5× (1-core) and 83.4× (16-core).

To confirm that speedups generalise to paired tumor/normal somatic VCFs, we benchmarked mafsmith on five additional datasets: MuTect2 and Strelka2 VCFs from the GIAB HG008 somatic benchmark (NYGC pipeline, GRCh38; HG008-T / HG008-N), and MuTect2 and Strelka VCFs from the SEQC2 WGS somatic dataset (HCC1395 / HCC1395BL). All five VCFs carried paired tumor and normal sample columns.

**Table 3. Somatic tumor/normal benchmark: mafsmith vs vcf2maf.pl.**

| Dataset | Caller | Variants | mafsmith 1-core (s) | mafsmith 16-core (s) | vcf2maf.pl (s) | Speedup (1-core) | Speedup (16-core) |
|---------|--------|----------|---------------------|----------------------|----------------|------------------|-------------------|
| GIAB HG008 | MuTect2 | 277,645 | 1.294 | 0.779 | 42.397 | 32.8× | 54.4× |
| GIAB HG008 | Strelka2 SNV | 1,562,847 | 5.273 | 2.960 | 245.899 | 46.6× | 83.1× |
| GIAB HG008 | Strelka2 INDEL | 293,719 | 1.239 | 0.684 | 44.676 | 36.1× | 65.3× |
| SEQC2 HCC1395 | MuTect2 | 271,945 | 1.178 | 0.631 | 38.759 | 32.9× | 61.4× |
| SEQC2 HCC1395 | Strelka | 2,191,720 | 7.301 | 4.147 | 345.556 | 47.3× | 83.3× |
| **Mean** | | | | | | **39.1×** | **69.5×** |

mafsmith achieved a mean single-core speedup of **39.1×** (range 32.8–47.3×) and a 16-core speedup of **69.5×** (range 54.4–83.3×) across paired tumor/normal VCFs. The lower bound of the range reflects MuTect2 VCFs, which carry larger per-variant INFO fields (TLOD, NLOD, per-allele annotations) that increase per-line parsing cost for both tools. Note also that vcf2maf.pl does not accept gzip-compressed input and requires decompression before processing; mafsmith reads gzip natively, a further practical advantage not captured in these timings.

To quantify the full annotation pipeline speedup, we re-ran all five datasets with annotation enabled. We used two vcf2maf.pl configurations: the default `--vep-forks 4` and a symmetric `--vep-forks 16` matching fastVEP's 16-core Rayon parallelism. Table 4 reports the symmetric comparison.

**Table 4. Full annotated pipeline benchmark: mafsmith + fastVEP vs. vcf2maf.pl + VEP 115 (symmetric 16-core comparison).**

| Dataset | Caller | Variants | mafsmith 1-core (s) | mafsmith 16-core (s) | vcf2maf.pl + VEP (--vep-forks 16) (s) | Speedup (1-core) | Speedup (16-core) |
|---------|--------|----------|---------------------|----------------------|----------------------------------------|------------------|-------------------|
| GIAB HG008 | MuTect2 | 277,645 | 25.755 | 11.615 | 459.675 | 17.8× | 39.6× |
| GIAB HG008 | Strelka2 SNV | 1,562,847 | 93.963 | 30.860 | 2851.569 | 30.3× | 92.4× |
| GIAB HG008 | Strelka2 INDEL | 293,719 | 25.720 | 11.577 | 613.111 | 23.8× | 53.0× |
| SEQC2 HCC1395 | MuTect2 | 271,945 | 20.049 | 10.888 | 445.663 | 22.2× | 40.9× |
| SEQC2 HCC1395 | Strelka | 2,191,720 | 128.527 | 40.961 | 4166.796 | 32.4× | 101.7× |
| **Mean** | | | | | | **25.3×** | **65.5×** |

Even in the symmetric 16-core comparison, the full mafsmith + fastVEP pipeline achieved a mean single-core speedup of **25.3×** (range 17.8–32.4×) and a 16-core speedup of **65.5×** (range 39.6–101.7×). With VEP at its default `--vep-forks 4`, speedups increase to a mean of **31.5×** (1-core) and **83.4×** (16-core), reflecting the asymmetry between fastVEP's Rayon parallelism and VEP's fork-based model at low fork counts. The MuTect2 datasets show the lowest speedups due to their larger per-variant INFO fields, consistent with the conversion-only results.

fastVEP and VEP 115 annotate using the same Ensembl GFF3 format but may use different gene model releases depending on the GFF3 version installed. In a comparison of 1,000 variants from the HG008 MuTect2 VCF, we observed transcript-level discordance between fastVEP (using `~/.mafsmith/hg38/genes.gff3`) and VEP 115 for approximately 26% of variants, arising from canonical transcript assignment differences between gene model versions. Consequence-level agreement was higher (87.6% concordance), with most discordances stemming from different gene model versions rather than algorithmic differences in consequence classification. Reconciling gene model versions between fastVEP and VEP 115 is planned for a future release.

---

## Discussion

The primary motivation for mafsmith is throughput: converting thousands of VCFs in large cancer genomics cohorts places substantial demands on compute infrastructure when using vcf2maf.pl, and the per-sample annotation step is a major bottleneck. By combining fastVEP's compiled annotation engine with mafsmith's Rust conversion logic, the full pipeline achieves substantially faster conversion. The MAF conversion logic (field population, allele assignment, consequence mapping) maintains field-for-field agreement with vcf2maf.pl, while annotation-level concordance with VEP 115 depends on the Ensembl gene model version used by fastVEP. The single-core benchmark (47.6× speedup) confirms that the performance advantage is primarily algorithmic — arising from compiled native code, efficient I/O, and zero-copy parsing — with Rayon parallelism providing an additional 1.67× on top on a 16-core instance.

Achieving full concordance with vcf2maf.pl required careful reverse-engineering of a number of non-obvious behaviours accumulated over years of production use. Key examples include: the treatment of absent normal samples vs. no-call normal GTs in homozygous-alt inference; the VAF-based Allele1 override and its interaction with single-sample vs. paired VCF configurations; the handling of truncated AD arrays from GATK multi-allelic calling; consequence severity ranking for multi-consequence VEP annotations; and the representation of structural variant secondary breakpoint rows.

Several known differences remain. In default mode (without `--strict`), mafsmith extracts partial depth counts when AD arrays are truncated, which we regard as more informative than outputting missing values. mafsmith correctly populates SV secondary breakpoint rows with the partner chromosome and position, whereas vcf2maf.pl leaves these fields blank (a bug in the reference implementation). These intentional differences are documented and `--strict` mode is provided for workflows requiring exact vcf2maf.pl compatibility.

mafsmith currently uses a prefix/suffix-stripping approach for allele normalisation rather than FASTA-guided left-alignment. For most variant types this produces identical results to vcf2maf.pl, but complex indels at homopolymer runs or short tandem repeats may differ. Full FASTA-based normalisation is planned for a future release.

---

## Conclusion

mafsmith is a fast, self-contained reimplementation of the vcf2maf.pl conversion pipeline in Rust. It produces field-for-field identical MAF output across thirteen distinct caller types and VCF formats — spanning germline, somatic, structural variant, and annotation-database VCFs — while achieving **47.6-fold** faster conversion on a single core (range 46.6–48.6×) and **79.4-fold** faster on 16 cores (range 74.3–84.1×) across seven GIAB reference samples totalling 27.5 million variants, making large-cohort VCF-to-MAF conversion practical on standard compute infrastructure. mafsmith is open source (Apache 2.0) and available at https://github.com/nf-osi/mafsmith.

---

## Data availability

mafsmith source code and release binaries are available at https://github.com/nf-osi/mafsmith. Validation VCF datasets are available from the sources listed below.

Validation and benchmark VCF datasets used in this work:

- **NF-OSI Synapse project (syn16858331):** Validation VCF files from multiple callers (syn31624545, syn31624535, syn31624939, syn68172710, syn21296193, syn6840402, syn6039268, syn31624525, syn31624637). Available at: https://synapse.org
- **GIAB HG008 somatic (NYGC pipeline, GRCh38):** Paired tumor/normal (HG008-T / HG008-N) MuTect2 and Strelka2 VCFs. Available at: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/NYGC-somatic-pipeline_20240412/GRCh38-GIABv3/
- **SEQC2 WGS somatic (HCC1395 / HCC1395BL, GRCh38):** Paired tumor/normal MuTect2, Strelka2, and SomaticSniper VCFs from the FDA Sequencing Quality Control Phase II (SEQC2) Somatic Mutation Working Group. Available at: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/
- **COLO829 somatic SV truth set (hg38):** Somatic structural variant truth set for the COLO829 melanoma cell line. Lift-over to GRCh38: `truthset_somaticSVs_COLO829_hg38lifted.vcf`. Zenodo: https://zenodo.org/records/7515830
- **GIAB germline benchmarks (HG001–HG007, GRCh38):** NIST GIAB v4.2.1 benchmark VCFs used for conversion-speed benchmarking and validation. Available at: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/
- **COSMIC v103 (GRCh38):** Genome Screens Mutant (Normal) and Non-Coding Variants VCFs used for validation of annotation-database VCF format handling. Available from COSMIC under the COSMIC licence: https://cancer.sanger.ac.uk/cosmic/download
- **ICGC PCAWG cell-line VCFs (GRCh37):** DKFZ SNV/MNV somatic VCFs for HCC1143 and HCC1954 cell lines from the ICGC Pan-Cancer Analysis of Whole Genomes (PCAWG) open data release. Available via the ICGC-ARGO open-access S3 bucket at https://object.genomeinformatics.org (bucket: `icgc25k-open`; no sign request required for open-tier data).

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
