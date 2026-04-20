# Somatic tumor/normal benchmark: mafsmith vs vcf2maf.pl

**Date:** 2026-04-20  
**Instance:** AWS c6a.4xlarge (AMD EPYC 7R13, 16 vCPU, 30 GiB RAM, us-east-1)  
**mafsmith:** v0.1.0 (`--skip-annotation`); timed at 1 core (`RAYON_NUM_THREADS=1`, 3 runs, mean reported) and 16 cores (default Rayon, 3 runs, mean reported)  
**vcf2maf.pl:** bioconda ensembl-vep/vcf2maf (`--inhibit-vep`, single-threaded, 3 runs, mean reported)  
**Reference:** Broad GATK GRCh38 bundle (`Homo_sapiens_assembly38.fasta`, chr-prefixed)  
**Method:** Gzipped VCFs decompressed to plain text before vcf2maf.pl (which does not accept gzip input); mafsmith reads gzip natively. Conversion only — no VEP annotation.

## Datasets

| Dataset | Specimen (T/N) | Caller | Variants | Source |
|---------|---------------|--------|----------|--------|
| GIAB HG008 | HG008-T / HG008-N | MuTect2 | 277,645 | NYGC somatic pipeline, GRCh38-GIABv3 |
| GIAB HG008 | HG008-T / HG008-N | Strelka2 SNV | 1,562,847 | NYGC somatic pipeline, GRCh38-GIABv3 |
| GIAB HG008 | HG008-T / HG008-N | Strelka2 INDEL | 293,719 | NYGC somatic pipeline, GRCh38-GIABv3 |
| SEQC2 | HCC1395 / HCC1395BL | MuTect2 | 271,945 | SEQC2 WGS, bwa, replicate 1 |
| SEQC2 | HCC1395 / HCC1395BL | Strelka | 2,191,720 | SEQC2 WGS, bwa, replicate 1 |

## Results

| Dataset | Caller | Variants | mafsmith 1-core (s) | mafsmith 16-core (s) | vcf2maf.pl (s) | Speedup (1-core) | Speedup (16-core) |
|---------|--------|----------|---------------------|----------------------|----------------|------------------|-------------------|
| GIAB HG008 | MuTect2 | 277,645 | 1.294 | 0.779 | 42.397 | 32.8× | 54.4× |
| GIAB HG008 | Strelka2 SNV | 1,562,847 | 5.273 | 2.960 | 245.899 | 46.6× | 83.1× |
| GIAB HG008 | Strelka2 INDEL | 293,719 | 1.239 | 0.684 | 44.676 | 36.1× | 65.3× |
| SEQC2 | MuTect2 | 271,945 | 1.178 | 0.631 | 38.759 | 32.9× | 61.4× |
| SEQC2 | Strelka | 2,191,720 | 7.301 | 4.147 | 345.556 | 47.3× | 83.3× |
| **Mean** | | | | | | **39.1×** | **69.5×** |

**Throughput — mafsmith (1-core):** 255,807 variants/s (mean)  
**Throughput — mafsmith (16-core):** 454,661 variants/s (mean)  
**Throughput — vcf2maf.pl:** 6,568 variants/s (mean)

## Notes

- The MuTect2 speedup (32.8–32.9×) is lower than Strelka2/Strelka (36.1–47.3×) because MuTect2 VCFs carry substantially larger INFO fields (TLOD, NLOD, per-allele annotations), making per-line parsing more expensive for both tools — but Perl pays a proportionally higher cost.
- vcf2maf.pl does not accept gzip-compressed VCF input; files were decompressed to plain text before timing. mafsmith reads gzip natively, so no decompression step is needed in practice — a further practical advantage not reflected in these timings.
- All timings are means of 3 runs (mafsmith and vcf2maf.pl).
