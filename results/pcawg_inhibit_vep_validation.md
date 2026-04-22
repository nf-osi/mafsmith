# PCAWG consensus SNV/MNV validation: mafsmith vs vcf2maf.pl (--inhibit-vep)

**Date:** 2026-04-22
**Instance:** AWS (16 vCPU, 30 GiB RAM); peak RAM used 9.3 GB (30%), no swap
**mafsmith:** `--skip-annotation --strict` (conversion only, no VEP)
**vcf2maf.pl:** `--inhibit-vep --strict` (conversion only, no VEP)
**Reference:** Ensembl GRCh37 release-112 (`Homo_sapiens.GRCh37.dna.toplevel.fa.gz`, bgzip+faidx)
**Comparison fields:** End_Position, Matched_Norm_Sample_Barcode, Reference_Allele, Start_Position, Tumor_Sample_Barcode, Tumor_Seq_Allele1, Tumor_Seq_Allele2, Variant_Type, n_alt_count, n_depth, n_ref_count, t_alt_count, t_depth, t_ref_count (`Variant_Classification` excluded — requires annotation)

## Dataset

ICGC PCAWG consensus SNV/MNV callset (open-access), downloaded from `s3://icgc25k-open/PCAWG/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public.tgz` via `--endpoint-url https://object.genomeinformatics.org`.

| Cohort | Samples | Total variants | Mean variants/sample |
|--------|---------|---------------|----------------------|
| ICGC PCAWG consensus SNV/MNV | 1,902 | 21,628,933 | 11,371 |

## Results

| Metric | Value |
|--------|-------|
| Samples processed | 1,902 / 1,902 |
| Samples with ≥1 conversion diff | **0** |
| Total conversion diffs | **0** |
| Peak memory used | 9.3 GB / 30.6 GB |
| Swap used | 0 MB |

**mafsmith and vcf2maf.pl produce byte-identical conversion output across all 21.6 M PCAWG variants** (allele normalization, coordinate conversion, depth extraction, Variant_Type) in `--inhibit-vep` mode.

## Script

`scripts/batch_vcf_diff_pcawg_inhibit_vep.sh` — runs `vcf_diff.py --inhibit-vep --strict` across all 1,902 VCFs at configurable parallelism (default 4 workers). Per-sample outputs written to `bench_vcfs/PCAWG_consensus/diffs_inhibit_vep/` (not tracked in git).
