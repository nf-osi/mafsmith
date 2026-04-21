# mafsmith benchmark results

**Date:** 2026-04-21  
**Instance:** ip-10-41-27-195.ec2.internal  
**mafsmith:** `/home/ssm-user/mafsmith/target/release/mafsmith`  
**vcf2maf.pl:** `/home/ssm-user/miniconda3/envs/vcf2maf-env/bin/vcf2maf.pl`  
**Reference:** `/home/ssm-user/.mafsmith/hg38/reference.fa`  


## Somatic T/N — GIAB HG008

| Dataset | Variants | ms 1-core (s) | ms 16-core (s) | vcf2maf.pl (s) | Speedup (1-core) | Speedup (16-core) |
|---------|----------|---------------|-----------------|----------------|-----------------|--------------|
| HG008 MuTect2 | 277,645 | 21.870 | 11.013 | 582.531 | 26.6× | 52.9× |
| HG008 Strelka2 SNV | 1,562,847 | 94.123 | 31.123 | 3591.814 | 38.2× | 115.4× |
| HG008 Strelka2 INDEL | 293,719 | 32.410 | 11.371 | 799.689 | 24.7× | 70.3× |
| **Mean** | | | | | **29.8×** | **79.5×** |

## Somatic T/N — SEQC2 HCC1395

| Dataset | Variants | ms 1-core (s) | ms 16-core (s) | vcf2maf.pl (s) | Speedup (1-core) | Speedup (16-core) |
|---------|----------|---------------|-----------------|----------------|-----------------|--------------|
| SEQC2 MuTect2 | 271,945 | 20.383 | 10.474 | 566.607 | 27.8× | 54.1× |
| SEQC2 Strelka | 2,191,720 | 128.165 | 41.725 | 5178.114 | 40.4× | 124.1× |
| **Mean** | | | | | **34.1×** | **89.1×** |
