# mafsmith benchmark results

**Date:** 2026-04-21  
**Instance:** ip-10-41-27-195.ec2.internal  
**mafsmith:** `/home/ssm-user/mafsmith/target/release/mafsmith`  
**vcf2maf.pl:** `/home/ssm-user/miniconda3/envs/vcf2maf-env/bin/vcf2maf.pl`  
**Reference:** `/home/ssm-user/.mafsmith/hg38/reference.fa`  


## Somatic T/N — GIAB HG008

| Dataset | Variants | ms 1-core (s) | ms 16-core (s) | vcf2maf.pl (s) | Speedup (1-core) | Speedup (16-core) |
|---------|----------|---------------|-----------------|----------------|-----------------|--------------|
| HG008 MuTect2 | 277,645 | 25.755 | 11.615 | 459.675 | 17.8× | 39.6× |
| HG008 Strelka2 SNV | 1,562,847 | 93.963 | 30.860 | 2851.569 | 30.3× | 92.4× |
| HG008 Strelka2 INDEL | 293,719 | 25.720 | 11.577 | 613.111 | 23.8× | 53.0× |
| **Mean** | | | | | **24.0×** | **61.6×** |

## Somatic T/N — SEQC2 HCC1395

| Dataset | Variants | ms 1-core (s) | ms 16-core (s) | vcf2maf.pl (s) | Speedup (1-core) | Speedup (16-core) |
|---------|----------|---------------|-----------------|----------------|-----------------|--------------|
| SEQC2 MuTect2 | 271,945 | 20.049 | 10.888 | 445.663 | 22.2× | 40.9× |
| SEQC2 Strelka | 2,191,720 | 128.527 | 40.961 | 4166.796 | 32.4× | 101.7× |
| **Mean** | | | | | **27.3×** | **71.3×** |
