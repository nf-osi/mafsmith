# mafsmith benchmark results

**Date:** 2026-04-20  
**Instance:** ip-10-41-27-195.ec2.internal  
**mafsmith:** `/home/ssm-user/mafsmith/target/release/mafsmith`  
**vcf2maf.pl:** `/home/ssm-user/miniconda3/envs/vcf2maf-env/bin/vcf2maf.pl`  
**Reference:** `/home/ssm-user/.mafsmith/hg38/reference.fa`  


## Somatic T/N — GIAB HG008

| Dataset | Variants | ms 1-core (s) | ms 16-core (s) | vcf2maf.pl (s) | Speedup (1-core) | Speedup (16-core) |
|---------|----------|---------------|-----------------|----------------|-----------------|--------------|
| HG008 MuTect2 | 277,645 | 1.294 | 0.779 | 42.397 | 32.8× | 54.4× |
| HG008 Strelka2 SNV | 1,562,847 | 5.273 | 2.960 | 245.899 | 46.6× | 83.1× |
| HG008 Strelka2 INDEL | 293,719 | 1.239 | 0.684 | 44.676 | 36.1× | 65.3× |
| **Mean** | | | | | **38.5×** | **67.6×** |

## Somatic T/N — SEQC2 HCC1395

| Dataset | Variants | ms 1-core (s) | ms 16-core (s) | vcf2maf.pl (s) | Speedup (1-core) | Speedup (16-core) |
|---------|----------|---------------|-----------------|----------------|-----------------|--------------|
| SEQC2 MuTect2 | 271,945 | 1.178 | 0.631 | 38.759 | 32.9× | 61.4× |
| SEQC2 Strelka | 2,191,720 | 7.301 | 4.147 | 345.556 | 47.3× | 83.3× |
| **Mean** | | | | | **40.1×** | **72.3×** |
