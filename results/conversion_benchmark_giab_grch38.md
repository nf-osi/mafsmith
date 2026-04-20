# Conversion-only benchmark: mafsmith vs vcf2maf.pl

**Date:** 2026-04-20  
**Instance:** AWS c6a.4xlarge (AMD EPYC 7R13, 16 vCPU, 30 GiB RAM, us-east-1)  
**mafsmith:** v0.1.0 (`--skip-annotation`, Rayon parallel, all 16 cores)  
**vcf2maf.pl:** bioconda ensembl-vep/vcf2maf (`--inhibit-vep`, single-threaded)  
**Input:** GIAB NIST v4.2.1 benchmark VCFs, GRCh38 chr1–22, HG001–HG007  
**Source:** https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/  
**Method:** Raw VCF passed directly to each tool (no pre-annotation); one run per file.

## Results

| Sample | Variants | mafsmith (s) | vcf2maf.pl (s) | Speedup |
|--------|----------|-------------|----------------|---------|
| HG001 (NA12878) | 3,893,341 | 7.398 | 549.988 | 74.3× |
| HG002 (NA24385) | 4,048,342 | 7.399 | 579.662 | 78.3× |
| HG003 (NA24149) | 4,000,097 | 7.194 | 573.716 | 79.7× |
| HG004 (NA24143) | 4,031,346 | 7.054 | 573.687 | 81.3× |
| HG005 (NA24631) | 3,856,856 | 7.210 | 558.302 | 77.4× |
| HG006 (NA24694) | 3,839,315 | 6.669 | 537.916 | 80.7× |
| HG007 (NA24695) | 3,859,704 | 6.424 | 540.243 | 84.1× |
| **Mean ± SD** | **3,932,714 ± 93,001** | **7.050 ± 0.371s** | **559.073 ± 17.016s** | **79.4× ± 3.1×** |

**Total variants processed:** 27,529,001  
**Throughput — mafsmith:** 558,914 variants/s (mean)  
**Throughput — vcf2maf.pl:** 7,036 variants/s (mean)

## Cost and carbon savings

### Methodology

Costs use AWS on-demand pricing for c6a.4xlarge in us-east-1 ($0.612/hr).

Power is estimated using the [Cloud Carbon Footprint (CCF) methodology](https://www.cloudcarbonfootprint.org/docs/methodology/):

```
instance_watts = (cpu_watts + memory_watts) × PUE
cpu_watts      = min_cpu_w + utilisation × (max_cpu_w − min_cpu_w)
min_cpu_w      = n_vcpus × 1.04   # AMD EPYC Milan (3rd Gen) idle coefficient
max_cpu_w      = n_vcpus × 5.12   # AMD EPYC Milan peak coefficient
memory_watts   = mem_gb × 0.392   # CCF DRAM coefficient (W/GB)
PUE            = 1.2              # AWS data-centre PUE
```

CPU utilisation:
- **mafsmith**: 100% (Rayon uses all 16 cores)
- **vcf2maf.pl**: 6.25% (single-threaded; 1/16 cores)

Carbon intensity: **0.386 kg CO₂e/kWh** — EPA eGRID 2022, SRVC subregion (Virginia),
location-based. AWS market-based rate is lower due to renewable energy certificates
but the location-based figure reflects actual grid emissions.

### Reproducing these calculations

```python
# results/compute_savings.py
import statistics

# --- Raw benchmark data (wall-clock seconds, one run per GIAB sample) ---
data = [
    # (sample,    n_variants, mafsmith_s, vcf2maf_s)
    ("HG001",  3_893_341,  7.398,  549.988),
    ("HG002",  4_048_342,  7.399,  579.662),
    ("HG003",  4_000_097,  7.194,  573.716),
    ("HG004",  4_031_346,  7.054,  573.687),
    ("HG005",  3_856_856,  7.210,  558.302),
    ("HG006",  3_839_315,  6.669,  537.916),
    ("HG007",  3_859_704,  6.424,  540.243),
]

# --- Instance parameters ---
VCPUS            = 16
HOURLY_RATE_USD  = 0.612    # c6a.4xlarge on-demand, us-east-1 (2026-04)
PUE              = 1.2      # AWS data-centre power usage effectiveness
MEM_GB           = 32
MEM_W_PER_GB     = 0.392    # CCF DRAM coefficient (W/GB)

# CCF AMD EPYC 3rd Gen (Milan) CPU coefficients (W/vCPU)
MIN_W_PER_VCPU   = 1.04
MAX_W_PER_VCPU   = 5.12

# Carbon intensity: EPA eGRID 2022 SRVC (Virginia), location-based
CARBON_KG_PER_KWH = 0.386

# --- Power model ---
min_cpu_w = VCPUS * MIN_W_PER_VCPU   # 16.64 W at idle
max_cpu_w = VCPUS * MAX_W_PER_VCPU   # 81.92 W at peak
mem_w     = MEM_GB * MEM_W_PER_GB    # 12.54 W

def instance_watts(cpu_util: float) -> float:
    cpu_w = min_cpu_w + cpu_util * (max_cpu_w - min_cpu_w)
    return (cpu_w + mem_w) * PUE

# mafsmith: all cores active; vcf2maf.pl: single-threaded (1/VCPUS)
pwr_ms = instance_watts(1.0)
pwr_vc = instance_watts(1.0 / VCPUS)

# --- Per-sample calculations ---
ms_times = [d[2] for d in data]
vc_times = [d[3] for d in data]
ms_mean  = statistics.mean(ms_times)
vc_mean  = statistics.mean(vc_times)

cost_ms      = (ms_mean / 3600) * HOURLY_RATE_USD
cost_vc      = (vc_mean / 3600) * HOURLY_RATE_USD
cost_saved   = cost_vc - cost_ms

energy_ms    = pwr_ms * (ms_mean / 3600) / 1000   # kWh
energy_vc    = pwr_vc * (vc_mean / 3600) / 1000   # kWh
energy_saved = energy_vc - energy_ms

carbon_saved_g = energy_saved * CARBON_KG_PER_KWH * 1000   # grams

print(f"Mean wall-clock — mafsmith: {ms_mean:.3f}s  vcf2maf.pl: {vc_mean:.3f}s")
print(f"Mean speedup: {statistics.mean(d[3]/d[2] for d in data):.1f}x")
print(f"Cost saved per sample:   ${cost_saved:.5f}")
print(f"Energy saved per sample: {energy_saved*1000:.4f} Wh")
print(f"CO2e saved per sample:   {carbon_saved_g:.3f} g")

print("\nAt scale:")
print(f"{'Samples':>10}  {'Cost (USD)':>12}  {'CO2e (kg)':>12}  {'Time (hrs)':>12}")
for n in [100, 1_000, 10_000, 100_000, 1_000_000]:
    print(f"{n:>10,}  ${n*cost_saved:>11.2f}  "
          f"{n*carbon_saved_g/1000:>12.3f}  "
          f"{n*(vc_mean-ms_mean)/3600:>12.1f}")
```

### Summary (c6a.4xlarge on-demand, mean across HG001–HG007)

| Metric | mafsmith | vcf2maf.pl | Saved per sample |
|--------|----------|------------|-----------------|
| Wall-clock time | 7.050s | 559.073s | 552.0s (9.2 min) |
| Instance cost (USD) | $0.00120 | $0.09504 | **$0.09384** |
| Est. power draw | 113.4 W | 39.9 W | — |
| Energy | 0.000222 kWh | 0.006199 kWh | **5.977 Wh** |
| CO₂e (location-based) | 0.086 g | 2.393 g | **2.307 g** |

### At scale

| Samples | Cost saved (USD) | CO₂e saved (kg) | Time saved (hrs) |
|---------|-----------------|-----------------|-----------------|
| 100 | $9.38 | 0.23 | 15.3 |
| 1,000 | $93.84 | 2.31 | 153.3 |
| 10,000 | $938.44 | 23.1 | 1,533 |
| 100,000 | $9,384 | 231 | 15,334 |
| 1,000,000 | $93,844 | 2,307 | 153,340 |
