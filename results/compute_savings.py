#!/usr/bin/env python3
"""
Reproduce cost and carbon savings estimates from the GIAB conversion benchmark.

Data: mafsmith v0.1.0 vs vcf2maf.pl (--inhibit-vep) on GIAB NIST v4.2.1 GRCh38,
      HG001-HG007, AWS c6a.4xlarge (us-east-1), 2026-04-20.

Power model: Cloud Carbon Footprint methodology
  https://www.cloudcarbonfootprint.org/docs/methodology/
Carbon intensity: EPA eGRID 2022 SRVC (Virginia), location-based
Instance price: AWS on-demand c6a.4xlarge us-east-1, April 2026
"""

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
min_cpu_w = VCPUS * MIN_W_PER_VCPU
max_cpu_w = VCPUS * MAX_W_PER_VCPU
mem_w     = MEM_GB * MEM_W_PER_GB

def instance_watts(cpu_util: float) -> float:
    cpu_w = min_cpu_w + cpu_util * (max_cpu_w - min_cpu_w)
    return (cpu_w + mem_w) * PUE

# mafsmith: all cores active; vcf2maf.pl: single-threaded (1/VCPUS)
pwr_ms = instance_watts(1.0)
pwr_vc = instance_watts(1.0 / VCPUS)

# --- Per-sample statistics ---
ms_times  = [d[2] for d in data]
vc_times  = [d[3] for d in data]
speedups  = [d[3] / d[2] for d in data]
ms_mean   = statistics.mean(ms_times)
vc_mean   = statistics.mean(vc_times)

cost_ms      = (ms_mean / 3600) * HOURLY_RATE_USD
cost_vc      = (vc_mean / 3600) * HOURLY_RATE_USD
cost_saved   = cost_vc - cost_ms

energy_ms    = pwr_ms * (ms_mean / 3600) / 1000
energy_vc    = pwr_vc * (vc_mean / 3600) / 1000
energy_saved = energy_vc - energy_ms

carbon_ms_g    = energy_ms    * CARBON_KG_PER_KWH * 1000
carbon_vc_g    = energy_vc    * CARBON_KG_PER_KWH * 1000
carbon_saved_g = energy_saved * CARBON_KG_PER_KWH * 1000

ms_vps = statistics.mean(d[1] / d[2] for d in data)
vc_vps = statistics.mean(d[1] / d[3] for d in data)

print("=" * 62)
print(f" Instance:  c6a.4xlarge  {VCPUS} vCPU, {MEM_GB} GiB, AMD EPYC 7R13")
print(f" Rate:      ${HOURLY_RATE_USD}/hr on-demand (us-east-1)")
print("=" * 62)

print(f"\n{'':32s} {'mafsmith':>12} {'vcf2maf.pl':>12}")
print(f" {'Mean wall-clock (s)':32s} {ms_mean:>12.3f} {vc_mean:>12.3f}")
print(f" {'SD (s)':32s} {statistics.stdev(ms_times):>12.3f} {statistics.stdev(vc_times):>12.3f}")
print(f" {'CPU utilisation':32s} {'100%':>12} {f'1/{VCPUS} cores':>12}")
print(f" {'Power draw (W)':32s} {pwr_ms:>12.1f} {pwr_vc:>12.1f}")
print(f" {'Energy (kWh)':32s} {energy_ms:>12.6f} {energy_vc:>12.6f}")
print(f" {'CO2e (g)':32s} {carbon_ms_g:>12.3f} {carbon_vc_g:>12.3f}")
print(f" {'Cost (USD)':32s} ${cost_ms:>11.5f} ${cost_vc:>11.5f}")
print(f" {'Throughput (variants/s)':32s} {ms_vps:>12,.0f} {vc_vps:>12,.0f}")

print(f"\n{'Speedup':32s} {statistics.mean(speedups):.1f}x "
      f"(SD {statistics.stdev(speedups):.1f}x, "
      f"range {min(speedups):.1f}–{max(speedups):.1f}x)")

print("\n" + "=" * 62)
print(" SAVINGS PER SAMPLE (mean)")
print("=" * 62)
print(f"  Cost saved:    ${cost_saved:.5f} USD")
print(f"  Energy saved:  {energy_saved * 1000:.4f} Wh")
print(f"  CO2e saved:    {carbon_saved_g:.3f} g")
print(f"  Time saved:    {vc_mean - ms_mean:.1f} s ({(vc_mean - ms_mean) / 60:.2f} min)")

print("\n" + "=" * 62)
print(" AT SCALE")
print("=" * 62)
print(f"  {'Samples':>10}  {'Cost (USD)':>12}  {'CO2e (kg)':>12}  {'Time (hrs)':>12}")
for n in [100, 1_000, 10_000, 100_000, 1_000_000]:
    print(f"  {n:>10,}  ${n * cost_saved:>11.2f}  "
          f"{n * carbon_saved_g / 1000:>12.3f}  "
          f"{n * (vc_mean - ms_mean) / 3600:>12.1f}")

print(f"""
Notes:
  Power model:      CCF methodology (AMD EPYC Milan coefficients)
                    min {MIN_W_PER_VCPU} W/vCPU, max {MAX_W_PER_VCPU} W/vCPU, PUE {PUE}
                    mafsmith {pwr_ms:.1f}W (100% util), vcf2maf.pl {pwr_vc:.1f}W (1/{VCPUS} util)
  Carbon intensity: {CARBON_KG_PER_KWH} kg CO2e/kWh (EPA eGRID 2022, location-based)
                    AWS market-based rate lower due to renewable energy credits
  Spot pricing:     ~$0.18/hr scales cost figures by ~0.29x
""")
