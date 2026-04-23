/// Allele depth information extracted from a VCF FORMAT/sample column.
#[derive(Debug, Default, Clone)]
pub struct AlleleDepth {
    pub ref_count: Option<u32>,
    pub alt_count: Option<u32>,
    pub total_depth: Option<u32>,
}

impl AlleleDepth {
    pub fn depth(&self) -> Option<u32> {
        self.total_depth
            .or_else(|| self.ref_count.zip(self.alt_count).map(|(r, a)| r + a))
    }
}

/// Extract allele depths from a VCF sample column given FORMAT keys.
///
/// Handles the following caller conventions:
/// - Standard (GATK/FreeBayes): AD=ref,alt
/// - VarScan 2: RD=ref, AD=alt (single value, not comma list)
/// - Strelka SNVs: AU,CU,GU,TU counts per base
/// - Strelka indels: TAR=ref_tier1,ref_tier2; TIR=alt_tier1,alt_tier2
/// - SomaticSniper: DP4 for depth + BCOUNT=A,C,G,T for per-allele ref/alt counts
/// - Ion Torrent: RO=ref, AO=alt
/// - Total depth fallback: DP
///
/// `alt_vcf_idx`: 1-based VCF allele index of the selected ALT (1 = first ALT, 2 = second, …).
/// This is used as the AD array index for the alt count in multi-allelic records.
///
/// `F` is generic so callers can pass `&[String]` (from VcfRecord) without first collecting
/// to `Vec<&str>`, avoiding a per-record allocation.
pub fn extract_depth<F: AsRef<str>>(
    format_keys: &[F],
    sample_values: &[&str],
    ref_allele: &str,
    alt_allele: &str,
    alt_vcf_idx: usize,
    strict: bool,
) -> AlleleDepth {
    let get = |key: &str| -> Option<&str> {
        format_keys
            .iter()
            .position(|k| k.as_ref() == key)
            .and_then(|i| sample_values.get(i).copied())
            .filter(|v| !v.is_empty() && *v != ".")
    };

    let first_num = |s: &str| s.split(',').next().and_then(|v| v.parse::<u32>().ok());

    // DP (total depth) — used as fallback
    let dp = get("DP").and_then(|v| v.parse::<u32>().ok());

    // 1. Standard AD field: ref,alt[,alt2,...]
    //    Values may be "." (unknown) — handle per-field, not by filtering the whole split.
    if let Some(ad_str) = get("AD") {
        let raw: Vec<&str> = ad_str.split(',').collect();
        if raw.len() >= 2 {
            let ref_count: Option<u32> = raw[0].parse().ok();
            // Use alt_vcf_idx to pick the correct AD field for multi-allelic records.
            let alt_count: Option<u32> = raw.get(alt_vcf_idx).and_then(|v| v.parse().ok());
            // Accept even if one side is unknown ("."), as long as we parsed something
            if ref_count.is_some() || alt_count.is_some() {
                // vcf2maf.pl uses max(DP, sum(AD)) for total depth: when reads that
                // pass allele-depth filters outnumber DP (e.g. GIAB ADALL-rich sites),
                // sum(AD) wins; when many reads are not allele-counted (e.g. ambiguous/
                // supplementary), DP wins. Mirrors the same max() used for Strelka TAR/TIR.
                let ad_total: u32 = raw.iter().filter_map(|v| v.parse::<u32>().ok()).sum();
                let total = match (dp, ad_total) {
                    (Some(d), a) if a > 0 => Some(d.max(a)),
                    (Some(d), _) => Some(d),
                    (None, a) if a > 0 => Some(a),
                    (None, _) => None,
                };
                return AlleleDepth {
                    ref_count,
                    alt_count,
                    total_depth: total,
                };
            }
        }
    }

    // 2. VarScan: RD (ref depth) + AD (single alt depth)
    if let (Some(rd), Some(ad)) = (get("RD"), get("AD")) {
        if let (Some(r), Some(a)) = (rd.parse::<u32>().ok(), ad.parse::<u32>().ok()) {
            return AlleleDepth {
                ref_count: Some(r),
                alt_count: Some(a),
                total_depth: dp.or(Some(r + a)),
            };
        }
    }

    // 3. Strelka SNVs: AU, CU, GU, TU (first value = tier1)
    let base_keys = ["AU", "CU", "GU", "TU"];
    if base_keys.iter().any(|k| get(k).is_some()) {
        let ref_base = ref_allele.chars().next().unwrap_or('N').to_ascii_uppercase();
        let alt_base = alt_allele.chars().next().unwrap_or('N').to_ascii_uppercase();
        let base_key = |b: char| match b {
            'A' => "AU",
            'C' => "CU",
            'G' => "GU",
            'T' => "TU",
            _ => "AU",
        };
        let ref_count = get(base_key(ref_base)).and_then(first_num);
        let alt_count = get(base_key(alt_base)).and_then(first_num);
        return AlleleDepth {
            ref_count,
            alt_count,
            total_depth: dp,
        };
    }

    // 4. Strelka indels: TAR (ref), TIR (alt)
    // vcf2maf.pl uses max(DP, TAR[0]+TIR[0]) for total depth, because TAR/TIR count
    // reads at a specific indel tier that can exceed the DP-field count.
    if get("TAR").is_some() || get("TIR").is_some() {
        let ref_count = get("TAR").and_then(first_num);
        let alt_count = get("TIR").and_then(first_num);
        let sum_tar_tir = ref_count.zip(alt_count).map(|(r, a)| r + a);
        let total = match (dp, sum_tar_tir) {
            (Some(d), Some(s)) => Some(d.max(s)),
            (Some(d), None)    => Some(d),
            (None, Some(s))    => Some(s),
            (None, None)       => None,
        };
        return AlleleDepth {
            ref_count,
            alt_count,
            total_depth: total,
        };
    }

    // 5. SomaticSniper: has both DP4 and BCOUNT.
    // When BCOUNT is present alongside DP4, use BCOUNT for per-allele ref/alt counts
    // (matches vcf2maf.pl behavior). DP4 alt_fwd+alt_rev counts ALL non-ref reads, not
    // just the specific alt allele — wrong for multi-allelic sites and single-allelic
    // sites where multiple alt bases are observed.
    if get("DP4").is_some() {
        if let Some(bc) = get("BCOUNT") {
            let bc_parts: Vec<u32> = bc.split(',').filter_map(|v| v.parse().ok()).collect();
            if bc_parts.len() == 4 {
                let idx_for = |b: char| match b {
                    'A' => 0usize, 'C' => 1, 'G' => 2, 'T' => 3, _ => 0,
                };
                let ref_char = ref_allele.chars().next().unwrap_or('N');
                let alt_base = alt_allele.chars().next().unwrap_or('N').to_ascii_uppercase();
                // In strict mode, vcf2maf.pl does a case-sensitive hash lookup over {A,C,G,T}.
                // Lowercase ref alleles (e.g. on decoy contigs) fail that lookup → ref_count='.'.
                let ref_count = if strict && !matches!(ref_char, 'A' | 'C' | 'G' | 'T') {
                    None
                } else {
                    bc_parts.get(idx_for(ref_char.to_ascii_uppercase())).copied()
                };
                let alt_count = bc_parts.get(idx_for(alt_base)).copied();
                let dp4_total = get("DP4").and_then(|dp4| {
                    let p: Vec<u32> = dp4.split(',').filter_map(|v| v.parse().ok()).collect();
                    if p.len() == 4 { Some(p.iter().sum()) } else { None }
                });
                return AlleleDepth {
                    ref_count,
                    alt_count,
                    total_depth: dp.or(dp4_total).or(Some(bc_parts.iter().sum())),
                };
            }
        }
        // DP4 only (no BCOUNT): ref_fwd+ref_rev, alt_fwd+alt_rev
        if let Some(dp4) = get("DP4") {
            let parts: Vec<u32> = dp4.split(',').filter_map(|v| v.parse().ok()).collect();
            if parts.len() == 4 {
                let ref_count = parts[0] + parts[1];
                let alt_count = parts[2] + parts[3];
                return AlleleDepth {
                    ref_count: Some(ref_count),
                    alt_count: Some(alt_count),
                    total_depth: dp.or(Some(ref_count + alt_count)),
                };
            }
        }
    }

    // 6. Ion Torrent: RO (ref obs), AO (alt obs).
    // AO is alt-only (no REF entry), so the correct element for a given VCF allele
    // index is AO[alt_vcf_idx - 1].  Older code always took AO[0] which was wrong
    // for multi-allelic sites where GT = 0/2, 0/3, etc.
    if let (Some(ro), Some(ao)) = (get("RO"), get("AO")) {
        let ref_count = ro.parse::<u32>().ok();
        let alt_count = ao.split(',')
            .nth(alt_vcf_idx.saturating_sub(1))
            .and_then(|v| v.parse().ok());
        return AlleleDepth {
            ref_count,
            alt_count,
            total_depth: dp,
        };
    }

    // 7. BCOUNT only (no DP4): allele counts as A,C,G,T
    if let Some(bc) = get("BCOUNT") {
        let parts: Vec<u32> = bc.split(',').filter_map(|v| v.parse().ok()).collect();
        if parts.len() == 4 {
            let idx_for = |b: char| match b {
                'A' => 0usize,
                'C' => 1,
                'G' => 2,
                'T' => 3,
                _ => 0,
            };
            let ref_char = ref_allele.chars().next().unwrap_or('N');
            let alt_base = alt_allele.chars().next().unwrap_or('N').to_ascii_uppercase();
            let ref_count = if strict && !matches!(ref_char, 'A' | 'C' | 'G' | 'T') {
                None
            } else {
                parts.get(idx_for(ref_char.to_ascii_uppercase())).copied()
            };
            let alt_count = parts.get(idx_for(alt_base)).copied();
            return AlleleDepth {
                ref_count,
                alt_count,
                total_depth: dp.or(Some(parts.iter().sum())),
            };
        }
    }

    // Depth-only fallback
    AlleleDepth {
        ref_count: None,
        alt_count: None,
        total_depth: dp,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_ad() {
        let keys = vec!["GT", "AD", "DP"];
        let vals = vec!["0/1", "45,12", "57"];
        let d = extract_depth(&keys, &vals, "A", "T", 1, false);
        assert_eq!(d.ref_count, Some(45));
        assert_eq!(d.alt_count, Some(12));
        assert_eq!(d.total_depth, Some(57));
    }

    #[test]
    fn ad_max_dp_adsum() {
        // GIAB NIST v4.2.1: sum(AD)=174 > DP=157 → max wins → 174.
        let keys = vec!["GT", "PS", "DP", "ADALL", "AD", "GQ"];
        let vals = vec!["1/1", ".", "157", "0,70", "87,87", "237"];
        let d = extract_depth(&keys, &vals, "T", "C", 1, false);
        assert_eq!(d.ref_count, Some(87));
        assert_eq!(d.alt_count, Some(87));
        assert_eq!(d.total_depth, Some(174)); // max(157, 87+87) = 174

        // DP=652 > sum(AD)=82 → max wins → 652.
        let keys2 = vec!["GT", "PS", "DP", "ADALL", "AD", "GQ"];
        let vals2 = vec!["1/1", ".", "652", "16,234", "0,82", "312"];
        let d2 = extract_depth(&keys2, &vals2, "A", "G", 1, false);
        assert_eq!(d2.ref_count, Some(0));
        assert_eq!(d2.alt_count, Some(82));
        assert_eq!(d2.total_depth, Some(652)); // max(652, 0+82) = 652
    }

    #[test]
    fn multiallelic_ad_second_alt() {
        let keys = vec!["GT", "AD", "DP"];
        let vals = vec!["0/2", "10,5,8", "23"];
        let d = extract_depth(&keys, &vals, "A", "G", 2, false);
        assert_eq!(d.ref_count, Some(10));
        assert_eq!(d.alt_count, Some(8));
        assert_eq!(d.total_depth, Some(23));
    }

    #[test]
    fn varscan() {
        let keys = vec!["GT", "RD", "AD", "DP"];
        let vals = vec!["0/1", "30", "10", "40"];
        let d = extract_depth(&keys, &vals, "A", "T", 1, false);
        assert_eq!(d.ref_count, Some(30));
        assert_eq!(d.alt_count, Some(10));
    }

    #[test]
    fn dp4_only() {
        let keys = vec!["GT", "DP4"];
        let vals = vec!["0/1", "10,20,5,3"];
        let d = extract_depth(&keys, &vals, "A", "T", 1, false);
        assert_eq!(d.ref_count, Some(30));
        assert_eq!(d.alt_count, Some(8));
    }

    #[test]
    fn somaticsniper_dp4_bcount() {
        // SomaticSniper T>C record: DP4 alt total=8 (5 fwd+3 rev), but BCOUNT[C]=7 (one
        // read supports G not C). vcf2maf.pl uses BCOUNT, so mafsmith must too.
        let keys = vec!["GT", "DP", "DP4", "BCOUNT"];
        let vals = vec!["0/1", "44", "17,19,5,3", "0,7,1,36"];
        let d = extract_depth(&keys, &vals, "T", "C", 1, false);
        assert_eq!(d.ref_count, Some(36)); // BCOUNT[T]
        assert_eq!(d.alt_count, Some(7));  // BCOUNT[C], not DP4 alt(8)
        assert_eq!(d.total_depth, Some(44)); // DP wins over DP4 sum
    }

    #[test]
    fn bcount_lowercase_ref_strict() {
        // SomaticSniper on decoy contig: ref allele is lowercase 't', alt uppercase 'C'.
        // vcf2maf.pl does a case-sensitive BCOUNT hash lookup → ref_count='.'.
        // In strict mode mafsmith must match; non-strict mode should still extract correctly.
        let keys = vec!["GT", "DP", "DP4", "BCOUNT"];
        let vals = vec!["0/1", "173", "69,88,4,12", "0,16,0,157"];
        let d_strict = extract_depth(&keys, &vals, "t", "C", 1, true);
        assert_eq!(d_strict.ref_count, None);   // lowercase ref → None in strict mode
        assert_eq!(d_strict.alt_count, Some(16)); // C uppercase → correct
        let d_normal = extract_depth(&keys, &vals, "t", "C", 1, false);
        assert_eq!(d_normal.ref_count, Some(157)); // converted to uppercase, finds T
        assert_eq!(d_normal.alt_count, Some(16));
    }
}
