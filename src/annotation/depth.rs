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
/// - SomaticSniper: DP4=ref_fwd,ref_rev,alt_fwd,alt_rev
/// - Ion Torrent: RO=ref, AO=alt
/// - Total depth fallback: DP
///
/// `alt_vcf_idx`: 1-based VCF allele index of the selected ALT (1 = first ALT, 2 = second, …).
/// This is used as the AD array index for the alt count in multi-allelic records.
pub fn extract_depth(
    format_keys: &[&str],
    sample_values: &[&str],
    ref_allele: &str,
    alt_allele: &str,
    alt_vcf_idx: usize,
) -> AlleleDepth {
    let get = |key: &str| -> Option<&str> {
        format_keys
            .iter()
            .position(|k| *k == key)
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
                // For multi-allelic AD (ref,alt1,alt2,...), sum all values for total depth.
                // This matches vcf2maf behavior where all allele counts contribute to depth.
                let ad_total: u32 = raw.iter().filter_map(|v| v.parse::<u32>().ok()).sum();
                let total = dp.or(if ad_total > 0 { Some(ad_total) } else { None });
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
    if get("TAR").is_some() || get("TIR").is_some() {
        let ref_count = get("TAR").and_then(first_num);
        let alt_count = get("TIR").and_then(first_num);
        return AlleleDepth {
            ref_count,
            alt_count,
            total_depth: dp,
        };
    }

    // 5. SomaticSniper: DP4=ref_fwd,ref_rev,alt_fwd,alt_rev
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

    // 6. Ion Torrent: RO (ref obs), AO (alt obs)
    if let (Some(ro), Some(ao)) = (get("RO"), get("AO")) {
        let ref_count = ro.parse::<u32>().ok();
        let alt_count = ao.split(',').next().and_then(|v| v.parse().ok());
        return AlleleDepth {
            ref_count,
            alt_count,
            total_depth: dp,
        };
    }

    // 7. BCOUNT: allele counts as A,C,G,T
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
            let ref_base = ref_allele.chars().next().unwrap_or('N').to_ascii_uppercase();
            let alt_base = alt_allele.chars().next().unwrap_or('N').to_ascii_uppercase();
            let ref_count = parts.get(idx_for(ref_base)).copied();
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
        let d = extract_depth(&keys, &vals, "A", "T", 1);
        assert_eq!(d.ref_count, Some(45));
        assert_eq!(d.alt_count, Some(12));
        assert_eq!(d.total_depth, Some(57));
    }

    #[test]
    fn multiallelic_ad_second_alt() {
        let keys = vec!["GT", "AD", "DP"];
        let vals = vec!["0/2", "10,5,8", "23"];
        let d = extract_depth(&keys, &vals, "A", "G", 2);
        assert_eq!(d.ref_count, Some(10));
        assert_eq!(d.alt_count, Some(8));
        assert_eq!(d.total_depth, Some(23));
    }

    #[test]
    fn varscan() {
        let keys = vec!["GT", "RD", "AD", "DP"];
        let vals = vec!["0/1", "30", "10", "40"];
        let d = extract_depth(&keys, &vals, "A", "T", 1);
        assert_eq!(d.ref_count, Some(30));
        assert_eq!(d.alt_count, Some(10));
    }

    #[test]
    fn dp4() {
        let keys = vec!["GT", "DP4"];
        let vals = vec!["0/1", "10,20,5,3"];
        let d = extract_depth(&keys, &vals, "A", "T", 1);
        assert_eq!(d.ref_count, Some(30));
        assert_eq!(d.alt_count, Some(8));
    }
}
