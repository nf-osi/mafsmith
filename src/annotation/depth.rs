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

    /// True when at least one field was populated (i.e. extract_depth found a depth source).
    /// Used to distinguish "extracted but None" (→ ".") from "no depth info at all" (→ "").
    pub fn has_data(&self) -> bool {
        self.ref_count.is_some() || self.alt_count.is_some() || self.total_depth.is_some()
    }
}

/// Reconcile the total depth with the per-allele depths, mirroring the sanity
/// check at the end of vcf2maf.pl's `FixAlleleDepths`:
///
/// - DP present, but a REF/ALT depth exceeds it (or the two together do): DP is
///   counting a different read set than the allele depths (filtered or
///   supplementary reads), so DP is replaced by the sum of all allele depths.
/// - DP absent: DP becomes the sum of all allele depths.
/// - Otherwise DP is reported as given.
///
/// `depths_sum` is the sum over *every* allele depth (REF plus all ALTs), which
/// is what vcf2maf.pl sums — not just the selected REF/ALT pair. `slots_defined`
/// says whether the caller filled both the REF and selected-ALT slots at all;
/// vcf2maf.pl only derives a missing DP when both slots exist.
///
/// The zero-guards reproduce Perl truthiness: a depth of 0 does not take part in
/// the comparisons, so e.g. AD="0,82" with DP=652 leaves DP alone.
fn fix_total_depth(
    ref_count: Option<u32>,
    alt_count: Option<u32>,
    dp: Option<u32>,
    depths_sum: u32,
    slots_defined: bool,
) -> Option<u32> {
    let r = ref_count.filter(|&v| v != 0);
    let a = alt_count.filter(|&v| v != 0);
    match dp {
        Some(d) => {
            let inconsistent = r.is_some_and(|v| v > d)
                || a.is_some_and(|v| v > d)
                || matches!((r, a), (Some(x), Some(y)) if x + y > d);
            Some(if inconsistent { depths_sum } else { d })
        }
        None => slots_defined.then_some(depths_sum),
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
                let ad_total: u32 = raw.iter().filter_map(|v| v.parse::<u32>().ok()).sum();
                // Both slots exist as long as the AD array reaches the selected ALT;
                // "." entries count as present (vcf2maf.pl tests `defined`, not numeric).
                let slots_defined = alt_vcf_idx < raw.len();
                return AlleleDepth {
                    ref_count,
                    alt_count,
                    total_depth: fix_total_depth(ref_count, alt_count, dp, ad_total, slots_defined),
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
                total_depth: fix_total_depth(Some(r), Some(a), dp, r + a, true),
            };
        }
    }

    // 3. Strelka SNVs: AU, CU, GU, TU (first value = tier1)
    let base_keys = ["AU", "CU", "GU", "TU"];
    if base_keys.iter().any(|k| get(k).is_some()) {
        let ref_base = ref_allele
            .chars()
            .next()
            .unwrap_or('N')
            .to_ascii_uppercase();
        let alt_base = alt_allele
            .chars()
            .next()
            .unwrap_or('N')
            .to_ascii_uppercase();
        let base_key = |b: char| match b {
            'A' => "AU",
            'C' => "CU",
            'G' => "GU",
            'T' => "TU",
            _ => "AU",
        };
        let ref_count = get(base_key(ref_base)).and_then(first_num);
        let alt_count = get(base_key(alt_base)).and_then(first_num);
        // vcf2maf.pl sums the per-allele base counts across all alleles; for the
        // bi-allelic SNVs Strelka emits that is exactly ref+alt.
        let sum = ref_count.unwrap_or(0) + alt_count.unwrap_or(0);
        return AlleleDepth {
            ref_count,
            alt_count,
            total_depth: fix_total_depth(ref_count, alt_count, dp, sum, true),
        };
    }

    // 4. Strelka indels: TAR (ref), TIR (alt).
    // TAR/TIR count reads at a specific indel tier, which can exceed the DP field —
    // the sanity check then rebuilds the total from TAR[0]+TIR[0].
    if get("TAR").is_some() || get("TIR").is_some() {
        let ref_count = get("TAR").and_then(first_num);
        let alt_count = get("TIR").and_then(first_num);
        let sum = ref_count.unwrap_or(0) + alt_count.unwrap_or(0);
        return AlleleDepth {
            ref_count,
            alt_count,
            total_depth: fix_total_depth(ref_count, alt_count, dp, sum, true),
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
                    'A' => 0usize,
                    'C' => 1,
                    'G' => 2,
                    'T' => 3,
                    _ => 0,
                };
                let ref_char = ref_allele.chars().next().unwrap_or('N');
                let alt_base = alt_allele
                    .chars()
                    .next()
                    .unwrap_or('N')
                    .to_ascii_uppercase();
                // In strict mode, vcf2maf.pl does a case-sensitive hash lookup over {A,C,G,T}.
                // Lowercase ref alleles (e.g. on decoy contigs) fail that lookup → ref_count='.'.
                let ref_count = if strict && !matches!(ref_char, 'A' | 'C' | 'G' | 'T') {
                    None
                } else {
                    bc_parts
                        .get(idx_for(ref_char.to_ascii_uppercase()))
                        .copied()
                };
                let alt_count = bc_parts.get(idx_for(alt_base)).copied();
                let sum = ref_count.unwrap_or(0) + alt_count.unwrap_or(0);
                return AlleleDepth {
                    ref_count,
                    alt_count,
                    total_depth: fix_total_depth(ref_count, alt_count, dp, sum, true),
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
                    total_depth: fix_total_depth(
                        Some(ref_count),
                        Some(alt_count),
                        dp,
                        ref_count + alt_count,
                        true,
                    ),
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
        let alt_count = ao
            .split(',')
            .nth(alt_vcf_idx.saturating_sub(1))
            .and_then(|v| v.parse().ok());
        // vcf2maf.pl's depth array here is (RO, AO...), so the sum spans every ALT.
        let ao_total: u32 = ao.split(',').filter_map(|v| v.parse::<u32>().ok()).sum();
        let slots_defined = alt_vcf_idx.saturating_sub(1) < ao.split(',').count();
        return AlleleDepth {
            ref_count,
            alt_count,
            total_depth: fix_total_depth(
                ref_count,
                alt_count,
                dp,
                ref_count.unwrap_or(0) + ao_total,
                slots_defined,
            ),
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
            let alt_base = alt_allele
                .chars()
                .next()
                .unwrap_or('N')
                .to_ascii_uppercase();
            let ref_count = if strict && !matches!(ref_char, 'A' | 'C' | 'G' | 'T') {
                None
            } else {
                parts.get(idx_for(ref_char.to_ascii_uppercase())).copied()
            };
            let alt_count = parts.get(idx_for(alt_base)).copied();
            let sum = ref_count.unwrap_or(0) + alt_count.unwrap_or(0);
            return AlleleDepth {
                ref_count,
                alt_count,
                total_depth: fix_total_depth(ref_count, alt_count, dp, sum, true),
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
    fn ad_sum_replaces_inconsistent_dp() {
        // GIAB NIST v4.2.1: ref+alt = 174 > DP = 157, so DP is rebuilt from sum(AD).
        let keys = vec!["GT", "PS", "DP", "ADALL", "AD", "GQ"];
        let vals = vec!["1/1", ".", "157", "0,70", "87,87", "237"];
        let d = extract_depth(&keys, &vals, "T", "C", 1, false);
        assert_eq!(d.ref_count, Some(87));
        assert_eq!(d.alt_count, Some(87));
        assert_eq!(d.total_depth, Some(174));

        // ref+alt = 82 fits inside DP = 652 → DP is reported as-is. The ref depth of 0
        // must not participate in the comparisons (Perl truthiness).
        let vals2 = vec!["1/1", ".", "652", "16,234", "0,82", "312"];
        let d2 = extract_depth(&keys, &vals2, "A", "G", 1, false);
        assert_eq!(d2.ref_count, Some(0));
        assert_eq!(d2.alt_count, Some(82));
        assert_eq!(d2.total_depth, Some(652));
    }

    #[test]
    fn multiallelic_other_alt_depth_does_not_inflate_dp() {
        // The selected REF/ALT pair (10+5) fits inside DP=100, so DP stands even though
        // the untouched third allele pushes sum(AD) to 215. A plain max(DP, sum(AD))
        // would wrongly report 215 here.
        let keys = vec!["GT", "AD", "DP"];
        let vals = vec!["0/1", "10,5,200", "100"];
        let d = extract_depth(&keys, &vals, "A", "G", 1, false);
        assert_eq!(d.total_depth, Some(100));

        // Selecting the deep third allele does trip the check, and then the total comes
        // from the whole AD array, as in vcf2maf.pl.
        let d2 = extract_depth(&keys, &vals, "A", "T", 2, false);
        assert_eq!(d2.alt_count, Some(200));
        assert_eq!(d2.total_depth, Some(215));
    }

    #[test]
    fn missing_dp_falls_back_to_ad_sum() {
        let keys = vec!["GT", "AD", "DP"];
        let vals = vec!["0/1", "45,12", "."];
        let d = extract_depth(&keys, &vals, "A", "T", 1, false);
        assert_eq!(d.total_depth, Some(57));
    }

    #[test]
    fn varscan_dp_rebuilt_when_below_allele_depths() {
        // VarScan DP can exclude reads that RD/AD count; vcf2maf.pl then uses RD+AD.
        let keys = vec!["GT", "RD", "AD", "DP"];
        let vals = vec!["0/1", "30", "20", "40"];
        let d = extract_depth(&keys, &vals, "A", "T", 1, false);
        assert_eq!(d.total_depth, Some(50));
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
        assert_eq!(d.alt_count, Some(7)); // BCOUNT[C], not DP4 alt(8)
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
        assert_eq!(d_strict.ref_count, None); // lowercase ref → None in strict mode
        assert_eq!(d_strict.alt_count, Some(16)); // C uppercase → correct
        let d_normal = extract_depth(&keys, &vals, "t", "C", 1, false);
        assert_eq!(d_normal.ref_count, Some(157)); // converted to uppercase, finds T
        assert_eq!(d_normal.alt_count, Some(16));
    }
}
