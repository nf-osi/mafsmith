use super::{consequence::consequence_severity, csq::CsqEntry};
use std::collections::HashSet;

/// Biotype priority (lower = preferred). Matches vcf2maf.pl ordering.
fn biotype_rank(biotype: &str) -> u8 {
    match biotype {
        "protein_coding" => 0,
        "LRG_gene" => 1,
        "IG_C_gene" | "IG_D_gene" | "IG_J_gene" | "IG_LV_gene" | "IG_V_gene" => 2,
        "TR_C_gene" | "TR_D_gene" | "TR_J_gene" | "TR_V_gene" => 3,
        bt if bt.ends_with("_pseudogene") || bt == "pseudogene" => 90,
        _ => 50,
    }
}

/// Select the single best transcript annotation for a variant.
///
/// Cascade (highest priority first):
/// 1. User-supplied custom ENST list
/// 2. VEP canonical flag
/// 3. Biotype (protein_coding preferred)
/// 4. Consequence severity (highest impact preferred)
/// 5. Transcript length (longest preferred as tiebreaker)
pub fn select_transcript<'a>(
    entries: &'a [CsqEntry],
    custom_enst: Option<&HashSet<String>>,
) -> Option<&'a CsqEntry> {
    if entries.is_empty() {
        return None;
    }

    // 1. Custom ENST
    if let Some(enst_set) = custom_enst {
        if let Some(e) = entries.iter().find(|e| enst_set.contains(&e.feature)) {
            return Some(e);
        }
    }

    // 2. VEP canonical
    if let Some(e) = entries.iter().find(|e| e.canonical) {
        return Some(e);
    }

    // 3–5. Score-based selection
    entries.iter().min_by_key(|e| {
        let csq_refs: Vec<&str> = e.consequences.iter().map(|s| s.as_str()).collect();
        (
            biotype_rank(&e.biotype),
            consequence_severity(&csq_refs),
            // negate length so longest sorts first under min_by_key
            -(e.transcript_length as i64),
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_entry(feature: &str, canonical: bool, biotype: &str, csq: &str) -> CsqEntry {
        CsqEntry {
            feature: feature.to_owned(),
            canonical,
            biotype: biotype.to_owned(),
            consequences: vec![csq.to_owned()],
            transcript_length: 1000,
            ..Default::default()
        }
    }

    #[test]
    fn canonical_wins_over_biotype() {
        let entries = vec![
            make_entry("ENST1", false, "protein_coding", "intron_variant"),
            make_entry("ENST2", true, "lincRNA", "missense_variant"),
        ];
        let sel = select_transcript(&entries, None).unwrap();
        assert_eq!(sel.feature, "ENST2");
    }

    #[test]
    fn custom_enst_wins_over_canonical() {
        let entries = vec![
            make_entry("ENST1", true, "protein_coding", "synonymous_variant"),
            make_entry("ENST2", false, "protein_coding", "missense_variant"),
        ];
        let mut custom = HashSet::new();
        custom.insert("ENST2".to_owned());
        let sel = select_transcript(&entries, Some(&custom)).unwrap();
        assert_eq!(sel.feature, "ENST2");
    }

    #[test]
    fn severity_tiebreaks_biotype() {
        let entries = vec![
            make_entry("ENST1", false, "protein_coding", "intron_variant"),
            make_entry("ENST2", false, "protein_coding", "stop_gained"),
        ];
        let sel = select_transcript(&entries, None).unwrap();
        assert_eq!(sel.feature, "ENST2");
    }
}
