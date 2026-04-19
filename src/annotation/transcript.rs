use super::{consequence::consequence_severity, csq::CsqEntry};
use std::collections::HashSet;

/// Biotype priority (lower = preferred). Matches vcf2maf.pl GetBiotypePriority exactly.
fn biotype_rank(biotype: &str) -> u8 {
    match biotype {
        "protein_coding" => 1,
        "LRG_gene" => 2,
        "IG_C_gene" | "IG_D_gene" | "IG_J_gene" | "IG_LV_gene" | "IG_V_gene" => 2,
        "TR_C_gene" | "TR_D_gene" | "TR_J_gene" | "TR_V_gene" => 2,
        "miRNA" | "snRNA" | "snoRNA" | "ribozyme" | "tRNA" | "sRNA" | "scaRNA" | "rRNA"
        | "scRNA" | "lincRNA" | "lncRNA" | "bidirectional_promoter_lncrna"
        | "bidirectional_promoter_lncRNA" => 3,
        "known_ncrna" | "vaultRNA" | "macro_lncRNA" | "Mt_tRNA" | "Mt_rRNA" => 4,
        "antisense" | "antisense_RNA" | "sense_intronic" | "sense_overlapping"
        | "3prime_overlapping_ncrna" | "3prime_overlapping_ncRNA" | "misc_RNA"
        | "non_coding" => 5,
        "regulatory_region" | "disrupted_domain" | "processed_transcript"
        | "protein_coding_CDS_not_defined" | "TEC" => 6,
        "TF_binding_site" | "CTCF_binding_site" | "promoter_flanking_region" | "enhancer"
        | "promoter" | "open_chromatin_region" | "retained_intron"
        | "nonsense_mediated_decay" | "non_stop_decay" | "ambiguous_orf" => 7,
        bt if bt.ends_with("_pseudogene") || bt == "pseudogene" => 8,
        "artifact" => 9,
        _ => 10,
    }
}

/// Select the single best transcript annotation for a variant.
///
/// Matches vcf2maf.pl logic:
/// 1. Sort all transcripts by biotype → consequence severity → transcript length (longest first).
/// 2. Identify the "worst affected gene": first sorted entry with a non-empty gene symbol.
/// 3. Pick the user-preferred (custom ENST) isoform of that gene.
/// 4. Else pick the canonical isoform of that gene.
/// 5. Else pick any user-preferred isoform with a gene symbol.
/// 6. Else pick any canonical isoform with a gene symbol.
/// 7. Else fall back to the first sorted entry.
///
/// Only Transcript-type features are considered (regulatory/TFBS entries are excluded).
pub fn select_transcript<'a>(
    entries: &'a [CsqEntry],
    custom_enst: Option<&HashSet<String>>,
) -> Option<&'a CsqEntry> {
    if entries.is_empty() {
        return None;
    }

    // Filter to Transcript features only (vcf2maf.pl skips non-Transcript entries).
    // Fall back to all entries if none are typed as Transcript.
    let candidates: Vec<&'a CsqEntry> = {
        let ts: Vec<&'a CsqEntry> = entries
            .iter()
            .filter(|e| e.feature_type.is_empty() || e.feature_type == "Transcript")
            .collect();
        if ts.is_empty() { entries.iter().collect() } else { ts }
    };

    // Sort by biotype priority → consequence severity → transcript length (longest wins).
    let mut sorted: Vec<&'a CsqEntry> = candidates;
    sorted.sort_by_key(|e| (
        biotype_rank(&e.biotype),
        consequence_severity(&e.consequences),  // generic: &[String] works directly
        -(e.transcript_length as i64),
    ));

    // Find the "worst affected gene": first sorted entry with a non-empty gene symbol.
    let maf_gene: Option<&str> = sorted
        .iter()
        .copied()
        .find(|e| !e.symbol.is_empty())
        .map(|e| e.symbol.as_str());

    // Priority 1: user-preferred isoform of the target gene.
    if let (Some(gene), Some(enst_set)) = (maf_gene, custom_enst) {
        if let Some(e) = sorted
            .iter()
            .copied()
            .find(|e| e.symbol == gene && enst_set.contains(&e.feature))
        {
            return Some(e);
        }
    }

    // Priority 2: canonical isoform of the target gene.
    if let Some(gene) = maf_gene {
        if let Some(e) = sorted
            .iter()
            .copied()
            .find(|e| e.symbol == gene && e.canonical)
        {
            return Some(e);
        }
    }

    // Priority 3: user-preferred isoform for any gene with a symbol.
    if let Some(enst_set) = custom_enst {
        if let Some(e) = sorted
            .iter()
            .copied()
            .find(|e| !e.symbol.is_empty() && enst_set.contains(&e.feature))
        {
            return Some(e);
        }
    }

    // Priority 4: canonical isoform for any gene with a symbol.
    if let Some(e) = sorted
        .iter()
        .copied()
        .find(|e| !e.symbol.is_empty() && e.canonical)
    {
        return Some(e);
    }

    // Priority 5: first sorted entry regardless of symbol.
    sorted.first().copied()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_entry_sym(feature: &str, canonical: bool, biotype: &str, csq: &str, symbol: &str) -> CsqEntry {
        CsqEntry {
            feature: feature.to_owned(),
            canonical,
            biotype: biotype.to_owned(),
            consequences: vec![csq.to_owned()],
            transcript_length: 1000,
            symbol: symbol.to_owned(),
            feature_type: String::from("Transcript"),
            ..Default::default()
        }
    }

    fn make_entry(feature: &str, canonical: bool, biotype: &str, csq: &str) -> CsqEntry {
        make_entry_sym(feature, canonical, biotype, csq, "")
    }

    #[test]
    fn protein_coding_beats_canonical_lncrna_no_symbol() {
        // When neither has a symbol, biotype sort determines winner.
        let entries = vec![
            make_entry("ENST1", false, "protein_coding", "intron_variant"),
            make_entry("ENST2", true, "lincRNA", "missense_variant"),
        ];
        let sel = select_transcript(&entries, None).unwrap();
        assert_eq!(sel.feature, "ENST1");
    }

    #[test]
    fn custom_enst_wins_for_target_gene() {
        let entries = vec![
            make_entry_sym("ENST1", true, "protein_coding", "synonymous_variant", "GENE"),
            make_entry_sym("ENST2", false, "protein_coding", "missense_variant", "GENE"),
        ];
        let mut custom = HashSet::new();
        custom.insert("ENST2".to_owned());
        let sel = select_transcript(&entries, Some(&custom)).unwrap();
        assert_eq!(sel.feature, "ENST2");
    }

    #[test]
    fn biotype_tiebreaks_among_canonicals() {
        let entries = vec![
            make_entry("ENST_RNA", true, "lncRNA", "non_coding_transcript_exon_variant"),
            make_entry("ENST_PC", true, "protein_coding", "missense_variant"),
        ];
        let sel = select_transcript(&entries, None).unwrap();
        assert_eq!(sel.feature, "ENST_PC");
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

    #[test]
    fn canonical_of_target_gene_wins() {
        // Gene A has non-canonical protein_coding transcripts; gene B has canonical lncRNA.
        // vcf2maf picks gene A (best biotype) then its canonical isoform.
        // If no canonical exists for A, falls back to first sorted (non-canonical A).
        let entries = vec![
            make_entry_sym("ENST_A1", false, "protein_coding", "intron_variant", "GENE_A"),
            make_entry_sym("ENST_A2", true,  "protein_coding", "intron_variant", "GENE_A"),
            make_entry_sym("ENST_B",  true,  "lincRNA",        "missense_variant", "GENE_B"),
        ];
        let sel = select_transcript(&entries, None).unwrap();
        // GENE_A is target gene (protein_coding wins biotype sort);
        // ENST_A2 is its canonical isoform.
        assert_eq!(sel.feature, "ENST_A2");
    }

    #[test]
    fn fallback_to_canonical_any_gene_when_target_gene_has_none() {
        // Gene A has only non-canonical protein_coding; gene B has canonical lncRNA.
        // vcf2maf: target=GENE_A (protein_coding), no canonical for A
        //   → priority 4: any canonical with symbol → picks ENST_B (GENE_B canonical).
        let entries = vec![
            make_entry_sym("ENST_A", false, "protein_coding", "intron_variant", "GENE_A"),
            make_entry_sym("ENST_B", true,  "lincRNA",        "missense_variant", "GENE_B"),
        ];
        let sel = select_transcript(&entries, None).unwrap();
        assert_eq!(sel.feature, "ENST_B");
    }

    #[test]
    fn fallback_to_first_sorted_when_nothing_canonical() {
        // No canonical entries anywhere → falls back to first sorted (protein_coding).
        let entries = vec![
            make_entry_sym("ENST_A", false, "protein_coding", "intron_variant", "GENE_A"),
            make_entry_sym("ENST_B", false, "lincRNA",        "missense_variant", "GENE_B"),
        ];
        let sel = select_transcript(&entries, None).unwrap();
        assert_eq!(sel.feature, "ENST_A");
    }
}
