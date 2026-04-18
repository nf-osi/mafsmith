/// Maps Sequence Ontology consequence terms to MAF Variant_Classification.
/// Priority order matches vcf2maf.pl for multi-consequence resolution.
pub fn so_to_variant_classification(
    consequences: &[&str],
    ref_allele: &str,
    alt_allele: &str,
) -> &'static str {
    // Walk consequences in SO priority order; first match wins.
    for csq in consequences {
        let cls = match *csq {
            "splice_acceptor_variant" => "Splice_Site",
            "splice_donor_variant" => "Splice_Site",
            "transcript_ablation" => "Splice_Site",
            "exon_loss_variant" => "Splice_Site",
            "stop_gained" => "Nonsense_Mutation",
            "frameshift_variant" => {
                let ref_len = ref_allele.len();
                let alt_len = alt_allele.len();
                if alt_len > ref_len {
                    "Frame_Shift_Ins"
                } else {
                    "Frame_Shift_Del"
                }
            }
            "stop_lost" => "Nonstop_Mutation",
            "initiator_codon_variant" | "start_lost" => "Translation_Start_Site",
            "missense_variant" | "conservative_missense_variant" | "rare_amino_acid_variant" => {
                "Missense_Mutation"
            }
            "protein_altering_variant" => "Missense_Mutation",
            "inframe_insertion" | "disruptive_inframe_insertion" => "In_Frame_Ins",
            "inframe_deletion" | "disruptive_inframe_deletion" => "In_Frame_Del",
            "synonymous_variant"
            | "stop_retained_variant"
            | "start_retained_variant"
            | "incomplete_terminal_codon_variant" => "Silent",
            "coding_sequence_variant" | "sequence_variant" => "Missense_Mutation",
            "transcript_amplification" | "transcript_variant" => "Intron",
            "5_prime_UTR_variant"
            | "5_prime_UTR_premature_start_codon_gain_variant"
            | "upstream_gene_variant" => "5'UTR",
            "3_prime_UTR_variant" | "downstream_gene_variant" => "3'UTR",
            "intron_variant" | "INTRAGENIC" | "intragenic_variant" => "Intron",
            "splice_region_variant"
            | "splice_donor_5th_base_variant"
            | "splice_donor_region_variant"
            | "splice_polypyrimidine_tract_variant" => "Splice_Region",
            "non_coding_transcript_exon_variant"
            | "non_coding_transcript_variant"
            | "mature_miRNA_variant"
            | "NMD_transcript_variant" => "RNA",
            "feature_truncation" => "Splice_Site",
            "feature_elongation" => "Intron",
            "regulatory_region_variant"
            | "regulatory_region_amplification"
            | "regulatory_region_ablation"
            | "TF_binding_site_variant"
            | "TFBS_amplification"
            | "TFBS_ablation" => "IGR",
            "intergenic_variant" | "intergenic_region" => "IGR",
            _ => continue,
        };
        return cls;
    }
    "IGR"
}

/// Consequence priority for transcript selection tiebreaking (lower = higher priority).
pub fn consequence_severity(consequences: &[&str]) -> u8 {
    for csq in consequences {
        let rank: u8 = match *csq {
            "transcript_ablation" => 1,
            "splice_acceptor_variant" | "splice_donor_variant" => 2,
            "stop_gained" => 3,
            "frameshift_variant" => 4,
            "stop_lost" => 5,
            "start_lost" | "initiator_codon_variant" => 6,
            "transcript_amplification" => 7,
            "inframe_insertion" | "disruptive_inframe_insertion" => 8,
            "inframe_deletion" | "disruptive_inframe_deletion" => 9,
            "missense_variant" | "protein_altering_variant" => 10,
            "splice_region_variant" => 11,
            "splice_donor_5th_base_variant" | "splice_donor_region_variant" => 12,
            "synonymous_variant" => 13,
            "stop_retained_variant" | "start_retained_variant" => 14,
            "coding_sequence_variant" | "sequence_variant" => 15,
            "5_prime_UTR_variant" => 16,
            "3_prime_UTR_variant" => 17,
            "non_coding_transcript_exon_variant" => 18,
            "intron_variant" => 19,
            "NMD_transcript_variant" => 20,
            "non_coding_transcript_variant" => 21,
            "upstream_gene_variant" => 22,
            "downstream_gene_variant" => 23,
            "TFBS_ablation" | "TFBS_amplification" => 24,
            "TF_binding_site_variant" => 25,
            "regulatory_region_variant"
            | "regulatory_region_amplification"
            | "regulatory_region_ablation" => 26,
            "intergenic_variant" | "intergenic_region" | "INTRAGENIC" => 27,
            _ => 28,
        };
        return rank;
    }
    28
}

/// Determine MAF Variant_Type from ref and alt alleles.
pub fn variant_type(ref_allele: &str, alt_allele: &str) -> &'static str {
    let r = ref_allele.len();
    let a = alt_allele.len();
    if r == 1 && a == 1 {
        "SNP"
    } else if r < a {
        "INS"
    } else if r > a {
        "DEL"
    } else if r == 2 {
        "DNP"
    } else if r == 3 {
        "TNP"
    } else {
        "ONP"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snp_classification() {
        assert_eq!(
            so_to_variant_classification(&["missense_variant"], "A", "T"),
            "Missense_Mutation"
        );
    }

    #[test]
    fn frameshift_del() {
        assert_eq!(
            so_to_variant_classification(&["frameshift_variant"], "AT", "A"),
            "Frame_Shift_Del"
        );
    }

    #[test]
    fn frameshift_ins() {
        assert_eq!(
            so_to_variant_classification(&["frameshift_variant"], "A", "AT"),
            "Frame_Shift_Ins"
        );
    }

    #[test]
    fn variant_types() {
        assert_eq!(variant_type("A", "T"), "SNP");
        assert_eq!(variant_type("A", "AT"), "INS");
        assert_eq!(variant_type("AT", "A"), "DEL");
        assert_eq!(variant_type("AT", "GC"), "DNP");
    }
}
