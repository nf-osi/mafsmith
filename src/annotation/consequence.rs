/// Maps Sequence Ontology consequence terms to MAF Variant_Classification.
/// Priority order matches vcf2maf.pl for multi-consequence resolution.
///
/// Generic over the slice element type so callers can pass `&[String]` or `&[&str]`
/// without an intermediate Vec allocation.
pub fn so_to_variant_classification<S: AsRef<str>>(
    consequences: &[S],
    ref_allele: &str,
    alt_allele: &str,
) -> &'static str {
    let inframe = {
        let r = ref_allele.len();
        let a = alt_allele.len();
        r != a && (r.abs_diff(a)) % 3 == 0
    };
    // Sort consequences by severity so the most impactful term is classified first.
    // fastVEP (and VEP) can emit '&'-joined multi-consequence strings where the order
    // does not reflect severity (e.g. "feature_truncation&splice_acceptor_variant").
    // vcf2maf.pl sorts by GetEffectPriority before classifying, so we must do the same.
    let mut sorted: Vec<&str> = consequences.iter().map(|s| s.as_ref()).collect();
    sorted.sort_by_key(|&c| consequence_severity::<&str>(&[c]));
    for csq in &sorted {
        let cls = match *csq {
            "splice_acceptor_variant" | "splice_donor_variant"
            | "transcript_ablation" | "exon_loss_variant" => "Splice_Site",
            "stop_gained" => "Nonsense_Mutation",
            "frameshift_variant" => {
                // Match vcf2maf.pl: Frame_Shift_{Del,Ins} only for DEL/INS var_type;
                // SNP/ONP allele mismatch (multiallelic edge case) → Targeted_Region.
                if ref_allele == "-" || alt_allele.len() > ref_allele.len() {
                    "Frame_Shift_Ins"
                } else if alt_allele == "-" || ref_allele.len() > alt_allele.len() {
                    "Frame_Shift_Del"
                } else {
                    "Targeted_Region"
                }
            }
            "protein_altering_variant" => {
                if inframe {
                    if alt_allele.len() > ref_allele.len() { "In_Frame_Ins" } else { "In_Frame_Del" }
                } else if ref_allele == "-" || alt_allele.len() > ref_allele.len() {
                    "Frame_Shift_Ins"
                } else if alt_allele == "-" || ref_allele.len() > alt_allele.len() {
                    "Frame_Shift_Del"
                } else {
                    "Targeted_Region"
                }
            }
            "stop_lost" => "Nonstop_Mutation",
            "initiator_codon_variant" | "start_lost" => "Translation_Start_Site",
            "disruptive_inframe_insertion"
            | "conservative_inframe_insertion"
            | "inframe_insertion" => "In_Frame_Ins",
            "disruptive_inframe_deletion"
            | "conservative_inframe_deletion"
            | "inframe_deletion" => "In_Frame_Del",
            "missense_variant" | "conservative_missense_variant" | "rare_amino_acid_variant" => {
                "Missense_Mutation"
            }
            "transcript_amplification" => "Intron",
            "splice_region_variant" => "Splice_Region",
            // vcf2maf.pl does not map these splice sub-types → they fall through to Targeted_Region.
            "splice_donor_5th_base_variant"
            | "splice_donor_region_variant"
            | "splice_polypyrimidine_tract_variant" => "Targeted_Region",
            "start_retained_variant"
            | "stop_retained_variant"
            | "synonymous_variant"
            | "incomplete_terminal_codon_variant" => "Silent",
            // NMD_transcript_variant → Silent in vcf2maf.pl.
            "NMD_transcript_variant" => "Silent",
            "coding_sequence_variant" => "Missense_Mutation",
            "mature_miRNA_variant"
            | "exon_variant"
            | "non_coding_exon_variant"
            | "non_coding_transcript_exon_variant"
            | "non_coding_transcript_variant"
            | "nc_transcript_variant" => "RNA",
            "5_prime_UTR_variant"
            | "5_prime_UTR_premature_start_codon_gain_variant" => "5'UTR",
            "3_prime_UTR_variant" => "3'UTR",
            "intron_variant" | "INTRAGENIC" | "intragenic_variant" => "Intron",
            "upstream_gene_variant" => "5'Flank",
            "downstream_gene_variant" => "3'Flank",
            "TF_binding_site_variant"
            | "regulatory_region_variant"
            | "regulatory_region"
            | "intergenic_variant"
            | "intergenic_region" => "IGR",
            // These are in vcf2maf.pl's effect priority table but not in GetVariantClassification
            // → they fall through to the "Targeted_Region" catch-all.
            "TFBS_ablation"
            | "TFBS_amplification"
            | "regulatory_region_ablation"
            | "regulatory_region_amplification"
            | "feature_elongation"
            | "feature_truncation"
            | "sequence_variant"
            | "transcript_variant" => "Targeted_Region",
            _ => continue,
        };
        return cls;
    }
    // vcf2maf.pl returns "Targeted_Region" for unrecognized/empty consequences.
    "Targeted_Region"
}

/// Consequence priority for transcript selection tiebreaking.
/// Matches vcf2maf.pl GetEffectPriority (lower number = higher priority).
///
/// Generic over the slice element type so callers can pass `&[String]` or `&[&str]`
/// without an intermediate Vec allocation.
pub fn consequence_severity<S: AsRef<str>>(consequences: &[S]) -> u8 {
    for csq in consequences {
        let rank: u8 = match csq.as_ref() {
            "transcript_ablation" | "exon_loss_variant" => 1,
            "splice_donor_variant" | "splice_acceptor_variant" => 2,
            "stop_gained" | "frameshift_variant" | "stop_lost" => 3,
            "start_lost" | "initiator_codon_variant" => 4,
            "disruptive_inframe_insertion"
            | "disruptive_inframe_deletion"
            | "conservative_inframe_insertion"
            | "conservative_inframe_deletion"
            | "inframe_insertion"
            | "inframe_deletion"
            | "protein_altering_variant" => 5,
            "missense_variant"
            | "conservative_missense_variant"
            | "rare_amino_acid_variant" => 6,
            "transcript_amplification" => 7,
            "splice_region_variant"
            | "splice_donor_5th_base_variant"
            | "splice_donor_region_variant"
            | "splice_polypyrimidine_tract_variant" => 8,
            "start_retained_variant" | "stop_retained_variant" | "synonymous_variant" => 9,
            "incomplete_terminal_codon_variant" => 10,
            "coding_sequence_variant" | "mature_miRNA_variant" | "exon_variant" => 11,
            "5_prime_UTR_variant" | "5_prime_UTR_premature_start_codon_gain_variant"
            | "3_prime_UTR_variant" => 12,
            "non_coding_exon_variant" | "non_coding_transcript_exon_variant" => 13,
            "non_coding_transcript_variant"
            | "nc_transcript_variant"
            | "intron_variant"
            | "intragenic_variant"
            | "INTRAGENIC" => 14,
            "NMD_transcript_variant" => 15,
            "upstream_gene_variant" | "downstream_gene_variant" => 16,
            "TFBS_ablation"
            | "TFBS_amplification"
            | "TF_binding_site_variant"
            | "regulatory_region_ablation"
            | "regulatory_region_amplification"
            | "regulatory_region_variant"
            | "regulatory_region" => 17,
            "feature_elongation" | "feature_truncation" => 18,
            "intergenic_variant" | "intergenic_region" => 19,
            _ => 20,
        };
        return rank;
    }
    20
}

/// Determine MAF Variant_Type from ref and alt alleles.
pub fn variant_type(ref_allele: &str, alt_allele: &str) -> &'static str {
    if ref_allele == "-" {
        return "INS";
    }
    if alt_allele == "-" {
        return "DEL";
    }
    let r = ref_allele.len();
    let a = alt_allele.len();
    if r < a {
        "INS"
    } else if r > a {
        "DEL"
    } else if r == 1 {
        "SNP"
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
            so_to_variant_classification(&["frameshift_variant"], "T", "-"),
            "Frame_Shift_Del"
        );
    }

    #[test]
    fn frameshift_ins() {
        assert_eq!(
            so_to_variant_classification(&["frameshift_variant"], "-", "T"),
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

    #[test]
    fn splice_polypyrimidine_is_targeted_region() {
        assert_eq!(
            so_to_variant_classification(&["splice_polypyrimidine_tract_variant", "intron_variant"], "A", "T"),
            "Targeted_Region"
        );
    }

    #[test]
    fn nmd_is_silent() {
        assert_eq!(
            so_to_variant_classification(&["NMD_transcript_variant"], "A", "T"),
            "Silent"
        );
    }

    #[test]
    fn unknown_consequence_is_targeted_region() {
        assert_eq!(
            so_to_variant_classification(&["some_unknown_variant"], "A", "T"),
            "Targeted_Region"
        );
    }
}
