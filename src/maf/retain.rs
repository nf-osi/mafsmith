//! `--retain-info` / `--retain-fmt`: copy VCF INFO and per-sample FORMAT values
//! into extra MAF columns, reproducing vcf2maf.pl's column naming, de-duplication
//! and value rules.
//!
//! The FORMAT values are not the raw VCF strings. vcf2maf.pl writes the retained
//! columns from its per-sample hashes *after* `FixAlleleDepths` and its genotype
//! fix-ups have rewritten them, so e.g. `t_AD` is the normalized REF,ALT depth
//! list (built from TAR/TIR, BCOUNT, DP4, … when the caller has no AD) and `t_GT`
//! is `./.`, `0/1` or `1/1` when the tumor genotype was missing or unusable. This
//! FORMAT values therefore come from `vcf2maf_pl::PerlSamples`, which ports that
//! Perl logic so the retained columns match vcf2maf.pl byte for byte.

use crate::annotation::vcf2maf_pl::{parse_info, split_sv_preprocess, PerlSamples};
use crate::vcf::VcfRecord;
use std::collections::HashSet;
use tracing::warn;

/// vcf2maf.pl's MAF header up to the retained columns: its standard MAF columns
/// plus every VEP annotation column it emits. A retained INFO/FORMAT column whose
/// name matches one of these (case-insensitively) is dropped, as vcf2maf.pl does.
/// mafsmith writes only a subset of these, but the de-duplication has to use the
/// full list or the two tools would disagree on which columns exist.
#[rustfmt::skip]
const VCF2MAF_PL_HEADER: &[&str] = &[
    "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position",
    "End_Position", "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele",
    "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status",
    "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1",
    "Match_Norm_Seq_Allele2", "Tumor_Validation_Allele1", "Tumor_Validation_Allele2",
    "Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele2", "Verification_Status",
    "Validation_Status", "Mutation_Status", "Sequencing_Phase", "Sequence_Source",
    "Validation_Method", "Score", "BAM_File", "Sequencer", "Tumor_Sample_UUID",
    "Matched_Norm_Sample_UUID", "HGVSc", "HGVSp", "HGVSp_Short", "Transcript_ID", "Exon_Number",
    "t_depth", "t_ref_count", "t_alt_count", "n_depth", "n_ref_count", "n_alt_count",
    "all_effects", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position",
    "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation",
    "ALLELE_NUM", "DISTANCE", "STRAND_VEP", "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "BIOTYPE",
    "CANONICAL", "CCDS", "ENSP", "SWISSPROT", "TREMBL", "UNIPARC", "RefSeq", "SIFT", "PolyPhen",
    "EXON", "INTRON", "DOMAINS", "AF", "AFR_AF", "AMR_AF", "ASN_AF", "EAS_AF", "EUR_AF", "SAS_AF",
    "AA_AF", "EA_AF", "CLIN_SIG", "SOMATIC", "PUBMED", "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE", "IMPACT", "PICK", "VARIANT_CLASS", "TSL", "HGVS_OFFSET", "PHENO",
    "MINIMISED", "GENE_PHENO", "FILTER", "flanking_bps", "vcf_id", "vcf_qual", "gnomADe_AF",
    "gnomADe_AFR_AF", "gnomADe_AMR_AF", "gnomADe_ASJ_AF", "gnomADe_EAS_AF", "gnomADe_FIN_AF",
    "gnomADe_NFE_AF", "gnomADe_OTH_AF", "gnomADe_SAS_AF",
];

/// The retained INFO and FORMAT columns, in output order.
#[derive(Debug, Default, Clone)]
pub struct RetainSpec {
    /// INFO keys; each column is named after its key.
    info: Vec<String>,
    /// FORMAT columns, named `t_<KEY>` / `n_<KEY>`.
    fmt: Vec<String>,
}

impl RetainSpec {
    /// Build the column list the way vcf2maf.pl does: INFO columns first, then a
    /// `t_`/`n_` pair per FORMAT key, skipping any name already in the header.
    pub fn new(retain_info: &[String], retain_fmt: &[String], retain_ann: &[String]) -> Self {
        let mut taken: HashSet<String> = VCF2MAF_PL_HEADER
            .iter()
            .map(|c| c.to_lowercase())
            .chain(retain_ann.iter().map(|c| c.to_lowercase()))
            .chain(std::iter::once("vcf_pos".to_owned()))
            .collect();

        // vcf2maf.pl checks INFO names against the header as it stood before any of
        // them were added, so a key repeated in --retain-info yields two columns.
        let info: Vec<String> = retain_info
            .iter()
            .filter(|key| {
                let free = !taken.contains(&key.to_lowercase());
                if !free {
                    warn!("--retain-info {key}: a MAF column with that name already exists; skipped (as vcf2maf.pl does)");
                }
                free
            })
            .cloned()
            .collect();
        taken.extend(info.iter().map(|c| c.to_lowercase()));

        // Likewise FORMAT names are checked against the header including the INFO
        // columns, but not against FORMAT columns added earlier in the same list.
        let mut fmt = Vec::new();
        for key in retain_fmt {
            for col in [format!("t_{key}"), format!("n_{key}")] {
                if taken.contains(&col.to_lowercase()) {
                    warn!("--retain-fmt {key}: a MAF column named {col} already exists; skipped (as vcf2maf.pl does)");
                } else {
                    fmt.push(col);
                }
            }
        }

        Self { info, fmt }
    }

    pub fn is_empty(&self) -> bool {
        self.info.is_empty() && self.fmt.is_empty()
    }

    /// Header names, in output order.
    pub fn columns(&self) -> Vec<String> {
        self.info.iter().chain(self.fmt.iter()).cloned().collect()
    }

    pub fn has_fmt(&self) -> bool {
        !self.fmt.is_empty()
    }

    /// Values for one VCF record, aligned with [`RetainSpec::columns`]. `samples`
    /// must be given whenever [`RetainSpec::has_fmt`] is true.
    pub fn values(&self, rec: &VcfRecord, samples: Option<&PerlSamples>) -> Vec<String> {
        if self.is_empty() {
            return Vec::new();
        }
        let mut out = Vec::with_capacity(self.info.len() + self.fmt.len());

        if !self.info.is_empty() {
            let (info_line, _) = split_sv_preprocess(&rec.info, &rec.ref_allele, &rec.all_alts);
            let info = parse_info(&info_line);
            out.extend(
                self.info
                    .iter()
                    .map(|k| info.get(k.as_str()).cloned().unwrap_or_default()),
            );
        }

        if let Some(samples) = samples.filter(|_| self.has_fmt()) {
            let tum = samples.tumor.as_ref().map(|s| &s.fmt);
            let nrm = samples.normal.as_ref().map(|s| &s.fmt);
            for col in &self.fmt {
                // vcf2maf.pl strips `t_` and then `n_` in turn, so a FORMAT key that
                // itself starts with `n_` ends up read from the normal sample.
                let mut key = col.as_str();
                let mut val: Option<&String> = None;
                if let Some(k) = key.strip_prefix("t_") {
                    key = k;
                    val = tum.and_then(|t| t.get(key));
                }
                if let Some(k) = key.strip_prefix("n_") {
                    key = k;
                    val = nrm.and_then(|n| n.get(key));
                }
                out.push(val.cloned().unwrap_or_default());
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(v: &[&str]) -> Vec<String> {
        v.iter().map(|x| x.to_string()).collect()
    }

    #[test]
    fn columns_follow_vcf2maf_pl_dedup() {
        let spec = RetainSpec::new(
            &s(&["DP", "AF", "filter", "TLOD"]),
            &s(&["AD", "depth"]),
            &[],
        );
        // AF/FILTER are vcf2maf.pl columns already; t_depth/n_depth too.
        assert_eq!(spec.columns(), s(&["DP", "TLOD", "t_AD", "n_AD"]));
    }
}
