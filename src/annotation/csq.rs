use anyhow::{Context, Result};
use std::collections::HashMap;

/// Parsed representation of one CSQ transcript entry.
#[derive(Debug, Clone, Default)]
pub struct CsqEntry {
    pub allele: String,
    pub consequences: Vec<String>,
    pub impact: String,
    pub symbol: String,
    pub gene: String,
    pub feature_type: String,
    pub feature: String,   // transcript ID (ENST...)
    pub biotype: String,
    pub exon: String,
    pub intron: String,
    pub hgvsc: String,
    pub hgvsp: String,
    pub cdna_position: String,
    pub cds_position: String,
    pub protein_position: String,
    pub amino_acids: String,
    pub codons: String,
    pub existing_variation: String,
    pub canonical: bool,
    pub mane_select: bool,
    pub mane_plus_clinical: bool,
    pub transcript_length: u64,
    pub sift: String,
    pub polyphen: String,
    pub strand: String,
    pub extra: HashMap<String, String>,
}

/// Extracted from the VCF header CSQ FORMAT string.
#[derive(Debug, Clone)]
pub struct CsqFormat {
    pub fields: Vec<String>,
    pub index: HashMap<String, usize>,
}

impl CsqFormat {
    /// Parse the CSQ FORMAT from a VCF INFO description line.
    ///
    /// VEP encodes the field list in the Description after "Format: "
    /// e.g. `##INFO=<ID=CSQ,...,Description="... Format: Allele|Consequence|IMPACT|...">`.
    pub fn from_header_description(description: &str) -> Result<Self> {
        let marker = "Format: ";
        let pos = description
            .find(marker)
            .context("CSQ header Description missing 'Format: '")?;
        let format_str = description[pos + marker.len()..].trim_end_matches('"');
        let fields: Vec<String> = format_str.split('|').map(|s| s.to_owned()).collect();
        let index: HashMap<String, usize> =
            fields.iter().enumerate().map(|(i, f)| (f.clone(), i)).collect();
        Ok(Self { fields, index })
    }

    fn get<'a>(&self, parts: &'a [&'a str], name: &str) -> &'a str {
        self.index
            .get(name)
            .and_then(|&i| parts.get(i))
            .copied()
            .unwrap_or("")
    }

    /// Parse a single pipe-delimited CSQ value string into a [`CsqEntry`].
    pub fn parse_entry(&self, raw: &str, retain_fields: &[String]) -> CsqEntry {
        let parts: Vec<&str> = raw.split('|').collect();
        let get = |name: &str| self.get(&parts, name).to_owned();

        let canonical_str = get("CANONICAL");
        let canonical = canonical_str == "YES";

        // vcf2maf.pl checks the MANE field for "MANE_Select" / "MANE_Plus_Clinical".
        // The MANE field holds the designation ("MANE_Select") while MANE_SELECT holds
        // the RefSeq accession. We must check MANE, not MANE_SELECT.
        let mane_field = get("MANE");
        let mane_select = mane_field.contains("MANE_Select");
        let mane_plus_clinical = mane_field.contains("MANE_Plus_Clinical");

        // Transcript length: use the denominator from cDNA_position (e.g. "1480/5971" → 5971).
        // This matches vcf2maf.pl which uses cds_length from the cDNA_position total for
        // tiebreaking when multiple transcripts have the same biotype and consequence severity.
        let cdna_pos = get("cDNA_position");
        let transcript_length: u64 = cdna_pos
            .split('/')
            .nth(1)
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        let consequences: Vec<String> = get("Consequence")
            .split('&')
            .map(|s| s.to_owned())
            .collect();

        let mut extra = HashMap::new();
        for field in retain_fields {
            let val = get(field);
            if !val.is_empty() {
                extra.insert(field.clone(), val);
            }
        }

        CsqEntry {
            allele: get("Allele"),
            consequences,
            impact: get("IMPACT"),
            symbol: get("SYMBOL"),
            gene: get("Gene"),
            feature_type: get("Feature_type"),
            feature: get("Feature"),
            biotype: get("BIOTYPE"),
            exon: get("EXON"),
            intron: get("INTRON"),
            hgvsc: strip_hgvs_prefix(&get("HGVSc")),
            hgvsp: strip_hgvs_prefix(&get("HGVSp")),
            cdna_position: get("cDNA_position"),
            cds_position: get("CDS_position"),
            protein_position: get("Protein_position"),
            amino_acids: get("Amino_acids"),
            codons: get("Codons"),
            existing_variation: get("Existing_variation"),
            canonical,
            mane_select,
            mane_plus_clinical,
            transcript_length,
            sift: get("SIFT"),
            polyphen: get("PolyPhen"),
            strand: get("STRAND"),
            extra,
        }
    }

    /// Parse all CSQ entries from the INFO CSQ value (comma-separated).
    pub fn parse_all(&self, csq_value: &str, retain_fields: &[String]) -> Vec<CsqEntry> {
        csq_value
            .split(',')
            .map(|raw| self.parse_entry(raw, retain_fields))
            .collect()
    }
}

/// Strip the transcript/protein accession prefix that some annotators (e.g. fastVEP) prepend
/// to HGVS notation, e.g. "ENST00000675419.1:c.1286T>G" → "c.1286T>G".
///
/// Uses the LAST colon to match vcf2maf.pl's greedy `s/^.*://` substitution, which correctly
/// handles BND HGVSc that contain a partner-breakpoint coordinate with an additional colon
/// (e.g. "ENST123:c.100ins]chr2:12345]" → "12345]").
fn strip_hgvs_prefix(s: &str) -> String {
    if let Some(pos) = s.rfind(':') {
        s[pos + 1..].to_owned()
    } else {
        s.to_owned()
    }
}

/// Shorten a 3-letter amino acid HGVSp notation to single-letter.
/// e.g. "p.Glu123Lys" → "p.E123K"
pub fn shorten_hgvsp(hgvsp: &str) -> String {
    let mut result = hgvsp.to_owned();
    for &(three, one) in THREE_TO_ONE {
        // contains() check avoids allocating a new String when there's no match,
        // since str::replace always allocates even on zero matches.
        if result.contains(three) {
            result = result.replace(three, one);
        }
    }
    result
}

/// For splice site variants where VEP leaves HGVSp empty, vcf2maf synthesizes a
/// pseudo-HGVSp_Short like "p.X170_splice" from the cDNA position in HGVSc.
/// Mirrors vcf2maf.pl behavior: extract first numeric position from "c.508-4_508+1dup" → aa 170.
pub fn splice_hgvsp_short(hgvsc: &str) -> String {
    let digits: String = hgvsc
        .trim_start_matches("c.")
        .chars()
        .take_while(|c| c.is_ascii_digit())
        .collect();
    if let Ok(cds_pos) = digits.parse::<u64>() {
        let aa_pos = (cds_pos + 2) / 3; // ceiling division: ceil(cds/3)
        format!("p.X{aa_pos}_splice")
    } else {
        String::new()
    }
}

const THREE_TO_ONE: &[(&str, &str)] = &[
    ("Ala", "A"),
    ("Arg", "R"),
    ("Asn", "N"),
    ("Asp", "D"),
    ("Cys", "C"),
    ("Gln", "Q"),
    ("Glu", "E"),
    ("Gly", "G"),
    ("His", "H"),
    ("Ile", "I"),
    ("Leu", "L"),
    ("Lys", "K"),
    ("Met", "M"),
    ("Phe", "F"),
    ("Pro", "P"),
    ("Ser", "S"),
    ("Thr", "T"),
    ("Trp", "W"),
    ("Tyr", "Y"),
    ("Val", "V"),
    ("Ter", "*"),
    ("Sec", "U"),
    ("Pyl", "O"),
    // vcf2maf.pl uses Xxx→X but not Xaa→X; Xaa is left as-is to match vcf2maf output.
    ("Xxx", "X"),
    ("Xle", "J"),
    ("Asx", "B"),
    ("Glx", "Z"),
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_format_from_description() {
        let desc = r#"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL""#;
        let fmt = CsqFormat::from_header_description(desc).unwrap();
        assert_eq!(fmt.index["Allele"], 0);
        assert_eq!(fmt.index["CANONICAL"], 23);
    }

    #[test]
    fn shorten_amino_acids() {
        assert_eq!(shorten_hgvsp("p.Glu123Lys"), "p.E123K");
        assert_eq!(shorten_hgvsp("p.Ter456Trp"), "p.*456W");
    }
}
