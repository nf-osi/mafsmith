use anyhow::{bail, Context, Result};
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

        let tsl = get("TSL");
        let transcript_length: u64 = tsl.parse().unwrap_or(0);

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
            hgvsc: get("HGVSc"),
            hgvsp: get("HGVSp"),
            cdna_position: get("cDNA_position"),
            cds_position: get("CDS_position"),
            protein_position: get("Protein_position"),
            amino_acids: get("Amino_acids"),
            codons: get("Codons"),
            existing_variation: get("Existing_variation"),
            canonical,
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

/// Shorten a 3-letter amino acid HGVSp notation to single-letter.
/// e.g. "p.Glu123Lys" → "p.E123K"
pub fn shorten_hgvsp(hgvsp: &str) -> String {
    let mut result = hgvsp.to_owned();
    for (three, one) in THREE_TO_ONE {
        result = result.replace(three, one);
    }
    result
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
    ("Xaa", "X"),
    ("Xle", "J"),
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
