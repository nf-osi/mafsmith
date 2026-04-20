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

/// Lightweight transcript info for selection (borrows from the raw CSQ string — no allocations).
pub struct CsqLightEntry<'a> {
    pub feature_type: &'a str,
    pub symbol: &'a str,
    pub feature: &'a str,
    pub biotype: &'a str,
    pub canonical: bool,
    pub transcript_length: u64,
    pub consequences: Vec<&'a str>,
    /// Raw HGVSp value (may have a transcript prefix "ENST…:p.X"); strip with rfind(':').
    pub hgvsp_raw: &'a str,
}

/// Field names extracted for lightweight transcript selection, in slot order.
/// Slots map 1-to-1 with the fields stored in `CsqLightEntry`.
const LIGHT_FIELD_NAMES: [&str; 8] = [
    "Feature_type",  // slot 0
    "SYMBOL",        // slot 1
    "Feature",       // slot 2
    "BIOTYPE",       // slot 3
    "CANONICAL",     // slot 4
    "cDNA_position", // slot 5
    "Consequence",   // slot 6
    "HGVSp",         // slot 7
];

/// Extracted from the VCF header CSQ FORMAT string.
#[derive(Debug, Clone)]
pub struct CsqFormat {
    pub fields: Vec<String>,
    pub index: HashMap<String, usize>,
    /// (csq_field_index, slot_index) pairs for the 8 light fields, sorted by csq_field_index.
    /// Used by parse_light for a single forward scan — no Vec<&str> parts allocation needed.
    light_scan: Vec<(usize, u8)>,
    /// Largest csq_field_index among the 8 light fields (scan can stop after this).
    light_max: usize,
}

impl CsqFormat {
    fn build_light_scan(index: &HashMap<String, usize>) -> (Vec<(usize, u8)>, usize) {
        let mut pairs: Vec<(usize, u8)> = LIGHT_FIELD_NAMES.iter()
            .enumerate()
            .filter_map(|(slot, name)| index.get(*name).map(|&fi| (fi, slot as u8)))
            .collect();
        pairs.sort_by_key(|(fi, _)| *fi);
        let max = pairs.iter().map(|(fi, _)| *fi).max().unwrap_or(0);
        (pairs, max)
    }

    /// An empty format with no fields — used when no CSQ header is present.
    /// parse_all() on any record returns an empty Vec, causing unannotated
    /// variants to be dropped (matching vcf2maf.pl --inhibit-vep behavior).
    pub fn empty() -> Self {
        Self { fields: Vec::new(), index: HashMap::new(), light_scan: Vec::new(), light_max: 0 }
    }

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
        let (light_scan, light_max) = Self::build_light_scan(&index);
        Ok(Self { fields, index, light_scan, light_max })
    }

    fn get<'a>(&self, parts: &'a [&'a str], name: &str) -> &'a str {
        self.index
            .get(name)
            .and_then(|&i| parts.get(i))
            .copied()
            .unwrap_or("")
    }

    /// Parse only the fields needed for transcript selection — no Vec<&str> parts allocation.
    /// Uses a single forward scan through the pipe-delimited entry guided by pre-sorted
    /// field indices, stopping as soon as all 8 target fields have been collected.
    /// All returned values borrow from `raw`.
    pub fn parse_light<'a>(&self, raw: &'a str) -> CsqLightEntry<'a> {
        let mut slots: [&'a str; 8] = [""; 8];
        let mut scan_pos = 0usize; // index into self.light_scan
        let total_targets = self.light_scan.len();

        for (field_idx, part) in raw.split('|').enumerate() {
            if field_idx > self.light_max || scan_pos >= total_targets {
                break;
            }
            if self.light_scan[scan_pos].0 == field_idx {
                slots[self.light_scan[scan_pos].1 as usize] = part;
                scan_pos += 1;
            }
        }

        // slot 5 = cDNA_position
        let transcript_length: u64 = slots[5]
            .split('/').nth(1)
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        CsqLightEntry {
            feature_type: slots[0],
            symbol: slots[1],
            feature: slots[2],
            biotype: slots[3],
            canonical: slots[4] == "YES",
            transcript_length,
            consequences: slots[6].split('&').collect(), // Vec<&'a str>, typically 1-3 items
            hgvsp_raw: slots[7],
        }
    }

    /// Parse all CSQ entries into lightweight selection keys paired with their raw strings.
    /// Returns `(light_entry, raw_csq_entry)` — call `parse_entry(raw, fields)` on the winner.
    pub fn parse_all_light<'a>(&self, csq_value: &'a str) -> Vec<(CsqLightEntry<'a>, &'a str)> {
        csq_value.split(',')
            .map(|raw| (self.parse_light(raw), raw))
            .collect()
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
