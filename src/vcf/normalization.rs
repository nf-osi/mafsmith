/// Normalized allele representation for MAF output.
#[derive(Debug, Clone, PartialEq)]
pub struct Normalized {
    /// 1-based position after stripping shared prefix.
    pub pos: u64,
    /// MAF reference allele. "-" for pure insertions.
    pub ref_allele: String,
    /// MAF alt allele. "-" for pure deletions.
    pub alt_allele: String,
}

/// Normalize a VCF REF/ALT pair to MAF convention.
///
/// VCF uses padding bases for indels (e.g. REF=A ALT=ACCTCTT for an insertion).
/// MAF uses "-" for the absent allele (REF=- ALT=CCTCTT) and positions that
/// flank the event rather than the VCF padding convention.
///
/// Algorithm (matches vcf2maf.pl):
/// 1. Strip common prefix characters, advancing position.
/// 2. Replace empty allele strings with "-".
/// Note: vcf2maf.pl does NOT strip common suffix, so ONPs like GGGT>TGGT stay as-is.
pub fn normalize(vcf_pos: u64, ref_allele: &str, alt_allele: &str) -> Normalized {
    let rb = ref_allele.as_bytes();
    let ab = alt_allele.as_bytes();

    // Strip common prefix only, advancing position.
    let prefix = rb.iter().zip(ab.iter()).take_while(|(r, a)| r == a).count();
    // Keep at least one base when alleles are identical (should not happen in valid VCFs).
    let prefix = if prefix == rb.len() && prefix == ab.len() {
        prefix.saturating_sub(1)
    } else {
        prefix
    };
    let pos = vcf_pos + prefix as u64;

    // VCF alleles are ASCII, so byte indices equal char indices — direct string slicing is safe.
    let ref_str = &ref_allele[prefix..];
    let alt_str = &alt_allele[prefix..];

    Normalized {
        pos,
        ref_allele: if ref_str.is_empty() { "-".to_owned() } else { ref_str.to_owned() },
        alt_allele: if alt_str.is_empty() { "-".to_owned() } else { alt_str.to_owned() },
    }
}

/// Compute MAF Start_Position and End_Position from normalized alleles.
///
/// Conventions (GDC/vcf2maf):
/// - Insertion  (ref="-"): Start = pos−1, End = pos   (flanking bases)
/// - Deletion   (alt="-"): Start = pos,   End = pos + len(ref) − 1
/// - SNP / MNP           : Start = pos,   End = pos + len(ref) − 1
pub fn maf_positions(norm: &Normalized) -> (u64, u64) {
    if norm.ref_allele == "-" {
        (norm.pos - 1, norm.pos)
    } else {
        let end = norm.pos + norm.ref_allele.len() as u64 - 1;
        (norm.pos, end)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snp_unchanged() {
        let n = normalize(100, "A", "T");
        assert_eq!(n.pos, 100);
        assert_eq!(n.ref_allele, "A");
        assert_eq!(n.alt_allele, "T");
        let (s, e) = maf_positions(&n);
        assert_eq!((s, e), (100, 100));
    }

    #[test]
    fn insertion_strip_prefix() {
        // VCF: POS=34880555 REF=A ALT=ACCTCTT
        let n = normalize(34880555, "A", "ACCTCTT");
        assert_eq!(n.pos, 34880556);
        assert_eq!(n.ref_allele, "-");
        assert_eq!(n.alt_allele, "CCTCTT");
        let (s, e) = maf_positions(&n);
        // MAF expected: Start=34880555, End=34880556
        assert_eq!((s, e), (34880555, 34880556));
    }

    #[test]
    fn single_base_insertion() {
        // VCF: POS=34880650 REF=T ALT=TG
        let n = normalize(34880650, "T", "TG");
        assert_eq!(n.pos, 34880651);
        assert_eq!(n.ref_allele, "-");
        assert_eq!(n.alt_allele, "G");
        let (s, e) = maf_positions(&n);
        assert_eq!((s, e), (34880650, 34880651));
    }

    #[test]
    fn deletion_strip_prefix() {
        // VCF: POS=41479182 REF=GT ALT=G
        let n = normalize(41479182, "GT", "G");
        assert_eq!(n.pos, 41479183);
        assert_eq!(n.ref_allele, "T");
        assert_eq!(n.alt_allele, "-");
        let (s, e) = maf_positions(&n);
        assert_eq!((s, e), (41479183, 41479183));
    }

    #[test]
    fn mnp_strip_prefix_only() {
        // REF=ACT ALT=AGT: shared prefix 'A' stripped; suffix NOT stripped (matches vcf2maf.pl).
        // Result: CT>GT (DNP) at pos+1, not C>G (SNP).
        let n = normalize(200, "ACT", "AGT");
        assert_eq!(n.pos, 201);
        assert_eq!(n.ref_allele, "CT");
        assert_eq!(n.alt_allele, "GT");
    }

    #[test]
    fn onp_no_suffix_strip() {
        // GGGT>TGGT: no common prefix; suffix GGT shared but NOT stripped (vcf2maf.pl behavior).
        let n = normalize(100, "GGGT", "TGGT");
        assert_eq!(n.pos, 100);
        assert_eq!(n.ref_allele, "GGGT");
        assert_eq!(n.alt_allele, "TGGT");
    }
}
