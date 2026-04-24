use crate::cli::Maf2vcfArgs;
use anyhow::{Context, Result};
use std::{
    collections::HashMap,
    fs,
    io::{BufRead, BufReader, Write},
    path::Path,
};
use tracing::info;

pub async fn run(args: Maf2vcfArgs) -> Result<()> {
    let data_dir = crate::commands::fetch::data_dir(None)?;
    let genome_str = args.genome.ncbi_build();

    let ref_fasta = args
        .ref_fasta
        .unwrap_or_else(|| data_dir.join(genome_str).join("reference.fa"));

    let maf = fs::File::open(&args.input_maf)
        .with_context(|| format!("Cannot open {}", args.input_maf.display()))?;
    let mut reader = BufReader::new(maf);

    // --- Pass 1: parse header and collect column indices ---
    let mut col_index: HashMap<String, usize> = HashMap::new();
    let mut line = String::new();
    loop {
        line.clear();
        reader.read_line(&mut line)?;
        let trimmed = line.trim();
        if trimmed.starts_with('#') {
            continue;
        }
        for (i, col) in trimmed.split('\t').enumerate() {
            col_index.insert(col.to_owned(), i);
        }
        break;
    }

    let req = |name: &str| -> Result<usize> {
        col_index
            .get(name)
            .copied()
            .with_context(|| format!("MAF missing required column '{name}'"))
    };
    let chrom_idx = req("Chromosome")?;
    let start_idx = req("Start_Position")?;
    let ref_idx = req("Reference_Allele")?;
    let alt_idx = req("Tumor_Seq_Allele2")?;
    let tsa1_idx = req("Tumor_Seq_Allele1")?;
    let tumor_idx = req("Tumor_Sample_Barcode")?;
    let normal_idx = req("Matched_Norm_Sample_Barcode")?;

    // Optional depth columns
    let t_depth_idx = col_index.get("t_depth").copied();
    let t_ref_idx = col_index.get("t_ref_count").copied();
    let t_alt_idx = col_index.get("t_alt_count").copied();
    let n_depth_idx = col_index.get("n_depth").copied();
    let n_ref_idx = col_index.get("n_ref_count").copied();
    let n_alt_idx = col_index.get("n_alt_count").copied();
    // Normal allele columns (for multi-allelic reconstruction when normal has a different alt)
    let mnsa1_idx = col_index.get("Match_Norm_Seq_Allele1").copied();
    let mnsa2_idx = col_index.get("Match_Norm_Seq_Allele2").copied();

    let mut records: Vec<Vec<String>> = Vec::new();
    let mut sample_pairs: indexmap::IndexSet<(String, String)> = indexmap::IndexSet::new();

    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<String> = trimmed.split('\t').map(|s| s.to_owned()).collect();
        let tumor = cols.get(tumor_idx).cloned().unwrap_or_default();
        let normal = cols.get(normal_idx).cloned().unwrap_or_default();
        sample_pairs.insert((tumor, normal));
        records.push(cols);
    }

    // --- Write combined VCF ---
    let samples: Vec<String> = sample_pairs
        .iter()
        .flat_map(|(t, n)| [t.clone(), n.clone()])
        .collect::<indexmap::IndexSet<String>>()
        .into_iter()
        .collect();

    let out = fs::File::create(&args.output_vcf)
        .with_context(|| format!("Cannot create {}", args.output_vcf.display()))?;
    let mut w = std::io::BufWriter::new(out);

    let has_depth = t_depth_idx.is_some()
        && t_ref_idx.is_some()
        && t_alt_idx.is_some()
        && n_depth_idx.is_some()
        && n_ref_idx.is_some()
        && n_alt_idx.is_some();

    write_vcf_header(&mut w, genome_str, &samples, has_depth)?;

    for row in &records {
        let chrom = row.get(chrom_idx).map(|s| s.as_str()).unwrap_or(".");
        let start: u64 = row.get(start_idx).and_then(|s| s.parse().ok()).unwrap_or(0);
        let ref_allele = row.get(ref_idx).map(|s| s.as_str()).unwrap_or(".");
        let alt_allele = row.get(alt_idx).map(|s| s.as_str()).unwrap_or(".");
        let tsa1 = row.get(tsa1_idx).map(|s| s.as_str()).unwrap_or(ref_allele);
        let tumor_id = row.get(tumor_idx).map(|s| s.as_str()).unwrap_or("");
        let normal_id = row.get(normal_idx).map(|s| s.as_str()).unwrap_or("");

        // Determine multi-allelic configuration from MAF allele columns.
        //
        // Case A — Normal has a different alt (e.g. SomaticSniper: normal has somatic mutation
        //   at the same site as tumor): TSA1=REF, TSA2=G (tumor), MNSA2=A (normal-specific alt)
        //   → ALT=G,A, tumor GT=0/1, normal GT=0/2
        //
        // Case B — Tumor is compound het: TSA1=G (secondary), TSA2=C (effective alt), MNSA2=REF
        //   → ALT=C,G, tumor GT=1/2, normal GT=0/0
        let normal_alt2 = mnsa2_idx
            .and_then(|i| row.get(i))
            .map(|s| s.as_str())
            .unwrap_or("");
        let normal_ref1 = mnsa1_idx
            .and_then(|i| row.get(i))
            .map(|s| s.as_str())
            .unwrap_or(ref_allele);

        // Normal MNSA2 is a different alt when it's non-empty, non-ref, non-"-", and differs
        // from the tumor's effective alt. We do NOT require MNSA1==REF so we also handle
        // hom-alt normal (MNSA1==MNSA2==G, both different from REF).
        let normal_has_diff_alt = !normal_alt2.is_empty()
            && normal_alt2 != "-"
            && normal_alt2 != ref_allele
            && normal_alt2 != alt_allele;

        let tumor_has_two_alts = tsa1 != ref_allele
            && tsa1 != "-"
            && tsa1 != alt_allele;

        let (vcf_pos, vcf_ref, vcf_alt) = if normal_has_diff_alt {
            if ref_allele != "-" && alt_allele != "-" && normal_alt2 != "-" {
                (start, ref_allele.to_owned(), format!("{alt_allele},{normal_alt2}"))
            } else {
                maf_alleles_to_vcf(&ref_fasta, chrom, start, ref_allele, alt_allele)
            }
        } else if tumor_has_two_alts {
            if ref_allele == "-" {
                // Both alts are insertions at the same anchor position
                let anchor = crate::fasta::fasta_fetch_base(&ref_fasta, chrom, start)
                    .map(|c| c.to_string())
                    .unwrap_or_else(|| "N".to_string());
                (start, anchor.clone(), format!("{anchor}{alt_allele},{anchor}{tsa1}"))
            } else if alt_allele != "-" && tsa1 != "-" {
                // Both alts are SNP/MNP
                (start, ref_allele.to_owned(), format!("{alt_allele},{tsa1}"))
            } else if alt_allele == "-" {
                // Effective alt is full deletion; TSA1 is a partial/different deletion allele.
                // Anchor at start-1; REF=anchor+ref, ALT1=anchor (full del), ALT2=anchor+TSA1 (partial del).
                let anchor = crate::fasta::fasta_fetch_base(&ref_fasta, chrom, start - 1)
                    .map(|c| c.to_string())
                    .unwrap_or_else(|| "N".to_string());
                let vcf_ref = format!("{anchor}{ref_allele}");
                let vcf_alt2 = format!("{anchor}{tsa1}");
                (start - 1, vcf_ref, format!("{anchor},{vcf_alt2}"))
            } else {
                maf_alleles_to_vcf(&ref_fasta, chrom, start, ref_allele, alt_allele)
            }
        } else {
            maf_alleles_to_vcf(&ref_fasta, chrom, start, ref_allele, alt_allele)
        };

        let is_multiallelic = normal_has_diff_alt || tumor_has_two_alts;

        let tumor_gt = if tumor_has_two_alts {
            "1/2"
        } else if tsa1 == alt_allele && tsa1 != ref_allele {
            // Both MAF alleles are the same alt → hom-alt (includes deletions where TSA1=TSA2="-")
            "1/1"
        } else {
            "0/1"
        };
        let normal_gt = if normal_has_diff_alt {
            // Hom-alt normal (MNSA1==MNSA2, both ≠ REF) → 2/2; otherwise het 0/2
            if normal_ref1 != ref_allele && normal_ref1 == normal_alt2 {
                "2/2"
            } else {
                "0/2"
            }
        } else {
            "0/0"
        };

        // Build FORMAT and per-sample genotype strings
        let (fmt, tumor_fmt, normal_fmt) = if has_depth {
            let t_ref = row.get(t_ref_idx.unwrap()).map(|s| s.as_str()).unwrap_or(".");
            let t_alt = row.get(t_alt_idx.unwrap()).map(|s| s.as_str()).unwrap_or(".");
            let t_dep = row.get(t_depth_idx.unwrap()).map(|s| s.as_str()).unwrap_or(".");
            let n_ref = row.get(n_ref_idx.unwrap()).map(|s| s.as_str()).unwrap_or(".");
            let n_alt = row.get(n_alt_idx.unwrap()).map(|s| s.as_str()).unwrap_or(".");
            let n_dep = row.get(n_depth_idx.unwrap()).map(|s| s.as_str()).unwrap_or(".");
            let (t_ad, n_ad) = if is_multiallelic {
                // Tri-allelic AD: REF,ALT1,ALT2 where the "other" allele gets "."
                (
                    if t_ref == "." && t_alt == "." { ".,.,.".to_string() } else { format!("{t_ref},{t_alt},.") },
                    if n_ref == "." && n_alt == "." { ".,.,.".to_string() } else { format!("{n_ref},{n_alt},.") },
                )
            } else {
                (
                    if t_ref == "." && t_alt == "." { ".,.".to_string() } else { format!("{t_ref},{t_alt}") },
                    if n_ref == "." && n_alt == "." { ".,.".to_string() } else { format!("{n_ref},{n_alt}") },
                )
            };
            (
                "GT:AD:DP",
                format!("{tumor_gt}:{t_ad}:{t_dep}"),
                format!("{normal_gt}:{n_ad}:{n_dep}"),
            )
        } else {
            ("GT", tumor_gt.to_string(), normal_gt.to_string())
        };

        write!(
            w,
            "{chrom}\t{vcf_pos}\t.\t{vcf_ref}\t{vcf_alt}\t.\t.\t.\t{fmt}"
        )?;
        for s in &samples {
            if s == tumor_id {
                write!(w, "\t{tumor_fmt}")?;
            } else if s == normal_id {
                write!(w, "\t{normal_fmt}")?;
            } else {
                write!(w, "\t./.")?;
            }
        }
        writeln!(w)?;
    }

    info!(
        "Wrote {} variants to {}",
        records.len(),
        args.output_vcf.display()
    );
    Ok(())
}

/// Convert MAF allele representation to VCF representation.
///
/// For indels, looks up the anchor base from the FASTA at position (start-1).
/// Falls back to "N" if the FASTA is unavailable or the position is out of range.
fn maf_alleles_to_vcf(
    fasta_path: &Path,
    chrom: &str,
    start: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> (u64, String, String) {
    if ref_allele == "-" {
        // Insertion: MAF Start_Position equals the anchor VCF POS.
        let anchor = crate::fasta::fasta_fetch_base(fasta_path, chrom, start)
            .map(|c| c.to_string())
            .unwrap_or_else(|| "N".to_string());
        (start, anchor.clone(), format!("{anchor}{alt_allele}"))
    } else if alt_allele == "-" {
        // Deletion: MAF Start_Position is the first deleted base; anchor is one before.
        let anchor = crate::fasta::fasta_fetch_base(fasta_path, chrom, start - 1)
            .map(|c| c.to_string())
            .unwrap_or_else(|| "N".to_string());
        (start - 1, format!("{anchor}{ref_allele}"), anchor)
    } else {
        (start, ref_allele.to_owned(), alt_allele.to_owned())
    }
}

fn write_vcf_header(
    w: &mut impl Write,
    genome: &str,
    samples: &[String],
    has_depth: bool,
) -> Result<()> {
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##reference={genome}")?;
    writeln!(
        w,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    if has_depth {
        writeln!(w, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">")?;
        writeln!(w, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">")?;
    }
    write!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for s in samples {
        write!(w, "\t{s}")?;
    }
    writeln!(w)?;
    Ok(())
}
