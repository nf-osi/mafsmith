use crate::cli::Maf2vcfArgs;
use anyhow::{Context, Result};
use std::{
    collections::HashMap,
    fs,
    io::{BufRead, BufReader, Write},
};
use tracing::info;

pub async fn run(args: Maf2vcfArgs) -> Result<()> {
    let data_dir = crate::commands::fetch::data_dir(None)?;
    let genome_str = args.genome.ncbi_build();

    let _ref_fasta = args.ref_fasta.unwrap_or_else(|| {
        data_dir.join(genome_str).join("reference.fa")
    });

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
        // First non-comment line is the column header
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
    let tumor_idx = req("Tumor_Sample_Barcode")?;
    let normal_idx = req("Matched_Norm_Sample_Barcode")?;

    // Collect all unique sample pairs for per-pair VCFs
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

    write_vcf_header(&mut w, genome_str, &samples)?;

    for row in &records {
        let chrom = row.get(chrom_idx).map(|s| s.as_str()).unwrap_or(".");
        let start: u64 = row
            .get(start_idx)
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let ref_allele = row.get(ref_idx).map(|s| s.as_str()).unwrap_or(".");
        let alt_allele = row.get(alt_idx).map(|s| s.as_str()).unwrap_or(".");

        // MAF uses 1-based closed coords for SNPs; VCF is also 1-based.
        // For insertions (ref_allele == "-"), add padding base — requires FASTA.
        // For deletions (alt_allele == "-"), same.
        let (vcf_pos, vcf_ref, vcf_alt) =
            maf_alleles_to_vcf(start, ref_allele, alt_allele);

        writeln!(
            w,
            "{chrom}\t{vcf_pos}\t.\t{vcf_ref}\t{vcf_alt}\t.\t.\t.\tGT",
        )?;
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
/// MAF uses "-" for the absent allele in indels and doesn't include padding bases.
/// VCF requires a shared prefix base for indels.
/// Without a FASTA, we use a placeholder "N" padding base when needed.
fn maf_alleles_to_vcf(start: u64, ref_allele: &str, alt_allele: &str) -> (u64, String, String) {
    if ref_allele == "-" {
        // Insertion in MAF → VCF INS with "N" padding
        (start - 1, String::from("N"), format!("N{alt_allele}"))
    } else if alt_allele == "-" {
        // Deletion in MAF → VCF DEL with "N" padding
        (start - 1, format!("N{ref_allele}"), String::from("N"))
    } else {
        (start, ref_allele.to_owned(), alt_allele.to_owned())
    }
}

fn write_vcf_header(w: &mut impl Write, genome: &str, samples: &[String]) -> Result<()> {
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##reference={genome}")?;
    writeln!(w, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    write!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for s in samples {
        write!(w, "\t{s}")?;
    }
    writeln!(w)?;
    Ok(())
}
