use crate::{
    annotation::{
        consequence::{end_position, so_to_variant_classification, variant_type},
        csq::{shorten_hgvsp, CsqFormat},
        depth::extract_depth,
        transcript::select_transcript,
    },
    cli::Vcf2mafArgs,
    maf::{record::MafRecord, writer::MafWriter},
    vcf::VcfReader,
};
use anyhow::{bail, Context, Result};
use std::{
    collections::HashSet,
    fs,
    io::{BufRead, BufReader, BufWriter},
    path::{Path, PathBuf},
    process::{Command, Stdio},
};
use tempfile::NamedTempFile;
use tracing::{info, warn};

pub async fn run(args: Vcf2mafArgs) -> Result<()> {
    let data_dir = crate::commands::fetch::data_dir(None)?;
    let genome_str = args.genome.ncbi_build();

    let ref_fasta = resolve_ref_fasta(&args.ref_fasta, &data_dir, genome_str)?;
    let gff3 = resolve_gff3(&args.gff3, &data_dir, genome_str)?;
    let fastvep = resolve_fastvep(&args.fastvep_path, &data_dir)?;

    let custom_enst = load_custom_enst(args.custom_enst.as_deref())?;

    // Run fastVEP → annotated VCF in a temp file
    info!("Annotating variants with fastVEP...");
    let annotated_vcf = run_fastvep(&fastvep, &args.input_vcf, &ref_fasta, &gff3)?;

    // Parse and convert
    info!("Converting annotated VCF to MAF...");
    let out_file = fs::File::create(&args.output_maf)
        .with_context(|| format!("Cannot create output MAF: {}", args.output_maf.display()))?;
    let mut writer = MafWriter::new(BufWriter::new(out_file), args.retain_ann.clone())?;

    let in_file = fs::File::open(annotated_vcf.path())
        .context("Cannot open annotated VCF temp file")?;
    let mut vcf = VcfReader::new(BufReader::new(in_file));

    // Pull CSQ FORMAT from header before iterating records
    let mut csq_format: Option<CsqFormat> = None;

    // We need to read header first by initializing the reader
    // The reader inits lazily on first next_record call, so prime it:
    let first_rec = vcf.next_record()?;

    // Parse CSQ format from header
    for hline in &vcf.header_lines {
        if hline.contains("ID=CSQ") {
            csq_format = CsqFormat::from_header_description(hline).ok();
            break;
        }
    }
    let csq_format = csq_format.context(
        "fastVEP output did not contain a CSQ INFO header. \
         Was fastVEP run with annotation enabled?",
    )?;

    let sample_names = vcf.sample_names.clone();
    let tumor_col = resolve_sample_col(
        args.vcf_tumor_id.as_deref(),
        &sample_names,
        "tumor",
        true,
    )?;
    let normal_col = resolve_sample_col(
        args.vcf_normal_id.as_deref(),
        &sample_names,
        "normal",
        false,
    )?;

    let tumor_barcode = args.tumor_id.as_deref().unwrap_or(tumor_col.as_deref().unwrap_or("TUMOR"));
    let normal_barcode = args.normal_id.as_deref().unwrap_or(normal_col.as_deref().unwrap_or("NORMAL"));

    let mut process_record = |rec: &crate::vcf::VcfRecord| -> Result<()> {
        let csq_val = match rec.csq_value() {
            Some(v) => v.to_owned(),
            None => return Ok(()), // unannotated variant — skip
        };

        let entries = csq_format.parse_all(&csq_val, &args.retain_ann);
        let transcript = match select_transcript(&entries, custom_enst.as_ref()) {
            Some(t) => t,
            None => return Ok(()),
        };

        // Allele depths
        let fk: Vec<&str> = rec.format_keys.iter().map(|s| s.as_str()).collect();
        let tumor_depth = tumor_col
            .as_deref()
            .and_then(|col| rec.sample_values(col))
            .map(|vals| {
                let fkr: Vec<&str> = fk.clone();
                extract_depth(&fkr, &vals, &rec.ref_allele, &rec.alt_allele)
            });
        let normal_depth = normal_col
            .as_deref()
            .and_then(|col| rec.sample_values(col))
            .map(|vals| extract_depth(&fk, &vals, &rec.ref_allele, &rec.alt_allele));

        let csq_refs: Vec<&str> = transcript.consequences.iter().map(|s| s.as_str()).collect();
        let var_class =
            so_to_variant_classification(&csq_refs, &rec.ref_allele, &rec.alt_allele);
        let var_type = variant_type(&rec.ref_allele, &rec.alt_allele);
        let end_pos = end_position(rec.pos, &rec.ref_allele, &rec.alt_allele);

        // dbSNP: pull from Existing_variation or ID column
        let dbsnp = if rec.id != "." {
            rec.id.clone()
        } else if transcript.existing_variation.starts_with("rs") {
            transcript.existing_variation.clone()
        } else {
            String::new()
        };

        // all_effects: semicolon-separated list for every annotated transcript
        let all_effects = entries
            .iter()
            .filter(|e| !e.feature.is_empty())
            .map(|e| {
                let csq_r: Vec<&str> = e.consequences.iter().map(|s| s.as_str()).collect();
                format!(
                    "{},{},{},{}",
                    e.symbol,
                    e.consequences.join("&"),
                    so_to_variant_classification(&csq_r, &rec.ref_allele, &rec.alt_allele),
                    e.hgvsp
                )
            })
            .collect::<Vec<_>>()
            .join(";");

        let hgvsp_short = shorten_hgvsp(&transcript.hgvsp);

        let fmt_depth = |d: Option<u32>| d.map(|v| v.to_string()).unwrap_or_default();

        let mut extra = indexmap::IndexMap::new();
        for field in &args.retain_ann {
            let val = transcript.extra.get(field).cloned().unwrap_or_default();
            extra.insert(field.clone(), val);
        }

        let record = MafRecord {
            hugo_symbol: transcript.symbol.clone(),
            entrez_gene_id: String::from("0"),
            center: args.maf_center.clone(),
            ncbi_build: genome_str.to_owned(),
            chromosome: rec.chrom.clone(),
            start_position: rec.pos,
            end_position: end_pos,
            strand: String::from("+"),
            variant_classification: var_class.to_owned(),
            variant_type: var_type.to_owned(),
            reference_allele: rec.ref_allele.clone(),
            tumor_seq_allele1: rec.ref_allele.clone(),
            tumor_seq_allele2: rec.alt_allele.clone(),
            dbsnp_rs: dbsnp,
            dbsnp_val_status: String::new(),
            tumor_sample_barcode: tumor_barcode.to_owned(),
            matched_norm_sample_barcode: normal_barcode.to_owned(),
            match_norm_seq_allele1: rec.ref_allele.clone(),
            match_norm_seq_allele2: rec.ref_allele.clone(),
            mutation_status: String::from("Somatic"),
            hgvsc: transcript.hgvsc.clone(),
            hgvsp: transcript.hgvsp.clone(),
            hgvsp_short,
            transcript_id: transcript.feature.clone(),
            exon_number: transcript.exon.clone(),
            t_depth: fmt_depth(tumor_depth.as_ref().and_then(|d| d.depth())),
            t_ref_count: fmt_depth(tumor_depth.as_ref().and_then(|d| d.ref_count)),
            t_alt_count: fmt_depth(tumor_depth.as_ref().and_then(|d| d.alt_count)),
            n_depth: fmt_depth(normal_depth.as_ref().and_then(|d| d.depth())),
            n_ref_count: fmt_depth(normal_depth.as_ref().and_then(|d| d.ref_count)),
            n_alt_count: fmt_depth(normal_depth.as_ref().and_then(|d| d.alt_count)),
            all_effects,
            vep_allele: transcript.allele.clone(),
            vep_gene: transcript.gene.clone(),
            vep_feature: transcript.feature.clone(),
            vep_biotype: transcript.biotype.clone(),
            vep_canonical: if transcript.canonical {
                String::from("YES")
            } else {
                String::new()
            },
            vep_sift: transcript.sift.clone(),
            vep_polyphen: transcript.polyphen.clone(),
            extra,
            ..Default::default()
        };

        writer.write_record(&record)?;
        Ok(())
    };

    if let Some(rec) = first_rec {
        process_record(&rec)?;
    }
    while let Some(rec) = vcf.next_record()? {
        process_record(&rec)?;
    }

    writer.flush()?;
    info!("MAF written to {}", args.output_maf.display());
    Ok(())
}

fn run_fastvep(
    fastvep: &Path,
    input_vcf: &Path,
    ref_fasta: &Path,
    gff3: &Path,
) -> Result<NamedTempFile> {
    let tmp = NamedTempFile::new().context("Cannot create temp file for fastVEP output")?;
    let status = Command::new(fastvep)
        .args([
            "annotate",
            "-i",
            input_vcf.to_str().unwrap(),
            "-o",
            tmp.path().to_str().unwrap(),
            "--fasta",
            ref_fasta.to_str().unwrap(),
            "--gff3",
            gff3.to_str().unwrap(),
            "--hgvs",
            "--output-format",
            "vcf",
        ])
        .status()
        .context("Failed to run fastVEP. Run `mafsmith fetch` to install it.")?;

    if !status.success() {
        bail!("fastVEP exited with status {}", status);
    }
    Ok(tmp)
}

fn resolve_ref_fasta(
    cli: &Option<PathBuf>,
    data_dir: &Path,
    genome: &str,
) -> Result<PathBuf> {
    if let Some(p) = cli {
        return Ok(p.clone());
    }
    let p = data_dir.join(genome).join("reference.fa");
    if !p.exists() {
        bail!(
            "Reference FASTA not found at {}. Run `mafsmith fetch --genome {}` first.",
            p.display(),
            genome.to_lowercase()
        );
    }
    Ok(p)
}

fn resolve_gff3(cli: &Option<PathBuf>, data_dir: &Path, genome: &str) -> Result<PathBuf> {
    if let Some(p) = cli {
        return Ok(p.clone());
    }
    // Accept both compressed and uncompressed
    for name in &["genes.gff3.gz", "genes.gff3"] {
        let p = data_dir.join(genome).join(name);
        if p.exists() {
            return Ok(p);
        }
    }
    bail!(
        "GFF3 not found under {}. Run `mafsmith fetch --genome {}` first.",
        data_dir.join(genome).display(),
        genome.to_lowercase()
    );
}

fn resolve_fastvep(cli: &Option<PathBuf>, data_dir: &Path) -> Result<PathBuf> {
    if let Some(p) = cli {
        return Ok(p.clone());
    }
    // Check our managed bin dir first
    let managed = data_dir.join("bin").join("fastvep");
    if managed.exists() {
        return Ok(managed);
    }
    // Fall back to $PATH
    if let Ok(path) = which_fastvep() {
        return Ok(path);
    }
    bail!("fastVEP binary not found. Run `mafsmith fetch` to install it.");
}

fn which_fastvep() -> Result<PathBuf> {
    let out = Command::new("which").arg("fastvep").output()?;
    if out.status.success() {
        let path = String::from_utf8_lossy(&out.stdout).trim().to_owned();
        Ok(PathBuf::from(path))
    } else {
        bail!("fastvep not in PATH")
    }
}

fn load_custom_enst(path: Option<&Path>) -> Result<Option<HashSet<String>>> {
    let Some(p) = path else { return Ok(None) };
    let content = fs::read_to_string(p)
        .with_context(|| format!("Cannot read custom ENST file: {}", p.display()))?;
    let set: HashSet<String> = content
        .lines()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| l.split_whitespace().next().unwrap_or(l).to_owned())
        .collect();
    Ok(Some(set))
}

/// Determine which VCF sample column to use for tumor or normal.
/// Returns None if no samples present (single-sample VCF without tumor/normal designation).
fn resolve_sample_col(
    cli: Option<&str>,
    sample_names: &[String],
    role: &str,
    is_tumor: bool,
) -> Result<Option<String>> {
    if let Some(name) = cli {
        if sample_names.contains(&name.to_owned()) {
            return Ok(Some(name.to_owned()));
        }
        warn!(
            "Specified {} sample '{}' not found in VCF samples {:?}",
            role, name, sample_names
        );
        return Ok(None);
    }
    if sample_names.is_empty() {
        return Ok(None);
    }
    // Heuristic: if two samples, first = tumor, second = normal
    if sample_names.len() >= 2 {
        return Ok(Some(sample_names[if is_tumor { 0 } else { 1 }].clone()));
    }
    if is_tumor {
        Ok(Some(sample_names[0].clone()))
    } else {
        Ok(None)
    }
}
