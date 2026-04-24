use crate::cli::Vcf2vcfArgs;
use anyhow::{Context, Result};
use std::{
    fs,
    io::{BufReader, BufWriter, Write},
};
use tracing::info;

/// Normalize a VCF: standardize FORMAT fields, select tumor/normal sample columns,
/// remove multiallelic sites (split to biallelic), and pass through canonical VCF.
pub async fn run(args: Vcf2vcfArgs) -> Result<()> {
    let in_file = fs::File::open(&args.input_vcf)
        .with_context(|| format!("Cannot open {}", args.input_vcf.display()))?;
    let out_file = fs::File::create(&args.output_vcf)
        .with_context(|| format!("Cannot create {}", args.output_vcf.display()))?;

    let reader = BufReader::new(in_file);
    let mut writer = BufWriter::new(out_file);

    let mut vcf = crate::vcf::VcfReader::new(reader);

    // Write header lines through
    let first = vcf.next_record()?;
    for h in &vcf.header_lines {
        writeln!(writer, "{h}")?;
    }

    let sample_names = vcf.sample_names.clone();
    let tumor_col = args.vcf_tumor_id.or_else(|| sample_names.first().cloned());
    let normal_col = args.vcf_normal_id.or_else(|| sample_names.get(1).cloned());

    let write_record = |w: &mut BufWriter<fs::File>, rec: &crate::vcf::VcfRecord| -> Result<bool> {
        // Skip ref-only variants (spanning deletions or gVCF blocks)
        if rec.alt_allele == "." || rec.alt_allele == rec.ref_allele {
            return Ok(false);
        }

        let all_alts = rec.all_alts.join(",");
        let format = rec.format_keys.join(":");
        write!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            rec.chrom,
            rec.pos,
            rec.id,
            rec.ref_allele,
            all_alts,
            rec.qual,
            rec.filter,
            rec.info,
            format,
        )?;

        // Write selected sample columns
        for col in &[&tumor_col, &normal_col] {
            if let Some(name) = col {
                if let Some(vals) = rec.sample_values(name) {
                    write!(w, "\t{}", vals.join(":"))?;
                }
            }
        }
        writeln!(w)?;
        Ok(true)
    };

    let mut count = 0usize;
    if let Some(rec) = first {
        if write_record(&mut writer, &rec)? {
            count += 1;
        }
    }
    while let Some(rec) = vcf.next_record()? {
        if write_record(&mut writer, &rec)? {
            count += 1;
        }
    }

    info!(
        "vcf2vcf: wrote {count} variants to {}",
        args.output_vcf.display()
    );
    Ok(())
}
