use crate::cli::Maf2mafArgs;
use anyhow::Result;
use tempfile::NamedTempFile;
use tracing::info;

/// Reannotate a MAF by round-tripping through fastVEP.
///
/// Pipeline: input MAF → VCF (maf2vcf) → fastVEP annotation → MAF (vcf2maf).
pub async fn run(args: Maf2mafArgs) -> Result<()> {
    info!("maf2maf: converting input MAF to intermediate VCF...");
    let intermediate_vcf = NamedTempFile::new()?;
    let _intermediate_maf = NamedTempFile::new()?;

    crate::commands::maf2vcf::run(crate::cli::Maf2vcfArgs {
        input_maf: args.input_maf,
        output_vcf: intermediate_vcf.path().to_path_buf(),
        genome: args.genome.clone(),
        ref_fasta: None,
        per_tn_vcfs: false,
    })
    .await?;

    info!("maf2maf: annotating intermediate VCF with fastVEP...");
    crate::commands::vcf2maf::run(crate::cli::Vcf2mafArgs {
        input_vcf: intermediate_vcf.path().to_path_buf(),
        output_maf: args.output_maf,
        genome: args.genome,
        tumor_id: None,
        normal_id: None,
        vcf_tumor_id: None,
        vcf_normal_id: None,
        custom_enst: args.custom_enst,
        annotator: crate::cli::Annotator::Fastvep,
        fastvep_path: args.fastvep_path,
        vep_path: None,
        vep_data: None,
        vep_forks: 0,
        ref_fasta: None,
        gff3: None,
        maf_center: String::new(),
        min_hom_vaf: 0.0,
        retain_ann: vec![],
        skip_annotation: false,
        strict: false,
    })
    .await?;

    Ok(())
}
