use clap::{Args, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "mafsmith", version, about = "Fast, self-contained VCF↔MAF converter")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Subcommand)]
pub enum Command {
    /// Convert VCF to MAF with fastVEP annotation
    Vcf2maf(Vcf2mafArgs),
    /// Convert MAF to VCF
    Maf2vcf(Maf2vcfArgs),
    /// Reannotate an existing MAF (round-trip through fastVEP)
    Maf2maf(Maf2mafArgs),
    /// Normalize a VCF file
    Vcf2vcf(Vcf2vcfArgs),
    /// Download reference data and install fastVEP
    Fetch(FetchArgs),
}

#[derive(ValueEnum, Clone, Debug)]
pub enum Genome {
    Grch38,
    Grch37,
    Grcm39,
}

impl Genome {
    pub fn ncbi_build(&self) -> &'static str {
        match self {
            Genome::Grch38 => "GRCh38",
            Genome::Grch37 => "GRCh37",
            Genome::Grcm39 => "GRCm39",
        }
    }
}

#[derive(Args, Debug)]
pub struct Vcf2mafArgs {
    /// Input VCF file
    #[arg(short = 'i', long, value_name = "VCF")]
    pub input_vcf: PathBuf,

    /// Output MAF file
    #[arg(short = 'o', long, value_name = "MAF")]
    pub output_maf: PathBuf,

    /// Reference genome assembly
    #[arg(long, value_enum, default_value = "grch38")]
    pub genome: Genome,

    /// Tumor sample barcode (defaults to VCF sample name)
    #[arg(long)]
    pub tumor_id: Option<String>,

    /// Matched normal sample barcode
    #[arg(long)]
    pub normal_id: Option<String>,

    /// Tumor sample column name in VCF
    #[arg(long)]
    pub vcf_tumor_id: Option<String>,

    /// Normal sample column name in VCF
    #[arg(long)]
    pub vcf_normal_id: Option<String>,

    /// Custom Ensembl transcript IDs to prefer (one per line)
    #[arg(long)]
    pub custom_enst: Option<PathBuf>,

    /// Override path to fastVEP binary
    #[arg(long)]
    pub fastvep_path: Option<PathBuf>,

    /// Override reference FASTA (default: from mafsmith data dir)
    #[arg(long)]
    pub ref_fasta: Option<PathBuf>,

    /// Override gene GFF3 file (default: from mafsmith data dir)
    #[arg(long)]
    pub gff3: Option<PathBuf>,

    /// Sequencing center for MAF header
    #[arg(long, default_value = "")]
    pub maf_center: String,

    /// Minimum alt allele frequency to include somatic variants
    #[arg(long, default_value_t = 0.0)]
    pub min_hom_vaf: f64,

    /// Retain extra VEP annotation fields in output (comma-separated CSQ field names)
    #[arg(long, value_delimiter = ',')]
    pub retain_ann: Vec<String>,

    /// Skip fastVEP annotation (input VCF must already contain CSQ annotations)
    #[arg(long, default_value_t = false)]
    pub skip_annotation: bool,
}

#[derive(Args, Debug)]
pub struct Maf2vcfArgs {
    /// Input MAF file
    #[arg(short = 'i', long, value_name = "MAF")]
    pub input_maf: PathBuf,

    /// Output VCF file
    #[arg(short = 'o', long, value_name = "VCF")]
    pub output_vcf: PathBuf,

    /// Reference genome assembly
    #[arg(long, value_enum, default_value = "grch38")]
    pub genome: Genome,

    /// Reference FASTA (for allele validation)
    #[arg(long)]
    pub ref_fasta: Option<PathBuf>,

    /// Per-tumor-normal-pair VCFs output directory
    #[arg(long)]
    pub per_tn_vcfs: bool,
}

#[derive(Args, Debug)]
pub struct Maf2mafArgs {
    /// Input MAF file
    #[arg(short = 'i', long, value_name = "MAF")]
    pub input_maf: PathBuf,

    /// Output MAF file
    #[arg(short = 'o', long, value_name = "MAF")]
    pub output_maf: PathBuf,

    /// Reference genome assembly
    #[arg(long, value_enum, default_value = "grch38")]
    pub genome: Genome,

    /// Override fastVEP path
    #[arg(long)]
    pub fastvep_path: Option<PathBuf>,

    /// Custom transcript list
    #[arg(long)]
    pub custom_enst: Option<PathBuf>,
}

#[derive(Args, Debug)]
pub struct Vcf2vcfArgs {
    /// Input VCF file
    #[arg(short = 'i', long, value_name = "VCF")]
    pub input_vcf: PathBuf,

    /// Output VCF file
    #[arg(short = 'o', long, value_name = "VCF")]
    pub output_vcf: PathBuf,

    /// Reference genome assembly
    #[arg(long, value_enum, default_value = "grch38")]
    pub genome: Genome,

    /// Reference FASTA for normalization
    #[arg(long)]
    pub ref_fasta: Option<PathBuf>,

    /// Tumor sample column name in VCF
    #[arg(long)]
    pub vcf_tumor_id: Option<String>,

    /// Normal sample column name in VCF
    #[arg(long)]
    pub vcf_normal_id: Option<String>,
}

#[derive(Args, Debug)]
pub struct FetchArgs {
    /// Genome(s) to download reference data for
    #[arg(long, value_enum, num_args = 1.., value_delimiter = ',', default_values = ["grch38"])]
    pub genome: Vec<Genome>,

    /// Data directory (default: ~/.mafsmith)
    #[arg(long)]
    pub data_dir: Option<PathBuf>,

    /// Use pre-existing GFF3 file instead of downloading
    #[arg(long)]
    pub gff3: Option<PathBuf>,

    /// Use pre-existing reference FASTA instead of downloading
    #[arg(long)]
    pub ref_fasta: Option<PathBuf>,

    /// Skip installing fastVEP
    #[arg(long)]
    pub skip_fastvep: bool,

    /// Ensembl release to use
    #[arg(long, default_value_t = 113)]
    pub ensembl_release: u32,
}
