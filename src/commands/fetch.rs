use crate::cli::{FetchArgs, Genome};
use anyhow::{bail, Context, Result};
use futures_util::StreamExt;
use indicatif::{ProgressBar, ProgressStyle};
use std::{
    fs,
    io::Write,
    path::{Path, PathBuf},
    process::Command,
};
use tracing::info;

/// Default mafsmith data directory.
pub fn data_dir(override_path: Option<PathBuf>) -> Result<PathBuf> {
    if let Some(p) = override_path {
        return Ok(p);
    }
    let home = dirs_or_home()?;
    Ok(home.join(".mafsmith"))
}

fn dirs_or_home() -> Result<PathBuf> {
    // Try $HOME first
    if let Ok(h) = std::env::var("HOME") {
        return Ok(PathBuf::from(h));
    }
    bail!("Cannot determine home directory. Set $HOME or pass --data-dir.");
}

pub async fn run(args: FetchArgs) -> Result<()> {
    let data_dir = data_dir(args.data_dir.clone())?;
    fs::create_dir_all(&data_dir).context("Cannot create mafsmith data directory")?;

    if !args.skip_fastvep {
        install_fastvep(&data_dir).await?;
    }

    for genome in &args.genome {
        let genome_str = genome.ncbi_build();
        let genome_dir = data_dir.join(genome_str);
        fs::create_dir_all(&genome_dir)?;

        if let (Some(gff3_src), Some(fasta_src)) = (&args.gff3, &args.ref_fasta) {
            // User-provided files: symlink or copy
            link_or_copy(gff3_src, &genome_dir.join("genes.gff3.gz"))?;
            link_or_copy(fasta_src, &genome_dir.join("reference.fa"))?;
            info!("Linked user-provided reference files for {}", genome_str);
        } else {
            download_reference_data(genome, &genome_dir, args.ensembl_release).await?;
        }
    }

    info!("mafsmith fetch complete. Data at {}", data_dir.display());
    Ok(())
}

async fn install_fastvep(data_dir: &Path) -> Result<()> {
    let bin_dir = data_dir.join("bin");
    fs::create_dir_all(&bin_dir)?;
    let fastvep_bin = bin_dir.join("fastvep");

    if fastvep_bin.exists() {
        info!("fastVEP already installed at {}", fastvep_bin.display());
        return Ok(());
    }

    info!("Installing fastVEP from source (requires cargo in PATH)...");
    let status = Command::new("cargo")
        .args([
            "install",
            "--git",
            "https://github.com/Huang-lab/fastVEP",
            "--branch",
            "master",
            "fastvep-cli",
            "--root",
            data_dir.to_str().unwrap(),
        ])
        .status()
        .context("Failed to run cargo. Is Rust installed? See https://rustup.rs")?;

    if !status.success() {
        bail!("cargo install fastvep-cli failed with status {}", status);
    }
    info!("fastVEP installed to {}", fastvep_bin.display());
    Ok(())
}

async fn download_reference_data(genome: &Genome, genome_dir: &Path, release: u32) -> Result<()> {
    let (gff3_url, fasta_url) = ensembl_urls(genome, release);
    let genome_str = genome.ncbi_build();

    let gff3_dest = genome_dir.join("genes.gff3.gz");
    if gff3_dest.exists() {
        info!("{} GFF3 already present, skipping download", genome_str);
    } else {
        info!("Downloading {} gene annotations...", genome_str);
        download_file(&gff3_url, &gff3_dest).await?;
    }

    let fasta_dest = genome_dir.join("reference.fa.gz");
    if fasta_dest.exists() {
        info!("{} FASTA already present, skipping download", genome_str);
    } else {
        info!("Downloading {} reference genome (~1 GB)...", genome_str);
        download_file(&fasta_url, &fasta_dest).await?;
    }

    // Decompress FASTA (fastVEP needs uncompressed or bgzip with .fai)
    let fasta_uncompressed = genome_dir.join("reference.fa");
    if !fasta_uncompressed.exists() {
        info!("Decompressing FASTA...");
        decompress_gz(&fasta_dest, &fasta_uncompressed)?;
    }

    info!("{} reference data ready.", genome_str);
    Ok(())
}

fn ensembl_urls(genome: &Genome, release: u32) -> (String, String) {
    match genome {
        Genome::Grch38 => (
            format!(
                "https://ftp.ensembl.org/pub/release-{release}/gff3/homo_sapiens/\
                 Homo_sapiens.GRCh38.{release}.chr.gff3.gz"
            ),
            format!(
                "https://ftp.ensembl.org/pub/release-{release}/fasta/homo_sapiens/dna/\
                 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
            ),
        ),
        Genome::Grch37 => (
            // GRCh37 is archived on Ensembl; use the latest grch37 release path.
            // GFF3 filename keeps the 87 suffix (Ensembl annotation version); FASTA does not.
            "https://ftp.ensembl.org/pub/grch37/release-115/gff3/homo_sapiens/\
             Homo_sapiens.GRCh37.87.chr.gff3.gz"
                .to_owned(),
            "https://ftp.ensembl.org/pub/grch37/release-115/fasta/homo_sapiens/dna/\
             Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
                .to_owned(),
        ),
        Genome::Grcm39 => (
            format!(
                "https://ftp.ensembl.org/pub/release-{release}/gff3/mus_musculus/\
                 Mus_musculus.GRCm39.{release}.chr.gff3.gz"
            ),
            format!(
                "https://ftp.ensembl.org/pub/release-{release}/fasta/mus_musculus/dna/\
                 Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
            ),
        ),
    }
}

async fn download_file(url: &str, dest: &Path) -> Result<()> {
    let client = reqwest::Client::new();
    let response = client
        .get(url)
        .send()
        .await
        .with_context(|| format!("HTTP GET failed: {url}"))?
        .error_for_status()
        .with_context(|| format!("HTTP error for: {url}"))?;

    let total = response.content_length();
    let pb = ProgressBar::new(total.unwrap_or(0));
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} ({eta})",
        )
        .unwrap()
        .progress_chars("#>-"),
    );

    let tmp_dest = dest.with_extension("tmp");
    let mut file = fs::File::create(&tmp_dest)
        .with_context(|| format!("Cannot create {}", tmp_dest.display()))?;

    let mut stream = response.bytes_stream();
    while let Some(chunk) = stream.next().await {
        let chunk = chunk.context("Stream error during download")?;
        file.write_all(&chunk)?;
        pb.inc(chunk.len() as u64);
    }
    pb.finish_with_message("done");
    fs::rename(&tmp_dest, dest)?;
    Ok(())
}

fn decompress_gz(src: &Path, dest: &Path) -> Result<()> {
    use flate2::read::MultiGzDecoder;
    use std::io::BufReader;

    let file = fs::File::open(src).with_context(|| format!("Cannot open {}", src.display()))?;
    // MultiGzDecoder reads all concatenated gzip members — required for genome FASTAs which
    // are typically multi-member gzip (one per chromosome, concatenated).
    let mut decoder = MultiGzDecoder::new(BufReader::new(file));
    let mut out =
        fs::File::create(dest).with_context(|| format!("Cannot create {}", dest.display()))?;
    std::io::copy(&mut decoder, &mut out)?;
    Ok(())
}

fn link_or_copy(src: &Path, dest: &Path) -> Result<()> {
    if dest.exists() {
        return Ok(());
    }
    // Try symlink first; fall back to copy
    #[cfg(unix)]
    {
        let abs_src = src
            .canonicalize()
            .with_context(|| format!("Cannot resolve {}", src.display()))?;
        std::os::unix::fs::symlink(&abs_src, dest).with_context(|| {
            format!("Cannot symlink {} → {}", abs_src.display(), dest.display())
        })?;
    }
    #[cfg(not(unix))]
    {
        fs::copy(src, dest)?;
    }
    Ok(())
}
