use anyhow::Result;
use clap::Parser;
use mafsmith::cli::{Cli, Command};
use mafsmith::commands;

#[global_allocator]
static ALLOC: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[tokio::main]
async fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::from_default_env()
                .add_directive("mafsmith=info".parse().unwrap()),
        )
        .with_writer(std::io::stderr)
        .init();

    let cli = Cli::parse();
    match cli.command {
        Command::Vcf2maf(args) => commands::vcf2maf::run(args).await,
        Command::Maf2vcf(args) => commands::maf2vcf::run(args).await,
        Command::Maf2maf(args) => commands::maf2maf::run(args).await,
        Command::Vcf2vcf(args) => commands::vcf2vcf::run(args).await,
        Command::Fetch(args) => commands::fetch::run(args).await,
    }
}
