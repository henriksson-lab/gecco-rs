//! The `gecco convert` subcommand — format conversion.

use std::path::PathBuf;

use clap::{Args, Subcommand, ValueEnum};

#[derive(Args)]
pub struct ConvertArgs {
    #[command(subcommand)]
    pub command: ConvertCommand,
}

#[derive(Subcommand)]
pub enum ConvertCommand {
    /// Convert GenBank cluster files to other formats.
    Gbk(ConvertGbkArgs),
    /// Convert cluster TSV tables to other formats.
    Clusters(ConvertClustersArgs),
}

/// Output format for GenBank conversion.
#[derive(Debug, Clone, ValueEnum)]
pub enum GbkFormat {
    /// BiG-SLiCE feature table format.
    Bigslice,
    /// Nucleotide FASTA.
    Fna,
    /// Protein FASTA.
    Faa,
}

/// Output format for cluster table conversion.
#[derive(Debug, Clone, ValueEnum)]
pub enum ClusterFormat {
    /// GFF3 annotation format.
    Gff,
}

#[derive(Args)]
pub struct ConvertGbkArgs {
    /// Input directory containing GenBank files.
    #[arg(short, long)]
    pub input_dir: PathBuf,

    /// Output directory (defaults to input directory).
    #[arg(short, long)]
    pub output_dir: Option<PathBuf>,

    /// Output format.
    #[arg(short, long)]
    pub format: GbkFormat,
}

#[derive(Args)]
pub struct ConvertClustersArgs {
    /// Input directory containing cluster TSV files.
    #[arg(short, long)]
    pub input_dir: PathBuf,

    /// Output directory (defaults to input directory).
    #[arg(short, long)]
    pub output_dir: Option<PathBuf>,

    /// Output format.
    #[arg(short, long)]
    pub format: ClusterFormat,
}

impl ConvertArgs {
    pub fn execute(&self) -> anyhow::Result<()> {
        match &self.command {
            ConvertCommand::Gbk(args) => {
                let out = args
                    .output_dir
                    .as_ref()
                    .unwrap_or(&args.input_dir);
                log::info!(
                    "Converting GenBank files from {:?} to {:?} ({:?})",
                    args.input_dir,
                    out,
                    args.format
                );
                anyhow::bail!("convert gbk not yet implemented")
            }
            ConvertCommand::Clusters(args) => {
                let out = args
                    .output_dir
                    .as_ref()
                    .unwrap_or(&args.input_dir);
                log::info!(
                    "Converting cluster tables from {:?} to {:?} ({:?})",
                    args.input_dir,
                    out,
                    args.format
                );
                anyhow::bail!("convert clusters not yet implemented")
            }
        }
    }
}
