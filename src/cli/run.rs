//! The `gecco run` subcommand — full pipeline.

use std::path::PathBuf;

use clap::Args;

#[derive(Args)]
pub struct RunArgs {
    /// Path to the input genome (FASTA or GenBank).
    #[arg(short, long)]
    pub genome: PathBuf,

    /// Output directory.
    #[arg(short, long)]
    pub output: PathBuf,

    /// Number of threads (0 = auto-detect).
    #[arg(short, long, default_value = "0")]
    pub jobs: usize,

    /// Minimum number of consecutive CDS in a cluster.
    #[arg(long, default_value = "3")]
    pub cds: usize,

    /// Probability threshold for gene cluster membership.
    #[arg(long, default_value = "0.8")]
    pub threshold: f64,
}

impl RunArgs {
    pub fn execute(&self) -> anyhow::Result<()> {
        log::info!("Running GECCO pipeline on {:?}", self.genome);
        log::info!("Output directory: {:?}", self.output);

        // TODO: Wire up the full pipeline:
        // 1. Load sequences (FASTA/GenBank)
        // 2. Find genes (orphos-core / CDSFinder)
        // 3. Annotate domains (hmmer-pure-rs)
        // 4. Predict probabilities (CRF)
        // 5. Refine clusters
        // 6. Predict types (RandomForest)
        // 7. Write output (TSV tables + GenBank files)

        anyhow::bail!("Pipeline not yet fully wired — individual modules are implemented")
    }
}
