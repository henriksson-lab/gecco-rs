//! The `gecco train` subcommand — train a new CRF model.

use std::path::PathBuf;

use clap::Args;

#[derive(Args)]
pub struct TrainArgs {
    /// Gene table (TSV).
    #[arg(short, long)]
    pub genes: PathBuf,

    /// Domain annotation table(s) (TSV).
    #[arg(short, long, num_args = 1..)]
    pub features: Vec<PathBuf>,

    /// Cluster annotation table (TSV).
    #[arg(short, long)]
    pub clusters: PathBuf,

    /// Output directory for the trained model.
    #[arg(short, long, default_value = ".")]
    pub output_dir: PathBuf,

    /// E-value cutoff for protein domains.
    #[arg(short, long)]
    pub e_filter: Option<f64>,

    /// P-value cutoff for protein domains.
    #[arg(short, long, default_value = "1e-9")]
    pub p_filter: f64,

    /// Disable shuffling of training data.
    #[arg(long)]
    pub no_shuffle: bool,

    /// Random seed.
    #[arg(long, default_value = "42")]
    pub seed: u64,

    /// CRF sliding window length.
    #[arg(short = 'W', long, default_value = "5")]
    pub window_size: usize,

    /// CRF sliding window step.
    #[arg(long, default_value = "1")]
    pub window_step: usize,

    /// L1 regularization strength.
    #[arg(long, default_value = "0.15")]
    pub c1: f64,

    /// L2 regularization strength.
    #[arg(long, default_value = "0.15")]
    pub c2: f64,

    /// Feature extraction level: protein or domain.
    #[arg(long, default_value = "protein")]
    pub feature_type: String,

    /// Fraction of significant features to select (0.0-1.0).
    #[arg(long)]
    pub select: Option<f64>,
}

impl TrainArgs {
    pub fn execute(&self) -> anyhow::Result<()> {
        log::info!("Training CRF model");
        log::info!(
            "Input: genes={:?}, features={:?}, clusters={:?}",
            self.genes,
            self.features,
            self.clusters
        );
        log::info!("Output directory: {:?}", self.output_dir);

        // 1. Load gene + feature + cluster tables
        // 2. Filter domains
        // 3. Label genes from cluster annotations
        // 4. Fit CRF model
        // 5. Save model + domain list

        anyhow::bail!("train subcommand not yet fully wired")
    }
}
