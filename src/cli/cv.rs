//! The `gecco cv` subcommand — cross-validation.

use std::path::PathBuf;

use clap::Args;

#[derive(Args)]
pub struct CvArgs {
    /// Gene table (TSV).
    #[arg(short, long)]
    pub genes: PathBuf,

    /// Domain annotation table(s) (TSV).
    #[arg(short, long, num_args = 1..)]
    pub features: Vec<PathBuf>,

    /// Cluster annotation table (TSV).
    #[arg(short, long)]
    pub clusters: PathBuf,

    /// Output file for CV results.
    #[arg(short, long, default_value = "cv.tsv")]
    pub output: PathBuf,

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

    /// Use Leave-One-Type-Out cross-validation instead of K-fold.
    #[arg(long)]
    pub loto: bool,

    /// Number of folds for K-fold cross-validation.
    #[arg(long, default_value = "10")]
    pub splits: usize,
}

impl CvArgs {
    pub fn execute(&self) -> anyhow::Result<()> {
        log::info!("Running cross-validation");
        let mode = if self.loto {
            "Leave-One-Type-Out".to_string()
        } else {
            format!("{}-fold", self.splits)
        };
        log::info!("Mode: {}", mode);

        // 1. Load gene + feature + cluster tables
        // 2. Filter domains
        // 3. Label genes from cluster annotations
        // 4. Split data (K-fold or LOTO)
        // 5. For each fold: train CRF, predict, evaluate
        // 6. Write CV results

        anyhow::bail!("cv subcommand not yet fully wired")
    }
}
