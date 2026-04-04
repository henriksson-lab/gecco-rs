//! The `gecco predict` subcommand — predict from pre-annotated features.

use std::path::PathBuf;

use clap::Args;

#[derive(Args)]
pub struct PredictArgs {
    /// Path to the input genome (FASTA or GenBank).
    #[arg(long)]
    pub genome: PathBuf,

    /// Gene table (TSV).
    #[arg(short, long)]
    pub genes: PathBuf,

    /// Domain annotation table(s) (TSV).
    #[arg(short, long, num_args = 1..)]
    pub features: Vec<PathBuf>,

    /// Output directory.
    #[arg(short, long, default_value = ".")]
    pub output_dir: PathBuf,

    /// Number of threads (0 = auto-detect).
    #[arg(short, long, default_value = "0")]
    pub jobs: usize,

    /// E-value cutoff for protein domains.
    #[arg(short, long)]
    pub e_filter: Option<f64>,

    /// P-value cutoff for protein domains.
    #[arg(short, long, default_value = "1e-9")]
    pub p_filter: f64,

    /// Alternative CRF model directory.
    #[arg(long)]
    pub model: Option<PathBuf>,

    /// Disable padding of short gene sequences.
    #[arg(long)]
    pub no_pad: bool,

    /// Minimum number of consecutive CDS in a cluster.
    #[arg(short, long, default_value = "3")]
    pub cds: usize,

    /// Probability threshold for cluster membership.
    #[arg(short = 'm', long, default_value = "0.8")]
    pub threshold: f64,

    /// Minimum genes separating a cluster from sequence edge.
    #[arg(short = 'E', long, default_value = "0")]
    pub edge_distance: usize,

    /// Disable trimming of genes without domain annotations.
    #[arg(long)]
    pub no_trim: bool,

    /// Write empty TSV files when no results are found.
    #[arg(long)]
    pub force_tsv: bool,

    /// Output a single merged GenBank file.
    #[arg(long)]
    pub merge_gbk: bool,

    /// Write AntiSMASH v6 sideload JSON.
    #[arg(long)]
    pub antismash_sideload: bool,
}

impl PredictArgs {
    pub fn execute(&self) -> anyhow::Result<()> {
        log::info!("Predicting clusters from pre-annotated features");

        // 1. Load gene + feature tables
        // 2. Filter domains
        // 3. Predict probabilities (CRF)
        // 4. Refine clusters
        // 5. Predict types (RandomForest)
        // 6. Write output

        anyhow::bail!("predict subcommand not yet fully wired")
    }
}
