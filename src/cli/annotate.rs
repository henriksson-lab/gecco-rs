//! The `gecco annotate` subcommand — domain annotation only.

use std::path::PathBuf;

use clap::Args;

#[derive(Args)]
pub struct AnnotateArgs {
    /// Path to the input genome (FASTA or GenBank).
    #[arg(short, long)]
    pub genome: PathBuf,

    /// Output directory.
    #[arg(short, long, default_value = ".")]
    pub output_dir: PathBuf,

    /// Number of threads (0 = auto-detect).
    #[arg(short, long, default_value = "0")]
    pub jobs: usize,

    /// Enable masking of ambiguous nucleotides.
    #[arg(short = 'M', long)]
    pub mask: bool,

    /// Extract genes from annotated records using this feature type.
    #[arg(long)]
    pub cds_feature: Option<String>,

    /// Feature qualifier for gene names when extracting CDS.
    #[arg(long, default_value = "locus_tag")]
    pub locus_tag: String,

    /// Path to additional HMM file(s).
    #[arg(long)]
    pub hmm: Vec<PathBuf>,

    /// E-value cutoff for protein domains.
    #[arg(short, long)]
    pub e_filter: Option<f64>,

    /// P-value cutoff for protein domains.
    #[arg(short, long, default_value = "1e-9")]
    pub p_filter: f64,

    /// Disentangle overlapping domains.
    #[arg(long)]
    pub disentangle: bool,

    /// Write empty TSV files when no results are found.
    #[arg(long)]
    pub force_tsv: bool,
}

impl AnnotateArgs {
    pub fn execute(&self) -> anyhow::Result<()> {
        log::info!("Annotating domains in {:?}", self.genome);
        log::info!("Output directory: {:?}", self.output_dir);

        // 1. Load sequences
        // 2. Find genes (Prodigal or CDSFinder)
        // 3. Annotate domains (HMMER)
        // 4. Filter domains
        // 5. Write gene + feature tables

        anyhow::bail!("annotate subcommand not yet fully wired")
    }
}
