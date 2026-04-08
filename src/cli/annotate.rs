//! The `gecco annotate` subcommand — domain annotation only.

use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Args;
use log::info;

use crate::data_dir;
use crate::hmmer::{self, DomainAnnotator, PyHMMER};
use crate::io::genbank;
use crate::io::tables::{FeatureTable, GeneTable};
use crate::orf::{ProdigalFinder, ORFFinder};

#[derive(Args)]
pub struct AnnotateArgs {
    /// Path to the input genome (FASTA or GenBank).
    #[arg(short, long)]
    pub genome: PathBuf,

    /// Output directory.
    #[arg(short, long, default_value = ".")]
    pub output_dir: PathBuf,

    /// Data directory containing HMM and InterPro files.
    /// Defaults to gecco_data/ next to the binary, or GECCO_DATA_DIR env var.
    #[arg(long)]
    pub data_dir: Option<PathBuf>,

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
    pub fn execute(&self) -> Result<()> {
        let base = self
            .genome
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();

        std::fs::create_dir_all(&self.output_dir)?;

        // 1. Load sequences
        info!("Loading sequences from {:?}", self.genome);
        let records = genbank::read_sequences(&self.genome)
            .with_context(|| format!("loading sequences from {:?}", self.genome))?;
        info!("Loaded {} sequence(s)", records.len());

        // 2. Find genes
        info!("Finding genes");
        let finder = ProdigalFinder {
            metagenome: true,
            mask: self.mask,
            cpus: self.jobs,
            ..Default::default()
        };
        let mut genes = finder.find_genes(&records)?;
        info!("Found {} genes", genes.len());

        if genes.is_empty() {
            log::warn!("No genes found");
            if self.force_tsv {
                GeneTable::write_from_genes(
                    std::fs::File::create(
                        self.output_dir.join(format!("{}.genes.tsv", base)),
                    )?,
                    &[],
                )?;
                FeatureTable::write_from_genes(
                    std::fs::File::create(
                        self.output_dir.join(format!("{}.features.tsv", base)),
                    )?,
                    &[],
                )?;
            }
            return Ok(());
        }

        // 3. Annotate domains
        info!("Annotating protein domains");
        let data_dir = data_dir::resolve(self.data_dir.as_ref());
        let interpro = super::run::load_interpro(&data_dir)?;
        let hmms = super::run::load_hmm_configs(&self.hmm, &data_dir)?;
        for hmm_config in &hmms {
            let annotator = PyHMMER::new(hmm_config.clone());
            annotator.run(&mut genes, &interpro, None)?;
        }

        let domain_count: usize = genes.iter().map(|g| g.protein.domains.len()).sum();
        info!("Found {} domains", domain_count);

        // 4. Filter
        if self.disentangle {
            for gene in &mut genes {
                hmmer::disentangle(gene);
            }
        }
        if let Some(e) = self.e_filter {
            hmmer::filter_by_evalue(&mut genes, e);
        }
        hmmer::filter_by_pvalue(&mut genes, self.p_filter);

        // Sort
        genes.sort_by(|a, b| {
            a.source_id
                .cmp(&b.source_id)
                .then(a.start.cmp(&b.start))
        });

        // 5. Write output
        info!("Writing results to {:?}", self.output_dir);
        GeneTable::write_from_genes(
            std::fs::File::create(self.output_dir.join(format!("{}.genes.tsv", base)))?,
            &genes,
        )?;
        FeatureTable::write_from_genes(
            std::fs::File::create(
                self.output_dir.join(format!("{}.features.tsv", base)),
            )?,
            &genes,
        )?;

        info!("Done");
        Ok(())
    }
}
