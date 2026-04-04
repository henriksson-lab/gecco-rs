//! The `gecco predict` subcommand — predict from pre-annotated features.

use std::collections::BTreeMap;
use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use log::info;

use crate::crf::ClusterCRF;
use crate::hmmer;
use crate::io::genbank;
use crate::io::tables::{ClusterTable, FeatureTable, GeneTable};
use crate::refine::ClusterRefiner;
use crate::types::backend::SmartcoreRF;
use crate::types::TypeClassifier;

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

    /// Alternative CRF model file.
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
    pub fn execute(&self) -> Result<()> {
        let base = self
            .genome
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();

        std::fs::create_dir_all(&self.output_dir)?;

        // 1. Load gene + feature tables
        info!("Loading gene table from {:?}", self.genes);
        let mut genes =
            GeneTable::read_to_genes(std::fs::File::open(&self.genes)?)?;
        info!("Loaded {} genes", genes.len());

        for feat_path in &self.features {
            info!("Loading features from {:?}", feat_path);
            let feat_genes =
                FeatureTable::read_to_genes(std::fs::File::open(feat_path)?)?;
            // Merge domains into existing genes by protein_id
            let domain_map: BTreeMap<String, Vec<_>> = feat_genes
                .into_iter()
                .map(|g| (g.protein.id.clone(), g.protein.domains))
                .collect();
            for gene in &mut genes {
                if let Some(domains) = domain_map.get(&gene.protein.id) {
                    gene.protein.domains = domains.clone();
                }
            }
        }

        // Load sequences for GenBank output
        let records = genbank::read_sequences(&self.genome)?;
        let source_seqs: BTreeMap<String, String> = records
            .iter()
            .map(|r| (r.id.clone(), r.seq.clone()))
            .collect();

        // 2. Filter domains
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

        // 3. Predict probabilities
        info!("Predicting cluster probabilities");
        let crf_model = super::run::load_crf_model(&self.model)?;
        let mut crf = ClusterCRF::new("protein", 5, 1);
        crf.set_model(Box::new(crf_model));
        genes = crf.predict_probabilities(&genes, !self.no_pad, None)?;

        // Write tables
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

        // 4. Extract clusters
        info!("Extracting clusters");
        let refiner = ClusterRefiner {
            threshold: self.threshold,
            n_cds: self.cds,
            edge_distance: self.edge_distance,
            trim: !self.no_trim,
            ..Default::default()
        };
        let mut clusters = refiner.iter_clusters(&genes);

        if clusters.is_empty() {
            log::warn!("No gene clusters found");
            if self.force_tsv {
                ClusterTable::write_from_clusters(
                    std::fs::File::create(
                        self.output_dir.join(format!("{}.clusters.tsv", base)),
                    )?,
                    &[],
                )?;
            }
            return Ok(());
        }
        info!("Found {} cluster(s)", clusters.len());

        // 5. Predict types
        info!("Predicting cluster types");
        let mut classifier = TypeClassifier::new(vec![
            "Alkaloid".into(),
            "NRP".into(),
            "Polyketide".into(),
            "RiPP".into(),
            "Saccharide".into(),
            "Terpene".into(),
        ]);
        let rf = SmartcoreRF::new(6);
        classifier.set_model(Box::new(rf));
        let _ = classifier.predict_types(&mut clusters);

        // 6. Write output
        info!("Writing results to {:?}", self.output_dir);
        ClusterTable::write_from_clusters(
            std::fs::File::create(
                self.output_dir.join(format!("{}.clusters.tsv", base)),
            )?,
            &clusters,
        )?;

        if self.merge_gbk {
            let gbk_path = self.output_dir.join(format!("{}.clusters.gbk", base));
            genbank::write_clusters_merged(
                std::fs::File::create(&gbk_path)?,
                &clusters,
                &source_seqs,
                env!("CARGO_PKG_VERSION"),
            )?;
        } else {
            for cluster in &clusters {
                let gbk_path = self.output_dir.join(format!("{}.gbk", cluster.id));
                let source_seq = source_seqs.get(cluster.source_id()).map(|s| s.as_str());
                genbank::write_cluster_gbk(
                    std::fs::File::create(&gbk_path)?,
                    cluster,
                    source_seq,
                    env!("CARGO_PKG_VERSION"),
                )?;
            }
        }

        Ok(())
    }
}
