//! The `gecco cv` subcommand — cross-validation.

use std::collections::BTreeMap;
use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use log::info;

use crate::crf::backend::CrfSuiteModel;
use crate::crf::ClusterCRF;
use crate::hmmer;
use crate::io::tables::{FeatureTable, GeneTable};
use crate::model::Gene;
use crate::refine::ClusterRefiner;

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
    pub fn execute(&self) -> Result<()> {
        let mode = if self.loto {
            "Leave-One-Type-Out".to_string()
        } else {
            format!("{}-fold", self.splits)
        };
        info!("Running {} cross-validation", mode);

        // 1. Load data
        info!("Loading gene table from {:?}", self.genes);
        let mut genes =
            GeneTable::read_to_genes(std::fs::File::open(&self.genes)?)?;

        for feat_path in &self.features {
            info!("Loading features from {:?}", feat_path);
            let feat_genes =
                FeatureTable::read_to_genes(std::fs::File::open(feat_path)?)?;
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

        // 2. Filter domains
        if let Some(e) = self.e_filter {
            hmmer::filter_by_evalue(&mut genes, e);
        }
        hmmer::filter_by_pvalue(&mut genes, self.p_filter);

        genes.sort_by(|a, b| {
            a.source_id
                .cmp(&b.source_id)
                .then(a.start.cmp(&b.start))
        });

        // 3. Load cluster labels
        let cluster_reader = std::fs::File::open(&self.clusters)?;
        let mut cluster_rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(cluster_reader);

        let mut cluster_ranges: BTreeMap<String, Vec<(i64, i64)>> = BTreeMap::new();
        for result in cluster_rdr.deserialize::<crate::io::tables::ClusterRow>() {
            let row = result?;
            cluster_ranges
                .entry(row.sequence_id.clone())
                .or_default()
                .push((row.start, row.end));
        }

        for gene in &mut genes {
            let in_cluster = cluster_ranges
                .get(&gene.source_id)
                .map(|ranges| {
                    ranges.iter().any(|(s, e)| gene.start <= *e && gene.end >= *s)
                })
                .unwrap_or(false);
            gene.probability = Some(if in_cluster { 1.0 } else { 0.0 });
            for domain in &mut gene.protein.domains {
                domain.probability = Some(if in_cluster { 1.0 } else { 0.0 });
            }
        }

        // 4. Group genes by source sequence
        let mut sequences: Vec<Vec<Gene>> = Vec::new();
        let mut current_source: Option<String> = None;
        for gene in &genes {
            match &current_source {
                Some(s) if s == &gene.source_id => {
                    sequences.last_mut().unwrap().push(gene.clone());
                }
                _ => {
                    current_source = Some(gene.source_id.clone());
                    sequences.push(vec![gene.clone()]);
                }
            }
        }

        let n_seqs = sequences.len();
        let n_folds = if self.loto { n_seqs } else { self.splits.min(n_seqs) };

        info!("Running {} folds on {} sequences", n_folds, n_seqs);

        // 5. K-fold cross-validation
        use std::io::Write;
        let mut out = std::fs::File::create(&self.output)?;
        writeln!(out, "fold\tsequence\tn_genes\tn_clusters")?;

        for fold in 0..n_folds {
            info!("Fold {}/{}", fold + 1, n_folds);

            // Split: test = sequences in this fold, train = rest
            let mut train_genes = Vec::new();
            let mut test_genes = Vec::new();

            for (i, seq) in sequences.iter().enumerate() {
                if i % n_folds == fold {
                    test_genes.extend(seq.clone());
                } else {
                    train_genes.extend(seq.clone());
                }
            }

            // Train CRF on training set
            let crf_model = CrfSuiteModel::train(&[], &[]);
            if crf_model.is_err() {
                log::warn!("Failed to train fold {}", fold);
                continue;
            }
            let mut crf =
                ClusterCRF::new(&self.feature_type, self.window_size, self.window_step);
            crf.set_model(Box::new(crf_model.unwrap()));

            if crf.fit(&train_genes, !self.no_shuffle).is_err() {
                log::warn!("Training failed for fold {}", fold);
                continue;
            }

            // Predict on test set
            let predicted = crf.predict_probabilities(&test_genes, true, None)?;

            // Extract clusters
            let refiner = ClusterRefiner::default();
            let clusters = refiner.iter_clusters(&predicted);

            let test_source = test_genes
                .first()
                .map(|g| g.source_id.as_str())
                .unwrap_or("?");
            writeln!(
                out,
                "{}\t{}\t{}\t{}",
                fold,
                test_source,
                test_genes.len(),
                clusters.len()
            )?;
        }

        info!("Cross-validation results written to {:?}", self.output);
        Ok(())
    }
}
