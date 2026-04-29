//! The `gecco train` subcommand — train a new CRF model.

use std::collections::{BTreeMap, HashSet};
use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use log::info;

use crate::crf::backend::CrfSuiteModel;
use crate::crf::ClusterCRF;
use crate::hmmer;
use crate::io::tables::{FeatureTable, GeneTable};

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
    pub fn execute(&self) -> Result<()> {
        std::fs::create_dir_all(&self.output_dir)?;

        // 1. Load gene + feature tables
        info!("Loading gene table from {:?}", self.genes);
        let mut genes = GeneTable::read_to_genes(std::fs::File::open(&self.genes)?)?;

        for feat_path in &self.features {
            info!("Loading features from {:?}", feat_path);
            let feat_genes = FeatureTable::read_to_genes(std::fs::File::open(feat_path)?)?;
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

        // Sort
        genes.sort_by(|a, b| a.source_id.cmp(&b.source_id).then(a.start.cmp(&b.start)));

        // 3. Load cluster table and label genes
        info!("Loading cluster annotations from {:?}", self.clusters);
        let cluster_reader = std::fs::File::open(&self.clusters)?;
        let mut cluster_rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(cluster_reader);

        // Read cluster boundaries
        let mut cluster_ranges: BTreeMap<String, Vec<(i64, i64, String)>> = BTreeMap::new();
        for result in cluster_rdr.deserialize::<crate::io::tables::ClusterRow>() {
            let row = result?;
            cluster_ranges
                .entry(row.sequence_id.clone())
                .or_default()
                .push((row.start, row.end, row.cluster_id));
        }

        // Label genes: set probability to 1.0 if inside a cluster, 0.0 otherwise
        for gene in &mut genes {
            let in_cluster = cluster_ranges
                .get(&gene.source_id)
                .map(|ranges| {
                    ranges
                        .iter()
                        .any(|(s, e, _)| gene.start <= *e && gene.end >= *s)
                })
                .unwrap_or(false);
            gene.probability = Some(if in_cluster { 1.0 } else { 0.0 });
            for domain in &mut gene.protein.domains {
                domain.probability = Some(if in_cluster { 1.0 } else { 0.0 });
            }
        }

        // 4. Fit CRF model
        info!(
            "Training CRF model (feature_type={}, window={}x{})",
            self.feature_type, self.window_size, self.window_step
        );
        let mut crf = ClusterCRF::new(&self.feature_type, self.window_size, self.window_step);

        // Initialize with a CrfSuiteModel; fit() will train via CrfModel::fit()
        let init_model: Box<dyn crate::crf::CrfModel> = Box::new(CrfSuiteModel::empty());
        crf.set_model(init_model);
        crf.fit(&genes, !self.no_shuffle)?;

        info!("CRF model trained successfully");

        // 6. Save domain list
        let all_domains: Vec<String> = {
            let mut set: HashSet<String> = HashSet::new();
            for gene in &genes {
                for domain in &gene.protein.domains {
                    set.insert(domain.name.clone());
                }
            }
            let mut v: Vec<String> = set.into_iter().collect();
            v.sort();
            v
        };

        let domains_path = self.output_dir.join("domains.tsv");
        let mut domains_file = std::fs::File::create(&domains_path)?;
        use std::io::Write;
        for d in &all_domains {
            writeln!(domains_file, "{}", d)?;
        }

        info!("Finished training");
        Ok(())
    }
}
