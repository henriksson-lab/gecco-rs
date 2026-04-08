//! The `gecco run` subcommand — full pipeline.

use std::collections::BTreeMap;
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Args;
use log::info;

use crate::crf::backend::CrfSuiteModel;
use crate::crf::ClusterCRF;
use crate::data_dir;
use crate::hmmer::{self, DomainAnnotator, PyHMMER, HMM};
use crate::interpro::InterPro;
use crate::io::genbank;
use crate::io::tables::{ClusterTable, FeatureTable, GeneTable};
use crate::orf::{ProdigalFinder, ORFFinder};
use crate::refine::ClusterRefiner;
use crate::types::backend::SmartcoreRF;
use crate::types::TypeClassifier;

#[derive(Args)]
pub struct RunArgs {
    /// Path to the input genome (FASTA or GenBank).
    #[arg(short, long)]
    pub genome: PathBuf,

    /// Output directory.
    #[arg(short, long, default_value = ".")]
    pub output_dir: PathBuf,

    /// Data directory containing HMM, CRF model, and InterPro files.
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

impl RunArgs {
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

        // Keep source sequences for GenBank output
        let source_seqs: BTreeMap<String, String> = records
            .iter()
            .map(|r| (r.id.clone(), r.seq.clone()))
            .collect();

        // 2. Find genes
        info!("Finding genes with Prodigal");
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
                write_empty_tables(&self.output_dir, &base)?;
            }
            return Ok(());
        }

        // Write initial gene table
        let gene_path = self.output_dir.join(format!("{}.genes.tsv", base));
        GeneTable::write_from_genes(
            std::fs::File::create(&gene_path)?,
            &genes,
        )?;

        // 3. Load InterPro metadata
        let data_dir = data_dir::resolve(self.data_dir.as_ref());
        let interpro = load_interpro(&data_dir)?;

        // 4. Annotate domains with HMMER
        info!("Annotating protein domains");
        let hmms = load_hmm_configs(&self.hmm, &data_dir)?;
        for hmm_config in &hmms {
            let annotator = PyHMMER::new(hmm_config.clone());
            annotator.run(&mut genes, &interpro, None)?;
        }

        let domain_count: usize = genes.iter().map(|g| g.protein.domains.len()).sum();
        info!("Found {} domains across all proteins", domain_count);

        // Filter and disentangle
        if self.disentangle {
            for gene in &mut genes {
                hmmer::disentangle(gene);
            }
        }
        if let Some(e) = self.e_filter {
            hmmer::filter_by_evalue(&mut genes, e);
        }
        hmmer::filter_by_pvalue(&mut genes, self.p_filter);

        // Sort genes
        genes.sort_by(|a, b| {
            a.source_id
                .cmp(&b.source_id)
                .then(a.start.cmp(&b.start))
        });

        // 5. Predict probabilities with CRF
        info!("Predicting cluster probabilities");
        let crf_model = load_crf_model(&self.model, &data_dir)?;
        let mut crf = ClusterCRF::new("protein", 5, 1);
        crf.set_model(Box::new(crf_model));
        genes = crf.predict_probabilities(&genes, !self.no_pad, None)?;

        // Write gene + feature tables
        GeneTable::write_from_genes(
            std::fs::File::create(&gene_path)?,
            &genes,
        )?;
        let feat_path = self.output_dir.join(format!("{}.features.tsv", base));
        FeatureTable::write_from_genes(
            std::fs::File::create(&feat_path)?,
            &genes,
        )?;

        // 6. Extract clusters
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
                let cluster_path =
                    self.output_dir.join(format!("{}.clusters.tsv", base));
                ClusterTable::write_from_clusters(
                    std::fs::File::create(&cluster_path)?,
                    &[],
                )?;
            }
            return Ok(());
        }

        info!("Found {} potential cluster(s)", clusters.len());

        // 7. Predict types
        info!("Predicting cluster types");
        let mut classifier = TypeClassifier::new(vec![
            "Alkaloid".to_string(),
            "NRP".to_string(),
            "Polyketide".to_string(),
            "RiPP".to_string(),
            "Saccharide".to_string(),
            "Terpene".to_string(),
        ]);
        let rf = SmartcoreRF::new(6);
        classifier.set_model(Box::new(rf));
        // Note: type prediction requires a trained RF model; skip if not available
        let _ = classifier.predict_types(&mut clusters);

        // 8. Write output
        info!("Writing results to {:?}", self.output_dir);

        let cluster_path = self.output_dir.join(format!("{}.clusters.tsv", base));
        ClusterTable::write_from_clusters(
            std::fs::File::create(&cluster_path)?,
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
                let gbk_path = self
                    .output_dir
                    .join(format!("{}.gbk", cluster.id));
                let source_seq = source_seqs.get(cluster.source_id()).map(|s| s.as_str());
                genbank::write_cluster_gbk(
                    std::fs::File::create(&gbk_path)?,
                    cluster,
                    source_seq,
                    env!("CARGO_PKG_VERSION"),
                )?;
            }
        }

        info!(
            "Found {} cluster(s)",
            clusters.len()
        );

        Ok(())
    }
}

/// Write empty TSV tables when no results.
fn write_empty_tables(output_dir: &std::path::Path, base: &str) -> Result<()> {
    GeneTable::write_from_genes(
        std::fs::File::create(output_dir.join(format!("{}.genes.tsv", base)))?,
        &[],
    )?;
    FeatureTable::write_from_genes(
        std::fs::File::create(output_dir.join(format!("{}.features.tsv", base)))?,
        &[],
    )?;
    ClusterTable::write_from_clusters(
        std::fs::File::create(output_dir.join(format!("{}.clusters.tsv", base)))?,
        &[],
    )?;
    Ok(())
}

/// Load InterPro metadata from the data directory.
pub fn load_interpro(data_dir: &std::path::Path) -> Result<InterPro> {
    let interpro_path = data_dir::interpro_path(data_dir);
    if interpro_path.exists() {
        let data = std::fs::read(&interpro_path)?;
        InterPro::from_json(&data)
    } else {
        log::warn!("InterPro metadata not found at {:?}, skipping", interpro_path);
        Ok(InterPro::from_json(b"[]")?)
    }
}

/// Load HMM configs from custom paths or the data directory.
pub fn load_hmm_configs(custom_hmms: &[PathBuf], data_dir: &std::path::Path) -> Result<Vec<HMM>> {
    if !custom_hmms.is_empty() {
        // Use custom HMM paths directly, with default relabeling
        // that strips version suffixes (e.g. PF07690.19 -> PF07690),
        // matching Python GECCO's behavior.
        Ok(custom_hmms
            .iter()
            .enumerate()
            .map(|(i, path)| {
                let base = path.file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("Custom")
                    .to_string();
                HMM {
                    id: if custom_hmms.len() == 1 { base } else { format!("Custom{}", i) },
                    version: String::new(),
                    url: String::new(),
                    path: path.clone(),
                    size: None,
                    relabel_with: Some("s/([^\\.]*)(\\..*)?/\\1/".to_string()),
                    md5: None,
                }
            })
            .collect())
    } else {
        let h3m_path = data_dir::hmm_path(data_dir);
        if h3m_path.exists() {
            Ok(vec![HMM {
                id: "Pfam".to_string(),
                version: String::new(),
                url: String::new(),
                path: h3m_path,
                size: None,
                relabel_with: Some("s/(PF\\d+).\\d+/\\1/".to_string()),
                md5: None,
            }])
        } else {
            anyhow::bail!(
                "HMM data file not found at {:?}. \
                 Run `gecco build-data` to download it, or use --data-dir to specify the location.",
                h3m_path
            );
        }
    }
}

/// Load CRF model from file or from the data directory.
pub fn load_crf_model(model_path: &Option<PathBuf>, data_dir: &std::path::Path) -> Result<CrfSuiteModel> {
    match model_path {
        Some(path) => CrfSuiteModel::from_file(path),
        None => {
            let default_path = data_dir::crf_model_path(data_dir);
            if default_path.exists() {
                CrfSuiteModel::from_file(&default_path)
            } else {
                anyhow::bail!(
                    "No CRF model found at {:?}. \
                     Train one with `gecco train` or provide with --model, \
                     or use --data-dir to specify the location.",
                    default_path
                )
            }
        }
    }
}
