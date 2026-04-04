//! High-level pipeline API for library consumers.
//!
//! # Example
//!
//! ```no_run
//! use gecco::pipeline::Pipeline;
//! use gecco::orf::SeqRecord;
//!
//! let pipeline = Pipeline::builder()
//!     .crf_model("path/to/model.crfsuite")
//!     .hmm("path/to/Pfam.hmm")
//!     .threshold(0.8)
//!     .build()
//!     .unwrap();
//!
//! let records = vec![SeqRecord {
//!     id: "contig_1".into(),
//!     seq: std::fs::read_to_string("genome.fna").unwrap(),
//! }];
//!
//! let clusters = pipeline.scan(&records).unwrap();
//! for cluster in &clusters {
//!     println!("{}: {} genes, type={:?}",
//!         cluster.id, cluster.genes.len(), cluster.cluster_type);
//! }
//! ```

use std::path::PathBuf;

use anyhow::{Context, Result};
use log::info;

use crate::crf::backend::CrfSuiteModel;
use crate::crf::ClusterCRF;
use crate::hmmer::{self, DomainAnnotator, PyHMMER, HMM};
use crate::interpro::InterPro;
use crate::model::{Cluster, Gene};
use crate::orf::{ProdigalFinder, SeqRecord, ORFFinder};
use crate::refine::ClusterRefiner;

/// Builder for constructing a [`Pipeline`] with custom settings.
pub struct PipelineBuilder {
    crf_model_path: Option<PathBuf>,
    hmm_configs: Vec<HMM>,
    hmm_paths: Vec<PathBuf>,
    interpro_path: Option<PathBuf>,
    metagenome: bool,
    mask: bool,
    jobs: usize,
    p_filter: f64,
    e_filter: Option<f64>,
    disentangle: bool,
    threshold: f64,
    n_cds: usize,
    edge_distance: usize,
    trim: bool,
    pad: bool,
    window_size: usize,
    feature_type: String,
}

impl Default for PipelineBuilder {
    fn default() -> Self {
        Self {
            crf_model_path: None,
            hmm_configs: Vec::new(),
            hmm_paths: Vec::new(),
            interpro_path: None,
            metagenome: true,
            mask: false,
            jobs: 0,
            p_filter: 1e-9,
            e_filter: None,
            disentangle: false,
            threshold: 0.8,
            n_cds: 3,
            edge_distance: 0,
            trim: true,
            pad: true,
            window_size: 5,
            feature_type: "protein".to_string(),
        }
    }
}

impl PipelineBuilder {
    /// Path to a CRFsuite model file (required).
    pub fn crf_model(mut self, path: impl Into<PathBuf>) -> Self {
        self.crf_model_path = Some(path.into());
        self
    }

    /// Add an HMM database file path. The file is loaded as a custom HMM
    /// without relabeling. Call multiple times for multiple databases.
    pub fn hmm(mut self, path: impl Into<PathBuf>) -> Self {
        self.hmm_paths.push(path.into());
        self
    }

    /// Add a fully configured [`HMM`] descriptor (with relabeling, version, etc.).
    pub fn hmm_config(mut self, hmm: HMM) -> Self {
        self.hmm_configs.push(hmm);
        self
    }

    /// Path to InterPro JSON metadata. If not set, InterPro annotations are skipped.
    pub fn interpro(mut self, path: impl Into<PathBuf>) -> Self {
        self.interpro_path = Some(path.into());
        self
    }

    /// Use metagenomic mode for gene prediction (default: true).
    pub fn metagenome(mut self, v: bool) -> Self {
        self.metagenome = v;
        self
    }

    /// Mask ambiguous nucleotides during gene prediction.
    pub fn mask(mut self, v: bool) -> Self {
        self.mask = v;
        self
    }

    /// Number of threads for parallelizable steps (0 = auto-detect).
    pub fn jobs(mut self, v: usize) -> Self {
        self.jobs = v;
        self
    }

    /// P-value cutoff for domain filtering (default: 1e-9).
    pub fn p_filter(mut self, v: f64) -> Self {
        self.p_filter = v;
        self
    }

    /// E-value cutoff for domain filtering. If not set, no e-value filter is applied.
    pub fn e_filter(mut self, v: f64) -> Self {
        self.e_filter = Some(v);
        self
    }

    /// Remove overlapping domains, keeping the one with the lowest p-value.
    pub fn disentangle(mut self, v: bool) -> Self {
        self.disentangle = v;
        self
    }

    /// Minimum probability for a gene to be considered part of a BGC (default: 0.8).
    pub fn threshold(mut self, v: f64) -> Self {
        self.threshold = v;
        self
    }

    /// Minimum number of annotated CDS in a cluster (default: 3).
    pub fn n_cds(mut self, v: usize) -> Self {
        self.n_cds = v;
        self
    }

    /// Minimum genes separating a cluster from the sequence edge (default: 0).
    pub fn edge_distance(mut self, v: usize) -> Self {
        self.edge_distance = v;
        self
    }

    /// Trim unannotated genes from cluster edges (default: true).
    pub fn trim(mut self, v: bool) -> Self {
        self.trim = v;
        self
    }

    /// Pad short sequences to fit the CRF window (default: true).
    pub fn pad(mut self, v: bool) -> Self {
        self.pad = v;
        self
    }

    /// CRF sliding window size (default: 5).
    pub fn window_size(mut self, v: usize) -> Self {
        self.window_size = v;
        self
    }

    /// Feature extraction mode: "protein" or "domain" (default: "protein").
    pub fn feature_type(mut self, v: impl Into<String>) -> Self {
        self.feature_type = v.into();
        self
    }

    /// Build the pipeline, loading models from the configured paths.
    pub fn build(self) -> Result<Pipeline> {
        let crf_model_path = self
            .crf_model_path
            .ok_or_else(|| anyhow::anyhow!("CRF model path is required — call .crf_model()"))?;

        let crf_model = CrfSuiteModel::from_file(&crf_model_path)
            .with_context(|| format!("loading CRF model from {:?}", crf_model_path))?;

        let mut hmms = self.hmm_configs;
        for (i, path) in self.hmm_paths.iter().enumerate() {
            hmms.push(HMM {
                id: format!("HMM{}", i),
                version: String::new(),
                url: String::new(),
                path: path.clone(),
                size: None,
                relabel_with: None,
                md5: None,
            });
        }

        let interpro = match &self.interpro_path {
            Some(path) => {
                let data = std::fs::read(path)
                    .with_context(|| format!("reading InterPro JSON from {:?}", path))?;
                InterPro::from_json(&data)?
            }
            None => InterPro::from_json(b"[]")?,
        };

        Ok(Pipeline {
            crf_model,
            hmms,
            interpro,
            metagenome: self.metagenome,
            mask: self.mask,
            jobs: self.jobs,
            p_filter: self.p_filter,
            e_filter: self.e_filter,
            disentangle: self.disentangle,
            threshold: self.threshold,
            n_cds: self.n_cds,
            edge_distance: self.edge_distance,
            trim: self.trim,
            pad: self.pad,
            window_size: self.window_size,
            feature_type: self.feature_type,
        })
    }
}

/// A configured GECCO pipeline that can scan sequences for biosynthetic gene clusters.
///
/// Use [`Pipeline::builder()`] to construct one.
pub struct Pipeline {
    crf_model: CrfSuiteModel,
    hmms: Vec<HMM>,
    interpro: InterPro,
    metagenome: bool,
    mask: bool,
    jobs: usize,
    p_filter: f64,
    e_filter: Option<f64>,
    disentangle: bool,
    threshold: f64,
    n_cds: usize,
    edge_distance: usize,
    trim: bool,
    pad: bool,
    window_size: usize,
    feature_type: String,
}

impl Pipeline {
    /// Create a new [`PipelineBuilder`].
    pub fn builder() -> PipelineBuilder {
        PipelineBuilder::default()
    }

    /// Run the full pipeline on the given sequences, returning predicted gene clusters.
    ///
    /// Steps: gene finding → domain annotation → CRF prediction → cluster refinement.
    pub fn scan(&self, records: &[SeqRecord]) -> Result<Vec<Cluster>> {
        // 1. Gene finding
        let mut genes = self.find_genes(records)?;
        if genes.is_empty() {
            return Ok(Vec::new());
        }

        // 2. Domain annotation
        self.annotate_domains(&mut genes)?;

        // 3. CRF prediction
        let genes = self.predict_probabilities(&genes)?;

        // 4. Cluster refinement
        Ok(self.extract_clusters(&genes))
    }

    /// Run the full pipeline, also returning the annotated genes alongside clusters.
    ///
    /// Useful when you need both gene-level probabilities and cluster predictions.
    pub fn scan_detailed(&self, records: &[SeqRecord]) -> Result<(Vec<Gene>, Vec<Cluster>)> {
        let mut genes = self.find_genes(records)?;
        if genes.is_empty() {
            return Ok((Vec::new(), Vec::new()));
        }

        self.annotate_domains(&mut genes)?;
        let genes = self.predict_probabilities(&genes)?;
        let clusters = self.extract_clusters(&genes);
        Ok((genes, clusters))
    }

    /// Run the pipeline from pre-annotated genes (skipping gene finding and HMMER).
    ///
    /// Use this when you already have genes with domain annotations, e.g. loaded
    /// from a features TSV file.
    pub fn predict_from_genes(&self, genes: &[Gene]) -> Result<(Vec<Gene>, Vec<Cluster>)> {
        let genes = self.predict_probabilities(genes)?;
        let clusters = self.extract_clusters(&genes);
        Ok((genes, clusters))
    }

    // -- individual pipeline stages, exposed for advanced use --

    /// Stage 1: Predict ORFs from DNA sequences.
    pub fn find_genes(&self, records: &[SeqRecord]) -> Result<Vec<Gene>> {
        info!("Finding genes with Prodigal ({} sequences)", records.len());
        let finder = ProdigalFinder {
            metagenome: self.metagenome,
            mask: self.mask,
            cpus: self.jobs,
            ..Default::default()
        };
        let genes = finder.find_genes(records)?;
        info!("Found {} genes", genes.len());
        Ok(genes)
    }

    /// Stage 2: Annotate genes with protein domains from HMM databases.
    ///
    /// Domains are appended to each gene's `protein.domains` in place.
    pub fn annotate_domains(&self, genes: &mut Vec<Gene>) -> Result<()> {
        info!("Annotating protein domains ({} HMM databases)", self.hmms.len());
        for hmm_config in &self.hmms {
            let annotator = PyHMMER::new(hmm_config.clone());
            annotator.run(genes, &self.interpro, None)?;
        }

        if self.disentangle {
            for gene in genes.iter_mut() {
                hmmer::disentangle(gene);
            }
        }
        if let Some(e) = self.e_filter {
            hmmer::filter_by_evalue(genes, e);
        }
        hmmer::filter_by_pvalue(genes, self.p_filter);

        genes.sort_by(|a, b| {
            a.source_id
                .cmp(&b.source_id)
                .then(a.start.cmp(&b.start))
        });

        let domain_count: usize = genes.iter().map(|g| g.protein.domains.len()).sum();
        info!("Found {} domains across all proteins", domain_count);
        Ok(())
    }

    /// Stage 3: Predict per-gene cluster probabilities with the CRF model.
    pub fn predict_probabilities(&self, genes: &[Gene]) -> Result<Vec<Gene>> {
        info!("Predicting cluster probabilities");
        let mut crf = ClusterCRF::new(&self.feature_type, self.window_size, 1);
        // CrfSuiteModel needs to be reloaded per call since ClusterCRF takes ownership.
        // We clone via the stored model.
        crf.set_model(Box::new(self.crf_model.clone()));
        crf.predict_probabilities(genes, self.pad, None)
    }

    /// Stage 4: Extract clusters from genes with predicted probabilities.
    pub fn extract_clusters(&self, genes: &[Gene]) -> Vec<Cluster> {
        info!("Extracting clusters");
        let refiner = ClusterRefiner {
            threshold: self.threshold,
            n_cds: self.n_cds,
            edge_distance: self.edge_distance,
            trim: self.trim,
            ..Default::default()
        };
        let clusters = refiner.iter_clusters(genes);
        info!("Found {} cluster(s)", clusters.len());
        clusters
    }
}
