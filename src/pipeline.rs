//! High-level pipeline API for library consumers.
//!
//! # Example
//!
//! ```no_run
//! use gecco::pipeline::Gecco;
//! use gecco::orf::SeqRecord;
//!
//! let pipeline = Gecco::builder()
//!     .data_dir("gecco_data")
//!     .threshold(0.8)
//!     .build()
//!     .unwrap();
//!
//! let records = vec![SeqRecord {
//!     id: "contig_1".into(),
//!     seq: std::fs::read_to_string("genome.fna").unwrap(),
//! }];
//!
//! let results = pipeline.scan(&records).unwrap();
//! for cluster in &results.clusters {
//!     println!("{}: {} genes, type={:?}",
//!         cluster.id, cluster.genes.len(), cluster.cluster_type);
//! }
//! ```

use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result};
use log::info;

use crate::crf::backend::CrfSuiteModel;
use crate::crf::ClusterCRF;
use crate::data_dir;
use crate::hmmer::{self, LoadedHMM, PyHMMER, HMM};
use crate::interpro::InterPro;
use crate::io::{
    genbank,
    tables::{ClusterTable, FeatureTable, GeneTable},
};
use crate::model::{Cluster, Gene};
use crate::orf::{ProdigalFinder, SeqRecord};
use crate::output::{RunOutput, StdioOutput};
use crate::refine::ClusterRefiner;

/// Results from a GECCO pipeline run.
pub struct GeccoResults {
    /// Annotated genes with cluster probabilities.
    pub genes: Vec<Gene>,
    /// Predicted biosynthetic gene clusters.
    pub clusters: Vec<Cluster>,
    /// Source sequences (for GenBank output). Empty if not available.
    pub source_seqs: BTreeMap<String, String>,
}

impl GeccoResults {
    /// Write the gene table (TSV) to a writer.
    pub fn write_gene_table(&self, writer: impl Write) -> Result<()> {
        GeneTable::write_from_genes(writer, &self.genes)
    }

    /// Write the feature/domain table (TSV) to a writer.
    pub fn write_feature_table(&self, writer: impl Write) -> Result<()> {
        FeatureTable::write_from_genes(writer, &self.genes)
    }

    /// Write the cluster table (TSV) to a writer.
    pub fn write_cluster_table(&self, writer: impl Write) -> Result<()> {
        ClusterTable::write_from_clusters(writer, &self.clusters)
    }

    /// Write each cluster as a separate GenBank file in the given directory.
    pub fn write_cluster_gbks(&self, dir: &Path) -> Result<()> {
        for cluster in &self.clusters {
            let path = dir.join(format!("{}.gbk", cluster.id));
            let source_seq = self
                .source_seqs
                .get(cluster.source_id())
                .map(|s| s.as_str());
            genbank::write_cluster_gbk(std::fs::File::create(&path)?, cluster, source_seq)?;
        }
        Ok(())
    }

    /// Write all clusters into a single merged GenBank file.
    pub fn write_clusters_merged_gbk(&self, writer: impl Write) -> Result<()> {
        genbank::write_clusters_merged(writer, &self.clusters, &self.source_seqs)
    }
}

/// Builder for constructing a [`Gecco`] with custom settings.
pub struct GeccoBuilder {
    data_dir: Option<PathBuf>,
    crf_model_path: Option<PathBuf>,
    hmm_configs: Vec<HMM>,
    hmm_paths: Vec<PathBuf>,
    interpro_path: Option<PathBuf>,
    metagenome: bool,
    mask: bool,
    jobs: usize,
    thread_pool: Option<Arc<rayon::ThreadPool>>,
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

impl Default for GeccoBuilder {
    fn default() -> Self {
        Self {
            data_dir: None,
            crf_model_path: None,
            hmm_configs: Vec::new(),
            hmm_paths: Vec::new(),
            interpro_path: None,
            metagenome: true,
            mask: false,
            jobs: 0,
            thread_pool: None,
            p_filter: 1e-9,
            e_filter: None,
            disentangle: false,
            threshold: 0.8,
            n_cds: 3,
            edge_distance: 0,
            trim: true,
            pad: true,
            window_size: 20,
            feature_type: "protein".to_string(),
        }
    }
}

impl GeccoBuilder {
    /// Set the data directory containing default HMM, CRF model, and InterPro files.
    ///
    /// When set, the builder will automatically load default files from this directory
    /// for any paths not explicitly configured. If not set, defaults to `gecco_data/`
    /// next to the binary or the `GECCO_DATA_DIR` environment variable.
    pub fn data_dir(mut self, path: impl Into<PathBuf>) -> Self {
        self.data_dir = Some(path.into());
        self
    }

    /// Path to a CRFsuite model file (required unless data_dir contains model.crfsuite).
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

    /// Use a caller-owned Rayon thread pool for GECCO stages that use Rayon.
    pub fn thread_pool(mut self, pool: Arc<rayon::ThreadPool>) -> Self {
        self.thread_pool = Some(pool);
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

    /// CRF sliding window size (default: 20).
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
    ///
    /// If no CRF model, HMM, or InterPro path was explicitly set, the builder
    /// will look for default files in the data directory.
    pub fn build(self) -> Result<Gecco> {
        let resolved_data_dir = data_dir::resolve(self.data_dir.as_ref());

        // CRF model: explicit path > data_dir default
        let crf_model_path = self
            .crf_model_path
            .unwrap_or_else(|| data_dir::crf_model_path(&resolved_data_dir));
        let crf_model = if crf_model_path.exists() {
            CrfSuiteModel::from_file(&crf_model_path)
                .with_context(|| format!("loading CRF model from {:?}", crf_model_path))?
        } else {
            #[cfg(feature = "bundled-data")]
            {
                CrfSuiteModel::from_bytes(crate::bundled_data::crf_model().to_vec())?
            }

            #[cfg(not(feature = "bundled-data"))]
            {
                CrfSuiteModel::from_file(&crf_model_path)
                    .with_context(|| format!("loading CRF model from {:?}", crf_model_path))?
            }
        };

        // HMMs: explicit configs/paths > data_dir default
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
        if hmms.is_empty() {
            let default_hmm = data_dir::hmm_path(&resolved_data_dir);
            if default_hmm.exists() {
                hmms.push(HMM {
                    id: "Pfam".to_string(),
                    version: String::new(),
                    url: String::new(),
                    path: default_hmm,
                    size: None,
                    relabel_with: Some("s/(PF\\d+).\\d+/\\1/".to_string()),
                    md5: None,
                });
            } else {
                #[cfg(feature = "bundled-data")]
                hmms.push(HMM {
                    id: "Pfam".to_string(),
                    version: String::new(),
                    url: String::new(),
                    path: crate::bundled_data::hmm_path(),
                    size: None,
                    relabel_with: Some("s/(PF\\d+).\\d+/\\1/".to_string()),
                    md5: None,
                });
            }
        }

        // InterPro: explicit path > data_dir default > empty
        let interpro = match &self.interpro_path {
            Some(path) => {
                let data = std::fs::read(path)
                    .with_context(|| format!("reading InterPro JSON from {:?}", path))?;
                InterPro::from_json(&data)?
            }
            None => {
                let default_path = data_dir::interpro_path(&resolved_data_dir);
                if default_path.exists() {
                    let data = std::fs::read(&default_path).with_context(|| {
                        format!("reading InterPro JSON from {:?}", default_path)
                    })?;
                    InterPro::from_json(&data)?
                } else {
                    #[cfg(feature = "bundled-data")]
                    {
                        InterPro::from_json(crate::bundled_data::interpro_json())?
                    }

                    #[cfg(not(feature = "bundled-data"))]
                    InterPro::from_json(b"[]")?
                }
            }
        };

        // Load HMM profiles into memory
        let mut loaded_hmms = Vec::with_capacity(hmms.len());
        for hmm_config in hmms {
            let annotator = PyHMMER::new(hmm_config.clone());
            let profiles = annotator
                .load_hmms()
                .with_context(|| format!("loading HMM profiles from {:?}", hmm_config.path))?;
            loaded_hmms.push(LoadedHMM {
                meta: hmm_config,
                profiles,
            });
        }

        Ok(Gecco {
            crf_model,
            hmms: loaded_hmms,
            interpro,
            metagenome: self.metagenome,
            mask: self.mask,
            jobs: self.jobs,
            thread_pool: self.thread_pool,
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
/// Use [`Gecco::builder()`] to construct one.
pub struct Gecco {
    crf_model: CrfSuiteModel,
    hmms: Vec<LoadedHMM>,
    interpro: InterPro,
    metagenome: bool,
    mask: bool,
    jobs: usize,
    thread_pool: Option<Arc<rayon::ThreadPool>>,
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

impl Gecco {
    /// Create a new [`GeccoBuilder`].
    pub fn builder() -> GeccoBuilder {
        GeccoBuilder::default()
    }

    /// Run the full pipeline on the given sequences, returning a [`GeccoResults`].
    ///
    /// Steps: gene finding → domain annotation → CRF prediction → cluster refinement.
    pub fn scan(&self, records: &[SeqRecord]) -> Result<GeccoResults> {
        self.scan_with_output(records, &StdioOutput)
    }

    /// Run the full pipeline on the given sequences, writing run diagnostics to `output`.
    pub fn scan_with_output(
        &self,
        records: &[SeqRecord],
        output: &dyn RunOutput,
    ) -> Result<GeccoResults> {
        let source_seqs: BTreeMap<String, String> = records
            .iter()
            .map(|r| (r.id.clone(), r.seq.clone()))
            .collect();

        let mut genes = self.find_genes_with_output(records, output)?;
        if genes.is_empty() {
            return Ok(GeccoResults {
                genes: Vec::new(),
                clusters: Vec::new(),
                source_seqs,
            });
        }

        self.annotate_domains_with_output(&mut genes, output)?;
        let genes = self.predict_probabilities_with_output(&genes, output)?;
        let clusters = self.extract_clusters_with_output(&genes, output);
        Ok(GeccoResults {
            genes,
            clusters,
            source_seqs,
        })
    }

    /// Run the pipeline from pre-annotated genes (skipping gene finding and HMMER).
    ///
    /// Use this when you already have genes with domain annotations, e.g. loaded
    /// from a features TSV file.
    pub fn predict_from_genes(&self, genes: &[Gene]) -> Result<GeccoResults> {
        let genes = self.predict_probabilities(genes)?;
        let clusters = self.extract_clusters(&genes);
        Ok(GeccoResults {
            genes,
            clusters,
            source_seqs: BTreeMap::new(),
        })
    }

    // -- individual pipeline stages, exposed for advanced use --

    /// Stage 1: Predict ORFs from DNA sequences.
    pub fn find_genes(&self, records: &[SeqRecord]) -> Result<Vec<Gene>> {
        self.find_genes_with_output(records, &StdioOutput)
    }

    pub fn find_genes_with_output(
        &self,
        records: &[SeqRecord],
        output: &dyn RunOutput,
    ) -> Result<Vec<Gene>> {
        output.stderr(format_args!(
            "Finding genes with Prodigal ({} sequences)",
            records.len()
        ))?;
        let finder = ProdigalFinder {
            metagenome: self.metagenome,
            mask: self.mask,
            cpus: self.jobs,
            ..Default::default()
        };
        let genes = finder.find_genes_with_output(records, output)?;
        output.stderr(format_args!("Found {} genes", genes.len()))?;
        Ok(genes)
    }

    /// Stage 2: Annotate genes with protein domains from HMM databases.
    ///
    /// Domains are appended to each gene's `protein.domains` in place.
    pub fn annotate_domains(&self, genes: &mut [Gene]) -> Result<()> {
        info!(
            "Annotating protein domains ({} HMM databases)",
            self.hmms.len()
        );
        self.annotate_domains_inner(genes, None)?;
        Ok(())
    }

    pub fn annotate_domains_with_output(
        &self,
        genes: &mut [Gene],
        output: &dyn RunOutput,
    ) -> Result<()> {
        output.stderr(format_args!(
            "Annotating protein domains ({} HMM databases)",
            self.hmms.len()
        ))?;
        self.annotate_domains_inner(genes, Some(output))?;
        Ok(())
    }

    fn annotate_domains_inner(
        &self,
        genes: &mut [Gene],
        output: Option<&dyn RunOutput>,
    ) -> Result<()> {
        for loaded in &self.hmms {
            let mut annotator = PyHMMER::new(loaded.meta.clone()).with_cpus(self.jobs);
            if let Some(pool) = &self.thread_pool {
                annotator = annotator.with_thread_pool(Arc::clone(pool));
            }
            annotator.run_with_profiles(genes, &self.interpro, &loaded.profiles, None)?;
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

        genes.sort_by(|a, b| a.source_id.cmp(&b.source_id).then(a.start.cmp(&b.start)));

        let domain_count: usize = genes.iter().map(|g| g.protein.domains.len()).sum();
        if let Some(output) = output {
            output.stderr(format_args!(
                "Found {} domains across all proteins",
                domain_count
            ))?;
        } else {
            info!("Found {} domains across all proteins", domain_count);
        }
        Ok(())
    }

    /// Stage 3: Predict per-gene cluster probabilities with the CRF model.
    pub fn predict_probabilities(&self, genes: &[Gene]) -> Result<Vec<Gene>> {
        info!("Predicting cluster probabilities");
        self.predict_probabilities_inner(genes)
    }

    pub fn predict_probabilities_with_output(
        &self,
        genes: &[Gene],
        output: &dyn RunOutput,
    ) -> Result<Vec<Gene>> {
        output.stderr(format_args!("Predicting cluster probabilities"))?;
        self.predict_probabilities_inner(genes)
    }

    fn predict_probabilities_inner(&self, genes: &[Gene]) -> Result<Vec<Gene>> {
        let mut crf = ClusterCRF::new(&self.feature_type, self.window_size, 1);
        // CrfSuiteModel needs to be reloaded per call since ClusterCRF takes ownership.
        // We clone via the stored model.
        crf.set_model(Box::new(self.crf_model.clone()));
        crf.predict_probabilities(genes, self.pad, None)
    }

    /// Stage 4: Extract clusters from genes with predicted probabilities.
    pub fn extract_clusters(&self, genes: &[Gene]) -> Vec<Cluster> {
        info!("Extracting clusters");
        let clusters = self.extract_clusters_inner(genes);
        info!("Found {} cluster(s)", clusters.len());
        clusters
    }

    pub fn extract_clusters_with_output(
        &self,
        genes: &[Gene],
        output: &dyn RunOutput,
    ) -> Vec<Cluster> {
        let _ = output.stderr(format_args!("Extracting clusters"));
        let clusters = self.extract_clusters_inner(genes);
        let _ = output.stderr(format_args!("Found {} cluster(s)", clusters.len()));
        clusters
    }

    fn extract_clusters_inner(&self, genes: &[Gene]) -> Vec<Cluster> {
        let refiner = ClusterRefiner {
            threshold: self.threshold,
            n_cds: self.n_cds,
            edge_distance: self.edge_distance,
            trim: self.trim,
            ..Default::default()
        };
        refiner.iter_clusters(genes)
    }
}
