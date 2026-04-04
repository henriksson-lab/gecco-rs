//! CRF-based gene cluster prediction.

pub mod backend;
pub mod features;

use std::collections::HashMap;

use anyhow::Result;

use crate::model::Gene;
use crate::util::sliding_window;

use self::features::{
    annotate_probabilities_domain, annotate_probabilities_protein, extract_features_domain,
    extract_features_protein,
};

/// Trait for a CRF model backend.
pub trait CrfModel: Send + Sync {
    /// Predict marginal probabilities for a single sequence of feature dicts.
    /// Returns a Vec of HashMap<label, probability> for each position.
    fn predict_marginals_single(
        &self,
        features: &[HashMap<String, bool>],
    ) -> Vec<HashMap<String, f64>>;

    /// Get state feature weights: (feature_name, label) -> weight.
    fn state_features(&self) -> &HashMap<(String, String), f64>;

    /// Fit the model on training data.
    fn fit(
        &mut self,
        x: &[Vec<HashMap<String, bool>>],
        y: &[Vec<String>],
    ) -> Result<()>;
}

/// Wrapper for CRF-based cluster prediction.
pub struct ClusterCRF {
    pub feature_type: String,
    pub window_size: usize,
    pub window_step: usize,
    pub significant_features: Option<Vec<String>>,
    model: Option<Box<dyn CrfModel>>,
}

impl ClusterCRF {
    pub fn new(feature_type: &str, window_size: usize, window_step: usize) -> Self {
        Self {
            feature_type: feature_type.to_string(),
            window_size,
            window_step,
            significant_features: None,
            model: None,
        }
    }

    pub fn set_model(&mut self, model: Box<dyn CrfModel>) {
        self.model = Some(model);
    }

    /// Predict cluster probabilities for each gene.
    pub fn predict_probabilities(
        &self,
        genes: &[Gene],
        pad: bool,
        progress: Option<&dyn Fn(usize, usize)>,
    ) -> Result<Vec<Gene>> {
        let model = self
            .model
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("CRF model not fitted/loaded"))?;

        let progress = progress.unwrap_or(&|_, _| {});

        // Sort genes by source and start
        let mut genes = genes.to_vec();
        genes.sort_by(|a, b| {
            a.source_id
                .cmp(&b.source_id)
                .then(a.start.cmp(&b.start))
        });
        for gene in &mut genes {
            gene.protein
                .domains
                .sort_by(|a, b| a.start.cmp(&b.start));
        }

        // Group by contig
        let mut contigs: Vec<(String, Vec<usize>)> = Vec::new();
        for (i, gene) in genes.iter().enumerate() {
            match contigs.last_mut() {
                Some((id, indices)) if id == &gene.source_id => indices.push(i),
                _ => contigs.push((gene.source_id.clone(), vec![i])),
            }
        }

        // Extract features per contig
        let extract_features: fn(&[Gene]) -> Vec<HashMap<String, bool>> =
            match self.feature_type.as_str() {
                "protein" => extract_features_protein,
                "domain" => extract_features_domain,
                _ => return Err(anyhow::anyhow!("invalid feature type")),
            };

        struct ContigData {
            feats: Vec<HashMap<String, bool>>,
            delta: usize,
        }

        let mut contig_data: HashMap<String, ContigData> = HashMap::new();
        let mut total_windows = 0usize;

        for (contig_id, indices) in &contigs {
            let contig_genes: Vec<Gene> = indices.iter().map(|&i| genes[i].clone()).collect();
            let mut feats = extract_features(&contig_genes);
            let mut delta = 0usize;

            if feats.len() < self.window_size {
                if pad {
                    delta = self.window_size - feats.len();
                    let pad_left = delta / 2;
                    let pad_right = (delta + 1) / 2;
                    let mut padded = vec![HashMap::new(); pad_left];
                    padded.append(&mut feats);
                    padded.extend(vec![HashMap::new(); pad_right]);
                    feats = padded;
                } else {
                    continue;
                }
            }

            total_windows += feats.len() - self.window_size + 1;
            contig_data.insert(contig_id.clone(), ContigData { feats, delta });
        }

        progress(0, total_windows);
        let mut window_index = 0usize;

        // Predict probabilities
        let mut predicted = Vec::new();

        for (contig_id, indices) in &contigs {
            let contig_genes: Vec<Gene> = indices.iter().map(|&i| genes[i].clone()).collect();

            let data = match contig_data.get(contig_id) {
                Some(d) => d,
                None => {
                    predicted.extend(contig_genes);
                    continue;
                }
            };

            let n = contig_genes.len().max(self.window_size);
            let mut probabilities = vec![0.0f64; n];

            for win in sliding_window(data.feats.len(), self.window_size, self.window_step) {
                let marginals = model.predict_marginals_single(&data.feats[win.clone()]);
                for (j, m) in marginals.iter().enumerate() {
                    let p = m.get("1").copied().unwrap_or(0.0);
                    let idx = win.start + j;
                    if p > probabilities[idx] {
                        probabilities[idx] = p;
                    }
                }
                window_index += 1;
                progress(window_index, total_windows);
            }

            // Annotate genes with probabilities
            let offset = data.delta / 2;
            let probs_slice = &probabilities[offset..offset + contig_genes.len()];

            let annotated: Vec<Gene> = match self.feature_type.as_str() {
                "protein" => annotate_probabilities_protein(&contig_genes, probs_slice),
                "domain" => annotate_probabilities_domain(&contig_genes, probs_slice),
                _ => contig_genes,
            };
            predicted.extend(annotated);
        }

        // Assign cluster weights from CRF state features
        let state_features = model.state_features();
        let predicted = predicted
            .into_iter()
            .map(|gene| {
                let new_domains: Vec<_> = gene
                    .protein
                    .domains
                    .iter()
                    .map(|d| {
                        let weight = state_features
                            .get(&(d.name.clone(), "1".to_string()))
                            .copied();
                        d.with_cluster_weight(weight)
                    })
                    .collect();
                gene.with_protein(gene.protein.with_domains(new_domains))
            })
            .collect();

        Ok(predicted)
    }

    /// Fit the CRF model to training data.
    pub fn fit(&mut self, genes: &[Gene], shuffle: bool) -> Result<()> {
        let model = self
            .model
            .as_mut()
            .ok_or_else(|| anyhow::anyhow!("No CRF model backend set"))?;

        let extract_features: fn(&[Gene]) -> Vec<HashMap<String, bool>> =
            match self.feature_type.as_str() {
                "protein" => extract_features_protein,
                "domain" => extract_features_domain,
                _ => return Err(anyhow::anyhow!("invalid feature type")),
            };
        let extract_labels: fn(&[Gene]) -> Vec<String> = match self.feature_type.as_str() {
            "protein" => features::extract_labels_protein,
            "domain" => features::extract_labels_domain,
            _ => return Err(anyhow::anyhow!("invalid feature type")),
        };

        // Sort and group genes by sequence
        let mut genes = genes.to_vec();
        genes.sort_by(|a, b| a.source_id.cmp(&b.source_id));
        for gene in &mut genes {
            gene.protein
                .domains
                .sort_by(|a, b| a.start.cmp(&b.start));
        }

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

        if shuffle {
            use rand::seq::SliceRandom;
            use rand::thread_rng;
            let mut rng = thread_rng();
            sequences.shuffle(&mut rng);
        }

        // Build training instances with sliding window
        let mut training_features = Vec::new();
        let mut training_labels = Vec::new();

        for sequence in &sequences {
            let feats = extract_features(sequence);
            let labels = extract_labels(sequence);
            if feats.len() < self.window_size {
                continue;
            }
            for win in sliding_window(feats.len(), self.window_size, self.window_step) {
                training_features.push(feats[win.clone()].to_vec());
                training_labels.push(labels[win].to_vec());
            }
        }

        model.fit(&training_features, &training_labels)
    }
}
