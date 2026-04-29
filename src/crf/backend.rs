//! CRFsuite backend using the `crfsuite_compliant_rs` crate.
//!
//! Uses `Crf1dTagger::marginal_point()` for marginal inference and
//! `ModelReader` for model loading and feature extraction.

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};

use crfsuite_compliant_rs::crf1d::tag::Crf1dTagger;
use crfsuite_compliant_rs::model::ModelReader;
use crfsuite_compliant_rs::types::{Attribute, Instance, Item};

use super::CrfModel;

/// CRF model backed by the `crfsuite_compliant_rs` crate.
#[derive(Clone)]
pub struct CrfSuiteModel {
    /// Raw model file bytes (ModelReader borrows from this).
    model_data: Vec<u8>,
    /// Label names in order of their IDs.
    labels: Vec<String>,
    /// Cached state features as (feature_name, label_name) -> weight.
    state_features_cache: HashMap<(String, String), f64>,
}

impl CrfSuiteModel {
    /// Load a CRFsuite model from a file path.
    pub fn from_file(path: &Path) -> Result<Self> {
        let data = std::fs::read(path)
            .with_context(|| format!("reading CRF model from {}", path.display()))?;
        Self::from_bytes(data)
    }

    /// Load a CRFsuite model from raw bytes.
    pub fn from_bytes(data: Vec<u8>) -> Result<Self> {
        let model =
            ModelReader::open(&data).ok_or_else(|| anyhow::anyhow!("failed to parse CRF model"))?;

        let num_labels = model.num_labels();
        let num_features = model.num_features();

        // Extract label names
        let labels: Vec<String> = (0..num_labels as i32)
            .map(|i| model.to_label(i).unwrap_or("?").to_string())
            .collect();

        // Extract state features directly from model features.
        // ftype 0 = state feature (attr -> label), ftype 1 = transition (label -> label)
        let mut state_features_cache: HashMap<(String, String), f64> = HashMap::new();

        for fid in 0..num_features {
            if let Some(feat) = model.get_feature(fid) {
                if feat.ftype == 0 {
                    if let (Some(attr_name), Some(label_name)) = (
                        model.to_attr(feat.src as i32),
                        model.to_label(feat.dst as i32),
                    ) {
                        state_features_cache
                            .insert((attr_name.to_string(), label_name.to_string()), feat.weight);
                    }
                }
            }
        }

        Ok(Self {
            model_data: data,
            labels,
            state_features_cache,
        })
    }

    /// Create an empty model (placeholder for training).
    pub fn empty() -> Self {
        Self {
            model_data: Vec::new(),
            labels: Vec::new(),
            state_features_cache: HashMap::new(),
        }
    }

    /// Save the model to a file.
    pub fn save(&self, path: &Path) -> Result<()> {
        std::fs::write(path, &self.model_data)
            .with_context(|| format!("saving CRF model to {}", path.display()))
    }

    /// Number of labels.
    pub fn num_labels(&self) -> usize {
        self.labels.len()
    }
}

impl CrfModel for CrfSuiteModel {
    fn predict_marginals_single(&self, features: &[Vec<String>]) -> Vec<HashMap<String, f64>> {
        let model = ModelReader::open(&self.model_data).unwrap();
        let mut tagger = Crf1dTagger::new(&model);

        let num_labels = model.num_labels();

        // Build items from features, converting string attribute names to IDs
        let items: Vec<Item> = features
            .iter()
            .map(|feat_map| {
                let attrs: Vec<Attribute> = feat_map
                    .iter()
                    .filter_map(|name| model.to_aid(name).map(|aid| Attribute { aid, value: 1.0 }))
                    .collect();
                Item { contents: attrs }
            })
            .collect();

        let instance = Instance {
            items,
            labels: vec![0i32; features.len()],
            weight: 1.0,
            group: 0,
        };

        tagger.set(&instance);

        // Get marginal probabilities for each position and label
        let mut result = Vec::with_capacity(features.len());
        for t in 0..features.len() {
            let mut probs = HashMap::new();
            for lid in 0..num_labels as i32 {
                let label = model.to_label(lid).unwrap_or("?");
                let prob = tagger.marginal_point(lid, t as i32);
                probs.insert(label.to_string(), prob);
            }
            result.push(probs);
        }

        result
    }

    fn state_features(&self) -> &HashMap<(String, String), f64> {
        &self.state_features_cache
    }

    fn fit(&mut self, x: &[Vec<Vec<String>>], y: &[Vec<String>]) -> Result<()> {
        use crfsuite_compliant_rs::crf1d::encode::Crf1dEncoder;
        use crfsuite_compliant_rs::quark::Quark;
        use crfsuite_compliant_rs::train::lbfgs::train_lbfgs;

        // Build string-to-ID mappings
        let mut label_quark = Quark::new();
        let mut attr_quark = Quark::new();

        for (xseq, yseq) in x.iter().zip(y.iter()) {
            for label in yseq {
                label_quark.get(label);
            }
            for feat_map in xseq {
                for attr in feat_map {
                    attr_quark.get(attr);
                }
            }
        }

        // Build instances
        let instances: Vec<Instance> = x
            .iter()
            .zip(y.iter())
            .map(|(xseq, yseq)| {
                let items: Vec<Item> = xseq
                    .iter()
                    .map(|feat_map| {
                        let attrs: Vec<Attribute> = feat_map
                            .iter()
                            .map(|name| {
                                let aid = attr_quark.to_id(name).unwrap();
                                Attribute { aid, value: 1.0 }
                            })
                            .collect();
                        Item { contents: attrs }
                    })
                    .collect();
                let labels: Vec<i32> = yseq.iter().map(|l| label_quark.to_id(l).unwrap()).collect();
                Instance {
                    items,
                    labels,
                    weight: 1.0,
                    group: 0,
                }
            })
            .collect();

        let num_labels = label_quark.num();
        let num_attrs = attr_quark.num();

        // Create encoder (generates features internally)
        let mut encoder = Crf1dEncoder::new(
            &instances, num_labels, num_attrs, 0.0,   // min_freq
            false, // possible_states (connect_all_attrs)
            false, // possible_transitions (connect_all_edges)
        );

        let mut log_fn: crfsuite_compliant_rs::train::LogFn = Box::new(|_| {});
        let weights = train_lbfgs(
            &mut encoder,
            &instances,
            1.0,           // c1 (L1 regularization)
            1.0,           // c2 (L2 regularization)
            1000,          // max_iterations
            6,             // num_memories (L-BFGS)
            1e-5,          // epsilon
            10,            // period
            1e-5,          // delta
            "MoreThuente", // linesearch
            20,            // max_linesearch
            &mut log_fn,
        );

        // Get label and attribute strings
        let label_strings: Vec<String> = (0..num_labels as i32)
            .map(|i| label_quark.to_string(i).unwrap().to_string())
            .collect();
        let attr_strings: Vec<String> = (0..num_attrs as i32)
            .map(|i| attr_quark.to_string(i).unwrap().to_string())
            .collect();

        let model_data = encoder.save_model(&weights, &label_strings, &attr_strings);

        let new_model = Self::from_bytes(model_data)?;
        *self = new_model;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_model() {
        let model = CrfSuiteModel::empty();
        assert_eq!(model.num_labels(), 0);
        assert!(model.state_features().is_empty());
    }
}
