//! Supervised classifier to predict the type of a gene cluster.

pub mod backend;

use std::collections::BTreeSet;

use anyhow::Result;

use crate::model::{Cluster, ClusterType};

/// A binarizer that converts between ClusterType and binary vectors.
pub struct TypeBinarizer {
    pub classes: Vec<String>,
}

impl TypeBinarizer {
    pub fn new(classes: Vec<String>) -> Self {
        Self { classes }
    }

    /// Convert a list of ClusterType into a binary matrix.
    pub fn transform(&self, types: &[ClusterType]) -> Vec<Vec<f64>> {
        types
            .iter()
            .map(|ct| {
                self.classes
                    .iter()
                    .map(|cls| if ct.names.contains(cls) { 1.0 } else { 0.0 })
                    .collect()
            })
            .collect()
    }

    /// Convert a binary matrix back to ClusterType instances.
    pub fn inverse_transform(&self, matrix: &[Vec<bool>]) -> Vec<ClusterType> {
        matrix
            .iter()
            .map(|row| {
                let names: BTreeSet<String> = row
                    .iter()
                    .zip(self.classes.iter())
                    .filter(|(&flag, _)| flag)
                    .map(|(_, cls)| cls.clone())
                    .collect();
                ClusterType { names }
            })
            .collect()
    }
}

/// Trait for a random forest backend used for type classification.
pub trait RandomForestModel: Send + Sync {
    /// Fit the model on training data.
    fn fit(&mut self, x: &[Vec<f64>], y: &[Vec<f64>]) -> Result<()>;
    /// Predict class probabilities for each sample.
    /// Returns shape [n_samples][n_classes] probabilities.
    fn predict_proba(&self, x: &[Vec<f64>]) -> Result<Vec<Vec<f64>>>;
}

/// Classifier that predicts cluster biosynthetic types.
pub struct TypeClassifier {
    pub binarizer: TypeBinarizer,
    pub domains: Vec<String>,
    model: Option<Box<dyn RandomForestModel>>,
}

impl TypeClassifier {
    pub fn new(classes: Vec<String>) -> Self {
        Self {
            binarizer: TypeBinarizer::new(classes),
            domains: Vec::new(),
            model: None,
        }
    }

    pub fn set_model(&mut self, model: Box<dyn RandomForestModel>) {
        self.model = Some(model);
    }

    pub fn set_domains(&mut self, domains: Vec<String>) {
        self.domains = domains;
    }

    /// Predict types for given clusters, mutating them in place.
    pub fn predict_types(&self, clusters: &mut [Cluster]) -> Result<()> {
        let model = self
            .model
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("TypeClassifier model not set"))?;

        // Compute domain compositions
        let comps: Vec<Vec<f64>> = clusters
            .iter()
            .map(|c| c.domain_composition(Some(&self.domains), true, false, true))
            .collect();

        // Predict probabilities
        let probas = model.predict_proba(&comps)?;

        // Assign types to clusters
        for (cluster, proba) in clusters.iter_mut().zip(probas.iter()) {
            // Threshold at 0.5
            let type_flags: Vec<bool> = proba.iter().map(|&p| p > 0.5).collect();
            let types = self.binarizer.inverse_transform(&[type_flags]);
            cluster.cluster_type = types.into_iter().next();
            cluster.type_probabilities = self
                .binarizer
                .classes
                .iter()
                .zip(proba.iter())
                .map(|(cls, &p)| (cls.clone(), p))
                .collect();
        }

        Ok(())
    }
}
