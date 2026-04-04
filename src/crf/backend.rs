//! CRFsuite backend using the `crfs` crate.
//!
//! Provides training via `crfs::Trainer` and marginal inference via a
//! custom forward-backward implementation (since `crfs` only exposes Viterbi).

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};

use super::CrfModel;

/// CRF model backed by the `crfs` crate.
///
/// Stores the raw model bytes (crfs::Model borrows them) plus extracted
/// weights for our own forward-backward marginal computation.
pub struct CrfSuiteModel {
    /// Raw model file bytes (crfs::Model borrows from this).
    model_data: Vec<u8>,
    /// Transition weights: trans\[i * L + j\] = weight for label i → label j.
    trans: Vec<f64>,
    /// State feature weights: (attr_name, label_index) → weight.
    state_weights: HashMap<(String, usize), f64>,
    /// Label names in order of their IDs.
    labels: Vec<String>,
    /// Cached state features as (feature_name, label_name) → weight.
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
        let model = crfs::Model::new(&data)
            .map_err(|e| anyhow::anyhow!("parsing CRF model: {}", e))?;

        let num_labels = model.num_labels() as usize;
        let _num_attrs = model.num_attrs() as usize;

        // Extract label names
        let labels: Vec<String> = (0..num_labels)
            .map(|i| {
                model
                    .to_label(i as u32)
                    .unwrap_or("?")
                    .to_string()
            })
            .collect();

        // Extract weights by parsing dump output.
        // This is the only public way to access feature weights.
        let mut dump_buf = Vec::new();
        model
            .dump(&mut dump_buf)
            .map_err(|e| anyhow::anyhow!("dumping CRF model: {}", e))?;
        let dump_str = String::from_utf8_lossy(&dump_buf);

        let mut trans = vec![0.0f64; num_labels * num_labels];
        let mut state_weights: HashMap<(String, usize), f64> = HashMap::new();
        let mut state_features_cache: HashMap<(String, String), f64> = HashMap::new();

        let label_to_id: HashMap<&str, usize> = labels
            .iter()
            .enumerate()
            .map(|(i, l)| (l.as_str(), i))
            .collect();

        let mut section = "";
        for line in dump_str.lines() {
            let line = line.trim();
            if line.starts_with("TRANSITIONS") {
                section = "trans";
                continue;
            } else if line.starts_with("STATE_FEATURES") {
                section = "state";
                continue;
            } else if line == "}" || line.is_empty() || line.starts_with("FILEHEADER")
                || line.starts_with("LABELS") || line.starts_with("ATTRIBUTES")
                || !line.starts_with('(')
            {
                if line.starts_with("FILEHEADER") || line.starts_with("LABELS")
                    || line.starts_with("ATTRIBUTES")
                {
                    section = "";
                }
                continue;
            }

            // Parse lines like:  (1) sunny --> rainy: 0.008212
            // or:                 (0) walk --> sunny: 0.443627
            if let Some(rest) = line.strip_prefix('(') {
                if let Some(after_paren) = rest.find(')') {
                    let content = &rest[after_paren + 1..].trim();
                    // Parse "source --> target: weight"
                    if let Some(arrow_pos) = content.find(" --> ") {
                        let source = &content[..arrow_pos];
                        let after_arrow = &content[arrow_pos + 5..];
                        if let Some(colon_pos) = after_arrow.find(": ") {
                            let target = &after_arrow[..colon_pos];
                            let weight_str = &after_arrow[colon_pos + 2..];
                            if let Ok(weight) = weight_str.parse::<f64>() {
                                match section {
                                    "trans" => {
                                        if let (Some(&si), Some(&ti)) =
                                            (label_to_id.get(source), label_to_id.get(target))
                                        {
                                            trans[si * num_labels + ti] = weight;
                                        }
                                    }
                                    "state" => {
                                        if let Some(&ti) = label_to_id.get(target) {
                                            state_weights
                                                .insert((source.to_string(), ti), weight);
                                            state_features_cache.insert(
                                                (source.to_string(), target.to_string()),
                                                weight,
                                            );
                                        }
                                    }
                                    _ => {}
                                }
                            }
                        }
                    }
                }
            }
        }

        Ok(Self {
            model_data: data,
            trans,
            state_weights,
            labels,
            state_features_cache,
        })
    }

    /// Train a new model and return it.
    pub fn train(
        x: &[Vec<HashMap<String, bool>>],
        y: &[Vec<String>],
    ) -> Result<Self> {
        let mut trainer = crfs::Trainer::lbfgs();

        for (xseq, yseq) in x.iter().zip(y.iter()) {
            let attrs: Vec<Vec<crfs::Attribute>> = xseq
                .iter()
                .map(|feat_map| {
                    feat_map
                        .iter()
                        .filter(|(_, &v)| v)
                        .map(|(k, _)| crfs::Attribute::new(k.as_str(), 1.0))
                        .collect()
                })
                .collect();
            let labels: Vec<&str> = yseq.iter().map(|s| s.as_str()).collect();
            trainer
                .append(&attrs, &labels)
                .map_err(|e| anyhow::anyhow!("appending training instance: {}", e))?;
        }

        // Train to a temporary file
        let tmp = tempfile_path();
        trainer
            .train(Path::new(&tmp))
            .map_err(|e| anyhow::anyhow!("training CRF model: {}", e))?;

        let data = std::fs::read(&tmp)
            .with_context(|| format!("reading trained model from {}", tmp))?;
        let _ = std::fs::remove_file(&tmp);

        Self::from_bytes(data)
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

    /// Compute marginal probabilities using forward-backward algorithm.
    fn forward_backward(
        &self,
        state_scores: &[Vec<f64>],
    ) -> Vec<Vec<f64>> {
        let t_len = state_scores.len();
        let l = self.labels.len();

        if t_len == 0 || l == 0 {
            return Vec::new();
        }

        // Forward pass (in log-space for numerical stability)
        // alpha[t][j] = log sum_i exp(alpha[t-1][i] + trans[i][j] + state[t][j])
        let mut alpha = vec![vec![0.0f64; l]; t_len];

        // t=0: alpha[0][j] = state_scores[0][j]
        alpha[0].copy_from_slice(&state_scores[0]);

        for t in 1..t_len {
            for j in 0..l {
                // Collect: alpha[t-1][i] + trans[i][j] for all i
                let mut terms = Vec::with_capacity(l);
                for i in 0..l {
                    terms.push(alpha[t - 1][i] + self.trans[i * l + j]);
                }
                alpha[t][j] = log_sum_exp(&terms) + state_scores[t][j];
            }
        }

        // Backward pass
        // beta[t][i] = log sum_j exp(trans[i][j] + state[t+1][j] + beta[t+1][j])
        let mut beta = vec![vec![0.0f64; l]; t_len];
        // beta[T-1][i] = 0 (log(1))

        for t in (0..t_len - 1).rev() {
            for i in 0..l {
                let mut terms = Vec::with_capacity(l);
                for j in 0..l {
                    terms.push(
                        self.trans[i * l + j]
                            + state_scores[t + 1][j]
                            + beta[t + 1][j],
                    );
                }
                beta[t][i] = log_sum_exp(&terms);
            }
        }

        // Log-normalization constant: Z = log sum_j exp(alpha[T-1][j])
        let last_alpha: Vec<f64> = alpha[t_len - 1].clone();
        let log_z = log_sum_exp(&last_alpha);

        // Marginals: P(y_t = j) = exp(alpha[t][j] + beta[t][j] - log_z)
        let mut marginals = vec![vec![0.0f64; l]; t_len];
        for t in 0..t_len {
            for j in 0..l {
                marginals[t][j] = (alpha[t][j] + beta[t][j] - log_z).exp();
            }
            // Normalize to handle numerical imprecision
            let sum: f64 = marginals[t].iter().sum();
            if sum > 0.0 {
                for j in 0..l {
                    marginals[t][j] /= sum;
                }
            }
        }

        marginals
    }

    /// Compute state scores for a sequence of feature dicts.
    fn compute_state_scores(&self, features: &[HashMap<String, bool>]) -> Vec<Vec<f64>> {
        let l = self.labels.len();
        features
            .iter()
            .map(|feat_map| {
                let mut scores = vec![0.0f64; l];
                for (attr_name, &present) in feat_map {
                    if !present {
                        continue;
                    }
                    for label_idx in 0..l {
                        if let Some(&w) = self.state_weights.get(&(attr_name.clone(), label_idx)) {
                            scores[label_idx] += w;
                        }
                    }
                }
                scores
            })
            .collect()
    }
}

impl CrfModel for CrfSuiteModel {
    fn predict_marginals_single(
        &self,
        features: &[HashMap<String, bool>],
    ) -> Vec<HashMap<String, f64>> {
        let state_scores = self.compute_state_scores(features);
        let marginals = self.forward_backward(&state_scores);

        marginals
            .into_iter()
            .map(|probs| {
                self.labels
                    .iter()
                    .enumerate()
                    .map(|(i, label)| (label.clone(), probs[i]))
                    .collect()
            })
            .collect()
    }

    fn state_features(&self) -> &HashMap<(String, String), f64> {
        &self.state_features_cache
    }

    fn fit(
        &mut self,
        x: &[Vec<HashMap<String, bool>>],
        y: &[Vec<String>],
    ) -> Result<()> {
        let new_model = CrfSuiteModel::train(x, y)?;
        *self = new_model;
        Ok(())
    }
}

/// Numerically stable log-sum-exp.
fn log_sum_exp(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NEG_INFINITY;
    }
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if max == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    let sum: f64 = values.iter().map(|&v| (v - max).exp()).sum();
    max + sum.ln()
}

/// Generate a temporary file path for training output.
fn tempfile_path() -> String {
    use std::time::SystemTime;
    let ts = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    format!("/tmp/gecco_crf_{}.crfsuite", ts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_log_sum_exp() {
        // log(exp(1) + exp(2)) = log(e + e^2) ≈ 2.3133
        let result = log_sum_exp(&[1.0, 2.0]);
        assert!((result - 2.3133).abs() < 0.001);

        // Single value
        assert!((log_sum_exp(&[5.0]) - 5.0).abs() < 1e-10);

        // Empty
        assert!(log_sum_exp(&[]).is_infinite());
    }

    #[test]
    fn test_forward_backward_uniform() {
        // With zero weights everywhere, marginals should be uniform
        let model = CrfSuiteModel {
            model_data: Vec::new(),
            trans: vec![0.0; 4], // 2x2
            state_weights: HashMap::new(),
            labels: vec!["0".to_string(), "1".to_string()],
            state_features_cache: HashMap::new(),
        };

        let state_scores = vec![vec![0.0, 0.0]; 3]; // 3 positions, 2 labels
        let marginals = model.forward_backward(&state_scores);

        assert_eq!(marginals.len(), 3);
        for t in 0..3 {
            assert!((marginals[t][0] - 0.5).abs() < 0.01);
            assert!((marginals[t][1] - 0.5).abs() < 0.01);
        }
    }

    #[test]
    fn test_forward_backward_biased() {
        // Strong bias toward label 1 via state scores
        let model = CrfSuiteModel {
            model_data: Vec::new(),
            trans: vec![0.0; 4],
            state_weights: HashMap::new(),
            labels: vec!["0".to_string(), "1".to_string()],
            state_features_cache: HashMap::new(),
        };

        // State scores strongly favor label 1
        let state_scores = vec![vec![-10.0, 10.0]; 3];
        let marginals = model.forward_backward(&state_scores);

        for t in 0..3 {
            assert!(marginals[t][1] > 0.99);
            assert!(marginals[t][0] < 0.01);
        }
    }

    #[test]
    fn test_predict_marginals() {
        let mut state_weights = HashMap::new();
        // Feature "PF00001" strongly associated with label "1"
        state_weights.insert(("PF00001".to_string(), 1), 5.0);
        state_weights.insert(("PF00001".to_string(), 0), -5.0);

        let mut state_features_cache = HashMap::new();
        state_features_cache.insert(("PF00001".to_string(), "1".to_string()), 5.0);
        state_features_cache.insert(("PF00001".to_string(), "0".to_string()), -5.0);

        let model = CrfSuiteModel {
            model_data: Vec::new(),
            trans: vec![0.0; 4],
            state_weights,
            labels: vec!["0".to_string(), "1".to_string()],
            state_features_cache,
        };

        let features = vec![
            {
                let mut m = HashMap::new();
                m.insert("PF00001".to_string(), true);
                m
            },
            HashMap::new(), // unannotated gene
        ];

        let result = model.predict_marginals_single(&features);
        assert_eq!(result.len(), 2);
        // First position should strongly favor "1"
        assert!(result[0]["1"] > 0.9);
        // Second position (no features) should be closer to uniform
        // but influenced by transitions
    }
}
