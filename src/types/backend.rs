//! Random Forest backend for BGC type classification.
//!
//! Uses `smartcore_proba` which provides `predict_proba` on Random Forest.
//! GECCO's type classifier uses N binary classifiers (one per BGC type),
//! each predicting P(positive) for that type.

use anyhow::Result;

use smartcore_proba::ensemble::random_forest_classifier::{
    RandomForestClassifier, RandomForestClassifierParameters,
};
use smartcore_proba::linalg::basic::arrays::Array;
use smartcore_proba::linalg::basic::matrix::DenseMatrix;

use super::RandomForestModel;

/// Random Forest model for multi-label type classification.
///
/// Wraps N binary `RandomForestClassifier`s, one per BGC type class.
/// Each classifier predicts P(positive) for its class.
pub struct SmartcoreRF {
    /// One binary classifier per class. None if not yet trained.
    classifiers: Vec<Option<RandomForestClassifier<f64, i64, DenseMatrix<f64>, Vec<i64>>>>,
    /// Number of output classes.
    n_classes: usize,
    /// Number of trees per classifier.
    n_trees: u16,
    /// Random seed for reproducibility.
    seed: u64,
}

impl SmartcoreRF {
    /// Create a new Random Forest model for `n_classes` binary classifiers.
    pub fn new(n_classes: usize) -> Self {
        let classifiers = (0..n_classes).map(|_| None).collect();
        Self {
            classifiers,
            n_classes,
            n_trees: 100,
            seed: 0,
        }
    }

    pub fn with_n_trees(mut self, n_trees: u16) -> Self {
        self.n_trees = n_trees;
        self
    }

    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = seed;
        self
    }
}

impl RandomForestModel for SmartcoreRF {
    fn fit(&mut self, x: &[Vec<f64>], y: &[Vec<f64>]) -> Result<()> {
        if x.is_empty() {
            return Err(anyhow::anyhow!("empty training data"));
        }

        let n_samples = x.len();
        let _n_features = x[0].len();

        // y shape: [n_samples][n_classes]
        if y.is_empty() || y[0].len() != self.n_classes {
            return Err(anyhow::anyhow!(
                "y must have shape [n_samples][{}], got [{}][{}]",
                self.n_classes,
                y.len(),
                y.first().map_or(0, |v| v.len())
            ));
        }

        // Build feature matrix once
        let x_mat = DenseMatrix::from_2d_vec(&x.to_vec())
            .map_err(|e| anyhow::anyhow!("building feature matrix: {}", e))?;

        let params = RandomForestClassifierParameters::default()
            .with_n_trees(self.n_trees)
            .with_seed(self.seed);

        // Train one binary classifier per class
        self.classifiers = Vec::with_capacity(self.n_classes);
        for class_idx in 0..self.n_classes {
            // Extract binary labels for this class
            let labels: Vec<i64> = y.iter().map(|row| row[class_idx] as i64).collect();

            // Check if this class has both positive and negative samples
            let n_pos = labels.iter().filter(|&&l| l == 1).count();
            if n_pos == 0 || n_pos == n_samples {
                // All same class — no point training, store None
                self.classifiers.push(None);
                continue;
            }

            let clf = RandomForestClassifier::fit(&x_mat, &labels, params.clone())
                .map_err(|e| anyhow::anyhow!("training RF for class {}: {}", class_idx, e))?;
            self.classifiers.push(Some(clf));
        }

        Ok(())
    }

    fn predict_proba(&self, x: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
        if x.is_empty() {
            return Ok(Vec::new());
        }

        let n_samples = x.len();
        let x_mat = DenseMatrix::from_2d_vec(&x.to_vec())
            .map_err(|e| anyhow::anyhow!("building feature matrix: {}", e))?;

        // Initialize result: [n_samples][n_classes] all zeros
        let mut result = vec![vec![0.0f64; self.n_classes]; n_samples];

        for (class_idx, clf_opt) in self.classifiers.iter().enumerate() {
            match clf_opt {
                Some(clf) => {
                    let probas: DenseMatrix<f64> = clf
                        .predict_proba(&x_mat)
                        .map_err(|e| {
                            anyhow::anyhow!("predicting probas for class {}: {}", class_idx, e)
                        })?;

                    // probas shape: [n_samples][n_unique_classes]
                    // For binary: column 0 = P(negative), column 1 = P(positive)
                    // We want P(positive)
                    let (n_rows, n_cols) = probas.shape();
                    let positive_col = if n_cols >= 2 { 1 } else { 0 };

                    for i in 0..n_rows.min(n_samples) {
                        result[i][class_idx] = *probas.get((i, positive_col));
                    }
                }
                None => {
                    // No classifier trained (single-class data) — leave at 0.0
                }
            }
        }

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_smartcore_rf_fit_predict() {
        // Simple 2-class, 3-feature problem
        let x = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
            vec![1.0, 1.0, 0.0],
            vec![0.0, 1.0, 1.0],
            vec![1.0, 0.0, 1.0],
        ];
        // 2 classes: class 0 = feature[0] > 0, class 1 = feature[2] > 0
        let y = vec![
            vec![1.0, 0.0], // sample 0: class0=yes, class1=no
            vec![0.0, 0.0], // sample 1: both no
            vec![0.0, 1.0], // sample 2: class0=no, class1=yes
            vec![1.0, 0.0], // sample 3: class0=yes, class1=no
            vec![0.0, 1.0], // sample 4: class0=no, class1=yes
            vec![1.0, 1.0], // sample 5: both yes
        ];

        let mut model = SmartcoreRF::new(2).with_n_trees(10).with_seed(42);
        model.fit(&x, &y).unwrap();

        let probas = model.predict_proba(&x).unwrap();
        assert_eq!(probas.len(), 6);
        assert_eq!(probas[0].len(), 2);

        // Probabilities should be between 0 and 1
        for row in &probas {
            for &p in row {
                assert!(p >= 0.0 && p <= 1.0, "probability out of range: {}", p);
            }
        }
    }

    #[test]
    fn test_smartcore_rf_empty_input() {
        let model = SmartcoreRF::new(3);
        let result = model.predict_proba(&[]).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_smartcore_rf_single_class() {
        // When all samples are positive for a class, that classifier is skipped
        let x = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 1.0],
        ];
        let y = vec![
            vec![1.0, 1.0], // all positive for both classes
            vec![1.0, 0.0],
            vec![1.0, 1.0],
        ];

        let mut model = SmartcoreRF::new(2).with_n_trees(10).with_seed(42);
        model.fit(&x, &y).unwrap();

        let probas = model.predict_proba(&x).unwrap();
        // Class 0 has all-positive labels → no classifier → probability = 0.0
        for row in &probas {
            assert_eq!(row[0], 0.0);
        }
    }
}
