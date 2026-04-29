//! Evaluator and small trainer for scikit-learn-shaped random forests.
//!
//! The model representation mirrors a fitted
//! `sklearn.ensemble.RandomForestClassifier` exported as tree arrays:
//! `children_left`, `children_right`, `feature`, `threshold`, and `value`.

use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Errors returned while loading or evaluating an exported forest.
#[derive(Debug, Error)]
pub enum Error {
    #[error("invalid model: {0}")]
    InvalidModel(String),
    #[error("invalid training data: {0}")]
    InvalidTrainingData(String),
    #[error("expected {expected} features, got {actual}")]
    FeatureCount { expected: usize, actual: usize },
    #[error(transparent)]
    Json(#[from] serde_json::Error),
}

/// Result alias for sklearn RF operations.
pub type Result<T> = std::result::Result<T, Error>;

/// A prediction-only representation of a fitted sklearn RandomForestClassifier.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RandomForestClassifier {
    pub n_features_in: usize,
    pub outputs: Vec<OutputClasses>,
    pub estimators: Vec<DecisionTreeClassifier>,
}

/// Class labels for one sklearn output.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputClasses {
    pub classes: Vec<f64>,
}

/// A prediction-only representation of one sklearn DecisionTreeClassifier.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DecisionTreeClassifier {
    pub children_left: Vec<i64>,
    pub children_right: Vec<i64>,
    pub feature: Vec<i64>,
    pub threshold: Vec<f64>,
    /// Leaf values with shape `[node][output][class]`.
    ///
    /// In scikit-learn 1.5.2 this is the probability tensor returned by
    /// `tree_.predict`, so inference uses it directly instead of renormalizing.
    pub value: Vec<Vec<Vec<f64>>>,
}

/// Training parameters for the supported RandomForestClassifier subset.
#[derive(Debug, Clone)]
pub struct RandomForestClassifierParameters {
    pub n_estimators: usize,
    pub max_features: MaxFeatures,
    pub bootstrap: bool,
    pub min_samples_split: usize,
    pub min_samples_leaf: usize,
    pub max_depth: Option<usize>,
    pub random_state: u64,
}

impl Default for RandomForestClassifierParameters {
    fn default() -> Self {
        Self {
            n_estimators: 100,
            max_features: MaxFeatures::Sqrt,
            bootstrap: true,
            min_samples_split: 2,
            min_samples_leaf: 1,
            max_depth: None,
            random_state: 0,
        }
    }
}

/// Feature-subset strategy at each split.
#[derive(Debug, Clone, Copy)]
pub enum MaxFeatures {
    All,
    Sqrt,
    Count(usize),
}

impl RandomForestClassifier {
    /// Train a RandomForestClassifier-like model.
    ///
    /// This implements the GECCO-relevant sklearn subset: dense numeric input,
    /// multi-output classification, Gini splits, bootstrap sampling, and
    /// `max_features="sqrt"` by default.
    pub fn fit(
        x: &[Vec<f64>],
        y: &[Vec<f64>],
        params: RandomForestClassifierParameters,
    ) -> Result<Self> {
        validate_training_data(x, y)?;
        if params.n_estimators == 0 {
            return Err(Error::InvalidTrainingData(
                "n_estimators must be positive".to_string(),
            ));
        }

        let n_features_in = x[0].len();
        let outputs = infer_outputs(y);
        let encoded_y = encode_y(y, &outputs)?;
        let sorted_features = sorted_feature_samples(x);
        let mut forest_rng = Mt19937::new(params.random_state as u32);
        let tree_seeds = (0..params.n_estimators)
            .map(|_| forest_rng.gen_interval(i32::MAX as u32 - 1))
            .collect::<Vec<_>>();
        let estimators = tree_seeds
            .into_par_iter()
            .map(|tree_seed| {
                build_tree(
                    x,
                    &encoded_y,
                    &outputs,
                    &params,
                    &sorted_features,
                    n_features_in,
                    tree_seed,
                )
            })
            .collect::<Result<Vec<_>>>()?;

        Self::new(Self {
            n_features_in,
            outputs,
            estimators,
        })
    }

    /// Load an exported forest from JSON bytes.
    pub fn from_json_slice(bytes: &[u8]) -> Result<Self> {
        let model = serde_json::from_slice(bytes)?;
        Self::new(model)
    }

    /// Serialize the forest to compact JSON bytes.
    pub fn to_json_vec(&self) -> Result<Vec<u8>> {
        serde_json::to_vec(self).map_err(Error::from)
    }

    /// Validate a deserialized forest.
    pub fn new(model: Self) -> Result<Self> {
        if model.outputs.is_empty() {
            return Err(Error::InvalidModel("forest has no outputs".to_string()));
        }
        if model.estimators.is_empty() {
            return Err(Error::InvalidModel("forest has no estimators".to_string()));
        }
        for (output_idx, output) in model.outputs.iter().enumerate() {
            if output.classes.is_empty() {
                return Err(Error::InvalidModel(format!(
                    "output {output_idx} has no classes"
                )));
            }
        }
        for (tree_idx, tree) in model.estimators.iter().enumerate() {
            tree.validate(model.outputs.len())
                .map_err(|e| Error::InvalidModel(format!("estimator {tree_idx}: {e}")))?;
        }
        Ok(model)
    }

    /// Predict class probabilities.
    ///
    /// The returned shape is `[sample][output][class]`, matching sklearn's
    /// multi-output `predict_proba` content with the sample axis moved first.
    pub fn predict_proba(&self, samples: &[Vec<f64>]) -> Result<Vec<Vec<Vec<f64>>>> {
        let mut result = Vec::with_capacity(samples.len());
        for sample in samples {
            result.push(self.predict_one(sample)?);
        }
        Ok(result)
    }

    /// Predict probabilities for one sample, returning `[output][class]`.
    pub fn predict_one(&self, sample: &[f64]) -> Result<Vec<Vec<f64>>> {
        if sample.len() != self.n_features_in {
            return Err(Error::FeatureCount {
                expected: self.n_features_in,
                actual: sample.len(),
            });
        }

        let mut totals: Vec<Vec<f64>> = self
            .outputs
            .iter()
            .map(|output| vec![0.0; output.classes.len()])
            .collect();

        for tree in &self.estimators {
            let leaf = tree.leaf_index(sample)?;
            for (output_idx, output) in self.outputs.iter().enumerate() {
                let values = tree
                    .value
                    .get(leaf)
                    .and_then(|node| node.get(output_idx))
                    .ok_or_else(|| {
                        Error::InvalidModel(format!(
                            "leaf {leaf} missing output {output_idx} values"
                        ))
                    })?;
                add_prefix(&mut totals[output_idx], values, output.classes.len());
            }
        }

        let scale = self.estimators.len() as f64;
        for output in &mut totals {
            for p in output {
                *p /= scale;
            }
        }

        Ok(totals)
    }

    /// Predict the probability of a requested class label for each output.
    ///
    /// This is useful for sklearn multi-label classifiers, where each output is
    /// a binary problem and GECCO wants the positive class probability.
    pub fn predict_label_proba(&self, samples: &[Vec<f64>], label: f64) -> Result<Vec<Vec<f64>>> {
        let probas = self.predict_proba(samples)?;
        let label_indices: Vec<Option<usize>> = self
            .outputs
            .iter()
            .map(|output| output.classes.iter().position(|&class| class == label))
            .collect();

        Ok(probas
            .into_iter()
            .map(|sample| {
                sample
                    .into_iter()
                    .enumerate()
                    .map(|(output_idx, output)| {
                        label_indices[output_idx].map_or(0.0, |class_idx| output[class_idx])
                    })
                    .collect()
            })
            .collect())
    }
}

impl DecisionTreeClassifier {
    fn empty() -> Self {
        Self {
            children_left: Vec::new(),
            children_right: Vec::new(),
            feature: Vec::new(),
            threshold: Vec::new(),
            value: Vec::new(),
        }
    }

    fn validate(&self, n_outputs: usize) -> std::result::Result<(), String> {
        let n_nodes = self.children_left.len();
        if n_nodes == 0 {
            return Err("tree has no nodes".to_string());
        }
        if self.children_right.len() != n_nodes
            || self.feature.len() != n_nodes
            || self.threshold.len() != n_nodes
            || self.value.len() != n_nodes
        {
            return Err("tree arrays have different lengths".to_string());
        }
        for node_idx in 0..n_nodes {
            let left = self.children_left[node_idx];
            let right = self.children_right[node_idx];
            if left == -1 && right == -1 {
                continue;
            }
            if left < 0 || right < 0 || left as usize >= n_nodes || right as usize >= n_nodes {
                return Err(format!("node {node_idx} has invalid children"));
            }
        }
        for (node_idx, node_values) in self.value.iter().enumerate() {
            if node_values.len() != n_outputs {
                return Err(format!(
                    "node {node_idx} has {} outputs, expected {n_outputs}",
                    node_values.len()
                ));
            }
        }
        Ok(())
    }

    fn leaf_index(&self, sample: &[f64]) -> Result<usize> {
        let mut node = 0usize;
        loop {
            let left = self.children_left[node];
            let right = self.children_right[node];
            if left == -1 && right == -1 {
                return Ok(node);
            }

            let feature_idx = self.feature[node];
            if feature_idx < 0 {
                return Err(Error::InvalidModel(format!(
                    "split node {node} has invalid feature {feature_idx}"
                )));
            }
            let feature_idx = feature_idx as usize;
            let value = sample.get(feature_idx).ok_or_else(|| Error::FeatureCount {
                expected: feature_idx + 1,
                actual: sample.len(),
            })?;
            let value = *value as f32 as f64;

            node = if value <= self.threshold[node] {
                left as usize
            } else {
                right as usize
            };
        }
    }
}

fn add_prefix(total: &mut [f64], values: &[f64], n_classes: usize) {
    let usable = values.len().min(n_classes);
    for (dst, src) in total.iter_mut().zip(values.iter()).take(usable) {
        *dst += *src;
    }
}

fn validate_training_data(x: &[Vec<f64>], y: &[Vec<f64>]) -> Result<()> {
    if x.is_empty() {
        return Err(Error::InvalidTrainingData("x is empty".to_string()));
    }
    if y.len() != x.len() {
        return Err(Error::InvalidTrainingData(format!(
            "x/y row mismatch: {} vs {}",
            x.len(),
            y.len()
        )));
    }
    let n_features = x[0].len();
    if n_features == 0 {
        return Err(Error::InvalidTrainingData("x has no features".to_string()));
    }
    for (row_idx, row) in x.iter().enumerate() {
        if row.len() != n_features {
            return Err(Error::InvalidTrainingData(format!(
                "x row {row_idx} has {} features, expected {n_features}",
                row.len()
            )));
        }
    }
    let n_outputs = y[0].len();
    if n_outputs == 0 {
        return Err(Error::InvalidTrainingData("y has no outputs".to_string()));
    }
    for (row_idx, row) in y.iter().enumerate() {
        if row.len() != n_outputs {
            return Err(Error::InvalidTrainingData(format!(
                "y row {row_idx} has {} outputs, expected {n_outputs}",
                row.len()
            )));
        }
    }
    Ok(())
}

fn infer_outputs(y: &[Vec<f64>]) -> Vec<OutputClasses> {
    (0..y[0].len())
        .map(|output_idx| {
            let mut classes = y.iter().map(|row| row[output_idx]).collect::<Vec<_>>();
            classes.sort_by(|a, b| a.total_cmp(b));
            classes.dedup();
            OutputClasses { classes }
        })
        .collect()
}

fn encode_y(y: &[Vec<f64>], outputs: &[OutputClasses]) -> Result<Vec<Vec<usize>>> {
    y.iter()
        .enumerate()
        .map(|(row_idx, row)| {
            row.iter()
                .enumerate()
                .map(|(output_idx, value)| {
                    outputs[output_idx]
                        .classes
                        .iter()
                        .position(|class| class == value)
                        .ok_or_else(|| {
                            Error::InvalidTrainingData(format!(
                                "unknown class {value} at row {row_idx}, output {output_idx}"
                            ))
                        })
                })
                .collect()
        })
        .collect()
}

fn sorted_feature_samples(x: &[Vec<f64>]) -> Vec<Vec<usize>> {
    let n_features = x[0].len();
    (0..n_features)
        .map(|feature| {
            let mut samples = (0..x.len()).collect::<Vec<_>>();
            samples.sort_by(|&a, &b| {
                (x[a][feature] as f32 as f64)
                    .total_cmp(&(x[b][feature] as f32 as f64))
                    .then(a.cmp(&b))
            });
            samples
        })
        .collect()
}

struct TreeBuilder<'a> {
    x: &'a [Vec<f64>],
    y: &'a [Vec<usize>],
    outputs: &'a [OutputClasses],
    params: &'a RandomForestClassifierParameters,
    sorted_features: &'a [Vec<usize>],
    rng: Mt19937,
    n_features: usize,
    tree: DecisionTreeClassifier,
}

fn build_tree(
    x: &[Vec<f64>],
    y: &[Vec<usize>],
    outputs: &[OutputClasses],
    params: &RandomForestClassifierParameters,
    sorted_features: &[Vec<usize>],
    n_features: usize,
    tree_seed: u32,
) -> Result<DecisionTreeClassifier> {
    let mut rng = Mt19937::new(tree_seed);
    let samples = if params.bootstrap {
        (0..x.len())
            .map(|_| rng.gen_range(x.len()))
            .collect::<Vec<_>>()
    } else {
        (0..x.len()).collect::<Vec<_>>()
    };
    let mut builder = TreeBuilder {
        x,
        y,
        outputs,
        params,
        sorted_features,
        rng,
        n_features,
        tree: DecisionTreeClassifier::empty(),
    };
    let _ = builder.build_node(&samples, 0)?;
    Ok(builder.tree)
}

impl TreeBuilder<'_> {
    fn build_node(&mut self, samples: &[usize], depth: usize) -> Result<usize> {
        let node_idx = self.tree.children_left.len();
        self.tree.children_left.push(-1);
        self.tree.children_right.push(-1);
        self.tree.feature.push(-2);
        self.tree.threshold.push(-2.0);
        self.tree
            .value
            .push(node_value(samples, self.y, self.outputs));

        if samples.len() < self.params.min_samples_split
            || self.params.max_depth.is_some_and(|max| depth >= max)
            || is_pure(samples, self.y)
        {
            return Ok(node_idx);
        }

        let Some(split) = self.best_split(samples) else {
            return Ok(node_idx);
        };
        if split.left.len() < self.params.min_samples_leaf
            || split.right.len() < self.params.min_samples_leaf
        {
            return Ok(node_idx);
        }

        self.tree.feature[node_idx] = split.feature as i64;
        self.tree.threshold[node_idx] = split.threshold;
        let left_idx = self.build_node(&split.left, depth + 1)?;
        let right_idx = self.build_node(&split.right, depth + 1)?;
        self.tree.children_left[node_idx] = left_idx as i64;
        self.tree.children_right[node_idx] = right_idx as i64;
        Ok(node_idx)
    }

    fn best_split(&mut self, samples: &[usize]) -> Option<Split> {
        let parent_impurity = gini(samples, self.y, self.outputs);
        if parent_impurity <= 0.0 {
            return None;
        }

        let mut best: Option<Split> = None;
        let min_features = self.max_features_count();
        let features = self.permuted_features();
        let mut in_node = vec![0usize; self.x.len()];
        for &sample in samples {
            in_node[sample] += 1;
        }
        for (attempted, feature) in features.into_iter().enumerate() {
            if attempted >= min_features && best.is_some() {
                break;
            }
            let n_pairs = samples.len();
            let mut left_counts = zero_counts(self.outputs);
            let mut right_counts = class_counts(samples, self.y, self.outputs);
            let sorted = &self.sorted_features[feature];
            let mut left_len = 0usize;

            let mut scan = ValueGroupScan {
                sorted,
                in_node: &in_node,
                idx: 0,
                feature,
            };

            while let Some(value) = self.move_value_group(
                &mut scan,
                &mut left_counts,
                &mut right_counts,
                &mut left_len,
            ) {
                let Some(next_value) = self.next_value_in_node(sorted, &in_node, scan.idx, feature)
                else {
                    break;
                };
                if value == next_value {
                    continue;
                }

                let right_len = n_pairs - left_len;
                if left_len < self.params.min_samples_leaf
                    || right_len < self.params.min_samples_leaf
                {
                    continue;
                }

                let weighted = (left_len as f64
                    * gini_from_counts(&left_counts, left_len, self.outputs)
                    + right_len as f64 * gini_from_counts(&right_counts, right_len, self.outputs))
                    / n_pairs as f64;
                let improvement = parent_impurity - weighted;
                if improvement < 0.0 {
                    continue;
                }

                let threshold = midpoint(value, next_value);
                let replace = best
                    .as_ref()
                    .is_none_or(|current| improvement > current.improvement);
                if replace {
                    let (left, right) = self.partition_samples(feature, threshold, &in_node);
                    best = Some(Split {
                        feature,
                        threshold,
                        improvement,
                        left,
                        right,
                    });
                }
            }
        }
        best
    }

    fn move_value_group(
        &self,
        scan: &mut ValueGroupScan<'_>,
        left_counts: &mut [Vec<usize>],
        right_counts: &mut [Vec<usize>],
        left_len: &mut usize,
    ) -> Option<f64> {
        while scan.idx < scan.sorted.len() && scan.in_node[scan.sorted[scan.idx]] == 0 {
            scan.idx += 1;
        }
        if scan.idx == scan.sorted.len() {
            return None;
        }

        let value = self.x[scan.sorted[scan.idx]][scan.feature] as f32 as f64;
        while scan.idx < scan.sorted.len() {
            let sample = scan.sorted[scan.idx];
            let count = scan.in_node[sample];
            if count == 0 {
                scan.idx += 1;
                continue;
            }
            let sample_value = self.x[sample][scan.feature] as f32 as f64;
            if sample_value != value {
                break;
            }
            for output_idx in 0..self.outputs.len() {
                let class_idx = self.y[sample][output_idx];
                left_counts[output_idx][class_idx] += count;
                right_counts[output_idx][class_idx] -= count;
            }
            *left_len += count;
            scan.idx += 1;
        }
        Some(value)
    }

    fn next_value_in_node(
        &self,
        sorted: &[usize],
        in_node: &[usize],
        mut idx: usize,
        feature: usize,
    ) -> Option<f64> {
        while idx < sorted.len() {
            let sample = sorted[idx];
            if in_node[sample] != 0 {
                return Some(self.x[sample][feature] as f32 as f64);
            }
            idx += 1;
        }
        None
    }

    fn partition_samples(
        &self,
        feature: usize,
        threshold: f64,
        in_node: &[usize],
    ) -> (Vec<usize>, Vec<usize>) {
        let mut left = Vec::new();
        let mut right = Vec::new();
        for &sample in &self.sorted_features[feature] {
            let count = in_node[sample];
            if count == 0 {
                continue;
            }
            let dst = if (self.x[sample][feature] as f32 as f64) <= threshold {
                &mut left
            } else {
                &mut right
            };
            dst.extend(std::iter::repeat_n(sample, count));
        }
        (left, right)
    }

    fn max_features_count(&self) -> usize {
        match self.params.max_features {
            MaxFeatures::All => self.n_features,
            MaxFeatures::Sqrt => (self.n_features as f64).sqrt() as usize,
            MaxFeatures::Count(n) => n,
        }
        .clamp(1, self.n_features)
    }

    fn permuted_features(&mut self) -> Vec<usize> {
        let mut features = (0..self.n_features).collect::<Vec<_>>();
        for i in 0..self.n_features {
            let j = i + self.rng.gen_range(self.n_features - i);
            features.swap(i, j);
        }
        features
    }
}

struct Split {
    feature: usize,
    threshold: f64,
    improvement: f64,
    left: Vec<usize>,
    right: Vec<usize>,
}

struct ValueGroupScan<'a> {
    sorted: &'a [usize],
    in_node: &'a [usize],
    idx: usize,
    feature: usize,
}

fn node_value(samples: &[usize], y: &[Vec<usize>], outputs: &[OutputClasses]) -> Vec<Vec<f64>> {
    outputs
        .iter()
        .enumerate()
        .map(|(output_idx, output)| {
            let mut counts = vec![0.0; output.classes.len()];
            for &sample in samples {
                counts[y[sample][output_idx]] += 1.0;
            }
            for count in &mut counts {
                *count /= samples.len() as f64;
            }
            counts
        })
        .collect()
}

fn is_pure(samples: &[usize], y: &[Vec<usize>]) -> bool {
    if let Some((&first, rest)) = samples.split_first() {
        rest.iter().all(|&sample| y[sample] == y[first])
    } else {
        true
    }
}

fn gini(samples: &[usize], y: &[Vec<usize>], outputs: &[OutputClasses]) -> f64 {
    if samples.is_empty() {
        return 0.0;
    }
    let counts = class_counts(samples, y, outputs);
    gini_from_counts(&counts, samples.len(), outputs)
}

fn zero_counts(outputs: &[OutputClasses]) -> Vec<Vec<usize>> {
    outputs
        .iter()
        .map(|output| vec![0; output.classes.len()])
        .collect()
}

fn class_counts(samples: &[usize], y: &[Vec<usize>], outputs: &[OutputClasses]) -> Vec<Vec<usize>> {
    let mut counts = zero_counts(outputs);
    for &sample in samples {
        for (output_idx, output_counts) in counts.iter_mut().enumerate() {
            output_counts[y[sample][output_idx]] += 1;
        }
    }
    counts
}

fn gini_from_counts(counts: &[Vec<usize>], n_samples: usize, outputs: &[OutputClasses]) -> f64 {
    if n_samples == 0 {
        return 0.0;
    }
    let mut total = 0.0;
    for (output_counts, _) in counts.iter().zip(outputs.iter()) {
        let impurity = 1.0
            - output_counts
                .iter()
                .map(|count| {
                    let p = *count as f64 / n_samples as f64;
                    p * p
                })
                .sum::<f64>();
        total += impurity;
    }
    total / outputs.len() as f64
}

fn midpoint(left: f64, right: f64) -> f64 {
    let mid = (left + right) / 2.0;
    if mid <= left || mid == f64::INFINITY {
        left
    } else {
        mid
    }
}

#[derive(Debug, Clone, Copy)]
struct Mt19937 {
    mt: [u32; 624],
    index: usize,
}

impl Mt19937 {
    fn new(seed: u32) -> Self {
        let mut mt = [0u32; 624];
        mt[0] = seed;
        for i in 1..624 {
            mt[i] = 1812433253u32
                .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Self { mt, index: 624 }
    }

    fn next_u32(&mut self) -> u32 {
        if self.index >= 624 {
            self.twist();
        }

        let mut y = self.mt[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= y >> 18;
        y
    }

    fn twist(&mut self) {
        const UPPER_MASK: u32 = 0x80000000;
        const LOWER_MASK: u32 = 0x7fffffff;
        const MATRIX_A: u32 = 0x9908b0df;

        for i in 0..624 {
            let y = (self.mt[i] & UPPER_MASK) | (self.mt[(i + 1) % 624] & LOWER_MASK);
            let mut next = self.mt[(i + 397) % 624] ^ (y >> 1);
            if y & 1 != 0 {
                next ^= MATRIX_A;
            }
            self.mt[i] = next;
        }
        self.index = 0;
    }

    fn gen_interval(&mut self, max: u32) -> u32 {
        let mut mask = max;
        mask |= mask >> 1;
        mask |= mask >> 2;
        mask |= mask >> 4;
        mask |= mask >> 8;
        mask |= mask >> 16;
        loop {
            let value = self.next_u32() & mask;
            if value <= max {
                return value;
            }
        }
    }

    fn gen_range(&mut self, upper: usize) -> usize {
        if upper <= 1 {
            return 0;
        }
        self.gen_interval((upper - 1) as u32) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toy_forest() -> RandomForestClassifier {
        RandomForestClassifier::new(RandomForestClassifier {
            n_features_in: 2,
            outputs: vec![
                OutputClasses {
                    classes: vec![0.0, 1.0],
                },
                OutputClasses {
                    classes: vec![0.0, 1.0],
                },
            ],
            estimators: vec![
                DecisionTreeClassifier {
                    children_left: vec![1, -1, -1],
                    children_right: vec![2, -1, -1],
                    feature: vec![0, -2, -2],
                    threshold: vec![0.5, -2.0, -2.0],
                    value: vec![
                        vec![vec![2.0, 2.0], vec![2.0, 2.0]],
                        vec![vec![0.75, 0.25], vec![0.0, 1.0]],
                        vec![vec![0.25, 0.75], vec![1.0, 0.0]],
                    ],
                },
                DecisionTreeClassifier {
                    children_left: vec![1, -1, -1],
                    children_right: vec![2, -1, -1],
                    feature: vec![1, -2, -2],
                    threshold: vec![1.5, -2.0, -2.0],
                    value: vec![
                        vec![vec![2.0, 2.0], vec![2.0, 2.0]],
                        vec![vec![0.5, 0.5], vec![0.75, 0.25]],
                        vec![vec![0.0, 1.0], vec![0.25, 0.75]],
                    ],
                },
            ],
        })
        .unwrap()
    }

    #[test]
    fn predicts_sklearn_style_tree_probabilities() {
        let forest = toy_forest();
        let probas = forest.predict_label_proba(&[vec![0.25, 2.0]], 1.0).unwrap();

        assert_eq!(probas.len(), 1);
        assert_eq!(probas[0].len(), 2);
        assert_eq!(probas[0][0], 0.625);
        assert_eq!(probas[0][1], 0.875);
    }

    #[test]
    fn threshold_equality_goes_left_like_sklearn() {
        let forest = toy_forest();
        let proba = forest.predict_one(&[0.5, 1.5]).unwrap();

        assert_eq!(proba[0][1], 0.375);
        assert_eq!(proba[1][1], 0.625);
    }

    #[test]
    fn loads_from_json() {
        let forest = toy_forest();
        let json = serde_json::to_vec(&forest).unwrap();
        let loaded = RandomForestClassifier::from_json_slice(&json).unwrap();
        assert_eq!(
            loaded
                .predict_label_proba(&[vec![0.25, 2.0], vec![1.0, 1.0]], 1.0)
                .unwrap(),
            vec![vec![0.625, 0.875], vec![0.625, 0.125]]
        );
    }

    #[test]
    fn rejects_wrong_feature_count() {
        let forest = toy_forest();
        let err = forest.predict_one(&[1.0]).unwrap_err();
        assert!(matches!(err, Error::FeatureCount { .. }));
    }

    #[test]
    fn trains_simple_multi_output_forest() {
        let x = vec![
            vec![0.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![2.0, 0.0],
            vec![2.0, 1.0],
        ];
        let y = vec![
            vec![0.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
        ];

        let forest = RandomForestClassifier::fit(
            &x,
            &y,
            RandomForestClassifierParameters {
                n_estimators: 7,
                bootstrap: false,
                max_features: MaxFeatures::All,
                random_state: 0,
                ..Default::default()
            },
        )
        .unwrap();

        let probas = forest.predict_label_proba(&x, 1.0).unwrap();
        assert!(probas[0][0] < 0.5);
        assert!(probas[2][0] > 0.5);
        assert!(probas[1][1] > 0.5);
        assert!(probas[2][1] < 0.5);
    }

    #[test]
    fn training_matches_sklearn_single_tree_no_bootstrap_all_features() {
        // Generated with scikit-learn 1.5.2:
        // RandomForestClassifier(
        //   n_estimators=1, bootstrap=False, max_features=None, random_state=0
        // )
        let x = vec![
            vec![0.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![2.0, 0.0],
            vec![2.0, 1.0],
        ];
        let y = vec![
            vec![0.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
        ];

        let forest = RandomForestClassifier::fit(
            &x,
            &y,
            RandomForestClassifierParameters {
                n_estimators: 1,
                bootstrap: false,
                max_features: MaxFeatures::All,
                random_state: 0,
                ..Default::default()
            },
        )
        .unwrap();
        let tree = &forest.estimators[0];

        assert_eq!(tree.children_left, vec![1, 2, -1, -1, 5, -1, -1]);
        assert_eq!(tree.children_right, vec![4, 3, -1, -1, 6, -1, -1]);
        assert_eq!(tree.feature, vec![1, 0, -2, -2, 0, -2, -2]);
        assert_eq!(tree.threshold, vec![0.5, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0]);
        assert_eq!(
            tree.value,
            vec![
                vec![vec![1.0 / 3.0, 2.0 / 3.0], vec![0.5, 0.5]],
                vec![vec![1.0 / 3.0, 2.0 / 3.0], vec![1.0, 0.0]],
                vec![vec![1.0, 0.0], vec![1.0, 0.0]],
                vec![vec![0.0, 1.0], vec![1.0, 0.0]],
                vec![vec![1.0 / 3.0, 2.0 / 3.0], vec![0.0, 1.0]],
                vec![vec![1.0, 0.0], vec![0.0, 1.0]],
                vec![vec![0.0, 1.0], vec![0.0, 1.0]],
            ]
        );
    }

    #[test]
    fn numpy_random_state_seed_zero_first_tree_seed_matches_sklearn() {
        let mut rng = Mt19937::new(0);
        assert_eq!(rng.gen_interval(i32::MAX as u32 - 1), 209652396);
    }

    #[test]
    fn numpy_random_state_bounded_draws_match_sklearn() {
        let mut rng = Mt19937::new(209652396);
        let draws = (0..10).map(|_| rng.gen_range(2)).collect::<Vec<_>>();
        assert_eq!(draws, vec![0, 1, 0, 0, 0, 0, 0, 1, 1, 1]);
    }

    #[test]
    fn training_matches_sklearn_single_tree_bootstrap_all_features() {
        // Generated with scikit-learn 1.5.2:
        // RandomForestClassifier(
        //   n_estimators=1, bootstrap=True, max_features=None, random_state=0
        // )
        let x = vec![
            vec![0.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![2.0, 0.0],
            vec![2.0, 1.0],
        ];
        let y = vec![
            vec![0.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
        ];

        let forest = RandomForestClassifier::fit(
            &x,
            &y,
            RandomForestClassifierParameters {
                n_estimators: 1,
                bootstrap: true,
                max_features: MaxFeatures::All,
                random_state: 0,
                ..Default::default()
            },
        )
        .unwrap();
        let tree = &forest.estimators[0];

        assert_eq!(tree.children_left, vec![1, -1, 3, -1, -1]);
        assert_eq!(tree.children_right, vec![2, -1, 4, -1, -1]);
        assert_eq!(tree.feature, vec![0, -2, 1, -2, -2]);
        assert_eq!(tree.threshold, vec![0.5, -2.0, 0.5, -2.0, -2.0]);
        assert_eq!(
            tree.value,
            vec![
                vec![vec![1.0 / 3.0, 2.0 / 3.0], vec![5.0 / 6.0, 1.0 / 6.0]],
                vec![vec![1.0, 0.0], vec![1.0, 0.0]],
                vec![vec![0.0, 1.0], vec![0.75, 0.25]],
                vec![vec![0.0, 1.0], vec![1.0, 0.0]],
                vec![vec![0.0, 1.0], vec![0.0, 1.0]],
            ]
        );
    }

    #[test]
    fn training_matches_sklearn_single_tree_no_bootstrap_sqrt_features() {
        // Generated with scikit-learn 1.5.2:
        // RandomForestClassifier(
        //   n_estimators=1, bootstrap=False, max_features="sqrt", random_state=0
        // )
        let x = vec![
            vec![0.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![2.0, 0.0],
            vec![2.0, 1.0],
        ];
        let y = vec![
            vec![0.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
        ];

        let forest = RandomForestClassifier::fit(
            &x,
            &y,
            RandomForestClassifierParameters {
                n_estimators: 1,
                bootstrap: false,
                max_features: MaxFeatures::Sqrt,
                random_state: 0,
                ..Default::default()
            },
        )
        .unwrap();
        let tree = &forest.estimators[0];

        assert_eq!(
            tree.children_left,
            vec![1, 2, -1, -1, 5, 6, -1, -1, 9, -1, -1]
        );
        assert_eq!(
            tree.children_right,
            vec![4, 3, -1, -1, 8, 7, -1, -1, 10, -1, -1]
        );
        assert_eq!(tree.feature, vec![0, 1, -2, -2, 0, 1, -2, -2, 1, -2, -2]);
        assert_eq!(
            tree.threshold,
            vec![0.5, 0.5, -2.0, -2.0, 1.5, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0]
        );
        assert_eq!(
            tree.value,
            vec![
                vec![vec![1.0 / 3.0, 2.0 / 3.0], vec![0.5, 0.5]],
                vec![vec![1.0, 0.0], vec![0.5, 0.5]],
                vec![vec![1.0, 0.0], vec![1.0, 0.0]],
                vec![vec![1.0, 0.0], vec![0.0, 1.0]],
                vec![vec![0.0, 1.0], vec![0.5, 0.5]],
                vec![vec![0.0, 1.0], vec![0.5, 0.5]],
                vec![vec![0.0, 1.0], vec![1.0, 0.0]],
                vec![vec![0.0, 1.0], vec![0.0, 1.0]],
                vec![vec![0.0, 1.0], vec![0.5, 0.5]],
                vec![vec![0.0, 1.0], vec![1.0, 0.0]],
                vec![vec![0.0, 1.0], vec![0.0, 1.0]],
            ]
        );
    }

    #[test]
    fn training_matches_sklearn_three_tree_default_shape() {
        // Generated with scikit-learn 1.5.2:
        // RandomForestClassifier(n_estimators=3, random_state=0)
        let x = vec![
            vec![0.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![2.0, 0.0],
            vec![2.0, 1.0],
        ];
        let y = vec![
            vec![0.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
        ];

        let forest = RandomForestClassifier::fit(
            &x,
            &y,
            RandomForestClassifierParameters {
                n_estimators: 3,
                random_state: 0,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(forest.estimators[0].children_left, vec![1, -1, 3, -1, -1]);
        assert_eq!(forest.estimators[0].children_right, vec![2, -1, 4, -1, -1]);
        assert_eq!(forest.estimators[0].feature, vec![0, -2, 1, -2, -2]);
        assert_eq!(
            forest.estimators[0].threshold,
            vec![0.5, -2.0, 0.5, -2.0, -2.0]
        );
        assert_eq!(forest.estimators[1].children_left, vec![1, 2, -1, -1, -1]);
        assert_eq!(forest.estimators[1].children_right, vec![4, 3, -1, -1, -1]);
        assert_eq!(forest.estimators[1].feature, vec![1, 0, -2, -2, -2]);
        assert_eq!(
            forest.estimators[1].threshold,
            vec![0.5, 1.0, -2.0, -2.0, -2.0]
        );
        assert_eq!(forest.estimators[2].children_left, vec![1, -1, -1]);
        assert_eq!(forest.estimators[2].children_right, vec![2, -1, -1]);
        assert_eq!(forest.estimators[2].feature, vec![1, -2, -2]);
        assert_eq!(forest.estimators[2].threshold, vec![0.5, -2.0, -2.0]);

        let positive = forest.predict_label_proba(&x, 1.0).unwrap();
        assert_eq!(
            positive,
            vec![
                vec![1.0 / 3.0, 0.0],
                vec![2.0 / 3.0, 2.0 / 3.0],
                vec![2.0 / 3.0, 0.0],
                vec![1.0, 1.0],
                vec![1.0, 0.0],
                vec![1.0, 1.0],
            ]
        );
    }
}
