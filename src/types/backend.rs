//! Random Forest backend for BGC type classification.

use anyhow::Result;

use super::RandomForestModel;
use crate::sklearn_rf;

/// Prediction-only backend for an exported sklearn RandomForestClassifier.
pub struct SklearnRF {
    model: Option<sklearn_rf::RandomForestClassifier>,
    params: sklearn_rf::RandomForestClassifierParameters,
}

impl SklearnRF {
    pub fn new(model: sklearn_rf::RandomForestClassifier) -> Self {
        Self {
            model: Some(model),
            params: Default::default(),
        }
    }

    pub fn untrained() -> Self {
        Self {
            model: None,
            params: Default::default(),
        }
    }

    pub fn with_params(mut self, params: sklearn_rf::RandomForestClassifierParameters) -> Self {
        self.params = params;
        self
    }

    pub fn model(&self) -> Option<&sklearn_rf::RandomForestClassifier> {
        self.model.as_ref()
    }

    pub fn to_json_vec(&self) -> Result<Vec<u8>> {
        self.model
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("sklearn RF model is not fitted"))?
            .to_json_vec()
            .map_err(|e| anyhow::anyhow!("serializing sklearn RF: {e}"))
    }

    pub fn from_json_slice(bytes: &[u8]) -> Result<Self> {
        let model = sklearn_rf::RandomForestClassifier::from_json_slice(bytes)
            .map_err(|e| anyhow::anyhow!("loading sklearn RF: {e}"))?;
        Ok(Self::new(model))
    }
}

impl RandomForestModel for SklearnRF {
    fn fit(&mut self, x: &[Vec<f64>], y: &[Vec<f64>]) -> Result<()> {
        let model = sklearn_rf::RandomForestClassifier::fit(x, y, self.params.clone())
            .map_err(|e| anyhow::anyhow!("training sklearn RF: {e}"))?;
        self.model = Some(model);
        Ok(())
    }

    fn predict_proba(&self, x: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
        let probas = self
            .model
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("sklearn RF model is not fitted"))?
            .predict_proba(x)
            .map_err(|e| anyhow::anyhow!("predicting with sklearn RF: {e}"))?;

        Ok(probas
            .into_iter()
            .map(|sample| {
                sample
                    .into_iter()
                    .map(|output| output.first().map_or(0.0, |p0| 1.0 - p0))
                    .collect()
            })
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Cluster, Domain, Gene, Protein, Strand};
    use std::collections::BTreeMap;

    #[test]
    fn test_sklearn_rf_backend_predicts_gecco_types() {
        let model = sklearn_rf::RandomForestClassifier::new(sklearn_rf::RandomForestClassifier {
            n_features_in: 2,
            outputs: vec![
                sklearn_rf::OutputClasses {
                    classes: vec![0.0, 1.0],
                },
                sklearn_rf::OutputClasses {
                    classes: vec![0.0, 1.0],
                },
            ],
            estimators: vec![sklearn_rf::DecisionTreeClassifier {
                children_left: vec![1, -1, -1],
                children_right: vec![2, -1, -1],
                feature: vec![0, -2, -2],
                threshold: vec![0.5, -2.0, -2.0],
                value: vec![
                    vec![vec![0.5, 0.5], vec![0.5, 0.5]],
                    vec![vec![1.0, 0.0], vec![0.0, 1.0]],
                    vec![vec![0.0, 1.0], vec![1.0, 0.0]],
                ],
            }],
        })
        .unwrap();

        let mut classifier = super::super::TypeClassifier::new(vec!["Left".into(), "Right".into()]);
        classifier.set_domains(vec!["PF00001".into(), "PF00002".into()]);
        classifier.set_model(Box::new(SklearnRF::new(model)));

        let domain = Domain {
            name: "PF00002".into(),
            start: 1,
            end: 10,
            hmm: "Pfam".into(),
            i_evalue: 1e-20,
            pvalue: 1e-20,
            probability: Some(1.0),
            cluster_weight: None,
            go_terms: Vec::new(),
            go_functions: Vec::new(),
            qualifiers: BTreeMap::new(),
        };
        let gene = Gene {
            source_id: "seq".into(),
            start: 1,
            end: 30,
            strand: Strand::Coding,
            protein: Protein {
                id: "prot".into(),
                seq: "MA".into(),
                domains: vec![domain],
            },
            qualifiers: BTreeMap::new(),
            probability: None,
        };
        let mut clusters = vec![Cluster::new("seq_cluster_1", vec![gene])];

        classifier.predict_types(&mut clusters).unwrap();

        let cluster_type = clusters[0].cluster_type.as_ref().unwrap();
        assert!(cluster_type.names.contains("Right"));
        assert!(!cluster_type.names.contains("Left"));
        assert_eq!(clusters[0].type_probabilities["Left"], 0.0);
        assert_eq!(clusters[0].type_probabilities["Right"], 1.0);
    }

    #[test]
    fn test_sklearn_rf_backend_fit_predict() {
        let x = vec![
            vec![0.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
        ];
        let y = vec![
            vec![0.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
        ];

        let mut model =
            SklearnRF::untrained().with_params(sklearn_rf::RandomForestClassifierParameters {
                n_estimators: 3,
                bootstrap: false,
                max_features: sklearn_rf::MaxFeatures::All,
                ..Default::default()
            });
        model.fit(&x, &y).unwrap();

        let probas = model.predict_proba(&x).unwrap();
        assert!(probas[0][0] < 0.5);
        assert!(probas[2][0] > 0.5);
        assert!(probas[1][1] > 0.5);
        assert!(probas[2][1] < 0.5);
    }

    #[test]
    fn test_type_classifier_fit_compositions_with_sklearn_rf() {
        let mut classifier = super::super::TypeClassifier::new(vec!["Left".into(), "Right".into()]);
        classifier.set_model(Box::new(SklearnRF::untrained().with_params(
            sklearn_rf::RandomForestClassifierParameters {
                n_estimators: 3,
                bootstrap: false,
                max_features: sklearn_rf::MaxFeatures::All,
                ..Default::default()
            },
        )));

        let x = vec![
            vec![0.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
        ];
        let y = vec![
            crate::model::ClusterType::new(["Left"]),
            crate::model::ClusterType::new(["Right"]),
            crate::model::ClusterType::new(["Left"]),
            crate::model::ClusterType::new(["Left", "Right"]),
        ];

        classifier.fit_compositions(&x, &y).unwrap();
        let probas = classifier
            .model
            .as_ref()
            .unwrap()
            .predict_proba(&x)
            .unwrap();
        assert!(probas[0][0] > 0.5);
        assert!(probas[1][1] > 0.5);
        assert!(probas[2][0] > 0.5);
        assert!(probas[3][0] > 0.5);
        assert!(probas[3][1] > 0.5);
    }
}
