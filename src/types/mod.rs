//! Supervised classifier to predict the type of a gene cluster.

pub mod backend;

use std::collections::BTreeSet;
use std::io::{Read, Seek};
use std::path::Path;

use anyhow::{Context, Result};

use crate::model::{Cluster, ClusterType};
use crate::sklearn_rf;

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

    /// Train a GECCO type classifier from Python GECCO's type-training assets.
    pub fn train_sklearn_rf_from_path(
        path: &Path,
        params: sklearn_rf::RandomForestClassifierParameters,
    ) -> Result<Self> {
        let data = load_training_data(path)?;
        let classes = classes_from_types(&data.types);
        let mut classifier = Self::new(classes);
        classifier.set_domains(data.domains);
        classifier.set_model(Box::new(
            backend::SklearnRF::untrained().with_params(params),
        ));
        classifier.fit_compositions(&data.compositions, &data.types)?;
        Ok(classifier)
    }

    /// Fit the configured random forest backend on precomputed compositions.
    pub fn fit_compositions(
        &mut self,
        compositions: &[Vec<f64>],
        types: &[ClusterType],
    ) -> Result<()> {
        let model = self
            .model
            .as_mut()
            .ok_or_else(|| anyhow::anyhow!("TypeClassifier model not set"))?;
        if compositions.len() != types.len() {
            anyhow::bail!(
                "composition/type row mismatch: {} vs {}",
                compositions.len(),
                types.len()
            );
        }
        let labels = self.binarizer.transform(types);
        model.fit(compositions, &labels)
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

fn classes_from_types(types: &[ClusterType]) -> Vec<String> {
    types
        .iter()
        .flat_map(|cluster_type| cluster_type.names.iter().cloned())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

/// GECCO type-classifier training data.
pub struct TypeTrainingData {
    pub domains: Vec<String>,
    pub compositions: Vec<Vec<f64>>,
    pub types: Vec<ClusterType>,
}

/// Load Python GECCO's embedded type-classifier training data.
pub fn load_training_data(path: &Path) -> Result<TypeTrainingData> {
    let domains = std::fs::read_to_string(path.join("domains.tsv"))
        .with_context(|| format!("reading {:?}", path.join("domains.tsv")))?
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(str::to_string)
        .collect::<Vec<_>>();

    let types = std::fs::read_to_string(path.join("types.tsv"))
        .with_context(|| format!("reading {:?}", path.join("types.tsv")))?
        .lines()
        .map(parse_type_line)
        .collect::<Result<Vec<_>>>()?;

    let compositions = load_scipy_coo_npz(&path.join("compositions.npz"))
        .with_context(|| format!("reading {:?}", path.join("compositions.npz")))?;
    if compositions.len() != types.len() {
        anyhow::bail!(
            "composition/type row mismatch: {} vs {}",
            compositions.len(),
            types.len()
        );
    }
    if compositions.first().map_or(0, Vec::len) != domains.len() {
        anyhow::bail!(
            "composition/domain column mismatch: {} vs {}",
            compositions.first().map_or(0, Vec::len),
            domains.len()
        );
    }

    Ok(TypeTrainingData {
        domains,
        compositions,
        types,
    })
}

fn parse_type_line(line: &str) -> Result<ClusterType> {
    let (_, labels) = line
        .split_once('\t')
        .ok_or_else(|| anyhow::anyhow!("invalid type row without tab: {line:?}"))?;
    Ok(ClusterType::new(
        labels
            .trim()
            .split(';')
            .filter(|label| !label.is_empty())
            .map(str::to_string),
    ))
}

fn load_scipy_coo_npz(path: &Path) -> Result<Vec<Vec<f64>>> {
    let file = std::fs::File::open(path)?;
    let mut zip = zip::ZipArchive::new(file)?;
    let row = read_npy_i32(&mut zip, "row.npy")?;
    let col = read_npy_i32(&mut zip, "col.npy")?;
    let data = read_npy_f64(&mut zip, "data.npy")?;
    let shape = read_npy_i64(&mut zip, "shape.npy")?;

    if shape.len() != 2 {
        anyhow::bail!("expected 2D sparse shape, got {:?}", shape);
    }
    if row.len() != col.len() || row.len() != data.len() {
        anyhow::bail!(
            "COO arrays have different lengths: row={}, col={}, data={}",
            row.len(),
            col.len(),
            data.len()
        );
    }

    let n_rows = usize::try_from(shape[0]).context("negative row count")?;
    let n_cols = usize::try_from(shape[1]).context("negative column count")?;
    let mut dense = vec![vec![0.0; n_cols]; n_rows];
    for ((&r, &c), &v) in row.iter().zip(col.iter()).zip(data.iter()) {
        let r = usize::try_from(r).context("negative COO row")?;
        let c = usize::try_from(c).context("negative COO column")?;
        let cell = dense
            .get_mut(r)
            .and_then(|dense_row| dense_row.get_mut(c))
            .ok_or_else(|| anyhow::anyhow!("COO coordinate out of bounds: ({r}, {c})"))?;
        *cell += v;
    }
    Ok(dense)
}

fn read_zip_member<R: Read + Seek>(zip: &mut zip::ZipArchive<R>, name: &str) -> Result<Vec<u8>> {
    let mut file = zip.by_name(name)?;
    let mut bytes = Vec::with_capacity(file.size() as usize);
    file.read_to_end(&mut bytes)?;
    Ok(bytes)
}

fn read_npy_i32<R: Read + Seek>(zip: &mut zip::ZipArchive<R>, name: &str) -> Result<Vec<i32>> {
    let bytes = read_zip_member(zip, name)?;
    let (descr, shape, data) = parse_npy_header(&bytes)?;
    if descr != "<i4" || shape.len() != 1 {
        anyhow::bail!("{name}: expected <i4 1D array, got {descr} {:?}", shape);
    }
    Ok(data
        .chunks_exact(4)
        .map(|chunk| i32::from_le_bytes(chunk.try_into().unwrap()))
        .collect())
}

fn read_npy_i64<R: Read + Seek>(zip: &mut zip::ZipArchive<R>, name: &str) -> Result<Vec<i64>> {
    let bytes = read_zip_member(zip, name)?;
    let (descr, shape, data) = parse_npy_header(&bytes)?;
    if descr != "<i8" || shape.len() != 1 {
        anyhow::bail!("{name}: expected <i8 1D array, got {descr} {:?}", shape);
    }
    Ok(data
        .chunks_exact(8)
        .map(|chunk| i64::from_le_bytes(chunk.try_into().unwrap()))
        .collect())
}

fn read_npy_f64<R: Read + Seek>(zip: &mut zip::ZipArchive<R>, name: &str) -> Result<Vec<f64>> {
    let bytes = read_zip_member(zip, name)?;
    let (descr, shape, data) = parse_npy_header(&bytes)?;
    if descr != "<f8" || shape.len() != 1 {
        anyhow::bail!("{name}: expected <f8 1D array, got {descr} {:?}", shape);
    }
    Ok(data
        .chunks_exact(8)
        .map(|chunk| f64::from_le_bytes(chunk.try_into().unwrap()))
        .collect())
}

fn parse_npy_header(bytes: &[u8]) -> Result<(String, Vec<usize>, &[u8])> {
    if bytes.len() < 10 || &bytes[..6] != b"\x93NUMPY" {
        anyhow::bail!("invalid NPY magic");
    }
    let major = bytes[6];
    let header_len_offset = 8usize;
    let (header_len, data_offset) = match major {
        1 => (
            u16::from_le_bytes(bytes[header_len_offset..header_len_offset + 2].try_into()?)
                as usize,
            10usize,
        ),
        2 | 3 => (
            u32::from_le_bytes(bytes[header_len_offset..header_len_offset + 4].try_into()?)
                as usize,
            12usize,
        ),
        _ => anyhow::bail!("unsupported NPY version {major}"),
    };
    let header_end = data_offset + header_len;
    let header = std::str::from_utf8(&bytes[data_offset..header_end])?;
    if header.contains("'fortran_order': True") {
        anyhow::bail!("Fortran-order NPY arrays are not supported");
    }
    let descr = extract_header_string(header, "'descr':")?;
    let shape = extract_header_shape(header)?;
    Ok((descr, shape, &bytes[header_end..]))
}

fn extract_header_string(header: &str, key: &str) -> Result<String> {
    let rest = header
        .split_once(key)
        .ok_or_else(|| anyhow::anyhow!("NPY header missing {key}"))?
        .1
        .trim_start();
    let quote = rest
        .chars()
        .next()
        .ok_or_else(|| anyhow::anyhow!("NPY header missing string value"))?;
    let rest = &rest[quote.len_utf8()..];
    let end = rest
        .find(quote)
        .ok_or_else(|| anyhow::anyhow!("unterminated NPY header string"))?;
    Ok(rest[..end].to_string())
}

fn extract_header_shape(header: &str) -> Result<Vec<usize>> {
    let rest = header
        .split_once("'shape':")
        .ok_or_else(|| anyhow::anyhow!("NPY header missing shape"))?
        .1;
    let start = rest
        .find('(')
        .ok_or_else(|| anyhow::anyhow!("NPY header shape missing '('"))?;
    let end = rest
        .find(')')
        .ok_or_else(|| anyhow::anyhow!("NPY header shape missing ')'"))?;
    rest[start + 1..end]
        .split(',')
        .map(str::trim)
        .filter(|part| !part.is_empty())
        .map(|part| {
            part.parse::<usize>()
                .with_context(|| format!("parsing NPY shape dimension {part:?}"))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use crate::sklearn_rf;

    #[test]
    fn loads_original_gecco_type_training_data() {
        let data = super::load_training_data(std::path::Path::new("GECCO/gecco/types")).unwrap();
        assert_eq!(data.domains.len(), 2766);
        assert_eq!(data.compositions.len(), 1870);
        assert_eq!(data.compositions[0].len(), 2766);
        assert_eq!(data.types.len(), 1870);
        assert!((data.compositions[0][1] - 0.025000000000593624).abs() < f64::EPSILON);
        assert!(data.types[0].names.contains("Polyketide"));
    }

    #[test]
    #[ignore = "trains on the full GECCO type matrix; run explicitly when checking trainer speed"]
    fn trains_one_tree_from_original_gecco_type_data() {
        let classifier = super::TypeClassifier::train_sklearn_rf_from_path(
            std::path::Path::new("GECCO/gecco/types"),
            sklearn_rf::RandomForestClassifierParameters {
                n_estimators: 1,
                random_state: 0,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            classifier.binarizer.classes,
            vec![
                "Alkaloid".to_string(),
                "NRP".to_string(),
                "Polyketide".to_string(),
                "RiPP".to_string(),
                "Saccharide".to_string(),
                "Terpene".to_string(),
            ]
        );
    }

    #[test]
    #[ignore = "trains the full 100-tree GECCO type forest; run explicitly when checking trainer speed"]
    fn trains_default_forest_from_original_gecco_type_data() {
        let classifier = super::TypeClassifier::train_sklearn_rf_from_path(
            std::path::Path::new("GECCO/gecco/types"),
            sklearn_rf::RandomForestClassifierParameters::default(),
        )
        .unwrap();
        assert_eq!(classifier.binarizer.classes.len(), 6);
    }
}
