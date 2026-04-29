//! Feature extraction and probability annotation for the CRF.

use crate::model::Gene;

/// Extract features at the gene (protein) level.
/// Each gene produces one ordered feature list; duplicate domains are merged.
pub fn extract_features_protein(genes: &[Gene]) -> Vec<Vec<String>> {
    genes
        .iter()
        .map(|gene| {
            let mut features = Vec::new();
            for domain in &gene.protein.domains {
                if !features.contains(&domain.name) {
                    features.push(domain.name.clone());
                }
            }
            features
        })
        .collect()
}

/// Extract features at the domain level.
/// Each domain produces its own feature list; unannotated genes get an empty list.
pub fn extract_features_domain(genes: &[Gene]) -> Vec<Vec<String>> {
    let mut features = Vec::new();
    for gene in genes {
        if gene.protein.domains.is_empty() {
            features.push(Vec::new());
        } else {
            for domain in &gene.protein.domains {
                features.push(vec![domain.name.clone()]);
            }
        }
    }
    features
}

/// Extract labels at the protein level for training.
pub fn extract_labels_protein(genes: &[Gene]) -> Vec<String> {
    genes
        .iter()
        .map(|gene| {
            if gene.average_probability().unwrap_or(0.0) > 0.5 {
                "1".to_string()
            } else {
                "0".to_string()
            }
        })
        .collect()
}

/// Extract labels at the domain level for training.
pub fn extract_labels_domain(genes: &[Gene]) -> Vec<String> {
    let mut labels = Vec::new();
    for gene in genes {
        if gene.protein.domains.is_empty() {
            let label = if gene.average_probability().unwrap_or(0.0) > 0.5 {
                "1"
            } else {
                "0"
            };
            labels.push(label.to_string());
        } else {
            for domain in &gene.protein.domains {
                let label = if domain.probability.unwrap_or(0.0) > 0.5 {
                    "1"
                } else {
                    "0"
                };
                labels.push(label.to_string());
            }
        }
    }
    labels
}

/// Annotate genes with CRF probabilities at the protein level.
pub fn annotate_probabilities_protein(genes: &[Gene], probabilities: &[f64]) -> Vec<Gene> {
    genes
        .iter()
        .zip(probabilities.iter())
        .map(|(gene, &p)| gene.with_probability(p))
        .collect()
}

/// Annotate genes with CRF probabilities at the domain level.
pub fn annotate_probabilities_domain(genes: &[Gene], probabilities: &[f64]) -> Vec<Gene> {
    let mut prob_iter = probabilities.iter();
    genes
        .iter()
        .map(|gene| {
            if gene.protein.domains.is_empty() {
                let p = *prob_iter.next().unwrap_or(&0.0);
                gene.with_probability(p)
            } else {
                let new_domains: Vec<_> = gene
                    .protein
                    .domains
                    .iter()
                    .map(|d| {
                        let p = *prob_iter.next().unwrap_or(&0.0);
                        d.with_probability(Some(p))
                    })
                    .collect();
                gene.with_protein(gene.protein.with_domains(new_domains))
            }
        })
        .collect()
}
