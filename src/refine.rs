//! Cluster refinement: extract contiguous gene clusters from CRF predictions.

use std::collections::HashSet;

use crate::model::{Cluster, Gene};

/// Pfam domains considered 'biosynthetic' by AntiSMASH.
pub static BIO_PFAMS: &[&str] = &[
    "PF00109", "PF02801", "PF08659", "PF00378", "PF08541", "PF08545",
    "PF02803", "PF00108", "PF02706", "PF03364", "PF08990", "PF00501",
    "PF00668", "PF08415", "PF00975", "PF03061", "PF00432", "PF00494",
    "PF03936", "PF01397", "PF04275", "PF00348", "PF02401",
    "PF04551", "PF00368", "PF00534", "PF00535", "PF02922", "PF01041",
    "PF00128", "PF00908", "PF02719", "PF04321", "PF01943", "PF02806",
    "PF02350", "PF02397", "PF04932", "PF01075", "PF00953", "PF01050",
    "PF03033", "PF01501", "PF05159", "PF04101", "PF02563", "PF08437",
    "PF02585", "PF01721", "PF02052", "PF02674", "PF03515", "PF04369",
    "PF08109", "PF08129", "PF09221", "PF09683", "PF10439", "PF11420",
    "PF11632", "PF11758", "PF12173", "PF04738", "PF04737", "PF04604",
    "PF05147", "PF08130", "PF00155", "PF00202",
    "PF00702", "PF06339", "PF04183", "PF10331", "PF03756", "PF00106",
    "PF01370", "PF00107", "PF08240", "PF00441", "PF02770", "PF02771",
    "PF08028", "PF01408", "PF02894", "PF00984", "PF00725", "PF03720",
    "PF03721", "PF07993", "PF02737", "PF00903", "PF00037", "PF04055",
    "PF00171", "PF00067", "PF01266", "PF01118", "PF02668", "PF00248",
    "PF01494", "PF01593", "PF03992", "PF00355", "PF01243", "PF00384",
    "PF01488", "PF00857", "PF04879", "PF08241", "PF08242", "PF00698",
    "PF00483", "PF00561", "PF00583", "PF01636", "PF01039", "PF00288",
    "PF00289", "PF02786", "PF01757", "PF02785", "PF02409", "PF01553",
    "PF02348", "PF00891", "PF01596", "PF04820", "PF02522", "PF08484",
    "PF08421",
];

fn bio_pfams_set() -> HashSet<&'static str> {
    BIO_PFAMS.iter().copied().collect()
}

/// Post-processor to extract contiguous clusters from CRF predictions.
pub struct ClusterRefiner {
    pub threshold: f64,
    pub criterion: String,
    pub n_cds: usize,
    pub n_biopfams: usize,
    pub average_threshold: f64,
    pub edge_distance: usize,
    pub trim: bool,
}

impl Default for ClusterRefiner {
    fn default() -> Self {
        Self {
            threshold: 0.8,
            criterion: "gecco".to_string(),
            n_cds: 5,
            n_biopfams: 5,
            average_threshold: 0.6,
            edge_distance: 0,
            trim: true,
        }
    }
}

impl ClusterRefiner {
    /// Find all valid clusters in a list of genes with probability annotations.
    pub fn iter_clusters(&self, genes: &[Gene]) -> Vec<Cluster> {
        let mut results = Vec::new();
        for (seq_genes, cluster) in self.raw_clusters(genes) {
            let cluster = if self.trim {
                self.trim_cluster(cluster)
            } else {
                cluster
            };
            if self.validate_cluster(&seq_genes, &cluster) {
                results.push(cluster);
            }
        }
        results
    }

    fn raw_clusters(&self, genes: &[Gene]) -> Vec<(Vec<Gene>, Cluster)> {
        // Group genes by source sequence
        let mut by_source: std::collections::BTreeMap<&str, Vec<&Gene>> =
            std::collections::BTreeMap::new();
        for gene in genes {
            by_source
                .entry(&gene.source_id)
                .or_default()
                .push(gene);
        }

        let mut results = Vec::new();

        for (seq_id, mut seq_genes) in by_source {
            // Sort by coordinate
            seq_genes.sort_by_key(|g| (g.start, g.end));

            // Group contiguous genes above/below threshold
            let mut in_cluster = false;
            let mut cluster_idx = 0usize;
            let mut current_genes: Vec<Gene> = Vec::new();
            let seq_cloned: Vec<Gene> = seq_genes.iter().map(|g| (*g).clone()).collect();

            for gene in &seq_genes {
                let new_state = gene
                    .average_probability()
                    .map(|p| p > self.threshold)
                    .unwrap_or(false);

                if new_state != in_cluster {
                    if in_cluster && !current_genes.is_empty() {
                        cluster_idx += 1;
                        let id = format!("{}_cluster_{}", seq_id, cluster_idx);
                        results.push((
                            seq_cloned.clone(),
                            Cluster::new(id, std::mem::take(&mut current_genes)),
                        ));
                    }
                    in_cluster = new_state;
                    current_genes.clear();
                }

                if in_cluster {
                    current_genes.push((*gene).clone());
                }
            }

            // Flush last cluster
            if in_cluster && !current_genes.is_empty() {
                cluster_idx += 1;
                let id = format!("{}_cluster_{}", seq_id, cluster_idx);
                results.push((
                    seq_cloned,
                    Cluster::new(id, current_genes),
                ));
            }
        }

        results
    }

    fn trim_cluster(&self, mut cluster: Cluster) -> Cluster {
        while !cluster.genes.is_empty() && cluster.genes[0].protein.domains.is_empty() {
            cluster.genes.remove(0);
        }
        while !cluster.genes.is_empty()
            && cluster.genes.last().unwrap().protein.domains.is_empty()
        {
            cluster.genes.pop();
        }
        cluster
    }

    fn validate_cluster(&self, seq: &[Gene], cluster: &Cluster) -> bool {
        if cluster.genes.is_empty() {
            return false;
        }

        match self.criterion.as_str() {
            "gecco" => {
                let annotated: Vec<&Gene> = cluster
                    .genes
                    .iter()
                    .filter(|g| !g.protein.domains.is_empty())
                    .collect();
                let cds_crit = annotated.len() >= self.n_cds;

                let edge_genes: HashSet<&str> = if self.edge_distance > 0 {
                    let annotated_ids: Vec<&str> = seq
                        .iter()
                        .filter(|g| !g.protein.domains.is_empty())
                        .map(|g| g.id())
                        .collect();
                    let mut edge = HashSet::new();
                    for id in annotated_ids.iter().take(self.edge_distance) {
                        edge.insert(*id);
                    }
                    for id in annotated_ids.iter().rev().take(self.edge_distance) {
                        edge.insert(*id);
                    }
                    edge
                } else {
                    HashSet::new()
                };

                let non_edge_count = cluster
                    .genes
                    .iter()
                    .filter(|g| !edge_genes.contains(g.id()))
                    .count();
                let edge_crit = non_edge_count >= self.n_cds;

                cds_crit && edge_crit
            }
            "antismash" => {
                let bio_set = bio_pfams_set();
                let domains: HashSet<&str> = cluster
                    .genes
                    .iter()
                    .flat_map(|g| g.protein.domains.iter().map(|d| d.name.as_str()))
                    .collect();

                let avg_prob: f64 = cluster
                    .genes
                    .iter()
                    .filter_map(|g| g.average_probability())
                    .sum::<f64>()
                    / cluster.genes.len() as f64;

                let p_crit = avg_prob >= self.average_threshold;
                let bio_crit = domains.iter().filter(|d| bio_set.contains(**d)).count()
                    >= self.n_biopfams;
                let cds_crit = cluster.genes.len() >= self.n_cds;

                p_crit && bio_crit && cds_crit
            }
            _ => false,
        }
    }
}
