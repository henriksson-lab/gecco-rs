//! TSV table I/O for FeatureTable, GeneTable, ClusterTable.

use std::collections::BTreeMap;
use std::io::{Read, Write};

use anyhow::Result;
use serde::{Deserialize, Serialize};

use crate::model::*;

// ---------------------------------------------------------------------------
// FeatureTable
// ---------------------------------------------------------------------------

#[derive(Debug, Serialize, Deserialize)]
pub struct FeatureRow {
    pub sequence_id: String,
    pub protein_id: String,
    pub start: i64,
    pub end: i64,
    pub strand: String,
    pub domain: String,
    pub hmm: String,
    pub i_evalue: f64,
    pub pvalue: f64,
    pub domain_start: i64,
    pub domain_end: i64,
    #[serde(default)]
    pub cluster_probability: Option<f64>,
}

pub struct FeatureTable;

impl FeatureTable {
    /// Write a feature table from genes.
    pub fn write_from_genes(writer: impl Write, genes: &[Gene]) -> Result<()> {
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(writer);
        for gene in genes {
            for domain in &gene.protein.domains {
                wtr.serialize(FeatureRow {
                    sequence_id: gene.source_id.clone(),
                    protein_id: gene.protein.id.clone(),
                    start: gene.start,
                    end: gene.end,
                    strand: gene.strand.sign().to_string(),
                    domain: domain.name.clone(),
                    hmm: domain.hmm.clone(),
                    i_evalue: domain.i_evalue,
                    pvalue: domain.pvalue,
                    domain_start: domain.start,
                    domain_end: domain.end,
                    cluster_probability: domain.probability,
                })?;
            }
        }
        wtr.flush()?;
        Ok(())
    }

    /// Load genes from a feature table TSV.
    pub fn read_to_genes(reader: impl Read) -> Result<Vec<Gene>> {
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(reader);

        // Group rows by protein_id
        let mut rows_by_protein: Vec<(String, Vec<FeatureRow>)> = Vec::new();
        let mut current_protein: Option<String> = None;
        for result in rdr.deserialize() {
            let row: FeatureRow = result?;
            match &current_protein {
                Some(pid) if pid == &row.protein_id => {
                    rows_by_protein.last_mut().unwrap().1.push(row);
                }
                _ => {
                    current_protein = Some(row.protein_id.clone());
                    rows_by_protein.push((row.protein_id.clone(), vec![row]));
                }
            }
        }

        let mut genes = Vec::new();
        for (_pid, rows) in rows_by_protein {
            let first = &rows[0];
            let domains: Vec<Domain> = rows
                .iter()
                .map(|r| Domain {
                    name: r.domain.clone(),
                    start: r.domain_start,
                    end: r.domain_end,
                    hmm: r.hmm.clone(),
                    i_evalue: r.i_evalue,
                    pvalue: r.pvalue,
                    probability: r.cluster_probability,
                    cluster_weight: None,
                    go_terms: Vec::new(),
                    go_functions: Vec::new(),
                    qualifiers: BTreeMap::new(),
                })
                .collect();
            let protein = Protein {
                id: first.protein_id.clone(),
                seq: String::new(),
                domains,
            };
            let strand = Strand::from_sign(&first.strand).unwrap_or(Strand::Coding);
            genes.push(Gene {
                source_id: first.sequence_id.clone(),
                start: first.start,
                end: first.end,
                strand,
                protein,
                qualifiers: BTreeMap::new(),
                probability: None,
            });
        }
        Ok(genes)
    }
}

// ---------------------------------------------------------------------------
// GeneTable
// ---------------------------------------------------------------------------

#[derive(Debug, Serialize, Deserialize)]
pub struct GeneRow {
    pub sequence_id: String,
    pub protein_id: String,
    pub start: i64,
    pub end: i64,
    pub strand: String,
    #[serde(default)]
    pub average_p: Option<f64>,
    #[serde(default)]
    pub max_p: Option<f64>,
}

pub struct GeneTable;

impl GeneTable {
    pub fn write_from_genes(writer: impl Write, genes: &[Gene]) -> Result<()> {
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(writer);
        for gene in genes {
            wtr.serialize(GeneRow {
                sequence_id: gene.source_id.clone(),
                protein_id: gene.protein.id.clone(),
                start: gene.start,
                end: gene.end,
                strand: gene.strand.sign().to_string(),
                average_p: gene.average_probability(),
                max_p: gene.maximum_probability(),
            })?;
        }
        wtr.flush()?;
        Ok(())
    }

    pub fn read_to_genes(reader: impl Read) -> Result<Vec<Gene>> {
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(reader);
        let mut genes = Vec::new();
        for result in rdr.deserialize() {
            let row: GeneRow = result?;
            let strand = Strand::from_sign(&row.strand).unwrap_or(Strand::Coding);
            let protein = Protein::new(&row.protein_id, "");
            genes.push(Gene {
                source_id: row.sequence_id,
                start: row.start,
                end: row.end,
                strand,
                protein,
                qualifiers: BTreeMap::new(),
                probability: row.average_p,
            });
        }
        Ok(genes)
    }
}

// ---------------------------------------------------------------------------
// ClusterTable
// ---------------------------------------------------------------------------

#[derive(Debug, Serialize, Deserialize)]
pub struct ClusterRow {
    pub sequence_id: String,
    pub cluster_id: String,
    pub start: i64,
    pub end: i64,
    #[serde(default)]
    pub average_p: Option<f64>,
    #[serde(default)]
    pub max_p: Option<f64>,
    #[serde(default = "default_type")]
    pub r#type: String,
    #[serde(default)]
    pub alkaloid_probability: f64,
    #[serde(default)]
    pub nrp_probability: f64,
    #[serde(default)]
    pub polyketide_probability: f64,
    #[serde(default)]
    pub ripp_probability: f64,
    #[serde(default)]
    pub saccharide_probability: f64,
    #[serde(default)]
    pub terpene_probability: f64,
    #[serde(default)]
    pub proteins: String,
    #[serde(default)]
    pub domains: String,
}

fn default_type() -> String {
    "Unknown".to_string()
}

pub struct ClusterTable;

impl ClusterTable {
    pub fn write_from_clusters(writer: impl Write, clusters: &[Cluster]) -> Result<()> {
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(writer);
        for cluster in clusters {
            let mut proteins: Vec<&str> = cluster
                .genes
                .iter()
                .map(|g| g.protein.id.as_str())
                .collect();
            proteins.sort();
            let mut domains: Vec<&str> = cluster
                .genes
                .iter()
                .flat_map(|g| g.protein.domains.iter().map(|d| d.name.as_str()))
                .collect();
            domains.sort();

            wtr.serialize(ClusterRow {
                sequence_id: cluster.source_id().to_string(),
                cluster_id: cluster.id.clone(),
                start: cluster.start(),
                end: cluster.end(),
                average_p: cluster.average_probability(),
                max_p: cluster.maximum_probability(),
                r#type: cluster
                    .cluster_type
                    .as_ref()
                    .map(|t| t.to_string())
                    .unwrap_or_else(|| "Unknown".to_string()),
                alkaloid_probability: type_probability(cluster, "Alkaloid"),
                nrp_probability: type_probability(cluster, "NRP"),
                polyketide_probability: type_probability(cluster, "Polyketide"),
                ripp_probability: type_probability(cluster, "RiPP"),
                saccharide_probability: type_probability(cluster, "Saccharide"),
                terpene_probability: type_probability(cluster, "Terpene"),
                proteins: proteins.join(";"),
                domains: domains.join(";"),
            })?;
        }
        wtr.flush()?;
        Ok(())
    }
}

fn type_probability(cluster: &Cluster, key: &str) -> f64 {
    cluster.type_probabilities.get(key).copied().unwrap_or(0.0)
}
