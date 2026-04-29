//! Data layer types for gene cluster detection.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fmt;

use serde::{Deserialize, Serialize};

use crate::interpro::GOTerm;

// ---------------------------------------------------------------------------
// ClusterType
// ---------------------------------------------------------------------------

/// An immutable storage for the type of a gene cluster.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ClusterType {
    pub names: BTreeSet<String>,
}

impl ClusterType {
    pub fn new(names: impl IntoIterator<Item = impl Into<String>>) -> Self {
        Self {
            names: names.into_iter().map(Into::into).collect(),
        }
    }

    pub fn unknown() -> Self {
        Self {
            names: BTreeSet::new(),
        }
    }

    pub fn is_unknown(&self) -> bool {
        self.names.is_empty()
    }

    /// Unpack a composite type into individual single types.
    pub fn unpack(&self) -> Vec<ClusterType> {
        self.names
            .iter()
            .map(|n| ClusterType::new(std::iter::once(n.clone())))
            .collect()
    }
}

impl fmt::Display for ClusterType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.names.is_empty() {
            write!(f, "Unknown")
        } else {
            let mut first = true;
            for name in &self.names {
                if !first {
                    write!(f, ";")?;
                }
                write!(f, "{}", name)?;
                first = false;
            }
            Ok(())
        }
    }
}

// ---------------------------------------------------------------------------
// Strand
// ---------------------------------------------------------------------------

/// DNA strand on which a gene is located.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Strand {
    Coding,
    Reverse,
}

impl Strand {
    pub fn sign(&self) -> &'static str {
        match self {
            Strand::Coding => "+",
            Strand::Reverse => "-",
        }
    }

    pub fn from_sign(s: &str) -> Option<Self> {
        match s {
            "+" => Some(Strand::Coding),
            "-" => Some(Strand::Reverse),
            _ => None,
        }
    }

    pub fn as_int(&self) -> i8 {
        match self {
            Strand::Coding => 1,
            Strand::Reverse => -1,
        }
    }
}

// ---------------------------------------------------------------------------
// Domain
// ---------------------------------------------------------------------------

/// A conserved region within a protein.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Domain {
    pub name: String,
    /// Start coordinate within protein (1-based, inclusive).
    pub start: i64,
    /// End coordinate within protein (inclusive).
    pub end: i64,
    /// HMM library name (e.g. "Pfam").
    pub hmm: String,
    /// Independent e-value from hmmsearch.
    pub i_evalue: f64,
    /// P-value from hmmsearch.
    pub pvalue: f64,
    /// Probability of being part of a gene cluster (None if not yet predicted).
    pub probability: Option<f64>,
    /// Cluster weight from CRF state features.
    pub cluster_weight: Option<f64>,
    pub go_terms: Vec<GOTerm>,
    pub go_functions: Vec<GOTerm>,
    pub qualifiers: BTreeMap<String, Vec<String>>,
}

impl Domain {
    pub fn with_probability(&self, probability: Option<f64>) -> Domain {
        let mut d = self.clone();
        d.probability = probability;
        d
    }

    pub fn with_cluster_weight(&self, cluster_weight: Option<f64>) -> Domain {
        let mut d = self.clone();
        d.cluster_weight = cluster_weight;
        d
    }
}

// ---------------------------------------------------------------------------
// Protein
// ---------------------------------------------------------------------------

/// A sequence of amino-acids translated from a gene.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Protein {
    pub id: String,
    pub seq: String,
    pub domains: Vec<Domain>,
}

impl Protein {
    pub fn new(id: impl Into<String>, seq: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            seq: seq.into(),
            domains: Vec::new(),
        }
    }

    pub fn with_seq(&self, seq: String) -> Protein {
        Protein {
            id: self.id.clone(),
            seq,
            domains: self.domains.clone(),
        }
    }

    pub fn with_domains(&self, domains: Vec<Domain>) -> Protein {
        Protein {
            id: self.id.clone(),
            seq: self.seq.clone(),
            domains,
        }
    }
}

// ---------------------------------------------------------------------------
// Gene
// ---------------------------------------------------------------------------

/// A nucleotide sequence coding a protein.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gene {
    /// Source sequence identifier.
    pub source_id: String,
    /// Start coordinate within source (1-based, inclusive).
    pub start: i64,
    /// End coordinate within source (inclusive).
    pub end: i64,
    pub strand: Strand,
    pub protein: Protein,
    pub qualifiers: BTreeMap<String, Vec<String>>,
    /// Override probability (takes precedence over domain-derived probability).
    pub probability: Option<f64>,
}

impl Gene {
    pub fn id(&self) -> &str {
        &self.protein.id
    }

    /// Average of domain probabilities (or override if set).
    pub fn average_probability(&self) -> Option<f64> {
        if let Some(p) = self.probability {
            return Some(p);
        }
        let probas: Vec<f64> = self
            .protein
            .domains
            .iter()
            .filter_map(|d| d.probability)
            .collect();
        if probas.is_empty() {
            None
        } else {
            Some(statistics_mean(&probas))
        }
    }

    /// Highest domain probability (or override if set).
    pub fn maximum_probability(&self) -> Option<f64> {
        if let Some(p) = self.probability {
            return Some(p);
        }
        self.protein
            .domains
            .iter()
            .filter_map(|d| d.probability)
            .fold(None, |acc, p| Some(acc.map_or(p, |a: f64| a.max(p))))
    }

    /// Gene functions from GO term annotations.
    pub fn functions(&self) -> HashSet<String> {
        let mut fns: HashSet<String> = self
            .protein
            .domains
            .iter()
            .flat_map(|d| d.go_functions.iter().map(|t| t.name.clone()))
            .collect();
        if fns.is_empty() {
            fns.insert("unknown".to_string());
        }
        fns
    }

    pub fn with_protein(&self, protein: Protein) -> Gene {
        Gene {
            source_id: self.source_id.clone(),
            start: self.start,
            end: self.end,
            strand: self.strand,
            protein,
            qualifiers: self.qualifiers.clone(),
            probability: self.probability,
        }
    }

    pub fn with_probability(&self, probability: f64) -> Gene {
        let new_domains: Vec<Domain> = self
            .protein
            .domains
            .iter()
            .map(|d| d.with_probability(Some(probability)))
            .collect();
        Gene {
            source_id: self.source_id.clone(),
            start: self.start,
            end: self.end,
            strand: self.strand,
            protein: self.protein.with_domains(new_domains),
            qualifiers: self.qualifiers.clone(),
            probability: Some(probability),
        }
    }
}

fn statistics_mean(values: &[f64]) -> f64 {
    let mut partials = Vec::<f64>::new();
    for &value in values {
        let mut x = value;
        let mut next = Vec::with_capacity(partials.len() + 1);
        for &partial in &partials {
            let (hi, lo) = match x.abs().partial_cmp(&partial.abs()) {
                Some(Ordering::Less) => (partial, x),
                _ => (x, partial),
            };
            x = hi + lo;
            let yr = x - hi;
            let lo = lo - yr;
            if lo != 0.0 {
                next.push(lo);
            }
        }
        next.push(x);
        partials = next;
    }
    partials.iter().sum::<f64>() / values.len() as f64
}

// ---------------------------------------------------------------------------
// Cluster
// ---------------------------------------------------------------------------

/// A sequence of contiguous genes forming a putative BGC.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Cluster {
    pub id: String,
    pub genes: Vec<Gene>,
    pub cluster_type: Option<ClusterType>,
    pub type_probabilities: BTreeMap<String, f64>,
}

impl Cluster {
    pub fn new(id: impl Into<String>, genes: Vec<Gene>) -> Self {
        Self {
            id: id.into(),
            genes,
            cluster_type: None,
            type_probabilities: BTreeMap::new(),
        }
    }

    pub fn source_id(&self) -> &str {
        &self.genes[0].source_id
    }

    pub fn start(&self) -> i64 {
        self.genes.iter().map(|g| g.start).min().unwrap_or(0)
    }

    pub fn end(&self) -> i64 {
        self.genes.iter().map(|g| g.end).max().unwrap_or(0)
    }

    pub fn average_probability(&self) -> Option<f64> {
        let probas: Vec<f64> = self
            .genes
            .iter()
            .filter_map(|g| g.average_probability())
            .collect();
        if probas.is_empty() {
            None
        } else {
            Some(statistics_mean(&probas))
        }
    }

    pub fn maximum_probability(&self) -> Option<f64> {
        self.genes
            .iter()
            .filter_map(|g| g.maximum_probability())
            .fold(None, |acc, p| Some(acc.map_or(p, |a: f64| a.max(p))))
    }

    /// Compute weighted domain composition vector.
    pub fn domain_composition(
        &self,
        all_possible: Option<&[String]>,
        normalize: bool,
        minlog_weights: bool,
        use_pvalue: bool,
    ) -> Vec<f64> {
        // Collect all domains from cluster genes
        let domains: Vec<&Domain> = self
            .genes
            .iter()
            .flat_map(|g| g.protein.domains.iter())
            .collect();

        let names: Vec<&str> = domains.iter().map(|d| d.name.as_str()).collect();
        let weights: Vec<f64> = domains
            .iter()
            .map(|d| {
                let v = if use_pvalue { d.pvalue } else { d.i_evalue };
                if minlog_weights {
                    -v.log10()
                } else {
                    1.0 - v
                }
            })
            .collect();

        let unique_names: HashSet<&str> = names.iter().copied().collect();

        // Determine the set of all possible domain names
        let default_possible: Vec<String>;
        let all_possible_ref: &[String] = match all_possible {
            Some(ap) => ap,
            None => {
                let mut u: Vec<String> = unique_names.iter().map(|s| s.to_string()).collect();
                u.sort();
                u.dedup();
                default_possible = u;
                &default_possible
            }
        };

        let mut composition = vec![0.0f64; all_possible_ref.len()];
        for (i, dom) in all_possible_ref.iter().enumerate() {
            if unique_names.contains(dom.as_str()) {
                for (j, name) in names.iter().enumerate() {
                    if *name == dom.as_str() {
                        composition[i] += weights[j];
                    }
                }
            }
        }

        if normalize {
            let total: f64 = composition.iter().sum();
            if total > 0.0 {
                for v in &mut composition {
                    *v /= total;
                }
            }
        }

        composition
    }
}
