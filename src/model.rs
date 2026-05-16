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
    /// Create a new product type from one or more base types.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// let t1 = ClusterType::unknown();                              // unknown type
    /// let t2 = ClusterType::new(["Polyketide"]);                    // single type
    /// let t3 = ClusterType::new(["Polyketide", "NRP"]);             // multiple types
    /// ```
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

    /// Unpack a composite `ClusterType` into a list of individual types.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// let ty = ClusterType::new(["Polyketide", "Saccharide"]);
    /// assert_eq!(
    ///     ty.unpack(),
    ///     vec![ClusterType::new(["Polyketide"]), ClusterType::new(["Saccharide"])],
    /// );
    /// ```
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

/// A flag to declare on which DNA strand a gene is located.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Strand {
    Coding,
    Reverse,
}

impl Strand {
    /// The strand as a single sign (`+` or `-`).
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
    /// The accession of the protein domain in the source HMM.
    pub name: String,
    /// The start coordinate of the domain within the protein sequence
    /// (first amino-acid at 1).
    pub start: i64,
    /// The end coordinate of the domain within the protein sequence
    /// (inclusive).
    pub end: i64,
    /// The name of the HMM library this domain belongs to (e.g. `Pfam`,
    /// `Panther`).
    pub hmm: String,
    /// The independent e-value reported by `hmmsearch`, measuring how
    /// reliable the domain annotation is.
    pub i_evalue: f64,
    /// The p-value reported by `hmmsearch`, measuring how likely the
    /// domain score is.
    pub pvalue: f64,
    /// The probability that this domain is part of a gene cluster, or
    /// `None` if no prediction has been made yet.
    pub probability: Option<f64>,
    /// The weight for this domain, measuring its importance as inferred
    /// from the training clusters by the CRF model.
    pub cluster_weight: Option<f64>,
    /// The Gene Ontology terms for this particular domain.
    pub go_terms: Vec<GOTerm>,
    /// The Gene Ontology term families for this particular domain.
    ///
    /// Term families are extracted by taking the highest superclasses
    /// (excluding the root) of each Gene Ontology term in the
    /// `molecular_function` namespace associated with this domain.
    pub go_functions: Vec<GOTerm>,
    /// Feature qualifiers added to the GenBank `misc_feature` built from
    /// this `Domain`.
    pub qualifiers: BTreeMap<String, Vec<String>>,
}

impl Domain {
    /// Copy the current domain and assign it a cluster probability.
    pub fn with_probability(&self, probability: Option<f64>) -> Domain {
        let mut d = self.clone();
        d.probability = probability;
        d
    }

    /// Copy the current domain and assign it a cluster weight.
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
    /// The identifier of the protein.
    pub id: String,
    /// The sequence of amino-acids of this protein.
    pub seq: String,
    /// A list of domains found in the protein sequence.
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

    /// Copy the current protein and assign it a new sequence.
    pub fn with_seq(&self, seq: String) -> Protein {
        Protein {
            id: self.id.clone(),
            seq,
            domains: self.domains.clone(),
        }
    }

    /// Copy the current protein and assign it new domains.
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
    /// The identifier of the DNA sequence this gene was found in.
    pub source_id: String,
    /// The index of the leftmost nucleotide of the gene within the source
    /// sequence, independent of the strandedness.
    pub start: i64,
    /// The index of the rightmost nucleotide of the gene within the source
    /// sequence.
    pub end: i64,
    /// The strand where the gene is located.
    pub strand: Strand,
    /// The protein translated from this gene.
    pub protein: Protein,
    /// Feature qualifiers added to the GenBank `CDS` feature built from
    /// this `Gene`.
    pub qualifiers: BTreeMap<String, Vec<String>>,
    /// Per-gene cluster probability override; if set, takes precedence over
    /// per-domain probabilities when computing `average_probability` and
    /// `maximum_probability`.
    pub probability: Option<f64>,
}

impl Gene {
    /// The identifier of the gene (same as the protein identifier).
    pub fn id(&self) -> &str {
        &self.protein.id
    }

    /// The average of domain probabilities of being in a cluster.
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

    /// The highest of domain probabilities of being in a cluster.
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

    /// Predict the function(s) of the gene from its domain annotations.
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

    /// Copy the current gene and assign it a different protein.
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

    /// Copy the current gene and assign it a different probability.
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

/// A sequence of contiguous genes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Cluster {
    /// The identifier of the gene cluster.
    pub id: String,
    /// A list of the genes belonging to this gene cluster.
    pub genes: Vec<Gene>,
    /// The putative type of this gene cluster, according to similarity in
    /// domain composition with curated clusters.
    pub cluster_type: Option<ClusterType>,
    /// The probability with which each cluster type was identified.
    pub type_probabilities: BTreeMap<String, f64>,
}

impl Cluster {
    /// Create a new cluster from a list of contiguous genes.
    pub fn new(id: impl Into<String>, genes: Vec<Gene>) -> Self {
        Self {
            id: id.into(),
            genes,
            cluster_type: None,
            type_probabilities: BTreeMap::new(),
        }
    }

    /// The identifier of the sequence this cluster was found in.
    pub fn source_id(&self) -> &str {
        &self.genes[0].source_id
    }

    /// The start of this cluster in the source sequence.
    pub fn start(&self) -> i64 {
        self.genes.iter().map(|g| g.start).min().unwrap_or(0)
    }

    /// The end of this cluster in the source sequence.
    pub fn end(&self) -> i64 {
        self.genes.iter().map(|g| g.end).max().unwrap_or(0)
    }

    /// The average of protein probabilities of being biosynthetic.
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

    /// The highest of protein probabilities of being biosynthetic.
    pub fn maximum_probability(&self) -> Option<f64> {
        self.genes
            .iter()
            .filter_map(|g| g.maximum_probability())
            .fold(None, |acc, p| Some(acc.map_or(p, |a: f64| a.max(p))))
    }

    /// Compute weighted domain composition with respect to `all_possible`.
    ///
    /// # Arguments
    ///
    /// * `all_possible` — A sequence containing all domain names to consider
    ///   when computing domain composition for the cluster. If `None`, only
    ///   domains within the cluster are taken into account.
    /// * `normalize` — Normalize the composition vector so that it sums to 1.
    /// * `minlog_weights` — Compute the weight for each domain as
    ///   `-log10(v)` (where `v` is either the p-value or the i-evalue,
    ///   depending on `use_pvalue`). Otherwise compute it as `1 - v`.
    /// * `use_pvalue` — Compute composition weights using the p-value of each
    ///   domain, instead of the i-evalue.
    ///
    /// # Returns
    ///
    /// A numerical vector containing the relative domain composition of the
    /// gene cluster.
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
