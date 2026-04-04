//! InterPro metadata for domain annotations.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

/// A single Gene Ontology term.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GOTerm {
    pub accession: String,
    pub name: String,
    #[serde(default)]
    pub namespace: String,
}

/// A single entry in the InterPro database.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InterProEntry {
    pub accession: String,
    pub members: Vec<String>,
    pub name: String,
    pub databases: Vec<String>,
    #[serde(rename = "type")]
    pub entry_type: String,
    pub go_terms: Vec<GOTerm>,
    pub go_functions: Vec<GOFunctionRef>,
}

/// GO function reference (without namespace, added at load time).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GOFunctionRef {
    pub accession: String,
    pub name: String,
}

/// A subset of the InterPro database.
pub struct InterPro {
    pub entries: Vec<InterProEntry>,
    /// Maps member accessions → InterPro entry index.
    pub by_accession: HashMap<String, usize>,
}

impl InterPro {
    /// Load InterPro metadata from a JSON byte slice.
    pub fn from_json(data: &[u8]) -> anyhow::Result<Self> {
        let entries: Vec<InterProEntry> = serde_json::from_slice(data)?;
        let mut by_accession = HashMap::new();
        for (idx, entry) in entries.iter().enumerate() {
            for member in &entry.members {
                by_accession.insert(member.clone(), idx);
            }
        }
        Ok(Self {
            entries,
            by_accession,
        })
    }

    /// Look up an entry by member domain accession.
    pub fn get(&self, accession: &str) -> Option<&InterProEntry> {
        self.by_accession
            .get(accession)
            .map(|&idx| &self.entries[idx])
    }

    /// Get GO terms for a domain accession.
    pub fn go_terms_for(&self, accession: &str) -> Vec<GOTerm> {
        match self.get(accession) {
            Some(entry) => entry.go_terms.clone(),
            None => Vec::new(),
        }
    }

    /// Get GO functions for a domain accession (with namespace filled in).
    pub fn go_functions_for(&self, accession: &str) -> Vec<GOTerm> {
        match self.get(accession) {
            Some(entry) => entry
                .go_functions
                .iter()
                .map(|gf| GOTerm {
                    accession: gf.accession.clone(),
                    name: gf.name.clone(),
                    namespace: "molecular_function".to_string(),
                })
                .collect(),
            None => Vec::new(),
        }
    }
}
