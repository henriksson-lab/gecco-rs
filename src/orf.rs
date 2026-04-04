//! ORF finding: gene prediction from DNA sequences.

use anyhow::Result;

use crate::model::Gene;

/// Record representing a DNA sequence.
#[derive(Debug, Clone)]
pub struct SeqRecord {
    pub id: String,
    pub seq: String,
}

/// Trait for ORF finders.
pub trait ORFFinder {
    fn find_genes(&self, records: &[SeqRecord]) -> Result<Vec<Gene>>;
}

/// ORFFinder that extracts existing CDS annotations (analogous to CDSFinder).
pub struct CDSFinder {
    pub feature: String,
    pub translation_table: u32,
    pub locus_tag: String,
}

impl Default for CDSFinder {
    fn default() -> Self {
        Self {
            feature: "CDS".to_string(),
            translation_table: 11,
            locus_tag: "locus_tag".to_string(),
        }
    }
}

// TODO: Implement PyrodigalFinder using orphos-core crate.
// The orphos-core integration will follow the same pattern:
// take SeqRecords, run gene finding, produce Gene objects.

/// Placeholder for Prodigal-based gene finder using orphos-core.
pub struct ProdigalFinder {
    pub metagenome: bool,
    pub mask: bool,
    pub cpus: usize,
}

impl Default for ProdigalFinder {
    fn default() -> Self {
        Self {
            metagenome: true,
            mask: false,
            cpus: 0,
        }
    }
}

impl ORFFinder for ProdigalFinder {
    fn find_genes(&self, _records: &[SeqRecord]) -> Result<Vec<Gene>> {
        // TODO: integrate orphos-core
        // For now, return an error indicating this needs the orphos-core crate
        anyhow::bail!("ProdigalFinder not yet implemented — requires orphos-core integration")
    }
}
