//! ORF finding: gene prediction from DNA sequences.
//!
//! Provides both Prodigal-based gene finding (via the `prodigal` crate) and
//! extraction of existing CDS annotations from GenBank records.

use std::collections::BTreeMap;

use anyhow::{Context, Result};
use log::{debug, info};

use crate::model::{Gene, Protein, Strand};

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

/// Prodigal-based gene finder using the `prodigal` crate.
pub struct ProdigalFinder {
    pub metagenome: bool,
    pub mask: bool,
    pub closed_ends: bool,
    pub translation_table: Option<u8>,
    pub cpus: usize,
}

impl Default for ProdigalFinder {
    fn default() -> Self {
        Self {
            metagenome: true,
            mask: false,
            closed_ends: false,
            translation_table: Some(11),
            cpus: 0,
        }
    }
}

impl ORFFinder for ProdigalFinder {
    fn find_genes(&self, records: &[SeqRecord]) -> Result<Vec<Gene>> {
        let config = prodigal_rs::ProdigalConfig {
            translation_table: self.translation_table.unwrap_or(11),
            closed_ends: self.closed_ends,
            mask_n_runs: self.mask,
            ..Default::default()
        };

        let tt = self.translation_table.unwrap_or(11) as u32;
        let mut all_genes = Vec::new();

        // Create a reusable predictor for metagenomic mode (caches models + buffers)
        let meta_predictor = if self.metagenome {
            Some(prodigal_rs::MetaPredictor::with_config(config.clone())
                .context("creating MetaPredictor")?)
        } else {
            None
        };

        for record in records {
            if record.seq.len() < 20 {
                debug!("Skipping short sequence {} ({} bp)", record.id, record.seq.len());
                continue;
            }

            let predicted = if let Some(ref predictor) = meta_predictor {
                predictor.predict(record.seq.as_bytes())
            } else {
                let training = prodigal_rs::train_with(record.seq.as_bytes(), &config)
                    .with_context(|| {
                        format!("training Prodigal on sequence {} ({} bp)", record.id, record.seq.len())
                    })?;
                prodigal_rs::predict_with(record.seq.as_bytes(), &training, &config)
            }
            .with_context(|| {
                format!(
                    "running gene prediction on sequence {} ({} bp)",
                    record.id,
                    record.seq.len()
                )
            })?;

            for (j, orf) in predicted.iter().enumerate() {
                let begin = orf.begin;
                let end = orf.end;
                let start = begin.min(end);
                let stop = begin.max(end);

                let strand = match orf.strand {
                    prodigal_rs::Strand::Forward => Strand::Coding,
                    prodigal_rs::Strand::Reverse => Strand::Reverse,
                };

                let nuc_seq = extract_subseq(&record.seq, start, stop, strand);
                let prot_seq = translate_dna(&nuc_seq, tt);

                let protein_id = format!("{}_{}", record.id, j + 1);
                let protein = Protein::new(&protein_id, prot_seq);

                let mut qualifiers = BTreeMap::new();
                qualifiers.insert(
                    "inference".to_string(),
                    vec!["ab initio prediction:Prodigal".to_string()],
                );
                qualifiers.insert(
                    "transl_table".to_string(),
                    vec![tt.to_string()],
                );

                all_genes.push(Gene {
                    source_id: record.id.clone(),
                    start: start as i64,
                    end: stop as i64,
                    strand,
                    protein,
                    qualifiers,
                    probability: None,
                });
            }

            info!(
                "Found {} genes in {} ({} bp)",
                predicted.len(),
                record.id,
                record.seq.len()
            );
        }

        Ok(all_genes)
    }
}

/// Extract a subsequence from DNA, reverse-complementing if on reverse strand.
fn extract_subseq(seq: &str, start: usize, end: usize, strand: Strand) -> Vec<u8> {
    // Coordinates are 1-based inclusive
    let s = (start.saturating_sub(1)).min(seq.len());
    let e = end.min(seq.len());
    let subseq = seq[s..e].as_bytes();

    match strand {
        Strand::Coding => subseq.to_vec(),
        Strand::Reverse => reverse_complement(subseq),
    }
}

/// Reverse complement a DNA sequence.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            other => other,
        })
        .collect()
}

/// Translate a DNA sequence to protein using the given translation table.
///
/// Supports table 11 (standard bacterial/archaeal) and table 4 (Mycoplasma).
/// Translates in reading frame 0, stopping at stop codons or end of sequence.
fn translate_dna(dna: &[u8], table: u32) -> String {
    let mut protein = String::with_capacity(dna.len() / 3);

    for codon in dna.chunks(3) {
        if codon.len() < 3 {
            break;
        }
        let aa = translate_codon(codon, table);
        if aa == '*' {
            break;
        }
        protein.push(aa);
    }

    protein
}

/// Translate a single codon to an amino acid.
fn translate_codon(codon: &[u8], table: u32) -> char {
    let c0 = codon[0].to_ascii_uppercase();
    let c1 = codon[1].to_ascii_uppercase();
    let c2 = codon[2].to_ascii_uppercase();

    // Standard genetic code (table 1/11 share the same codon→AA mapping;
    // table 11 only differs in alternative start codons, which Prodigal
    // already handles). Table 4 has UGA=Trp instead of Stop.
    match (c0, c1, c2) {
        // Phe
        (b'T', b'T', b'T') | (b'T', b'T', b'C') => 'F',
        // Leu
        (b'T', b'T', b'A') | (b'T', b'T', b'G') => 'L',
        (b'C', b'T', _) => 'L',
        // Ile
        (b'A', b'T', b'T') | (b'A', b'T', b'C') | (b'A', b'T', b'A') => 'I',
        // Met
        (b'A', b'T', b'G') => 'M',
        // Val
        (b'G', b'T', _) => 'V',
        // Ser
        (b'T', b'C', _) => 'S',
        (b'A', b'G', b'T') | (b'A', b'G', b'C') => 'S',
        // Pro
        (b'C', b'C', _) => 'P',
        // Thr
        (b'A', b'C', _) => 'T',
        // Ala
        (b'G', b'C', _) => 'A',
        // Tyr
        (b'T', b'A', b'T') | (b'T', b'A', b'C') => 'Y',
        // Stop (TAA, TAG)
        (b'T', b'A', b'A') | (b'T', b'A', b'G') => '*',
        // Stop or Trp (TGA) — table 4: Trp
        (b'T', b'G', b'A') => {
            if table == 4 {
                'W'
            } else {
                '*'
            }
        }
        // His
        (b'C', b'A', b'T') | (b'C', b'A', b'C') => 'H',
        // Gln
        (b'C', b'A', b'A') | (b'C', b'A', b'G') => 'Q',
        // Asn
        (b'A', b'A', b'T') | (b'A', b'A', b'C') => 'N',
        // Lys
        (b'A', b'A', b'A') | (b'A', b'A', b'G') => 'K',
        // Asp
        (b'G', b'A', b'T') | (b'G', b'A', b'C') => 'D',
        // Glu
        (b'G', b'A', b'A') | (b'G', b'A', b'G') => 'E',
        // Cys
        (b'T', b'G', b'T') | (b'T', b'G', b'C') => 'C',
        // Trp
        (b'T', b'G', b'G') => 'W',
        // Arg
        (b'C', b'G', _) => 'R',
        (b'A', b'G', b'A') | (b'A', b'G', b'G') => 'R',
        // Gly
        (b'G', b'G', _) => 'G',
        // Unknown
        _ => 'X',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(
            reverse_complement(b"ATGCNN"),
            b"NNGCAT".to_vec()
        );
        assert_eq!(reverse_complement(b""), b"".to_vec());
    }

    #[test]
    fn test_translate_simple() {
        // ATG=M, AAA=K, TAA=*
        assert_eq!(translate_dna(b"ATGAAATAA", 11), "MK");
        // With trailing partial codon
        assert_eq!(translate_dna(b"ATGAAATA", 11), "MK");
    }

    #[test]
    fn test_translate_all_amino_acids() {
        // Test representative codons for each amino acid
        assert_eq!(translate_codon(b"TTT", 11), 'F');
        assert_eq!(translate_codon(b"TTA", 11), 'L');
        assert_eq!(translate_codon(b"ATT", 11), 'I');
        assert_eq!(translate_codon(b"ATG", 11), 'M');
        assert_eq!(translate_codon(b"GTT", 11), 'V');
        assert_eq!(translate_codon(b"TCT", 11), 'S');
        assert_eq!(translate_codon(b"CCT", 11), 'P');
        assert_eq!(translate_codon(b"ACT", 11), 'T');
        assert_eq!(translate_codon(b"GCT", 11), 'A');
        assert_eq!(translate_codon(b"TAT", 11), 'Y');
        assert_eq!(translate_codon(b"CAT", 11), 'H');
        assert_eq!(translate_codon(b"CAA", 11), 'Q');
        assert_eq!(translate_codon(b"AAT", 11), 'N');
        assert_eq!(translate_codon(b"AAA", 11), 'K');
        assert_eq!(translate_codon(b"GAT", 11), 'D');
        assert_eq!(translate_codon(b"GAA", 11), 'E');
        assert_eq!(translate_codon(b"TGT", 11), 'C');
        assert_eq!(translate_codon(b"TGG", 11), 'W');
        assert_eq!(translate_codon(b"CGT", 11), 'R');
        assert_eq!(translate_codon(b"GGT", 11), 'G');
    }

    #[test]
    fn test_translate_table4_tga_is_trp() {
        assert_eq!(translate_codon(b"TGA", 11), '*');
        assert_eq!(translate_codon(b"TGA", 4), 'W');
    }

    #[test]
    fn test_extract_subseq_forward() {
        let seq = "ATGAAACCCGGGTTT";
        // 1-based coords: positions 1-6 = "ATGAAA"
        let sub = extract_subseq(seq, 1, 6, Strand::Coding);
        assert_eq!(sub, b"ATGAAA");
    }

    #[test]
    fn test_extract_subseq_reverse() {
        let seq = "ATGAAACCCGGGTTT";
        // 1-based coords: positions 1-6 = "ATGAAA", rev-comp = "TTTCAT"
        let sub = extract_subseq(seq, 1, 6, Strand::Reverse);
        assert_eq!(sub, b"TTTCAT");
    }
}
