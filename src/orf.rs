//! Generic protocol for ORF detection in DNA sequences.
//!
//! Provides both Prodigal-based gene finding (via the `prodigal` crate) and
//! extraction of existing CDS annotations from GenBank records.

use std::collections::BTreeMap;
use std::sync::Arc;

use anyhow::{Context, Result};
use log::debug;

use crate::model::{Gene, Protein, Strand};
use crate::output::{RunOutput, StdioOutput};

/// Record representing a DNA sequence.
#[derive(Debug, Clone)]
pub struct SeqRecord {
    pub id: String,
    pub seq: String,
}

/// An abstract base trait providing a generic ORF finder.
pub trait ORFFinder {
    /// Find all genes from a DNA sequence.
    fn find_genes(&self, records: &[SeqRecord]) -> Result<Vec<Gene>>;
}

/// An `ORFFinder` that simply extracts CDS annotations from records.
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

/// An `ORFFinder` that uses the `prodigal-rs` bindings to Prodigal.
///
/// Prodigal is a fast and reliable protein-coding gene prediction for
/// prokaryotic genomes, with support for draft genomes and metagenomes.
///
/// # References
///
/// Doug Hyatt, Gwo-Liang Chen, Philip F. LoCascio, Miriam L. Land,
/// Frank W. Larimer and Loren J. Hauser. "Prodigal: Prokaryotic Gene
/// Recognition and Translation Initiation Site Identification",
/// BMC Bioinformatics 11 (8 March 2010), p119.
/// <https://doi.org/10.1186/1471-2105-11-119>
pub struct ProdigalFinder {
    /// Whether or not to run Prodigal in metagenome mode (default `true`).
    pub metagenome: bool,
    /// Whether or not to mask genes running across regions containing
    /// unknown nucleotides (default `false`).
    pub mask: bool,
    pub closed_ends: bool,
    pub translation_table: Option<u8>,
    /// The number of threads to use to run Prodigal in parallel. Pass `0`
    /// to use the number of CPUs on the machine.
    pub cpus: usize,
    pub thread_pool: Option<Arc<rayon::ThreadPool>>,
}

impl Default for ProdigalFinder {
    fn default() -> Self {
        Self {
            metagenome: true,
            mask: false,
            closed_ends: false,
            translation_table: Some(11),
            cpus: 0,
            thread_pool: None,
        }
    }
}

impl ORFFinder for ProdigalFinder {
    /// Find all genes contained in a sequence of DNA records.
    fn find_genes(&self, records: &[SeqRecord]) -> Result<Vec<Gene>> {
        self.find_genes_with_output(records, &StdioOutput)
    }
}

impl ProdigalFinder {
    pub fn find_genes_with_output(
        &self,
        records: &[SeqRecord],
        output: &dyn RunOutput,
    ) -> Result<Vec<Gene>> {
        let config = prodigal_rs::ProdigalConfig {
            translation_table: self.translation_table.unwrap_or(11),
            closed_ends: self.closed_ends,
            mask_n_runs: self.mask,
            ..Default::default()
        };

        let thread_pool = if self.cpus > 0 {
            Some(Arc::new(
                rayon::ThreadPoolBuilder::new()
                    .num_threads(self.cpus)
                    .stack_size(prodigal_rs::META_PREDICTOR_STACK_SIZE)
                    .build()
                    .context("building Prodigal thread pool")?,
            ))
        } else if let Some(pool) = &self.thread_pool {
            Some(Arc::clone(pool))
        } else {
            None
        };

        // Create a reusable predictor for metagenomic mode (caches models + buffers).
        let meta_predictor = if self.metagenome {
            let predictor = if let Some(pool) = &thread_pool {
                prodigal_rs::MetaPredictor::with_config_and_thread_pool(
                    config.clone(),
                    Arc::clone(pool),
                )
            } else {
                prodigal_rs::MetaPredictor::with_config(config.clone())
            };
            Some(predictor.context("creating MetaPredictor")?)
        } else {
            None
        };

        let prediction_results = if let Some(predictor) = meta_predictor.as_ref() {
            predict_records_batch(records, predictor).with_context(|| {
                format!(
                    "running metagenomic gene prediction on {} sequence(s)",
                    records.len()
                )
            })?
        } else {
            records
                .iter()
                .map(|record| predict_record(record, &config))
                .collect::<Result<Vec<_>>>()?
        };

        let mut all_genes = Vec::new();
        for (record, predicted) in prediction_results {
            if record.seq.len() < 20 {
                debug!(
                    "Skipping short sequence {} ({} bp)",
                    record.id,
                    record.seq.len()
                );
                continue;
            }

            for (j, orf) in predicted.iter().enumerate() {
                let strand = match orf.strand {
                    prodigal_rs::Strand::Forward => Strand::Coding,
                    prodigal_rs::Strand::Reverse => Strand::Reverse,
                };

                let begin = orf.begin;
                let end = orf.end;
                let start = begin.min(end);
                let stop = begin.max(end);

                let orf_tt = orf.translation_table as u32;
                let nuc_seq = extract_subseq(&record.seq, start, stop, strand);
                let mut prot_seq = translate_dna(&nuc_seq, orf_tt);
                if orf.start_codon != prodigal_rs::StartCodon::Edge && !prot_seq.is_empty() {
                    prot_seq.replace_range(0..1, "M");
                }

                let protein_id = format!("{}_{}", record.id, j + 1);
                let protein = Protein::new(&protein_id, prot_seq);

                let mut qualifiers = BTreeMap::new();
                qualifiers.insert(
                    "inference".to_string(),
                    vec!["ab initio prediction:Prodigal".to_string()],
                );
                qualifiers.insert("transl_table".to_string(), vec![orf_tt.to_string()]);

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

            output.stderr(format_args!(
                "Found {} genes in {} ({} bp)",
                predicted.len(),
                record.id,
                record.seq.len()
            ))?;
        }

        Ok(all_genes)
    }
}

fn predict_record<'a>(
    record: &'a SeqRecord,
    config: &prodigal_rs::ProdigalConfig,
) -> Result<(&'a SeqRecord, Vec<prodigal_rs::PredictedGene>)> {
    if record.seq.len() < 20 {
        return Ok((record, Vec::new()));
    }

    let training = prodigal_rs::train_with(record.seq.as_bytes(), config).with_context(|| {
        format!(
            "training Prodigal on sequence {} ({} bp)",
            record.id,
            record.seq.len()
        )
    })?;
    let predicted = prodigal_rs::predict_with(record.seq.as_bytes(), &training, config)
        .with_context(|| {
            format!(
                "running gene prediction on sequence {} ({} bp)",
                record.id,
                record.seq.len()
            )
        })?;

    Ok((record, predicted))
}

fn predict_records_batch<'a>(
    records: &'a [SeqRecord],
    predictor: &prodigal_rs::MetaPredictor,
) -> Result<Vec<(&'a SeqRecord, Vec<prodigal_rs::PredictedGene>)>> {
    let predicted_records = records
        .iter()
        .filter(|record| record.seq.len() >= 20)
        .collect::<Vec<_>>();
    let seqs = predicted_records
        .iter()
        .map(|record| record.seq.as_bytes())
        .collect::<Vec<_>>();
    let mut predicted = predictor.predict_batch(&seqs)?.into_iter();

    records
        .iter()
        .map(|record| {
            if record.seq.len() < 20 {
                Ok((record, Vec::new()))
            } else {
                let genes = predicted.next().with_context(|| {
                    format!("missing Prodigal batch result for sequence {}", record.id)
                })?;
                Ok((record, genes))
            }
        })
        .collect()
}

/// Extract a subsequence from DNA, reverse-complementing if on reverse strand.
fn extract_subseq(seq: &str, start: usize, end: usize, strand: Strand) -> Vec<u8> {
    // Coordinates are 1-based inclusive
    let s = (start.saturating_sub(1)).min(seq.len());
    let e = end.min(seq.len());
    let subseq = &seq.as_bytes()[s..e];

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
/// Translates in reading frame 0, preserving the terminal stop codon.
fn translate_dna(dna: &[u8], table: u32) -> String {
    let mut protein = String::with_capacity(dna.len() / 3);

    for codon in dna.chunks(3) {
        if codon.len() < 3 {
            break;
        }
        let aa = translate_codon(codon, table);
        protein.push(aa);
        if aa == '*' {
            break;
        }
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
        assert_eq!(reverse_complement(b"ATGCNN"), b"NNGCAT".to_vec());
        assert_eq!(reverse_complement(b""), b"".to_vec());
    }

    #[test]
    fn test_translate_simple() {
        // ATG=M, AAA=K, TAA=*
        assert_eq!(translate_dna(b"ATGAAATAA", 11), "MK*");
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
