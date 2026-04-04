//! HMMER domain annotation module.
//!
//! Runs hmmsearch against HMM profile databases to annotate genes with
//! protein domains. Uses the pure-Rust hmmer crate as backend.

use std::collections::{BTreeMap, HashSet};
use std::io::BufRead;
use std::path::PathBuf;

use anyhow::{Context, Result};
use log::{debug, info};

use crate::interpro::InterPro;
use crate::model::{Domain, Gene};

/// Metadata for an HMM library (loaded from .ini config).
#[derive(Debug, Clone)]
pub struct HMM {
    pub id: String,
    pub version: String,
    pub url: String,
    pub path: PathBuf,
    pub size: Option<usize>,
    pub relabel_with: Option<String>,
    pub md5: Option<String>,
}

impl HMM {
    /// Load HMM metadata from an INI-style config string.
    pub fn from_ini(ini: &str, base_path: &std::path::Path) -> Result<Self> {
        let mut id = String::new();
        let mut version = String::new();
        let mut url = String::new();
        let mut size = None;
        let mut relabel_with = None;
        let mut md5 = None;

        for line in ini.lines() {
            let line = line.trim();
            if line.starts_with('[') || line.is_empty() {
                continue;
            }
            if let Some((key, val)) = line.split_once('=') {
                let key = key.trim();
                let val = val.trim();
                match key {
                    "id" => id = val.to_string(),
                    "version" => version = val.to_string(),
                    "url" => url = val.to_string(),
                    "size" => size = val.parse().ok(),
                    "relabel_with" => relabel_with = Some(val.to_string()),
                    "md5" => md5 = Some(val.to_string()),
                    _ => {}
                }
            }
        }

        // Look for the .hmm file next to the .ini
        let hmm_path = base_path.join(format!("{}.hmm", id));

        Ok(HMM {
            id,
            version,
            url,
            path: hmm_path,
            size,
            relabel_with,
            md5,
        })
    }

    /// Relabel a domain accession according to the relabeling rule.
    ///
    /// For example, Pfam uses `s/(PF\d+).\d+/\1/` to strip version suffixes.
    pub fn relabel(&self, domain: &str) -> String {
        match &self.relabel_with {
            Some(pattern) => {
                // Pattern format: s/before/after/
                if let Some(stripped) = pattern.strip_prefix("s/") {
                    if let Some(idx) = stripped[..stripped.len() - 1].find('/') {
                        let before = &stripped[..idx];
                        let after = &stripped[idx + 1..stripped.len() - 1];
                        // Convert sed-style \1 backrefs to regex-style $1
                        let after = after.replace("\\1", "$1")
                            .replace("\\2", "$2")
                            .replace("\\3", "$3");
                        if let Ok(re) = regex_lite::Regex::new(before) {
                            return re.replace(domain, after.as_str()).into_owned();
                        }
                    }
                }
                domain.to_string()
            }
            None => domain.to_string(),
        }
    }
}

/// Trait for domain annotation backends.
pub trait DomainAnnotator {
    /// Annotate genes with protein domains. Domains are appended to
    /// each gene's `protein.domains` list.
    fn run(
        &self,
        genes: &mut [Gene],
        interpro: &InterPro,
        progress: Option<&dyn Fn(usize, usize)>,
    ) -> Result<()>;
}

/// Domain annotator using the pure-Rust hmmer crate.
pub struct PyHMMER {
    pub hmm: HMM,
    pub cpus: Option<usize>,
    pub whitelist: Option<HashSet<String>>,
}

impl PyHMMER {
    pub fn new(hmm: HMM) -> Self {
        Self {
            hmm,
            cpus: None,
            whitelist: None,
        }
    }

    pub fn with_cpus(mut self, cpus: usize) -> Self {
        self.cpus = Some(cpus);
        self
    }

    pub fn with_whitelist(mut self, whitelist: HashSet<String>) -> Self {
        self.whitelist = Some(whitelist);
        self
    }

    /// Load and optionally filter HMM profiles from the database file.
    fn load_hmms(&self) -> Result<Vec<hmmer::core::Hmm>> {
        let file = std::fs::File::open(&self.hmm.path)
            .with_context(|| format!("opening HMM file: {}", self.hmm.path.display()))?;
        let mut reader = std::io::BufReader::new(file);

        // Try reading as gzip first (magic bytes)
        let hmms = self.read_hmms_maybe_compressed(&mut reader)?;

        // Filter by whitelist if set
        let filtered: Vec<hmmer::core::Hmm> = match &self.whitelist {
            Some(wl) => hmms
                .into_iter()
                .filter(|h| {
                    match &h.acc {
                        Some(acc) => {
                            let relabeled = self.hmm.relabel(acc);
                            wl.contains(&relabeled)
                        }
                        // Keep HMMs without accession (they pass by default)
                        None => true,
                    }
                })
                .collect(),
            None => hmms,
        };

        info!(
            "Loaded {} HMM profiles from {} ({})",
            filtered.len(),
            self.hmm.id,
            self.hmm.path.display()
        );
        Ok(filtered)
    }

    fn read_hmms_maybe_compressed(
        &self,
        reader: &mut std::io::BufReader<std::fs::File>,
    ) -> Result<Vec<hmmer::core::Hmm>> {
        // Check for gzip magic bytes
        let buf = reader.fill_buf()?;
        if buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b {
            let decoder = flate2::read::GzDecoder::new(reader);
            let mut gz_reader = std::io::BufReader::new(decoder);
            hmmer::io::hmm_file::read_all_hmms(&mut gz_reader)
                .map_err(|e| anyhow::anyhow!("reading gzipped HMM file: {}", e))
        } else {
            hmmer::io::hmm_file::read_all_hmms(reader)
                .map_err(|e| anyhow::anyhow!("reading HMM file: {}", e))
        }
    }
}

impl DomainAnnotator for PyHMMER {
    fn run(
        &self,
        genes: &mut [Gene],
        interpro: &InterPro,
        progress: Option<&dyn Fn(usize, usize)>,
    ) -> Result<()> {
        if genes.is_empty() {
            return Ok(());
        }

        // Step 1: Convert gene proteins to digital sequences.
        // Name each sequence by its index so we can map hits back to genes.
        let abc = hmmer::alphabet::Alphabet::amino();
        let mut seqs = Vec::with_capacity(genes.len());
        for (idx, gene) in genes.iter().enumerate() {
            let seq = hmmer::alphabet::DigitalSequence::new(
                &idx.to_string(),
                "",
                gene.protein.seq.as_bytes(),
                &abc,
            )
            .with_context(|| {
                format!(
                    "digitizing protein sequence for gene {}",
                    gene.protein.id
                )
            })?;
            seqs.push(seq);
        }

        // Step 2: Load HMM profiles
        let hmms = self.load_hmms()?;
        let total_hmms = hmms.len();

        if let Some(cb) = progress {
            cb(0, total_hmms);
        }

        // The Z parameter for E-value calculation.
        // GECCO uses hmm.size (total # of HMMs in database) rather than
        // number of sequences. We'll recompute E-values after search.
        let z_database = self.hmm.size.unwrap_or(total_hmms) as f64;
        let z_seqs = seqs.len() as f64;

        // Step 3: Run hmmsearch for each HMM profile against all sequences
        let config = hmmer::pipeline::PipelineConfig::default();

        for (hmm_idx, hmm_profile) in hmms.iter().enumerate() {
            let (hits, _stats) = hmmer::pipeline::hmmsearch(hmm_profile, &seqs, &config);

            // Get accession (relabeled) for this HMM
            let raw_acc = hmm_profile
                .acc
                .as_deref()
                .unwrap_or(&hmm_profile.name);
            let accession = self.hmm.relabel(raw_acc);

            // Look up InterPro metadata
            let go_terms = interpro.go_terms_for(&accession);
            let go_functions = interpro.go_functions_for(&accession);
            let entry = interpro.get(&accession);

            // Step 4: Map hits back to genes
            for hit in &hits {
                let gene_idx: usize = hit
                    .name
                    .parse()
                    .with_context(|| format!("parsing gene index from hit name '{}'", hit.name))?;

                if gene_idx >= genes.len() {
                    continue;
                }

                for dom in &hit.domains {
                    debug_assert!(dom.ienv <= dom.jenv);
                    debug_assert!(dom.ievalue >= 0.0);

                    // Recompute ievalue using database Z instead of sequence count Z.
                    // Original: ievalue = dom_pval * z_seqs
                    // Corrected: ievalue = dom_pval * z_database
                    let ievalue_corrected = if z_seqs > 0.0 {
                        dom.ievalue / z_seqs * z_database
                    } else {
                        dom.ievalue
                    };

                    // Compute p-value from ievalue: pvalue = ievalue / Z
                    let pvalue = ievalue_corrected / z_database;

                    // Build qualifiers
                    let mut qualifiers = BTreeMap::new();
                    qualifiers.insert(
                        "inference".to_string(),
                        vec!["protein motif".to_string()],
                    );

                    let mut db_xref = vec![format!("{}:{}", self.hmm.id, accession)];
                    if let Some(e) = entry {
                        db_xref.push(format!("InterPro:{}", e.accession));
                    }
                    qualifiers.insert("db_xref".to_string(), db_xref);

                    qualifiers.insert(
                        "note".to_string(),
                        vec![
                            format!("e-value: {:e}", ievalue_corrected),
                            format!("p-value: {:e}", pvalue),
                        ],
                    );

                    if let Some(e) = entry {
                        qualifiers
                            .insert("function".to_string(), vec![e.name.clone()]);
                    }

                    let domain = Domain {
                        name: accession.clone(),
                        start: dom.ienv as i64,
                        end: dom.jenv as i64,
                        hmm: self.hmm.id.clone(),
                        i_evalue: ievalue_corrected,
                        pvalue,
                        probability: None,
                        cluster_weight: None,
                        go_terms: go_terms.clone(),
                        go_functions: go_functions.clone(),
                        qualifiers,
                    };

                    genes[gene_idx].protein.domains.push(domain);
                }
            }

            if let Some(cb) = progress {
                cb(hmm_idx + 1, total_hmms);
            }
        }

        debug!(
            "Annotated {} genes with domains from {} ({} profiles)",
            genes.len(),
            self.hmm.id,
            total_hmms
        );

        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Post-processing utilities (used by CLI pipeline)
// ---------------------------------------------------------------------------

/// Remove overlapping domains from a gene, keeping the one with the
/// lowest p-value for each overlapping group.
pub fn disentangle(gene: &mut Gene) {
    if gene.protein.domains.len() <= 1 {
        return;
    }

    // Sort domains by start position
    gene.protein.domains.sort_by_key(|d| d.start);

    let mut keep = vec![true; gene.protein.domains.len()];

    for i in 0..gene.protein.domains.len() {
        if !keep[i] {
            continue;
        }
        for j in (i + 1)..gene.protein.domains.len() {
            if !keep[j] {
                continue;
            }
            let di = &gene.protein.domains[i];
            let dj = &gene.protein.domains[j];

            // Check overlap
            if dj.start <= di.end && di.start <= dj.end {
                // Remove the one with higher p-value
                if di.pvalue <= dj.pvalue {
                    keep[j] = false;
                } else {
                    keep[i] = false;
                    break;
                }
            }
        }
    }

    let mut idx = 0;
    gene.protein.domains.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });
}

/// Filter domains by e-value threshold. Removes domains with
/// `i_evalue >= e_filter`.
pub fn filter_by_evalue(genes: &mut [Gene], e_filter: f64) {
    for gene in genes.iter_mut() {
        gene.protein.domains.retain(|d| d.i_evalue < e_filter);
    }
}

/// Filter domains by p-value threshold. Removes domains with
/// `pvalue >= p_filter`.
pub fn filter_by_pvalue(genes: &mut [Gene], p_filter: f64) {
    for gene in genes.iter_mut() {
        gene.protein.domains.retain(|d| d.pvalue < p_filter);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_relabel_pfam() {
        let hmm = HMM {
            id: "Pfam".to_string(),
            version: "35.0".to_string(),
            url: String::new(),
            path: PathBuf::new(),
            size: Some(2766),
            relabel_with: Some("s/(PF\\d+).\\d+/\\1/".to_string()),
            md5: None,
        };

        assert_eq!(hmm.relabel("PF00001.21"), "PF00001");
        assert_eq!(hmm.relabel("PF12345.3"), "PF12345");
        // No match → unchanged
        assert_eq!(hmm.relabel("TIGR00001"), "TIGR00001");
    }

    #[test]
    fn test_relabel_none() {
        let hmm = HMM {
            id: "Tigrfam".to_string(),
            version: "1.0".to_string(),
            url: String::new(),
            path: PathBuf::new(),
            size: None,
            relabel_with: None,
            md5: None,
        };
        assert_eq!(hmm.relabel("TIGR00001"), "TIGR00001");
    }

    #[test]
    fn test_disentangle_no_overlap() {
        use crate::model::{Gene, Protein, Strand};

        let mut gene = Gene {
            source_id: "seq1".into(),
            start: 1,
            end: 300,
            strand: Strand::Coding,
            protein: Protein {
                id: "gene1".into(),
                seq: String::new(),
                domains: vec![
                    Domain {
                        name: "PF00001".into(),
                        start: 1,
                        end: 50,
                        hmm: "Pfam".into(),
                        i_evalue: 1e-10,
                        pvalue: 1e-14,
                        probability: None,
                        cluster_weight: None,
                        go_terms: vec![],
                        go_functions: vec![],
                        qualifiers: BTreeMap::new(),
                    },
                    Domain {
                        name: "PF00002".into(),
                        start: 60,
                        end: 100,
                        hmm: "Pfam".into(),
                        i_evalue: 1e-5,
                        pvalue: 1e-9,
                        probability: None,
                        cluster_weight: None,
                        go_terms: vec![],
                        go_functions: vec![],
                        qualifiers: BTreeMap::new(),
                    },
                ],
            },
            qualifiers: BTreeMap::new(),
            probability: None,
        };

        disentangle(&mut gene);
        assert_eq!(gene.protein.domains.len(), 2);
    }

    #[test]
    fn test_disentangle_overlap_keeps_best() {
        use crate::model::{Gene, Protein, Strand};

        let mut gene = Gene {
            source_id: "seq1".into(),
            start: 1,
            end: 300,
            strand: Strand::Coding,
            protein: Protein {
                id: "gene1".into(),
                seq: String::new(),
                domains: vec![
                    Domain {
                        name: "PF00001".into(),
                        start: 1,
                        end: 50,
                        hmm: "Pfam".into(),
                        i_evalue: 1e-5,
                        pvalue: 1e-9,
                        probability: None,
                        cluster_weight: None,
                        go_terms: vec![],
                        go_functions: vec![],
                        qualifiers: BTreeMap::new(),
                    },
                    Domain {
                        name: "PF00002".into(),
                        start: 30,
                        end: 80,
                        hmm: "Pfam".into(),
                        i_evalue: 1e-10,
                        pvalue: 1e-14,
                        probability: None,
                        cluster_weight: None,
                        go_terms: vec![],
                        go_functions: vec![],
                        qualifiers: BTreeMap::new(),
                    },
                ],
            },
            qualifiers: BTreeMap::new(),
            probability: None,
        };

        disentangle(&mut gene);
        assert_eq!(gene.protein.domains.len(), 1);
        // Should keep PF00002 (lower p-value)
        assert_eq!(gene.protein.domains[0].name, "PF00002");
    }

    #[test]
    fn test_hmm_from_ini() {
        let ini = r#"
[hmm]
id = Pfam
version = 35.0
size = 2766
url = ftp://example.com/Pfam-A.hmm.gz
relabel_with = s/(PF\d+).\d+/\1/
md5 = abc123
"#;
        let hmm = HMM::from_ini(ini, std::path::Path::new("/data/hmms")).unwrap();
        assert_eq!(hmm.id, "Pfam");
        assert_eq!(hmm.version, "35.0");
        assert_eq!(hmm.size, Some(2766));
        assert_eq!(hmm.path, PathBuf::from("/data/hmms/Pfam.hmm"));
        assert!(hmm.relabel_with.is_some());
    }
}
