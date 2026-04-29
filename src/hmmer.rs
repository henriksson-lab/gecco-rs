//! HMMER domain annotation module.
//!
//! Runs hmmsearch against HMM profile databases to annotate genes with
//! protein domains. Uses the pure-Rust hmmer crate as backend.

use std::collections::{BTreeMap, HashSet};
use std::io::BufRead;
#[cfg(feature = "bundled-data")]
use std::io::Cursor;
use std::path::PathBuf;
use std::sync::Arc;

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
                        let after = after
                            .replace("\\1", "$1")
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

/// A loaded HMM database: metadata + parsed profiles ready for search.
pub struct LoadedHMM {
    pub meta: HMM,
    pub profiles: Vec<hmmer::Hmm>,
}

/// Domain annotator using the pure-Rust hmmer crate.
pub struct PyHMMER {
    pub hmm: HMM,
    pub cpus: Option<usize>,
    pub whitelist: Option<HashSet<String>>,
    pub thread_pool: Option<Arc<rayon::ThreadPool>>,
}

impl PyHMMER {
    pub fn new(hmm: HMM) -> Self {
        Self {
            hmm,
            cpus: None,
            whitelist: None,
            thread_pool: None,
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

    pub fn with_thread_pool(mut self, thread_pool: Arc<rayon::ThreadPool>) -> Self {
        self.thread_pool = Some(thread_pool);
        self
    }

    /// Load and optionally filter HMM profiles from the database file.
    /// Supports text `.hmm`, binary `.h3m`, and gzipped variants.
    pub fn load_hmms(&self) -> Result<Vec<hmmer::Hmm>> {
        #[cfg(feature = "bundled-data")]
        if crate::bundled_data::is_hmm_path(&self.hmm.path) {
            let mut reader = Cursor::new(crate::bundled_data::pfam_h3m());
            let mut hmms = Vec::new();
            while let Some(hmm) = hmmer::hmmfile_binary::read_binary_hmm(&mut reader)
                .map_err(|e| anyhow::anyhow!("reading bundled binary HMM data: {}", e))?
            {
                hmms.push(hmm);
            }
            return self.filter_hmms(hmms);
        }

        let path_str = self.hmm.path.to_string_lossy();

        // Try binary .h3m format first
        if path_str.ends_with(".h3m") {
            let hmms = hmmer::hmmfile_binary::read_binary_hmm_file(&self.hmm.path)
                .map_err(|e| anyhow::anyhow!("reading binary HMM file: {}", e))?;
            return self.filter_hmms(hmms);
        }

        let hmms = self.read_hmms_maybe_compressed()?;

        self.filter_hmms(hmms)
    }

    fn filter_hmms(&self, hmms: Vec<hmmer::Hmm>) -> Result<Vec<hmmer::Hmm>> {
        let filtered: Vec<hmmer::Hmm> = match &self.whitelist {
            Some(wl) => hmms
                .into_iter()
                .filter(|h| match &h.acc {
                    Some(acc) => {
                        let relabeled = self.hmm.relabel(acc);
                        wl.contains(&relabeled)
                    }
                    None => true,
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

    fn read_hmms_maybe_compressed(&self) -> Result<Vec<hmmer::Hmm>> {
        let file = std::fs::File::open(&self.hmm.path)
            .with_context(|| format!("opening HMM file: {}", self.hmm.path.display()))?;
        let mut reader = std::io::BufReader::new(file);

        let buf = reader.fill_buf()?;
        if buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b {
            let file = std::fs::File::open(&self.hmm.path)
                .with_context(|| format!("opening HMM file: {}", self.hmm.path.display()))?;
            let decoder = flate2::read::GzDecoder::new(file);
            let gz_reader = std::io::BufReader::new(decoder);
            hmmer::hmmfile::read_hmms(gz_reader)
                .map_err(|e| anyhow::anyhow!("reading gzipped HMM file: {}", e))
        } else {
            hmmer::hmmfile::read_hmms(reader)
                .map_err(|e| anyhow::anyhow!("reading HMM file: {}", e))
        }
    }
}

impl PyHMMER {
    /// Run domain annotation using pre-loaded HMM profiles.
    pub fn run_with_profiles(
        &self,
        genes: &mut [Gene],
        interpro: &InterPro,
        profiles: &[hmmer::Hmm],
        progress: Option<&dyn Fn(usize, usize)>,
    ) -> Result<()> {
        if genes.is_empty() {
            return Ok(());
        }

        let seqs = digitize_genes(genes);

        let total_hmms = profiles.len();

        self.run_search(genes, interpro, profiles, total_hmms, &seqs, progress)
    }

    fn run_search(
        &self,
        genes: &mut [Gene],
        interpro: &InterPro,
        hmms: &[hmmer::Hmm],
        total_hmms: usize,
        seqs: &[hmmer::sequence::Sequence],
        progress: Option<&dyn Fn(usize, usize)>,
    ) -> Result<()> {
        if let Some(cb) = progress {
            cb(0, total_hmms);
        }

        let z_database = self.hmm.size.unwrap_or(total_hmms) as f64;

        let hmm_id = &self.hmm.id;
        let hmm_meta = &self.hmm;

        let process_profile = |hmm_profile: &hmmer::Hmm| {
            let hits = search_profile(hmm_profile, seqs, z_database);

            let raw_acc = hmm_profile.acc.as_deref().unwrap_or(&hmm_profile.name);
            let accession = hmm_meta.relabel(raw_acc);

            let go_terms = interpro.go_terms_for(&accession);
            let go_functions = interpro.go_functions_for(&accession);
            let entry = interpro.get(&accession);

            let mut results = Vec::new();
            for hit in &hits.hits {
                let gene_idx: usize = match hit.name.parse() {
                    Ok(idx) => idx,
                    Err(_) => continue,
                };

                for dom in hit.dcl.iter().filter(|dom| dom.is_reported) {
                    let pvalue = dom.lnp.exp();
                    let ievalue_corrected = pvalue * z_database;

                    let mut qualifiers = BTreeMap::new();
                    qualifiers.insert("inference".to_string(), vec!["protein motif".to_string()]);

                    let mut db_xref = vec![format!("{}:{}", hmm_id, accession)];
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
                        qualifiers.insert("function".to_string(), vec![e.name.clone()]);
                    }

                    results.push((
                        gene_idx,
                        Domain {
                            name: accession.clone(),
                            start: dom.iali,
                            end: dom.jali,
                            hmm: hmm_id.clone(),
                            i_evalue: ievalue_corrected,
                            pvalue,
                            probability: None,
                            cluster_weight: None,
                            go_terms: go_terms.clone(),
                            go_functions: go_functions.clone(),
                            qualifiers,
                        },
                    ));
                }
            }
            results
        };

        let search_serial = || hmms.iter().flat_map(process_profile).collect();

        let search_parallel = || {
            use rayon::prelude::*;
            hmms.par_iter().flat_map(&process_profile).collect()
        };

        let all_domains: Vec<(usize, Domain)> = if self.cpus == Some(1) {
            search_serial()
        } else if let Some(pool) = &self.thread_pool {
            pool.install(search_parallel)
        } else if let Some(cpus) = self.cpus.filter(|&cpus| cpus > 0) {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(cpus)
                .build()
                .context("building HMMER thread pool")?;
            pool.install(search_parallel)
        } else {
            search_parallel()
        };

        for (gene_idx, domain) in all_domains {
            if gene_idx < genes.len() {
                genes[gene_idx].protein.domains.push(domain);
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

        let seqs = digitize_genes(genes);

        let hmms = self.load_hmms()?;
        let total_hmms = hmms.len();

        self.run_search(genes, interpro, &hmms, total_hmms, &seqs, progress)
    }
}

fn digitize_genes(genes: &[Gene]) -> Vec<hmmer::sequence::Sequence> {
    let abc = hmmer::alphabet::Alphabet::amino();
    genes
        .iter()
        .enumerate()
        .map(|(idx, gene)| {
            let dsq = abc.digitize(gene.protein.seq.as_bytes());
            let n = dsq.len().saturating_sub(2);
            hmmer::sequence::Sequence {
                name: idx.to_string(),
                acc: String::new(),
                desc: gene.protein.id.clone(),
                dsq,
                n,
                l: n,
            }
        })
        .collect()
}

fn search_profile(
    hmm_profile: &hmmer::Hmm,
    seqs: &[hmmer::sequence::Sequence],
    z_database: f64,
) -> hmmer::TopHits {
    let abc = hmmer::Alphabet::new(hmm_profile.abc_type);
    let mut bg = hmmer::Bg::new(&abc);
    let mut gm = hmmer::Profile::new(hmm_profile.m, &abc);
    hmmer::profile::profile_config(hmm_profile, &bg, &mut gm, 100, hmmer::profile::P7_LOCAL);
    bg.set_filter(hmm_profile.m, &hmm_profile.compo);

    let mut om = hmmer::OProfile::convert(&gm);
    let mut pli = hmmer::Pipeline::new();
    pli.new_model(&gm);
    pli.z = z_database;
    pli.domz = z_database;
    pli.z_setby = hmmer::pipeline::ZSetBy::Option;
    pli.domz_setby = hmmer::pipeline::ZSetBy::Option;
    pli.do_alignment_display = false;

    let mut hits = hmmer::TopHits::new();
    for seq in seqs {
        pli.n_targets = 0;
        pli.n_past_msv = 0;
        pli.n_past_bias = 0;
        pli.n_past_vit = 0;
        pli.n_past_fwd = 0;
        let mut local_bg = bg.clone();
        local_bg.set_length(seq.n);
        let mut local_hits = hmmer::TopHits::new();
        if pli.run(
            &mut gm,
            &mut om,
            &local_bg,
            hmm_profile,
            seq,
            &mut local_hits,
        ) {
            hits.hits.extend(local_hits.hits.into_iter());
        }
    }
    hits.sort_by_sortkey();
    hits.threshold(&pli, z_database, z_database);
    hits
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
