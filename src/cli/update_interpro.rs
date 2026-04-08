//! The `gecco update-interpro` subcommand — rebuild interpro.json from upstream sources.
//!
//! Downloads the Gene Ontology (OBO format) and InterPro XML database,
//! then builds the interpro.json metadata file. Replaces the Python
//! `setup.py update_interpro` command that required the `pronto` library.

use std::collections::{HashMap, HashSet, VecDeque};
use std::io::{BufRead, BufReader, Read};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Args;
use log::info;

const GO_OBO_URL: &str = "http://purl.obolibrary.org/obo/go.obo";
const INTERPRO_XML_URL: &str =
    "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro.xml.gz";

/// Molecular function root term.
const GO_MOLECULAR_FUNCTION: &str = "GO:0003674";

// ---------------------------------------------------------------------------
// CLI args
// ---------------------------------------------------------------------------

#[derive(Args)]
pub struct UpdateInterProArgs {
    /// Output directory for interpro.json (default: ./gecco_data).
    #[arg(short, long)]
    pub output_dir: Option<PathBuf>,

    /// Force re-download even if files exist.
    #[arg(short, long)]
    pub force: bool,
}

impl UpdateInterProArgs {
    pub fn execute(&self) -> Result<()> {
        let output_dir = self
            .output_dir
            .clone()
            .unwrap_or_else(|| PathBuf::from("gecco_data"));
        std::fs::create_dir_all(&output_dir)?;

        let interpro_path = output_dir.join("interpro.json");
        if interpro_path.exists() && !self.force {
            info!(
                "{:?} already exists, skipping (use --force to re-download)",
                interpro_path
            );
            return Ok(());
        }

        // 1. Download and parse Gene Ontology
        info!("Downloading Gene Ontology from {}", GO_OBO_URL);
        let go_data = ureq_get(GO_OBO_URL)?;
        info!("Parsing OBO file ({} bytes)", go_data.len());
        let go = parse_obo(&go_data)?;
        info!("Loaded {} GO terms", go.terms.len());

        // Find top-level molecular functions (direct children of GO:0003674)
        let top_functions = go.children_of(GO_MOLECULAR_FUNCTION);
        info!(
            "Found {} top-level molecular function terms",
            top_functions.len()
        );

        // 2. Download and parse InterPro XML
        info!("Downloading InterPro XML from {}", INTERPRO_XML_URL);
        info!("(this file is ~1.5 GB compressed, download may take a while)");
        let xml_gz_data = ureq_get_large(INTERPRO_XML_URL)?;
        info!(
            "Downloaded {} MB, decompressing and parsing...",
            xml_gz_data.len() / (1024 * 1024)
        );

        let decoder = flate2::read::GzDecoder::new(xml_gz_data.as_slice());
        let entries = parse_interpro_xml(decoder, &go, &top_functions)?;
        info!("Parsed {} InterPro entries", entries.len());

        // 3. Write JSON
        let json = serde_json::to_string_pretty(&entries)?;
        std::fs::write(&interpro_path, &json)?;
        info!("Wrote interpro.json to {:?}", interpro_path);

        Ok(())
    }
}

// ---------------------------------------------------------------------------
// OBO parser — minimal parser for Gene Ontology .obo files
// ---------------------------------------------------------------------------

/// A parsed GO term.
#[derive(Debug, Clone)]
struct OboTerm {
    id: String,
    name: String,
    namespace: String,
    is_a: Vec<String>, // parent term IDs
}

/// The parsed Gene Ontology graph.
struct GeneOntology {
    terms: HashMap<String, OboTerm>,
    /// Reverse index: parent_id → set of direct child IDs.
    children: HashMap<String, HashSet<String>>,
}

impl GeneOntology {
    /// Get direct children of a term.
    fn children_of(&self, term_id: &str) -> HashSet<String> {
        self.children
            .get(term_id)
            .cloned()
            .unwrap_or_default()
    }

    /// Get all ancestor term IDs (transitive closure of is_a), including self.
    fn superclasses(&self, term_id: &str) -> HashSet<String> {
        let mut result = HashSet::new();
        let mut queue = VecDeque::new();
        queue.push_back(term_id.to_string());
        while let Some(current) = queue.pop_front() {
            if !result.insert(current.clone()) {
                continue;
            }
            if let Some(term) = self.terms.get(&current) {
                for parent in &term.is_a {
                    queue.push_back(parent.clone());
                }
            }
        }
        result
    }
}

/// Parse an OBO file into a GeneOntology.
fn parse_obo(data: &[u8]) -> Result<GeneOntology> {
    let reader = BufReader::new(data);
    let mut terms = HashMap::new();
    let mut children: HashMap<String, HashSet<String>> = HashMap::new();

    let mut in_term = false;
    let mut id = String::new();
    let mut name = String::new();
    let mut namespace = String::new();
    let mut is_a: Vec<String> = Vec::new();
    let mut is_obsolete = false;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim_end();

        if line == "[Term]" {
            // Save previous term if any
            if in_term && !id.is_empty() && !is_obsolete {
                for parent in &is_a {
                    children
                        .entry(parent.clone())
                        .or_default()
                        .insert(id.clone());
                }
                terms.insert(
                    id.clone(),
                    OboTerm {
                        id: id.clone(),
                        name: name.clone(),
                        namespace: namespace.clone(),
                        is_a: is_a.clone(),
                    },
                );
            }
            in_term = true;
            id.clear();
            name.clear();
            namespace.clear();
            is_a.clear();
            is_obsolete = false;
            continue;
        }

        if line.starts_with('[') {
            // End of [Term] block (e.g. [Typedef])
            if in_term && !id.is_empty() && !is_obsolete {
                for parent in &is_a {
                    children
                        .entry(parent.clone())
                        .or_default()
                        .insert(id.clone());
                }
                terms.insert(
                    id.clone(),
                    OboTerm {
                        id: id.clone(),
                        name: name.clone(),
                        namespace: namespace.clone(),
                        is_a: is_a.clone(),
                    },
                );
            }
            in_term = false;
            continue;
        }

        if !in_term {
            continue;
        }

        if let Some(val) = line.strip_prefix("id: ") {
            id = val.to_string();
        } else if let Some(val) = line.strip_prefix("name: ") {
            name = val.to_string();
        } else if let Some(val) = line.strip_prefix("namespace: ") {
            namespace = val.to_string();
        } else if let Some(val) = line.strip_prefix("is_a: ") {
            // "is_a: GO:0008150 ! biological_process" — take just the ID
            if let Some(parent_id) = val.split_whitespace().next() {
                is_a.push(parent_id.to_string());
            }
        } else if line == "is_obsolete: true" {
            is_obsolete = true;
        }
    }

    // Don't forget the last term
    if in_term && !id.is_empty() && !is_obsolete {
        for parent in &is_a {
            children
                .entry(parent.clone())
                .or_default()
                .insert(id.clone());
        }
        terms.insert(
            id.clone(),
            OboTerm {
                id,
                name,
                namespace,
                is_a,
                is_obsolete,
            },
        );
    }

    Ok(GeneOntology { terms, children })
}

// ---------------------------------------------------------------------------
// InterPro XML parser — streaming parser for interpro.xml
// ---------------------------------------------------------------------------

/// Output entry matching Python GECCO's interpro.json schema.
#[derive(Debug, serde::Serialize)]
struct InterProJsonEntry {
    accession: String,
    databases: Vec<String>,
    go_functions: Vec<GoRef>,
    go_terms: Vec<GoTermRef>,
    members: Vec<String>,
    name: String,
    #[serde(rename = "type")]
    entry_type: String,
}

#[derive(Debug, serde::Serialize)]
struct GoRef {
    accession: String,
    name: String,
}

#[derive(Debug, serde::Serialize)]
struct GoTermRef {
    accession: String,
    name: String,
    namespace: String,
}

/// Parse the InterPro XML and build JSON entries.
fn parse_interpro_xml<R: Read>(
    reader: R,
    go: &GeneOntology,
    top_functions: &HashSet<String>,
) -> Result<Vec<InterProJsonEntry>> {
    use quick_xml::events::Event;
    use quick_xml::Reader;

    let buf_reader = BufReader::new(reader);
    let mut xml_reader = Reader::from_reader(buf_reader);
    xml_reader.config_mut().trim_text(true);

    let mut entries = Vec::new();
    let mut buf = Vec::new();

    // State for current <interpro> element
    let mut in_interpro = false;
    let mut accession = String::new();
    let mut name = String::new();
    let mut entry_type = String::new();
    let mut members: HashSet<String> = HashSet::new();
    let mut databases: HashSet<String> = HashSet::new();
    let mut go_term_ids: HashSet<String> = HashSet::new();

    // Nested element tracking
    let mut in_name = false;
    let mut in_member_list = false;
    let mut in_class_list = false;
    let mut depth: u32 = 0;

    loop {
        match xml_reader.read_event_into(&mut buf) {
            Ok(Event::Eof) => break,
            Ok(Event::Start(ref e)) => {
                let local = e.local_name();
                match local.as_ref() {
                    b"interpro" => {
                        in_interpro = true;
                        depth = 0;
                        accession.clear();
                        name.clear();
                        entry_type.clear();
                        members.clear();
                        databases.clear();
                        go_term_ids.clear();
                        in_name = false;
                        in_member_list = false;
                        in_class_list = false;

                        for attr in e.attributes().flatten() {
                            match attr.key.as_ref() {
                                b"id" => {
                                    accession =
                                        String::from_utf8_lossy(&attr.value).to_string();
                                }
                                b"type" => {
                                    entry_type = String::from_utf8_lossy(&attr.value)
                                        .to_lowercase();
                                }
                                _ => {}
                            }
                        }
                    }
                    b"name" if in_interpro && depth == 0 => {
                        in_name = true;
                    }
                    b"member_list" if in_interpro => {
                        in_member_list = true;
                    }
                    b"db_xref" if in_interpro && in_member_list => {
                        let mut dbkey = String::new();
                        let mut db = String::new();
                        for attr in e.attributes().flatten() {
                            match attr.key.as_ref() {
                                b"dbkey" => {
                                    dbkey =
                                        String::from_utf8_lossy(&attr.value).to_string();
                                }
                                b"db" => {
                                    db = String::from_utf8_lossy(&attr.value).to_string();
                                }
                                _ => {}
                            }
                        }
                        if !dbkey.is_empty() {
                            members.insert(dbkey);
                        }
                        if !db.is_empty() {
                            databases.insert(db);
                        }
                    }
                    b"class_list" if in_interpro => {
                        in_class_list = true;
                    }
                    b"classification" if in_interpro && in_class_list => {
                        let mut class_type = String::new();
                        let mut class_id = String::new();
                        for attr in e.attributes().flatten() {
                            match attr.key.as_ref() {
                                b"class_type" => {
                                    class_type =
                                        String::from_utf8_lossy(&attr.value).to_string();
                                }
                                b"id" => {
                                    class_id =
                                        String::from_utf8_lossy(&attr.value).to_string();
                                }
                                _ => {}
                            }
                        }
                        if class_type == "GO" && !class_id.is_empty() {
                            go_term_ids.insert(class_id);
                        }
                    }
                    _ => {
                        if in_interpro {
                            depth += 1;
                        }
                    }
                }
            }
            Ok(Event::End(ref e)) => {
                let local = e.local_name();
                match local.as_ref() {
                    b"interpro" => {
                        if in_interpro {
                            // Build GO terms and functions
                            let mut go_terms: Vec<GoTermRef> = Vec::new();
                            let mut go_function_set: HashSet<String> = HashSet::new();

                            for go_id in &go_term_ids {
                                if let Some(term) = go.terms.get(go_id) {
                                    go_terms.push(GoTermRef {
                                        accession: term.id.clone(),
                                        name: term.name.clone(),
                                        namespace: term.namespace.clone(),
                                    });

                                    // Find top-level function ancestors
                                    let ancestors = go.superclasses(go_id);
                                    for ancestor_id in ancestors.intersection(top_functions) {
                                        go_function_set.insert(ancestor_id.clone());
                                    }
                                }
                            }

                            go_terms.sort_by(|a, b| a.accession.cmp(&b.accession));

                            let mut go_functions: Vec<GoRef> = go_function_set
                                .iter()
                                .filter_map(|fid| {
                                    go.terms.get(fid).map(|t| GoRef {
                                        accession: t.id.clone(),
                                        name: t.name.clone(),
                                    })
                                })
                                .collect();
                            go_functions.sort_by(|a, b| a.accession.cmp(&b.accession));

                            let mut members_sorted: Vec<String> =
                                members.iter().cloned().collect();
                            members_sorted.sort();

                            let mut databases_sorted: Vec<String> =
                                databases.iter().cloned().collect();
                            databases_sorted.sort();

                            entries.push(InterProJsonEntry {
                                accession: std::mem::take(&mut accession),
                                databases: databases_sorted,
                                go_functions,
                                go_terms,
                                members: members_sorted,
                                name: name.clone(),
                                entry_type: entry_type.clone(),
                            });

                            in_interpro = false;
                        }
                    }
                    b"name" => {
                        in_name = false;
                    }
                    b"member_list" => {
                        in_member_list = false;
                    }
                    b"class_list" => {
                        in_class_list = false;
                    }
                    _ => {
                        if in_interpro && depth > 0 {
                            depth -= 1;
                        }
                    }
                }
            }
            Ok(Event::Empty(ref e)) => {
                let local = e.local_name();
                // Handle self-closing tags like <db_xref .../> and <classification .../>
                match local.as_ref() {
                    b"db_xref" if in_interpro && in_member_list => {
                        let mut dbkey = String::new();
                        let mut db = String::new();
                        for attr in e.attributes().flatten() {
                            match attr.key.as_ref() {
                                b"dbkey" => {
                                    dbkey =
                                        String::from_utf8_lossy(&attr.value).to_string();
                                }
                                b"db" => {
                                    db = String::from_utf8_lossy(&attr.value).to_string();
                                }
                                _ => {}
                            }
                        }
                        if !dbkey.is_empty() {
                            members.insert(dbkey);
                        }
                        if !db.is_empty() {
                            databases.insert(db);
                        }
                    }
                    b"classification" if in_interpro && in_class_list => {
                        let mut class_type = String::new();
                        let mut class_id = String::new();
                        for attr in e.attributes().flatten() {
                            match attr.key.as_ref() {
                                b"class_type" => {
                                    class_type =
                                        String::from_utf8_lossy(&attr.value).to_string();
                                }
                                b"id" => {
                                    class_id =
                                        String::from_utf8_lossy(&attr.value).to_string();
                                }
                                _ => {}
                            }
                        }
                        if class_type == "GO" && !class_id.is_empty() {
                            go_term_ids.insert(class_id);
                        }
                    }
                    _ => {}
                }
            }
            Ok(Event::Text(ref e)) => {
                if in_name && in_interpro {
                    name = e.unescape().unwrap_or_default().to_string();
                    in_name = false;
                }
            }
            Err(e) => {
                return Err(anyhow::anyhow!(
                    "XML parse error at position {}: {}",
                    xml_reader.error_position(),
                    e
                ));
            }
            _ => {}
        }
        buf.clear();
    }

    // Sort by accession like Python does
    entries.sort_by(|a, b| a.accession.cmp(&b.accession));

    Ok(entries)
}

// ---------------------------------------------------------------------------
// HTTP helpers
// ---------------------------------------------------------------------------

fn ureq_get(url: &str) -> Result<Vec<u8>> {
    let agent = ureq::Agent::new_with_config(
        ureq::config::Config::builder()
            .timeout_global(Some(std::time::Duration::from_secs(300)))
            .build(),
    );
    let body = agent
        .get(url)
        .call()
        .with_context(|| format!("downloading {}", url))?
        .into_body()
        .with_config()
        .limit(200 * 1024 * 1024)
        .read_to_vec()
        .with_context(|| format!("reading response from {}", url))?;
    Ok(body)
}

/// Download a large file (up to 2 GB) — needed for interpro.xml.gz.
fn ureq_get_large(url: &str) -> Result<Vec<u8>> {
    let agent = ureq::Agent::new_with_config(
        ureq::config::Config::builder()
            .timeout_global(Some(std::time::Duration::from_secs(1800)))
            .build(),
    );
    let body = agent
        .get(url)
        .call()
        .with_context(|| format!("downloading {}", url))?
        .into_body()
        .with_config()
        .limit(2 * 1024 * 1024 * 1024) // 2 GB
        .read_to_vec()
        .with_context(|| format!("reading response from {}", url))?;
    Ok(body)
}
