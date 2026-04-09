//! GenBank I/O for reading input genomes and writing annotated cluster records.

use std::borrow::Cow;
use std::collections::BTreeMap;
use std::io::{Read, Write};

use anyhow::{Context, Result};
use gb_io::reader::SeqReader;
use gb_io::seq::{Date, Feature, Location, Reference, Seq, Source, Topology};

use crate::model::{Cluster, Domain, Gene, Strand};
use crate::orf::SeqRecord;

// ---------------------------------------------------------------------------
// Reading
// ---------------------------------------------------------------------------

/// Read FASTA sequences from a file, returning SeqRecords.
pub fn read_fasta(reader: impl Read) -> Result<Vec<SeqRecord>> {
    let mut records = Vec::new();
    let mut buf_reader = std::io::BufReader::new(reader);
    let mut content = String::new();
    buf_reader.read_to_string(&mut content)?;

    let mut current_id = String::new();
    let mut current_seq = String::new();

    for line in content.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if !current_id.is_empty() {
                records.push(SeqRecord {
                    id: current_id.clone(),
                    seq: current_seq.clone(),
                });
                current_seq.clear();
            }
            current_id = header.split_whitespace().next().unwrap_or("").to_string();
        } else {
            current_seq.push_str(line.trim());
        }
    }
    if !current_id.is_empty() {
        records.push(SeqRecord {
            id: current_id,
            seq: current_seq,
        });
    }

    Ok(records)
}

/// Read sequences from a GenBank file, returning SeqRecords with DNA sequences.
pub fn read_genbank(reader: impl Read) -> Result<Vec<SeqRecord>> {
    let mut records = Vec::new();
    for result in SeqReader::new(reader) {
        let seq = result.map_err(|e| anyhow::anyhow!("parsing GenBank: {}", e))?;
        let id = seq.name.unwrap_or_default();
        let dna = String::from_utf8_lossy(&seq.seq).to_uppercase();
        records.push(SeqRecord { id, seq: dna });
    }
    Ok(records)
}

/// Auto-detect format and read sequences.
pub fn read_sequences(path: &std::path::Path) -> Result<Vec<SeqRecord>> {
    let reader = crate::io::compression::zopen_path(path)?;
    let mut buf = Vec::new();
    let mut buf_reader = std::io::BufReader::new(reader);
    buf_reader.read_to_end(&mut buf)?;

    // Peek at first non-whitespace character to detect format
    let first_char = buf.iter().find(|&&b| !b.is_ascii_whitespace());
    match first_char {
        Some(b'>') => read_fasta(&buf[..]),
        Some(b'L') => read_genbank(&buf[..]), // LOCUS line
        _ => {
            // Try GenBank first, fall back to FASTA
            read_genbank(&buf[..]).or_else(|_| read_fasta(&buf[..]))
        }
    }
}

// ---------------------------------------------------------------------------
// Gene function colors
// ---------------------------------------------------------------------------

struct GeneColor {
    rgb: (u8, u8, u8),
}

impl GeneColor {
    fn rgb_str(&self) -> String {
        format!("{} {} {}", self.rgb.0, self.rgb.1, self.rgb.2)
    }

    fn hex_str(&self) -> String {
        format!("#{:02x}{:02x}{:02x}", self.rgb.0, self.rgb.1, self.rgb.2)
    }
}

fn gene_color(gene: &Gene) -> GeneColor {
    let funcs = gene.functions();
    let color = if funcs.contains("transporter activity") {
        (100, 149, 237)
    } else if funcs.iter().any(|f| {
        f.contains("regulation") || f.contains("regulatory") || f.contains("signaling")
    }) {
        (46, 139, 86)
    } else if funcs.contains("catalytic activity") || funcs.contains("binding") {
        (129, 14, 21)
    } else if funcs.contains("biosynthetic process") {
        (241, 109, 117)
    } else if funcs.contains("unknown") {
        (128, 128, 128)
    } else {
        (189, 183, 107)
    };
    GeneColor { rgb: color }
}

// ---------------------------------------------------------------------------
// Writing — Cluster to GenBank
// ---------------------------------------------------------------------------

/// Convert a Cluster to a gb-io Seq record for GenBank output.
///
/// The record contains:
/// - LOCUS, DEFINITION, ACCESSION, VERSION
/// - GECCO reference
/// - Structured comment with GECCO-Data
/// - CDS features for genes with qualifiers
/// - misc_feature entries for protein domains
/// - Full DNA sequence (if `source_seq` is provided)
pub fn cluster_to_seq(
    cluster: &Cluster,
    source_seq: Option<&str>,
) -> Seq {
    let version = env!("CARGO_PKG_VERSION");
    let now = current_date();

    // Extract cluster DNA if source sequence provided
    let cluster_dna = source_seq
        .map(|seq| {
            let start = (cluster.start() as usize).saturating_sub(1);
            let end = (cluster.end() as usize).min(seq.len());
            seq[start..end].as_bytes().to_vec()
        })
        .unwrap_or_default();

    let cluster_offset = cluster.start() - 1; // 0-based offset

    // Build features
    let mut features = Vec::new();

    for gene in &cluster.genes {
        // CDS feature
        let loc = gene_location(gene, cluster_offset);
        let mut qualifiers = Vec::new();

        // Standard qualifiers
        if let Some(inf) = gene.qualifiers.get("inference") {
            for v in inf {
                qualifiers.push((Cow::from("inference"), Some(v.clone())));
            }
        }
        if let Some(tt) = gene.qualifiers.get("transl_table") {
            for v in tt {
                qualifiers.push((
                    Cow::from("transl_table"),
                    Some(v.clone()),
                ));
            }
        }
        qualifiers.push((
            Cow::from("locus_tag"),
            Some(gene.protein.id.clone()),
        ));
        if !gene.protein.seq.is_empty() {
            qualifiers.push((
                Cow::from("translation"),
                Some(format!("{}*", gene.protein.seq)),
            ));
        }

        // Function
        let mut funcs: Vec<String> = gene.functions().into_iter().collect();
        funcs.sort();
        for f in &funcs {
            qualifiers.push((Cow::from("function"), Some(f.clone())));
        }

        // Colors
        let color = gene_color(gene);
        qualifiers.push((
            Cow::from("colour"),
            Some(color.rgb_str()),
        ));
        qualifiers.push((
            Cow::from("ApEinfo_fwdcolor"),
            Some(color.hex_str()),
        ));
        qualifiers.push((
            Cow::from("ApEinfo_revcolor"),
            Some(color.hex_str()),
        ));

        features.push(Feature {
            kind: Cow::from("CDS"),
            location: loc,
            qualifiers,
        });

        // misc_feature entries for domains
        for domain in &gene.protein.domains {
            let dom_loc = domain_location(gene, domain, cluster_offset);
            let mut dom_quals = Vec::new();

            // inference
            dom_quals.push((
                Cow::from("inference"),
                Some("protein motif".to_string()),
            ));

            // db_xref
            if let Some(xrefs) = domain.qualifiers.get("db_xref") {
                for xref in xrefs {
                    dom_quals.push((
                        Cow::from("db_xref"),
                        Some(xref.clone()),
                    ));
                }
            }
            for go_term in &domain.go_terms {
                dom_quals.push((
                    Cow::from("db_xref"),
                    Some(go_term.accession.clone()),
                ));
            }

            // note
            dom_quals.push((
                Cow::from("note"),
                Some(format!("e-value: {:e}", domain.i_evalue)),
            ));
            dom_quals.push((
                Cow::from("note"),
                Some(format!("p-value: {:e}", domain.pvalue)),
            ));

            // function
            if let Some(funcs) = domain.qualifiers.get("function") {
                for f in funcs {
                    dom_quals.push((
                        Cow::from("function"),
                        Some(f.clone()),
                    ));
                }
            }

            // standard_name
            dom_quals.push((
                Cow::from("standard_name"),
                Some(domain.name.clone()),
            ));

            features.push(Feature {
                kind: Cow::from("misc_feature"),
                location: dom_loc,
                qualifiers: dom_quals,
            });
        }
    }

    // Build structured comment
    let comment = build_gecco_comment(cluster, version);

    // Build reference
    let reference = Reference {
        description: "1".to_string(),
        authors: Some(
            "Laura M Carroll, Martin Larralde, Jonas Simon Fleck, Ruby Ponnudurai, \
             Alessio Milanese, Elisa Cappio Barazzone, Georg Zeller"
                .to_string(),
        ),
        consortium: None,
        title: "Accurate de novo identification of biosynthetic gene clusters with GECCO"
            .to_string(),
        journal: Some("bioRxiv (2021.05.03.442509)".to_string()),
        pubmed: None,
        remark: Some("doi:10.1101/2021.05.03.442509".to_string()),
    };

    Seq {
        name: Some(cluster.id.clone()),
        topology: Topology::Linear,
        date: Some(now),
        len: Some(cluster_dna.len()),
        molecule_type: Some("DNA".to_string()),
        division: "UNK".to_string(),
        definition: None,
        accession: Some(cluster.id.clone()),
        version: Some(cluster.id.clone()),
        source: Some(Source {
            source: ".".to_string(),
            organism: None,
        }),
        dblink: None,
        keywords: Some(".".to_string()),
        references: vec![reference],
        comments: vec![comment],
        seq: cluster_dna,
        contig: None,
        features,
    }
}

/// Write a cluster as a GenBank file.
pub fn write_cluster_gbk(
    writer: impl Write,
    cluster: &Cluster,
    source_seq: Option<&str>,
) -> Result<()> {
    let seq = cluster_to_seq(cluster, source_seq);
    seq.write(writer)
        .with_context(|| format!("writing GenBank for cluster {}", cluster.id))
}

/// Write multiple clusters to a single merged GenBank file.
pub fn write_clusters_merged(
    writer: impl Write,
    clusters: &[Cluster],
    source_seqs: &BTreeMap<String, String>,
) -> Result<()> {
    let mut w = std::io::BufWriter::new(writer);
    for cluster in clusters {
        let source_seq = source_seqs.get(cluster.source_id()).map(|s| s.as_str());
        let seq = cluster_to_seq(cluster, source_seq);
        seq.write(&mut w)
            .with_context(|| format!("writing GenBank for cluster {}", cluster.id))?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Build the GECCO-Data structured comment string.
fn build_gecco_comment(cluster: &Cluster, version: &str) -> String {
    let mut lines = Vec::new();
    lines.push("##GECCO-Data-START##".to_string());
    lines.push(format!("{:<27}:: GECCO v{}", "version", version));
    lines.push(format!(
        "{:<27}:: {}",
        "creation_date",
        chrono_like_now()
    ));

    if let Some(ref ct) = cluster.cluster_type {
        if !ct.is_unknown() {
            lines.push(format!("{:<27}:: {}", "cluster_type", ct));
        }
    }

    // Type probabilities
    let type_order = [
        "Alkaloid",
        "NRP",
        "Polyketide",
        "RiPP",
        "Saccharide",
        "Terpene",
    ];
    for type_name in &type_order {
        let key = format!("{}_probability", type_name.to_lowercase());
        let prob = cluster
            .type_probabilities
            .get(*type_name)
            .copied()
            .unwrap_or(0.0);
        lines.push(format!("{:<27}:: {:.3}", key, prob));
    }

    lines.push("##GECCO-Data-END##".to_string());
    lines.join("\n")
}

/// Convert gene genomic coordinates to a gb-io Location, offset by cluster start.
fn gene_location(gene: &Gene, offset: i64) -> Location {
    // gb-io uses 0-based exclusive-end coordinates
    let start = gene.start - 1 - offset;
    let end = gene.end - offset;
    let loc = Location::simple_range(start, end);
    match gene.strand {
        Strand::Reverse => Location::Complement(Box::new(loc)),
        Strand::Coding => loc,
    }
}

/// Convert domain coordinates (protein-relative, 1-based) to nucleotide
/// location in cluster coordinates.
fn domain_location(gene: &Gene, domain: &Domain, cluster_offset: i64) -> Location {
    let gene_start = gene.start - 1 - cluster_offset; // 0-based in cluster
    match gene.strand {
        Strand::Coding => {
            let nuc_start = gene_start + (domain.start - 1) * 3;
            let nuc_end = gene_start + domain.end * 3;
            Location::simple_range(nuc_start, nuc_end)
        }
        Strand::Reverse => {
            let gene_end = gene.end - cluster_offset; // exclusive, 0-based in cluster
            let nuc_end = gene_end - (domain.start - 1) * 3;
            let nuc_start = gene_end - domain.end * 3;
            Location::Complement(Box::new(Location::simple_range(nuc_start, nuc_end)))
        }
    }
}

/// Get current date as gb-io Date.
fn current_date() -> Date {
    // Use a simple approach without chrono dependency
    let since_epoch = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();

    // Approximate date calculation
    let days = since_epoch / 86400;
    let (year, month, day) = days_to_ymd(days as i64);
    Date::from_ymd(year as i32, month as u32, day as u32).unwrap_or_else(|_| {
        Date::from_ymd(2025, 1, 1).unwrap()
    })
}

/// Simple ISO-like timestamp for creation_date.
fn chrono_like_now() -> String {
    let since_epoch = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default();
    let secs = since_epoch.as_secs();
    let (year, month, day) = days_to_ymd((secs / 86400) as i64);
    let time_of_day = secs % 86400;
    let hour = time_of_day / 3600;
    let minute = (time_of_day % 3600) / 60;
    let second = time_of_day % 60;
    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}",
        year, month, day, hour, minute, second
    )
}

/// Convert days since Unix epoch to (year, month, day).
fn days_to_ymd(days: i64) -> (i64, i64, i64) {
    // Algorithm from http://howardhinnant.github.io/date_algorithms.html
    let z = days + 719468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    (y, m, d)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{ClusterType, Protein};

    fn make_test_cluster() -> Cluster {
        let domain = Domain {
            name: "PF00001".to_string(),
            start: 10,
            end: 50,
            hmm: "Pfam".to_string(),
            i_evalue: 1e-10,
            pvalue: 1e-14,
            probability: Some(0.95),
            cluster_weight: None,
            go_terms: vec![],
            go_functions: vec![],
            qualifiers: {
                let mut q = BTreeMap::new();
                q.insert("db_xref".to_string(), vec!["PFAM:PF00001".to_string()]);
                q.insert(
                    "function".to_string(),
                    vec!["Test function".to_string()],
                );
                q
            },
        };

        let protein = Protein {
            id: "gene_1".to_string(),
            seq: "MKVL".to_string(),
            domains: vec![domain],
        };

        let gene = Gene {
            source_id: "seq1".to_string(),
            start: 100,
            end: 400,
            strand: Strand::Coding,
            protein,
            qualifiers: {
                let mut q = BTreeMap::new();
                q.insert("transl_table".to_string(), vec!["11".to_string()]);
                q
            },
            probability: Some(0.9),
        };

        let mut type_probs = BTreeMap::new();
        type_probs.insert("Polyketide".to_string(), 0.96);
        type_probs.insert("NRP".to_string(), 0.14);
        type_probs.insert("Alkaloid".to_string(), 0.01);
        type_probs.insert("RiPP".to_string(), 0.0);
        type_probs.insert("Saccharide".to_string(), 0.0);
        type_probs.insert("Terpene".to_string(), 0.01);

        Cluster {
            id: "seq1_cluster_1".to_string(),
            genes: vec![gene],
            cluster_type: Some(ClusterType::new(["Polyketide"])),
            type_probabilities: type_probs,
        }
    }

    #[test]
    fn test_cluster_to_seq() {
        let cluster = make_test_cluster();
        let seq = cluster_to_seq(&cluster, None);

        assert_eq!(seq.name.as_deref(), Some("seq1_cluster_1"));
        assert_eq!(seq.topology, Topology::Linear);
        assert!(seq.molecule_type.as_deref() == Some("DNA"));

        // Should have 1 CDS + 1 misc_feature
        assert_eq!(seq.features.len(), 2);
        assert_eq!(
            seq.features[0].kind,
            Cow::from("CDS")
        );
        assert_eq!(
            seq.features[1].kind,
            Cow::from("misc_feature")
        );
    }

    #[test]
    fn test_gecco_comment() {
        let cluster = make_test_cluster();
        let comment = build_gecco_comment(&cluster, "0.1.0");

        assert!(comment.contains("##GECCO-Data-START##"));
        assert!(comment.contains("##GECCO-Data-END##"));
        assert!(comment.contains("GECCO v0.1.0"));
        assert!(comment.contains("cluster_type"));
        assert!(comment.contains("Polyketide"));
        assert!(comment.contains("polyketide_probability"));
        assert!(comment.contains("0.960"));
    }

    #[test]
    fn test_write_cluster_gbk() {
        let cluster = make_test_cluster();
        let source = "A".repeat(500);
        let mut buf = Vec::new();
        write_cluster_gbk(&mut buf, &cluster, Some(&source)).unwrap();
        let output = String::from_utf8(buf).unwrap();

        assert!(output.contains("LOCUS"));
        assert!(output.contains("seq1_cluster_1"));
        assert!(output.contains("CDS"));
        assert!(output.contains("misc_feature"));
        assert!(output.contains("PF00001"));
        assert!(output.contains("GECCO-Data"));
    }

    #[test]
    fn test_read_fasta() {
        let fasta = b">seq1 description\nATGCATGC\nAAAA\n>seq2\nTTTT\n";
        let records = read_fasta(&fasta[..]).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].seq, "ATGCATGCAAAA");
        assert_eq!(records[1].id, "seq2");
        assert_eq!(records[1].seq, "TTTT");
    }

    #[test]
    fn test_days_to_ymd() {
        // 2025-01-01 = day 20089 since epoch
        let (y, m, d) = days_to_ymd(20089);
        assert_eq!((y, m, d), (2025, 1, 1));
    }
}
