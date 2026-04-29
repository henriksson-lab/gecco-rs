//! Integration tests matching the Python test suite.

use std::collections::BTreeMap;
use std::io::Read;

use gecco::model::*;

// ---------------------------------------------------------------------------
// Test helpers
// ---------------------------------------------------------------------------

fn make_domain(name: &str, start: i64, end: i64, prob: Option<f64>) -> Domain {
    Domain {
        name: name.to_string(),
        start,
        end,
        hmm: "Pfam".to_string(),
        i_evalue: 1e-10,
        pvalue: 1e-14,
        probability: prob,
        cluster_weight: None,
        go_terms: vec![],
        go_functions: vec![],
        qualifiers: BTreeMap::new(),
    }
}

fn make_gene(
    source: &str,
    id: &str,
    start: i64,
    end: i64,
    strand: Strand,
    domains: Vec<Domain>,
    probability: Option<f64>,
) -> Gene {
    Gene {
        source_id: source.to_string(),
        start,
        end,
        strand,
        protein: Protein {
            id: id.to_string(),
            seq: String::new(),
            domains,
        },
        qualifiers: BTreeMap::new(),
        probability,
    }
}

// ---------------------------------------------------------------------------
// Model tests (matching test_model/)
// ---------------------------------------------------------------------------

#[test]
fn test_gene_average_probability_from_domains() {
    let gene = make_gene(
        "seq1",
        "p1",
        1,
        300,
        Strand::Coding,
        vec![
            make_domain("PF00106", 1, 20, Some(0.8)),
            make_domain("PF00107", 25, 50, Some(0.6)),
        ],
        None,
    );
    let avg = gene.average_probability().unwrap();
    assert!((avg - 0.7).abs() < 1e-10);
}

#[test]
fn test_gene_maximum_probability_from_domains() {
    let gene = make_gene(
        "seq1",
        "p1",
        1,
        300,
        Strand::Coding,
        vec![
            make_domain("PF00106", 1, 20, Some(0.8)),
            make_domain("PF00107", 25, 50, Some(0.6)),
        ],
        None,
    );
    let max = gene.maximum_probability().unwrap();
    assert!((max - 0.8).abs() < 1e-10);
}

#[test]
fn test_gene_probability_override() {
    let gene = make_gene(
        "seq1",
        "p1",
        1,
        300,
        Strand::Coding,
        vec![make_domain("PF00106", 1, 20, Some(0.8))],
        Some(0.5),
    );
    // Override takes precedence
    assert!((gene.average_probability().unwrap() - 0.5).abs() < 1e-10);
    assert!((gene.maximum_probability().unwrap() - 0.5).abs() < 1e-10);
}

#[test]
fn test_gene_no_probability() {
    let gene = make_gene("seq1", "p1", 1, 300, Strand::Coding, vec![], None);
    assert!(gene.average_probability().is_none());
    assert!(gene.maximum_probability().is_none());
}

#[test]
fn test_cluster_type_display() {
    let ct = ClusterType::new(["Polyketide", "NRP"]);
    let s = ct.to_string();
    assert!(s.contains("NRP"));
    assert!(s.contains("Polyketide"));
}

#[test]
fn test_cluster_type_unknown() {
    let ct = ClusterType::unknown();
    assert!(ct.is_unknown());
    assert_eq!(ct.to_string(), "Unknown");
}

#[test]
fn test_cluster_type_unpack() {
    let ct = ClusterType::new(["Polyketide", "NRP"]);
    let unpacked = ct.unpack();
    assert_eq!(unpacked.len(), 2);
}

#[test]
fn test_strand_sign() {
    assert_eq!(Strand::Coding.sign(), "+");
    assert_eq!(Strand::Reverse.sign(), "-");
}

#[test]
fn test_strand_from_sign() {
    assert_eq!(Strand::from_sign("+"), Some(Strand::Coding));
    assert_eq!(Strand::from_sign("-"), Some(Strand::Reverse));
    assert_eq!(Strand::from_sign("?"), None);
}

#[test]
fn test_cluster_coordinates() {
    let cluster = Cluster::new(
        "c1",
        vec![
            make_gene("seq1", "p1", 100, 200, Strand::Coding, vec![], Some(0.9)),
            make_gene("seq1", "p2", 250, 400, Strand::Coding, vec![], Some(0.8)),
        ],
    );
    assert_eq!(cluster.start(), 100);
    assert_eq!(cluster.end(), 400);
    assert_eq!(cluster.source_id(), "seq1");
}

#[test]
fn test_cluster_probabilities() {
    let cluster = Cluster::new(
        "c1",
        vec![
            make_gene("seq1", "p1", 1, 100, Strand::Coding, vec![], Some(0.9)),
            make_gene("seq1", "p2", 101, 200, Strand::Coding, vec![], Some(0.7)),
        ],
    );
    let avg = cluster.average_probability().unwrap();
    assert!((avg - 0.8).abs() < 1e-10);
    let max = cluster.maximum_probability().unwrap();
    assert!((max - 0.9).abs() < 1e-10);
}

#[test]
fn test_domain_composition() {
    let cluster = Cluster::new(
        "c1",
        vec![make_gene(
            "seq1",
            "p1",
            1,
            300,
            Strand::Coding,
            vec![
                make_domain("PF00001", 1, 20, None),
                make_domain("PF00002", 25, 50, None),
            ],
            None,
        )],
    );

    let all_domains = vec![
        "PF00001".to_string(),
        "PF00002".to_string(),
        "PF00003".to_string(),
    ];
    let comp = cluster.domain_composition(Some(&all_domains), false, false, true);
    assert_eq!(comp.len(), 3);
    // PF00001 and PF00002 present, PF00003 absent
    assert!(comp[0] > 0.0); // PF00001
    assert!(comp[1] > 0.0); // PF00002
    assert_eq!(comp[2], 0.0); // PF00003
}

// ---------------------------------------------------------------------------
// Table I/O tests (matching test_model/test_genetable.py)
// ---------------------------------------------------------------------------

#[test]
fn test_gene_table_roundtrip() {
    use gecco::io::tables::GeneTable;

    let genes = vec![
        make_gene(
            "seq1",
            "seq1_1",
            1,
            300,
            Strand::Coding,
            vec![make_domain("PF00106", 1, 20, Some(0.8))],
            Some(0.9),
        ),
        make_gene("seq1", "seq1_2", 400, 600, Strand::Reverse, vec![], None),
    ];

    // Write
    let mut buf = Vec::new();
    GeneTable::write_from_genes(&mut buf, &genes).unwrap();
    let tsv = String::from_utf8(buf.clone()).unwrap();

    // Should contain header and 2 data rows
    let lines: Vec<&str> = tsv.lines().collect();
    assert_eq!(lines.len(), 3); // header + 2 genes
    assert!(lines[0].contains("sequence_id"));
    assert!(lines[1].contains("seq1_1"));
    assert!(lines[2].contains("seq1_2"));

    // Read back
    let loaded = GeneTable::read_to_genes(&buf[..]).unwrap();
    assert_eq!(loaded.len(), 2);
    assert_eq!(loaded[0].source_id, "seq1");
    assert_eq!(loaded[0].protein.id, "seq1_1");
    assert_eq!(loaded[0].start, 1);
    assert_eq!(loaded[0].end, 300);
    assert_eq!(loaded[0].strand, Strand::Coding);
    assert_eq!(loaded[1].strand, Strand::Reverse);
}

#[test]
fn test_feature_table_roundtrip() {
    use gecco::io::tables::FeatureTable;

    let genes = vec![make_gene(
        "seq1",
        "seq1_1",
        1,
        300,
        Strand::Coding,
        vec![
            make_domain("PF00106", 1, 20, Some(0.8)),
            make_domain("PF00107", 25, 50, Some(0.6)),
        ],
        None,
    )];

    // Write
    let mut buf = Vec::new();
    FeatureTable::write_from_genes(&mut buf, &genes).unwrap();

    // Read back
    let loaded = FeatureTable::read_to_genes(&buf[..]).unwrap();
    assert_eq!(loaded.len(), 1);
    assert_eq!(loaded[0].protein.domains.len(), 2);
    assert_eq!(loaded[0].protein.domains[0].name, "PF00106");
    assert_eq!(loaded[0].protein.domains[1].name, "PF00107");
}

#[test]
fn test_cluster_table_write() {
    use gecco::io::tables::ClusterTable;

    let cluster = Cluster {
        id: "c1".to_string(),
        genes: vec![make_gene(
            "seq1",
            "p1",
            1,
            300,
            Strand::Coding,
            vec![make_domain("PF00106", 1, 20, None)],
            Some(0.9),
        )],
        cluster_type: Some(ClusterType::new(["Polyketide"])),
        type_probabilities: BTreeMap::new(),
    };

    let mut buf = Vec::new();
    ClusterTable::write_from_clusters(&mut buf, &[cluster]).unwrap();
    let tsv = String::from_utf8(buf).unwrap();

    assert!(tsv.contains("Polyketide"));
    assert!(tsv.contains("seq1"));
    assert!(tsv.contains("PF00106"));
}

// ---------------------------------------------------------------------------
// CRF feature extraction tests (matching test_crf/test_features.py)
// ---------------------------------------------------------------------------

#[test]
fn test_extract_features_domain() {
    use gecco::crf::features::extract_features_domain;

    let genes = vec![
        make_gene(
            "s1",
            "p1",
            0,
            1,
            Strand::Coding,
            vec![
                make_domain("A", 0, 1, Some(1.0)),
                make_domain("B", 0, 1, Some(1.0)),
            ],
            None,
        ),
        make_gene(
            "s1",
            "p2",
            1,
            2,
            Strand::Coding,
            vec![make_domain("C", 0, 1, Some(0.0))],
            None,
        ),
    ];

    let feats = extract_features_domain(&genes);
    // Domain mode: one ordered feature list per domain
    assert_eq!(feats.len(), 3);
    assert_eq!(feats[0], vec!["A"]);
    assert_eq!(feats[1], vec!["B"]);
    assert_eq!(feats[2], vec!["C"]);
}

#[test]
fn test_extract_features_protein() {
    use gecco::crf::features::extract_features_protein;

    let genes = vec![
        make_gene(
            "s1",
            "p1",
            0,
            1,
            Strand::Coding,
            vec![
                make_domain("A", 0, 1, Some(1.0)),
                make_domain("B", 0, 1, Some(1.0)),
            ],
            None,
        ),
        make_gene(
            "s1",
            "p2",
            1,
            2,
            Strand::Coding,
            vec![make_domain("C", 0, 1, Some(0.0))],
            None,
        ),
    ];

    let feats = extract_features_protein(&genes);
    // Protein mode: one ordered feature list per gene
    assert_eq!(feats.len(), 2);
    assert_eq!(feats[0], vec!["A", "B"]);
    assert_eq!(feats[1], vec!["C"]);
}

#[test]
fn test_extract_labels_protein() {
    use gecco::crf::features::extract_labels_protein;

    let genes = vec![
        make_gene("s1", "p1", 0, 1, Strand::Coding, vec![], Some(0.9)),
        make_gene("s1", "p2", 1, 2, Strand::Coding, vec![], Some(0.3)),
    ];

    let labels = extract_labels_protein(&genes);
    assert_eq!(labels, vec!["1", "0"]);
}

// ---------------------------------------------------------------------------
// Cluster refiner tests
// ---------------------------------------------------------------------------

#[test]
fn test_refiner_basic() {
    use gecco::refine::ClusterRefiner;

    let refiner = ClusterRefiner {
        threshold: 0.5,
        n_cds: 2,
        edge_distance: 0,
        trim: false,
        ..Default::default()
    };

    // 5 genes above threshold, 2 below
    let genes: Vec<Gene> = (0..7)
        .map(|i| {
            let prob = if i < 5 { 0.9 } else { 0.1 };
            make_gene(
                "seq1",
                &format!("p{}", i),
                i * 100 + 1,
                (i + 1) * 100,
                Strand::Coding,
                vec![make_domain("PF00001", 1, 20, Some(prob))],
                Some(prob),
            )
        })
        .collect();

    let clusters = refiner.iter_clusters(&genes);
    assert_eq!(clusters.len(), 1);
    assert_eq!(clusters[0].genes.len(), 5);
}

#[test]
fn test_refiner_trim() {
    use gecco::refine::ClusterRefiner;

    let refiner = ClusterRefiner {
        threshold: 0.5,
        n_cds: 1,
        edge_distance: 0,
        trim: true,
        ..Default::default()
    };

    // Gene with no domains at start → should be trimmed
    let genes = vec![
        make_gene("seq1", "p0", 1, 100, Strand::Coding, vec![], Some(0.9)),
        make_gene(
            "seq1",
            "p1",
            101,
            200,
            Strand::Coding,
            vec![make_domain("PF00001", 1, 20, Some(0.9))],
            Some(0.9),
        ),
        make_gene("seq1", "p2", 201, 300, Strand::Coding, vec![], Some(0.9)),
    ];

    let clusters = refiner.iter_clusters(&genes);
    assert_eq!(clusters.len(), 1);
    // Should trim p0 (no domains at start) and p2 (no domains at end)
    assert_eq!(clusters[0].genes.len(), 1);
    assert_eq!(clusters[0].genes[0].protein.id, "p1");
}

#[test]
fn test_refiner_no_clusters() {
    use gecco::refine::ClusterRefiner;

    let refiner = ClusterRefiner::default();

    // All below threshold
    let genes: Vec<Gene> = (0..3)
        .map(|i| {
            make_gene(
                "seq1",
                &format!("p{}", i),
                i * 100 + 1,
                (i + 1) * 100,
                Strand::Coding,
                vec![],
                Some(0.1),
            )
        })
        .collect();

    let clusters = refiner.iter_clusters(&genes);
    assert!(clusters.is_empty());
}

#[test]
fn test_refiner_multiple_contigs() {
    use gecco::refine::ClusterRefiner;

    let refiner = ClusterRefiner {
        threshold: 0.5,
        n_cds: 1,
        edge_distance: 0,
        trim: false,
        ..Default::default()
    };

    let mut genes = Vec::new();
    // 3 genes on seq1
    for i in 0..3 {
        genes.push(make_gene(
            "seq1",
            &format!("s1_p{}", i),
            i * 100 + 1,
            (i + 1) * 100,
            Strand::Coding,
            vec![make_domain("PF00001", 1, 20, Some(0.9))],
            Some(0.9),
        ));
    }
    // 2 genes on seq2
    for i in 0..2 {
        genes.push(make_gene(
            "seq2",
            &format!("s2_p{}", i),
            i * 100 + 1,
            (i + 1) * 100,
            Strand::Coding,
            vec![make_domain("PF00001", 1, 20, Some(0.9))],
            Some(0.9),
        ));
    }

    let clusters = refiner.iter_clusters(&genes);
    assert_eq!(clusters.len(), 2);
    assert_eq!(clusters[0].source_id(), "seq1");
    assert_eq!(clusters[1].source_id(), "seq2");
}

// ---------------------------------------------------------------------------
// Type binarizer tests
// ---------------------------------------------------------------------------

#[test]
fn test_type_binarizer_transform() {
    use gecco::types::TypeBinarizer;

    let binarizer = TypeBinarizer::new(vec![
        "Alkaloid".to_string(),
        "NRP".to_string(),
        "Polyketide".to_string(),
    ]);

    let types = vec![ClusterType::new(["Polyketide", "NRP"])];
    let matrix = binarizer.transform(&types);
    assert_eq!(matrix.len(), 1);
    assert_eq!(matrix[0], vec![0.0, 1.0, 1.0]);
}

#[test]
fn test_type_binarizer_inverse() {
    use gecco::types::TypeBinarizer;

    let binarizer = TypeBinarizer::new(vec![
        "Alkaloid".to_string(),
        "NRP".to_string(),
        "Polyketide".to_string(),
    ]);

    let flags = vec![vec![false, true, true]];
    let types = binarizer.inverse_transform(&flags);
    assert_eq!(types.len(), 1);
    assert!(types[0].names.contains("NRP"));
    assert!(types[0].names.contains("Polyketide"));
    assert!(!types[0].names.contains("Alkaloid"));
}

#[test]
fn test_type_binarizer_roundtrip() {
    use gecco::types::TypeBinarizer;

    let binarizer = TypeBinarizer::new(vec![
        "Alkaloid".to_string(),
        "NRP".to_string(),
        "Polyketide".to_string(),
    ]);

    let original = vec![ClusterType::new(["Polyketide"])];
    let matrix = binarizer.transform(&original);
    let flags: Vec<Vec<bool>> = matrix
        .iter()
        .map(|row| row.iter().map(|&v| v > 0.5).collect())
        .collect();
    let recovered = binarizer.inverse_transform(&flags);
    assert_eq!(recovered[0], original[0]);
}

// ---------------------------------------------------------------------------
// FASTA reading test
// ---------------------------------------------------------------------------

#[test]
fn test_read_fasta() {
    use gecco::io::genbank::read_fasta;

    let data = b">seq1 some description\nATGCATGC\nAAAA\n>seq2\nTTTTGGGG\n";
    let records = read_fasta(&data[..]).unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].id, "seq1");
    assert_eq!(records[0].seq, "ATGCATGCAAAA");
    assert_eq!(records[1].id, "seq2");
    assert_eq!(records[1].seq, "TTTTGGGG");
}

// ---------------------------------------------------------------------------
// GenBank output test
// ---------------------------------------------------------------------------

#[test]
fn test_genbank_output_structure() {
    use gecco::io::genbank::write_cluster_gbk;

    let cluster = Cluster {
        id: "seq1_cluster_1".to_string(),
        genes: vec![make_gene(
            "seq1",
            "seq1_1",
            100,
            400,
            Strand::Coding,
            vec![make_domain("PF00394", 10, 50, Some(0.95))],
            Some(0.9),
        )],
        cluster_type: Some(ClusterType::new(["Polyketide"])),
        type_probabilities: {
            let mut m = BTreeMap::new();
            m.insert("Polyketide".to_string(), 0.96);
            m.insert("NRP".to_string(), 0.14);
            m
        },
    };

    let source_dna = "A".repeat(500);
    let mut buf = Vec::new();
    write_cluster_gbk(&mut buf, &cluster, Some(&source_dna)).unwrap();
    let output = String::from_utf8(buf).unwrap();

    // Verify key sections present
    assert!(output.contains("LOCUS"));
    assert!(output.contains("seq1_cluster_1"));
    assert!(output.contains("GECCO-Data-START"));
    assert!(output.contains("GECCO-Data-END"));
    assert!(output.contains("Polyketide"));
    assert!(output.contains("polyketide_probability"));
    assert!(output.contains("CDS"));
    assert!(output.contains("misc_feature"));
    assert!(output.contains("PF00394"));
    assert!(output.contains(concat!("GECCO v", env!("CARGO_PKG_VERSION"))));
    assert!(output.contains("Laura M Carroll"));
}

// ---------------------------------------------------------------------------
// Sliding window test
// ---------------------------------------------------------------------------

#[test]
fn test_sliding_window() {
    use gecco::util::sliding_window;

    let windows = sliding_window(10, 3, 1);
    assert_eq!(windows.len(), 8);
    assert_eq!(windows[0], 0..3);
    assert_eq!(windows[7], 7..10);
}

#[test]
fn test_sliding_window_step2() {
    use gecco::util::sliding_window;

    let windows = sliding_window(10, 3, 2);
    assert_eq!(windows.len(), 4);
    assert_eq!(windows[0], 0..3);
    assert_eq!(windows[1], 2..5);
}

// ---------------------------------------------------------------------------
// Compression test
// ---------------------------------------------------------------------------

#[test]
fn test_zopen_plain() {
    use gecco::io::compression::zopen;

    let data = b"hello world";
    let mut reader = zopen(&data[..]).unwrap();
    let mut buf = String::new();
    reader.read_to_string(&mut buf).unwrap();
    assert_eq!(buf, "hello world");
}

#[test]
fn test_zopen_gzip() {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use gecco::io::compression::zopen_path;
    use std::io::Write;

    // Write gzipped data to a temp file and read back with zopen_path
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(b"compressed data").unwrap();
    let compressed = encoder.finish().unwrap();

    let tmp = std::env::temp_dir().join("gecco_test_gzip.gz");
    std::fs::write(&tmp, &compressed).unwrap();

    let mut reader = zopen_path(&tmp).unwrap();
    let mut buf = String::new();
    reader.read_to_string(&mut buf).unwrap();
    assert_eq!(buf, "compressed data");

    let _ = std::fs::remove_file(&tmp);
}

// ---------------------------------------------------------------------------
// HMMER helpers test
// ---------------------------------------------------------------------------

#[test]
fn test_hmm_relabel_pfam() {
    use gecco::hmmer::HMM;
    use std::path::PathBuf;

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
    assert_eq!(hmm.relabel("TIGR00001"), "TIGR00001");
}

#[test]
fn test_disentangle_keeps_best_pvalue() {
    use gecco::hmmer::disentangle;

    let mut gene = make_gene(
        "seq1",
        "p1",
        1,
        300,
        Strand::Coding,
        vec![
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
        None,
    );

    disentangle(&mut gene);
    assert_eq!(gene.protein.domains.len(), 1);
    assert_eq!(gene.protein.domains[0].name, "PF00002");
}

#[test]
fn test_filter_by_evalue() {
    use gecco::hmmer::filter_by_evalue;

    let mut genes = vec![make_gene(
        "seq1",
        "p1",
        1,
        300,
        Strand::Coding,
        vec![
            Domain {
                name: "good".into(),
                start: 1,
                end: 20,
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
                name: "bad".into(),
                start: 30,
                end: 50,
                hmm: "Pfam".into(),
                i_evalue: 1.0,
                pvalue: 0.5,
                probability: None,
                cluster_weight: None,
                go_terms: vec![],
                go_functions: vec![],
                qualifiers: BTreeMap::new(),
            },
        ],
        None,
    )];

    filter_by_evalue(&mut genes, 1e-5);
    assert_eq!(genes[0].protein.domains.len(), 1);
    assert_eq!(genes[0].protein.domains[0].name, "good");
}

// ---------------------------------------------------------------------------
// Merged GenBank output test
// ---------------------------------------------------------------------------

#[test]
fn test_genbank_merged_output() {
    use gecco::io::genbank::write_clusters_merged;

    // Two clusters on different source sequences.
    let cluster1 = Cluster {
        id: "seq1_cluster_1".to_string(),
        genes: vec![
            make_gene(
                "seq1",
                "seq1_1",
                100,
                400,
                Strand::Coding,
                vec![make_domain("PF00394", 10, 50, Some(0.95))],
                Some(0.9),
            ),
            make_gene(
                "seq1",
                "seq1_2",
                500,
                800,
                Strand::Coding,
                vec![make_domain("PF07690", 5, 30, Some(0.80))],
                Some(0.85),
            ),
        ],
        cluster_type: Some(ClusterType::new(["Polyketide"])),
        type_probabilities: {
            let mut m = BTreeMap::new();
            m.insert("Polyketide".to_string(), 0.96);
            m.insert("NRP".to_string(), 0.14);
            m
        },
    };

    let cluster2 = Cluster {
        id: "seq2_cluster_1".to_string(),
        genes: vec![make_gene(
            "seq2",
            "seq2_1",
            200,
            600,
            Strand::Reverse,
            vec![make_domain("PF00109", 15, 60, Some(0.88))],
            Some(0.92),
        )],
        cluster_type: Some(ClusterType::new(["NRP"])),
        type_probabilities: {
            let mut m = BTreeMap::new();
            m.insert("NRP".to_string(), 0.91);
            m.insert("Polyketide".to_string(), 0.05);
            m
        },
    };

    let mut source_seqs = BTreeMap::new();
    source_seqs.insert("seq1".to_string(), "A".repeat(1000));
    source_seqs.insert("seq2".to_string(), "C".repeat(1000));

    let mut buf = Vec::new();
    write_clusters_merged(&mut buf, &[cluster1, cluster2], &source_seqs).unwrap();
    let output = String::from_utf8(buf).unwrap();

    // Should contain two LOCUS entries (one per cluster).
    assert_eq!(
        output.matches("LOCUS").count(),
        2,
        "merged output should contain exactly 2 LOCUS entries"
    );

    // Both cluster IDs present.
    assert!(output.contains("seq1_cluster_1"));
    assert!(output.contains("seq2_cluster_1"));

    // Both have structured GECCO comments.
    assert_eq!(output.matches("GECCO-Data-START").count(), 2);
    assert_eq!(output.matches("GECCO-Data-END").count(), 2);

    // Cluster types present.
    assert!(output.contains("Polyketide"));
    assert!(output.contains("NRP"));

    // Both have CDS features.
    // cluster1 has 2 genes, cluster2 has 1 gene => 3 CDS total.
    assert_eq!(
        output.matches("\n     CDS ").count(),
        3,
        "expected 3 CDS features total"
    );

    // Both have misc_feature for domains.
    assert!(output.contains("PF00394"));
    assert!(output.contains("PF07690"));
    assert!(output.contains("PF00109"));

    // Version string present in both records.
    assert_eq!(
        output
            .matches(concat!("GECCO v", env!("CARGO_PKG_VERSION")))
            .count(),
        2
    );

    // Each record is terminated by //.
    assert_eq!(
        output.matches("\n//\n").count(),
        2,
        "each GenBank record should end with //"
    );
}

/// Run Python GECCO to produce merged GenBank output and compare against Rust.
///
/// This test requires Python GECCO to be installed (`pip install -e GECCO/`).
/// It is ignored by default; run with `cargo test -- --ignored` to execute.
#[test]
#[ignore]
fn test_genbank_merged_matches_python() {
    use gecco::io::genbank::write_clusters_merged;
    use std::path::Path;

    let script_path = Path::new("tests/data/gen_merged_gbk.py");
    if !script_path.exists() {
        panic!("Python helper script not found: {}", script_path.display());
    }

    // Run Python script to generate reference merged GenBank.
    let py_output = std::process::Command::new("python3")
        .arg(script_path)
        .output()
        .expect("failed to run python3");

    if !py_output.status.success() {
        let stderr = String::from_utf8_lossy(&py_output.stderr);
        panic!("Python script failed:\n{}", stderr);
    }

    let py_gbk = String::from_utf8(py_output.stdout).expect("Python output is not valid UTF-8");

    // Build equivalent clusters in Rust.
    let cluster1 = Cluster {
        id: "seq1_cluster_1".to_string(),
        genes: vec![make_gene(
            "seq1",
            "seq1_1",
            100,
            400,
            Strand::Coding,
            vec![make_domain("PF00394", 10, 50, None)],
            Some(0.9),
        )],
        cluster_type: Some(ClusterType::new(["Polyketide"])),
        type_probabilities: {
            let mut m = BTreeMap::new();
            m.insert("Polyketide".to_string(), 0.96);
            m
        },
    };

    let cluster2 = Cluster {
        id: "seq2_cluster_1".to_string(),
        genes: vec![make_gene(
            "seq2",
            "seq2_1",
            200,
            600,
            Strand::Reverse,
            vec![make_domain("PF00109", 15, 60, None)],
            Some(0.85),
        )],
        cluster_type: Some(ClusterType::new(["NRP"])),
        type_probabilities: {
            let mut m = BTreeMap::new();
            m.insert("NRP".to_string(), 0.91);
            m
        },
    };

    let mut source_seqs = BTreeMap::new();
    source_seqs.insert("seq1".to_string(), "A".repeat(1000));
    source_seqs.insert("seq2".to_string(), "C".repeat(1000));

    let mut buf = Vec::new();
    write_clusters_merged(&mut buf, &[cluster1, cluster2], &source_seqs).unwrap();
    let rs_gbk = String::from_utf8(buf).unwrap();

    // Extract sections from both outputs for structured comparison.
    // Rather than line-by-line (fragile due to formatting), compare key sections.
    let count_occurrences = |s: &str, pat: &str| -> usize { s.matches(pat).count() };

    // --- Structural checks: same number of records ---
    assert_eq!(
        count_occurrences(&py_gbk, "\n//\n"),
        count_occurrences(&rs_gbk, "\n//\n"),
        "different number of GenBank records"
    );

    // --- Both contain the same cluster IDs ---
    assert!(py_gbk.contains("seq1_cluster_1") && rs_gbk.contains("seq1_cluster_1"));
    assert!(py_gbk.contains("seq2_cluster_1") && rs_gbk.contains("seq2_cluster_1"));

    // --- Same number of CDS and misc_feature entries ---
    assert_eq!(
        count_occurrences(&py_gbk, "\n     CDS "),
        count_occurrences(&rs_gbk, "\n     CDS "),
        "different number of CDS features"
    );
    assert_eq!(
        count_occurrences(&py_gbk, "\n     misc_feature"),
        count_occurrences(&rs_gbk, "\n     misc_feature"),
        "different number of misc_feature entries"
    );

    // --- Same domain accessions present ---
    assert!(py_gbk.contains("PF00394") && rs_gbk.contains("PF00394"));
    assert!(py_gbk.contains("PF00109") && rs_gbk.contains("PF00109"));

    // --- Same cluster types in structured comments ---
    assert!(py_gbk.contains("Polyketide") && rs_gbk.contains("Polyketide"));
    assert!(py_gbk.contains("NRP") && rs_gbk.contains("NRP"));

    // --- Same non-zero probabilities ---
    assert!(py_gbk.contains("polyketide_probability") && rs_gbk.contains("polyketide_probability"));
    assert!(py_gbk.contains("nrp_probability") && rs_gbk.contains("nrp_probability"));

    // --- GECCO-Data structured comments in both records ---
    assert_eq!(count_occurrences(&py_gbk, "GECCO-Data-START"), 2);
    assert_eq!(count_occurrences(&rs_gbk, "GECCO-Data-START"), 2);

    // --- Same ORIGIN sequence content (case-insensitive) ---
    // Extract ORIGIN blocks and compare ignoring case and whitespace.
    let extract_origins = |s: &str| -> Vec<String> {
        let mut origins = Vec::new();
        let mut in_origin = false;
        let mut current = String::new();
        for line in s.lines() {
            if line.starts_with("ORIGIN") {
                in_origin = true;
                current.clear();
            } else if line.starts_with("//") {
                if in_origin {
                    origins.push(current.clone());
                    in_origin = false;
                }
            } else if in_origin {
                // Strip line numbers and spaces, keep only sequence chars.
                let seq: String = line.chars().filter(|c| c.is_alphabetic()).collect();
                current.push_str(&seq.to_lowercase());
            }
        }
        origins
    };

    let py_origins = extract_origins(&py_gbk);
    let rs_origins = extract_origins(&rs_gbk);
    assert_eq!(
        py_origins.len(),
        rs_origins.len(),
        "different number of ORIGIN blocks"
    );
    for (i, (py_seq, rs_seq)) in py_origins.iter().zip(rs_origins.iter()).enumerate() {
        assert_eq!(
            py_seq,
            rs_seq,
            "ORIGIN sequence differs for record {} (py len={}, rs len={})",
            i,
            py_seq.len(),
            rs_seq.len()
        );
    }

    // --- Report remaining formatting differences for visibility ---
    // This does not fail the test, just prints diffs for review.
    let py_lines: Vec<&str> = py_gbk.lines().collect();
    let rs_lines: Vec<&str> = rs_gbk.lines().collect();
    let max_len = py_lines.len().max(rs_lines.len());
    let mut diffs = Vec::new();
    for i in 0..max_len {
        let py_line = py_lines.get(i).copied().unwrap_or("<missing>");
        let rs_line = rs_lines.get(i).copied().unwrap_or("<missing>");
        if py_line != rs_line {
            diffs.push(format!("  line {:3}: py: {}", i + 1, py_line));
            diffs.push(format!("           rs: {}", rs_line));
        }
    }
    if !diffs.is_empty() {
        eprintln!(
            "\n=== Formatting differences ({} lines differ, not failing) ===\n{}",
            diffs.len() / 2,
            diffs.join("\n")
        );
    }
}
