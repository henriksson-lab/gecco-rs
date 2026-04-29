//! Smoke tests exercising the pipeline stages with real test data.
//!
//! These tests use the test data from the Python GECCO test suite.

use std::path::Path;

use gecco::orf::ORFFinder;

// ---------------------------------------------------------------------------
// 1. Sequence loading
// ---------------------------------------------------------------------------

#[test]
fn test_load_fasta() {
    let path = Path::new("GECCO/tests/test_cli/data/BGC0001866.fna");
    if !path.exists() {
        eprintln!("Skipping: test data not found at {:?}", path);
        return;
    }
    let records = gecco::io::genbank::read_sequences(path).unwrap();
    assert!(!records.is_empty(), "should load at least 1 sequence");
    assert!(
        records[0].seq.len() > 1000,
        "sequence should be >1000 bp, got {}",
        records[0].seq.len()
    );
    eprintln!(
        "Loaded {} sequence(s), first: {} ({} bp)",
        records.len(),
        records[0].id,
        records[0].seq.len()
    );
}

#[test]
fn test_load_genbank() {
    let path = Path::new("GECCO/tests/test_orf/data/BGC0001377.gbk");
    if !path.exists() {
        eprintln!("Skipping: test data not found at {:?}", path);
        return;
    }
    let records = gecco::io::genbank::read_sequences(path).unwrap();
    assert!(!records.is_empty());
    eprintln!(
        "Loaded GenBank: {} ({} bp)",
        records[0].id,
        records[0].seq.len()
    );
}

// ---------------------------------------------------------------------------
// 2. Gene finding (Prodigal)
// ---------------------------------------------------------------------------

#[test]
fn test_prodigal_gene_finding() {
    // orphos-core uses deep recursion; run in a thread with larger stack
    std::thread::Builder::new()
        .stack_size(16 * 1024 * 1024) // 16 MB stack
        .spawn(test_prodigal_gene_finding_inner)
        .unwrap()
        .join()
        .unwrap();
}

fn test_prodigal_gene_finding_inner() {
    let path = Path::new("GECCO/tests/test_orf/data/BGC0001737.fna");
    if !path.exists() {
        eprintln!("Skipping: test data not found");
        return;
    }
    let records = gecco::io::genbank::read_sequences(path).unwrap();

    let finder = gecco::orf::ProdigalFinder {
        metagenome: true,
        cpus: 1,
        ..Default::default()
    };
    let genes = finder.find_genes(&records).unwrap();

    assert!(!genes.is_empty(), "should find genes");
    eprintln!("Found {} genes", genes.len());

    // Validate basic properties
    for gene in &genes {
        assert!(gene.start < gene.end, "start < end");
        assert!(!gene.protein.seq.is_empty(), "protein should be translated");
        assert!(
            gene.protein
                .seq
                .char_indices()
                .all(|(i, c)| c.is_ascii_alphabetic()
                    || (c == '*' && i + 1 == gene.protein.seq.len())),
            "protein should be amino acids with an optional terminal stop"
        );
    }

    // Python test expects ~10 genes from this file
    assert!(
        genes.len() >= 5 && genes.len() <= 20,
        "expected ~10 genes, got {}",
        genes.len()
    );
}

// ---------------------------------------------------------------------------
// 3. HMMER domain annotation
// ---------------------------------------------------------------------------

#[test]
fn test_hmmer_annotation() {
    use gecco::hmmer::{DomainAnnotator, PyHMMER, HMM};
    use gecco::interpro::InterPro;
    use gecco::model::{Gene, Protein, Strand};
    use std::collections::BTreeMap;

    let hmm_path = Path::new("GECCO/tests/test_hmmer/data/minipfam.hmm");
    let prot_path = Path::new("GECCO/tests/test_hmmer/data/proteins.faa");
    if !hmm_path.exists() || !prot_path.exists() {
        eprintln!("Skipping: HMMER test data not found");
        return;
    }

    // Load test proteins as genes
    let prot_data = std::fs::read_to_string(prot_path).unwrap();
    let mut genes: Vec<Gene> = Vec::new();
    let mut current_id = String::new();
    let mut current_seq = String::new();

    for line in prot_data.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if !current_id.is_empty() {
                genes.push(Gene {
                    source_id: "test".to_string(),
                    start: 1,
                    end: (current_seq.len() * 3 + 1) as i64,
                    strand: Strand::Coding,
                    protein: Protein::new(&current_id, &current_seq),
                    qualifiers: BTreeMap::new(),
                    probability: None,
                });
                current_seq.clear();
            }
            current_id = header.split_whitespace().next().unwrap_or("").to_string();
        } else {
            current_seq.push_str(line.trim());
        }
    }
    if !current_id.is_empty() {
        genes.push(Gene {
            source_id: "test".to_string(),
            start: 1,
            end: (current_seq.len() * 3 + 1) as i64,
            strand: Strand::Coding,
            protein: Protein::new(&current_id, &current_seq),
            qualifiers: BTreeMap::new(),
            probability: None,
        });
    }

    assert!(!genes.is_empty(), "should load test proteins");
    eprintln!("Loaded {} test proteins", genes.len());

    // Run HMMER
    let hmm = HMM {
        id: "Pfam".to_string(),
        version: "test".to_string(),
        url: String::new(),
        path: hmm_path.to_path_buf(),
        size: Some(10),
        relabel_with: Some("s/(PF\\d+).\\d+/\\1/".to_string()),
        md5: None,
    };

    let interpro = InterPro::from_json(b"[]").unwrap();
    let annotator = PyHMMER::new(hmm);
    annotator.run(&mut genes, &interpro, None).unwrap();

    let total_domains: usize = genes.iter().map(|g| g.protein.domains.len()).sum();
    eprintln!(
        "Found {} domains across {} genes",
        total_domains,
        genes.len()
    );

    // Python test: at least some genes should get domains
    let annotated = genes
        .iter()
        .filter(|g| !g.protein.domains.is_empty())
        .count();
    assert!(annotated > 0, "at least some genes should have domains");

    // Validate domain properties
    for gene in &genes {
        for domain in &gene.protein.domains {
            assert!(!domain.name.is_empty());
            assert_eq!(domain.hmm, "Pfam");
            assert!(domain.i_evalue >= 0.0);
            assert!(domain.pvalue >= 0.0);
            assert!(domain.start <= domain.end);
        }
    }
}

// ---------------------------------------------------------------------------
// 4. Table I/O roundtrip with real data
// ---------------------------------------------------------------------------

#[test]
fn test_load_feature_table() {
    use gecco::io::tables::FeatureTable;

    let path = Path::new("GECCO/tests/test_cli/data/BGC0001866.features.tsv");
    if !path.exists() {
        eprintln!("Skipping: feature table not found");
        return;
    }

    let genes = FeatureTable::read_to_genes(std::fs::File::open(path).unwrap()).unwrap();
    assert!(!genes.is_empty());

    let total_domains: usize = genes.iter().map(|g| g.protein.domains.len()).sum();
    eprintln!(
        "Loaded {} genes with {} domains from feature table",
        genes.len(),
        total_domains
    );

    // Validate
    for gene in &genes {
        assert!(!gene.source_id.is_empty());
        assert!(!gene.protein.id.is_empty());
        assert!(gene.start < gene.end);
    }
}

#[test]
fn test_load_gene_table() {
    use gecco::io::tables::GeneTable;

    let path = Path::new("GECCO/tests/test_cli/data/BGC0001866.genes.tsv");
    if !path.exists() {
        eprintln!("Skipping: gene table not found");
        return;
    }

    let genes = GeneTable::read_to_genes(std::fs::File::open(path).unwrap()).unwrap();
    assert!(!genes.is_empty());
    eprintln!("Loaded {} genes from gene table", genes.len());

    // Python test expects 23 genes
    assert!(genes.len() >= 20, "expected ~23 genes, got {}", genes.len());
}

// ---------------------------------------------------------------------------
// 5. Cluster refinement with real-ish data
// ---------------------------------------------------------------------------

#[test]
fn test_refine_from_feature_table() {
    use gecco::io::tables::FeatureTable;
    let path = Path::new("GECCO/tests/test_cli/data/BGC0001866.features.tsv");
    if !path.exists() {
        eprintln!("Skipping: feature table not found");
        return;
    }

    let mut genes = FeatureTable::read_to_genes(std::fs::File::open(path).unwrap()).unwrap();

    // Simulate CRF predictions: set high probability for all genes
    for gene in &mut genes {
        gene.probability = Some(0.95);
    }

    let refiner = gecco::refine::ClusterRefiner {
        threshold: 0.8,
        n_cds: 3,
        edge_distance: 0,
        trim: true,
        ..Default::default()
    };

    let clusters = refiner.iter_clusters(&genes);
    eprintln!(
        "Found {} cluster(s) from simulated predictions",
        clusters.len()
    );

    // With all probabilities at 0.95, should find at least one cluster
    assert!(
        !clusters.is_empty(),
        "should find clusters with high probabilities"
    );

    for cluster in &clusters {
        assert!(!cluster.genes.is_empty());
        assert!(cluster.start() < cluster.end());
        eprintln!(
            "  Cluster {} ({} genes, {}..{})",
            cluster.id,
            cluster.genes.len(),
            cluster.start(),
            cluster.end()
        );
    }
}

// ---------------------------------------------------------------------------
// 6. GenBank output from real data
// ---------------------------------------------------------------------------

#[test]
fn test_genbank_roundtrip() {
    use gecco::io::tables::FeatureTable;
    use gecco::model::{Cluster, ClusterType};
    use std::collections::BTreeMap;

    let feat_path = Path::new("GECCO/tests/test_cli/data/BGC0001866.features.tsv");
    let fna_path = Path::new("GECCO/tests/test_cli/data/BGC0001866.fna");
    if !feat_path.exists() || !fna_path.exists() {
        eprintln!("Skipping: test data not found");
        return;
    }

    let mut genes = FeatureTable::read_to_genes(std::fs::File::open(feat_path).unwrap()).unwrap();
    for gene in &mut genes {
        gene.probability = Some(0.95);
    }

    let records = gecco::io::genbank::read_sequences(fna_path).unwrap();
    let source_seq = &records[0].seq;

    let mut type_probs = BTreeMap::new();
    type_probs.insert("Polyketide".to_string(), 0.96);

    let cluster = Cluster {
        id: "test_cluster_1".to_string(),
        genes: genes[..5.min(genes.len())].to_vec(),
        cluster_type: Some(ClusterType::new(["Polyketide"])),
        type_probabilities: type_probs,
    };

    let mut buf = Vec::new();
    gecco::io::genbank::write_cluster_gbk(&mut buf, &cluster, Some(source_seq)).unwrap();

    let output = String::from_utf8(buf).unwrap();
    assert!(output.contains("LOCUS"));
    assert!(output.contains("test_cluster_1"));
    assert!(output.contains("GECCO-Data-START"));
    assert!(output.contains("Polyketide"));
    assert!(output.contains("CDS"));

    eprintln!("Generated GenBank output: {} bytes", output.len());
}

// ---------------------------------------------------------------------------
// 7. CRF model loading and prediction
// ---------------------------------------------------------------------------

#[test]
fn test_load_crfsuite_model() {
    use gecco::crf::backend::CrfSuiteModel;
    use gecco::crf::CrfModel;

    let model_path = Path::new("GECCO/gecco/crf/model.crfsuite");
    if !model_path.exists() {
        eprintln!("Skipping: CRFsuite model not found (run Python export first)");
        return;
    }

    let model = CrfSuiteModel::from_file(model_path).unwrap();
    assert_eq!(model.num_labels(), 2);

    let state_feats = model.state_features();
    assert!(!state_feats.is_empty(), "should have state features");
    eprintln!("Loaded {} state features", state_feats.len());

    // Check a known feature exists
    if let Some(&w) = state_feats.get(&("PF00750".to_string(), "0".to_string())) {
        eprintln!("PF00750 -> 0: weight = {:.6}", w);
        assert!((w - 0.042199).abs() < 0.01, "weight should match Python");
    }

    // Test marginal prediction with some features
    let features = vec![
        vec!["PF00750".to_string()],
        Vec::new(),
        vec!["PF13471".to_string()],
    ];

    let marginals = model.predict_marginals_single(&features);
    assert_eq!(marginals.len(), 3);

    for (i, m) in marginals.iter().enumerate() {
        let p0 = m.get("0").copied().unwrap_or(0.0);
        let p1 = m.get("1").copied().unwrap_or(0.0);
        eprintln!("  Position {}: P(0)={:.4}, P(1)={:.4}", i, p0, p1);
        // Probabilities should sum to ~1
        assert!(
            (p0 + p1 - 1.0).abs() < 0.01,
            "probabilities should sum to 1"
        );
    }

    // PF13471 has negative weight for class 0, so position 2 should favor class 1
    let p1_pos2 = marginals[2].get("1").copied().unwrap_or(0.0);
    eprintln!("PF13471 position P(1) = {:.4}", p1_pos2);
    assert!(p1_pos2 > 0.5, "PF13471 should favor biosynthetic class");
}

#[test]
fn test_crf_predict_on_real_features() {
    use gecco::crf::backend::CrfSuiteModel;
    use gecco::crf::ClusterCRF;
    use gecco::io::tables::FeatureTable;

    let model_path = Path::new("GECCO/gecco/crf/model.crfsuite");
    let feat_path = Path::new("GECCO/tests/test_cli/data/BGC0001866.features.tsv");
    if !model_path.exists() || !feat_path.exists() {
        eprintln!("Skipping: model or feature data not found");
        return;
    }

    let model = CrfSuiteModel::from_file(model_path).unwrap();
    let genes = FeatureTable::read_to_genes(std::fs::File::open(feat_path).unwrap()).unwrap();

    eprintln!("Loaded {} genes with features", genes.len());

    // Use window=20, step=1 matching the Python model
    let mut crf = ClusterCRF::new("protein", 20, 1);
    crf.set_model(Box::new(model));

    let predicted = crf.predict_probabilities(&genes, true, None).unwrap();
    assert_eq!(predicted.len(), genes.len());

    let mut has_high = false;
    for gene in &predicted {
        if let Some(p) = gene.average_probability() {
            if p > 0.8 {
                has_high = true;
            }
        }
    }

    eprintln!("Predictions:");
    for (i, gene) in predicted.iter().enumerate() {
        let p = gene.average_probability().unwrap_or(0.0);
        let n_dom = gene.protein.domains.len();
        if n_dom > 0 || p > 0.5 {
            eprintln!("  gene {}: {} domains, P={:.4}", i, n_dom, p);
        }
    }

    // BGC0001866 is a known BGC — should have high probability genes
    assert!(has_high, "should find high-probability genes in BGC0001866");
}

// ---------------------------------------------------------------------------
// 8. Full pipeline (without HMMER — use pre-annotated features)
// ---------------------------------------------------------------------------

#[test]
fn test_full_pipeline_from_features() {
    use gecco::crf::backend::CrfSuiteModel;
    use gecco::crf::ClusterCRF;
    use gecco::io::tables::FeatureTable;
    use gecco::refine::ClusterRefiner;

    let model_path = Path::new("GECCO/gecco/crf/model.crfsuite");
    let feat_path = Path::new("GECCO/tests/test_cli/data/BGC0001866.features.tsv");
    if !model_path.exists() || !feat_path.exists() {
        eprintln!("Skipping: model or feature data not found");
        return;
    }

    // 1. Load features
    let genes = FeatureTable::read_to_genes(std::fs::File::open(feat_path).unwrap()).unwrap();
    eprintln!("Step 1: Loaded {} genes", genes.len());

    // 2. CRF prediction
    let model = CrfSuiteModel::from_file(model_path).unwrap();
    let mut crf = ClusterCRF::new("protein", 20, 1);
    crf.set_model(Box::new(model));
    let predicted = crf.predict_probabilities(&genes, true, None).unwrap();
    eprintln!(
        "Step 2: Predicted probabilities for {} genes",
        predicted.len()
    );

    // 3. Cluster extraction
    let refiner = ClusterRefiner {
        threshold: 0.8,
        n_cds: 3,
        edge_distance: 0,
        trim: true,
        ..Default::default()
    };
    let clusters = refiner.iter_clusters(&predicted);
    eprintln!("Step 3: Found {} cluster(s)", clusters.len());

    // BGC0001866 is a single-cluster BGC
    assert!(
        !clusters.is_empty(),
        "should find at least 1 cluster in BGC0001866"
    );

    for cluster in &clusters {
        eprintln!(
            "  Cluster: {} ({} genes, {}..{}, avg_p={:.4})",
            cluster.id,
            cluster.genes.len(),
            cluster.start(),
            cluster.end(),
            cluster.average_probability().unwrap_or(0.0)
        );
    }
}
