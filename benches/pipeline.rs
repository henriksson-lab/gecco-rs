//! Benchmark Rust GECCO pipeline stages.
//! Run with: cargo run --release --example bench_pipeline

use std::path::Path;
use std::time::Instant;

fn main() {
    // --- 1. Gene finding ---
    let genome = Path::new("GECCO/tests/test_orf/data/BGC0001737.fna");
    if !genome.exists() {
        eprintln!("Test data not found");
        return;
    }

    let records = gecco::io::genbank::read_sequences(genome).unwrap();

    let t0 = Instant::now();
    let finder = gecco::orf::ProdigalFinder {
        metagenome: true,
        cpus: 1,
        ..Default::default()
    };
    // Run in thread with larger stack for orphos-core
    let records_clone = records.clone();
    let genes = std::thread::Builder::new()
        .stack_size(16 * 1024 * 1024)
        .spawn(move || {
            use gecco::orf::ORFFinder;
            finder.find_genes(&records_clone).unwrap()
        })
        .unwrap()
        .join()
        .unwrap();
    let t1 = t0.elapsed();
    println!("Gene finding:  {:.3}s  ({} genes)", t1.as_secs_f64(), genes.len());

    // --- 2. HMMER annotation ---
    let hmm_path = Path::new("GECCO/tests/test_hmmer/data/minipfam.hmm");
    if hmm_path.exists() {
        use gecco::hmmer::{DomainAnnotator, PyHMMER, HMM};
        use gecco::interpro::InterPro;

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

        let mut genes_for_hmmer = genes.clone();
        let t0 = Instant::now();
        let annotator = PyHMMER::new(hmm);
        annotator.run(&mut genes_for_hmmer, &interpro, None).unwrap();
        let t1 = t0.elapsed();
        let n_dom: usize = genes_for_hmmer.iter().map(|g| g.protein.domains.len()).sum();
        println!("HMMER annot:   {:.3}s  ({} domains)", t1.as_secs_f64(), n_dom);
    }

    // --- 3. CRF prediction ---
    let model_path = Path::new("GECCO/gecco/crf/model.crfsuite");
    let feat_path = Path::new("GECCO/tests/test_cli/data/BGC0001866.features.tsv");
    if model_path.exists() && feat_path.exists() {
        use gecco::crf::backend::CrfSuiteModel;
        use gecco::crf::ClusterCRF;
        use gecco::io::tables::FeatureTable;

        let model = CrfSuiteModel::from_file(model_path).unwrap();
        let feat_genes =
            FeatureTable::read_to_genes(std::fs::File::open(feat_path).unwrap()).unwrap();

        let mut total = std::time::Duration::ZERO;
        let n_runs = 10;
        for _ in 0..n_runs {
            let mut crf = ClusterCRF::new("protein", 20, 1);
            crf.set_model(Box::new(
                CrfSuiteModel::from_file(model_path).unwrap(),
            ));
            let t0 = Instant::now();
            let _ = crf.predict_probabilities(&feat_genes, true, None).unwrap();
            total += t0.elapsed();
        }
        println!(
            "CRF predict:   {:.3}s  ({} runs avg, {} genes)",
            total.as_secs_f64() / n_runs as f64,
            n_runs,
            feat_genes.len()
        );

        // --- 4. Feature extraction ---
        use gecco::crf::features::extract_features_protein;

        let n_runs = 1000;
        let t0 = Instant::now();
        for _ in 0..n_runs {
            let _ = extract_features_protein(&feat_genes);
        }
        let t1 = t0.elapsed();
        println!(
            "Feature extract: {:.3}ms  ({} runs avg)",
            t1.as_secs_f64() / n_runs as f64 * 1000.0,
            n_runs
        );
    }
}
