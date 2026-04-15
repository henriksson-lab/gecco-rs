use std::path::Path;
use std::time::Instant;

fn main() {
    let genome_path = std::env::args().nth(1);
    let genome = match &genome_path {
        Some(p) => Path::new(p),
        None => Path::new("data/CP157504.1.fna"),
    };
    if !genome.exists() {
        eprintln!("Genome not found: {}", genome.display());
        eprintln!("Usage: bench_full [genome.fna]");
        eprintln!("Default: data/CP157504.1.fna");
        std::process::exit(1);
    }

    let cpus: usize = std::env::var("BENCH_CPUS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(4);
    println!("Genome: {} ({} bytes)", genome.display(), std::fs::metadata(genome).unwrap().len());
    println!("Threads: {}", cpus);
    println!();

    let t_total = Instant::now();

    // 1. Gene finding
    let records = gecco::io::genbank::read_sequences(genome).unwrap();
    let t0 = Instant::now();
    let records_clone = records.clone();
    let mut genes = std::thread::Builder::new()
        .stack_size(16 * 1024 * 1024)
        .spawn(move || {
            use gecco::orf::ORFFinder;
            let finder = gecco::orf::ProdigalFinder {
                metagenome: true,
                cpus,
                ..Default::default()
            };
            finder.find_genes(&records_clone).unwrap()
        })
        .unwrap()
        .join()
        .unwrap();
    println!("Gene finding:    {:>6.1}s  ({} genes)", t0.elapsed().as_secs_f64(), genes.len());

    // 2. HMM annotation
    let hmm_path = Path::new("GECCO/gecco/hmmer/Pfam.h3m");
    if hmm_path.exists() {
        use gecco::hmmer::{DomainAnnotator, PyHMMER, HMM};
        use gecco::interpro::InterPro;

        let hmm = HMM {
            id: "Pfam".to_string(),
            version: "35.0".to_string(),
            url: String::new(),
            path: hmm_path.to_path_buf(),
            size: Some(2766),
            relabel_with: Some("s/(PF\\d+).\\d+/\\1/".to_string()),
            md5: None,
        };
        let interpro_path = Path::new("GECCO/gecco/interpro/interpro.json");
        let interpro = if interpro_path.exists() {
            InterPro::from_json(&std::fs::read(interpro_path).unwrap()).unwrap()
        } else {
            InterPro::from_json(b"[]").unwrap()
        };

        let t0 = Instant::now();
        let annotator = PyHMMER::new(hmm);
        annotator.run(&mut genes, &interpro, None).unwrap();
        let n_dom: usize = genes.iter().map(|g| g.protein.domains.len()).sum();
        println!("HMM annotation:  {:>6.1}s  ({} domains)", t0.elapsed().as_secs_f64(), n_dom);
    } else {
        eprintln!("Warning: {} not found, skipping HMM annotation", hmm_path.display());
    }

    // 3. CRF prediction
    let model_path = Path::new("GECCO/gecco/crf/model.crfsuite");
    let model = gecco::crf::backend::CrfSuiteModel::from_file(model_path).unwrap();
    let mut crf = gecco::crf::ClusterCRF::new("protein", 20, 1);
    crf.set_model(Box::new(model));

    let t0 = Instant::now();
    let predicted = crf.predict_probabilities(&genes, true, None).unwrap();
    println!("CRF predict:     {:>6.1}s  ({} genes)", t0.elapsed().as_secs_f64(), predicted.len());

    // 4. Cluster refinement
    let t0 = Instant::now();
    let refiner = gecco::refine::ClusterRefiner {
        threshold: 0.8,
        n_cds: 3,
        edge_distance: 0,
        trim: true,
        ..Default::default()
    };
    let clusters = refiner.iter_clusters(&predicted);
    println!("Cluster refine:  {:>6.1}s  ({} clusters)", t0.elapsed().as_secs_f64(), clusters.len());

    println!();
    println!("Total:           {:>6.1}s", t_total.elapsed().as_secs_f64());
}
