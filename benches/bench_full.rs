use std::path::Path;
use std::time::Instant;

fn main() {
    let t_total = Instant::now();

    // 1. Gene finding (33kb genome)
    let genome = Path::new("GECCO/tests/test_cli/data/BGC0001866.fna");
    let records = gecco::io::genbank::read_sequences(genome).unwrap();

    let t0 = Instant::now();
    let records_clone = records.clone();
    let _genes = std::thread::Builder::new()
        .stack_size(16 * 1024 * 1024)
        .spawn(move || {
            use gecco::orf::ORFFinder;
            let finder = gecco::orf::ProdigalFinder {
                metagenome: true,
                cpus: 1,
                ..Default::default()
            };
            finder.find_genes(&records_clone).unwrap()
        })
        .unwrap()
        .join()
        .unwrap();
    println!(
        "Gene finding:    {:.3}s  ({} genes)",
        t0.elapsed().as_secs_f64(),
        _genes.len()
    );

    // 2. CRF prediction (from pre-annotated features)
    let model_path = Path::new("GECCO/gecco/crf/model.crfsuite");
    let feat_path = Path::new("GECCO/tests/test_cli/data/BGC0001866.features.tsv");

    let feat_genes = gecco::io::tables::FeatureTable::read_to_genes(
        std::fs::File::open(feat_path).unwrap(),
    )
    .unwrap();

    let model = gecco::crf::backend::CrfSuiteModel::from_file(model_path).unwrap();
    let mut crf = gecco::crf::ClusterCRF::new("protein", 20, 1);
    crf.set_model(Box::new(model));

    let t0 = Instant::now();
    let predicted = crf.predict_probabilities(&feat_genes, true, None).unwrap();
    println!(
        "CRF predict:     {:.3}s  ({} genes)",
        t0.elapsed().as_secs_f64(),
        predicted.len()
    );

    // 3. Cluster refinement
    let t0 = Instant::now();
    let refiner = gecco::refine::ClusterRefiner {
        threshold: 0.8,
        n_cds: 3,
        edge_distance: 0,
        trim: true,
        ..Default::default()
    };
    let clusters = refiner.iter_clusters(&predicted);
    println!(
        "Cluster refine:  {:.3}s  ({} clusters)",
        t0.elapsed().as_secs_f64(),
        clusters.len()
    );

    println!(
        "Total (no HMMER): {:.3}s",
        t_total.elapsed().as_secs_f64()
    );
}
