//! Bitwise conformance tests against original Python GECCO full-pipeline TSV output.
//!
//! Per-cluster GenBank records contain wall-clock creation metadata and
//! package-version strings, so they are not a stable byte-for-byte regression
//! target without a fixed-clock writer mode.
//!
//! These tests are ignored by default because they are intended as a strict
//! translation gate. Run with:
//!
//!     cargo test --test conformance_full_pipeline -- --ignored

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use gecco::cli::run::RunArgs;

const GENOME: &str = "data/CP157504.1.fna";
const HMM: &str = "data/Pfam.h3m";
const MODEL: &str = "GECCO/gecco/crf/model.crfsuite";
const INTERPRO_DATA_DIR: &str = "GECCO/gecco/interpro";
const PYTHON_GOLDEN_DIR: &str = "data/bench_python_full";
const CHECKED_IN_RUST_DIR: &str = "data/bench_rust";
const BASE: &str = "CP157504.1";

#[test]
#[ignore]
fn checked_in_full_pipeline_fixture_is_bitwise_identical_to_python_gecco() {
    assert_output_dir_bitwise_eq(Path::new(PYTHON_GOLDEN_DIR), Path::new(CHECKED_IN_RUST_DIR));
}

#[test]
#[ignore]
fn generated_full_pipeline_output_is_bitwise_identical_to_python_gecco() {
    let result = std::thread::Builder::new()
        .name("generated-full-pipeline-conformance".to_string())
        .stack_size(32 * 1024 * 1024)
        .spawn(|| {
            std::panic::catch_unwind(
                generated_full_pipeline_output_is_bitwise_identical_to_python_gecco_inner,
            )
        })
        .expect("spawning large-stack conformance test thread")
        .join()
        .expect("joining large-stack conformance test thread");

    if let Err(payload) = result {
        std::panic::resume_unwind(payload);
    }
}

fn generated_full_pipeline_output_is_bitwise_identical_to_python_gecco_inner() {
    require_file(GENOME);
    require_file(HMM);
    require_file(MODEL);
    require_file(&format!("{INTERPRO_DATA_DIR}/interpro.json"));

    let output_dir = fresh_output_dir("gecco-rs-full-pipeline-conformance");

    let args = RunArgs {
        genome: PathBuf::from(GENOME),
        output_dir: output_dir.clone(),
        data_dir: Some(PathBuf::from(INTERPRO_DATA_DIR)),
        jobs: 1,
        mask: false,
        cds_feature: None,
        locus_tag: "locus_tag".to_string(),
        hmm: vec![PathBuf::from(HMM)],
        e_filter: None,
        p_filter: 1e-9,
        disentangle: false,
        model: Some(PathBuf::from(MODEL)),
        no_pad: false,
        cds: 3,
        threshold: 0.8,
        edge_distance: 0,
        no_trim: false,
        force_tsv: false,
        merge_gbk: false,
        antismash_sideload: false,
    };

    args.execute().expect("running Rust GECCO full pipeline");

    assert_output_dir_bitwise_eq(Path::new(PYTHON_GOLDEN_DIR), &output_dir);
}

fn assert_output_dir_bitwise_eq(expected_dir: &Path, actual_dir: &Path) {
    require_dir(expected_dir);
    require_dir(actual_dir);

    let expected_files = output_files(expected_dir);
    let actual_files = output_files(actual_dir);

    assert_eq!(
        expected_files, actual_files,
        "output file set differs\nexpected: {:#?}\nactual: {:#?}",
        expected_files, actual_files
    );

    for rel in expected_files {
        assert_file_bitwise_eq(&expected_dir.join(&rel), &actual_dir.join(&rel));
    }
}

fn output_files(dir: &Path) -> BTreeSet<PathBuf> {
    std::fs::read_dir(dir)
        .unwrap_or_else(|err| panic!("reading output dir {}: {err}", dir.display()))
        .map(|entry| {
            let entry = entry.unwrap_or_else(|err| {
                panic!("reading entry from output dir {}: {err}", dir.display())
            });
            let file_type = entry.file_type().unwrap_or_else(|err| {
                panic!("reading file type for {}: {err}", entry.path().display())
            });
            assert!(
                file_type.is_file(),
                "unexpected non-file output: {}",
                entry.path().display()
            );
            entry.file_name().into()
        })
        .filter(|path: &PathBuf| {
            let name = path.to_string_lossy();
            name == format!("{BASE}.clusters.tsv")
                || name == format!("{BASE}.features.tsv")
                || name == format!("{BASE}.genes.tsv")
        })
        .collect()
}

fn assert_file_bitwise_eq(expected: &Path, actual: &Path) {
    let expected_bytes = std::fs::read(expected)
        .unwrap_or_else(|err| panic!("reading expected file {}: {err}", expected.display()));
    let actual_bytes = std::fs::read(actual)
        .unwrap_or_else(|err| panic!("reading actual file {}: {err}", actual.display()));

    if expected_bytes == actual_bytes {
        return;
    }

    let first_diff = expected_bytes
        .iter()
        .zip(actual_bytes.iter())
        .position(|(left, right)| left != right)
        .unwrap_or_else(|| expected_bytes.len().min(actual_bytes.len()));

    let expected_line = byte_line_number(&expected_bytes, first_diff);
    let actual_line = byte_line_number(&actual_bytes, first_diff);

    panic!(
        "bitwise mismatch: {} != {}\nexpected bytes: {}, actual bytes: {}\nfirst differing byte offset: {}\nexpected line: {}, actual line: {}",
        expected.display(),
        actual.display(),
        expected_bytes.len(),
        actual_bytes.len(),
        first_diff,
        expected_line,
        actual_line
    );
}

fn byte_line_number(bytes: &[u8], offset: usize) -> usize {
    bytes[..offset.min(bytes.len())]
        .iter()
        .filter(|&&byte| byte == b'\n')
        .count()
        + 1
}

fn fresh_output_dir(prefix: &str) -> PathBuf {
    let mut dir = std::env::temp_dir();
    dir.push(format!(
        "{}-{}-{}",
        prefix,
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("system time before UNIX_EPOCH")
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir)
        .unwrap_or_else(|err| panic!("creating output dir {}: {err}", dir.display()));
    dir
}

fn require_file(path: &str) {
    assert!(
        Path::new(path).is_file(),
        "required conformance fixture is missing: {path}"
    );
}

fn require_dir(path: &Path) {
    assert!(
        path.is_dir(),
        "required conformance directory is missing: {}",
        path.display()
    );
}
