//! Runtime materialization for data embedded with the `bundled-data` feature.

use std::io;
use std::path::{Path, PathBuf};

pub const HMM_SENTINEL: &str = "<bundled:Pfam.h3m>";

struct BundledFile {
    name: &'static str,
    bytes: &'static [u8],
}

const FILES: &[BundledFile] = &[
    BundledFile {
        name: "Pfam.h3m",
        bytes: include_bytes!(concat!(env!("OUT_DIR"), "/Pfam.h3m")),
    },
    BundledFile {
        name: "interpro.json",
        bytes: include_bytes!(concat!(env!("OUT_DIR"), "/interpro.json")),
    },
    BundledFile {
        name: "model.crfsuite",
        bytes: include_bytes!("../bundled-data/model.crfsuite"),
    },
    BundledFile {
        name: "domains.tsv",
        bytes: include_bytes!(concat!(env!("OUT_DIR"), "/domains.tsv")),
    },
    BundledFile {
        name: "type_classifier.rf.json",
        bytes: include_bytes!("../bundled-data/type_classifier.rf.json"),
    },
    BundledFile {
        name: "Pfam.ini",
        bytes: include_bytes!(concat!(env!("OUT_DIR"), "/Pfam.ini")),
    },
];

pub fn hmm_path() -> PathBuf {
    PathBuf::from(HMM_SENTINEL)
}

pub fn is_hmm_path(path: &Path) -> bool {
    path == Path::new(HMM_SENTINEL)
}

pub fn pfam_h3m() -> &'static [u8] {
    FILES[0].bytes
}

pub fn interpro_json() -> &'static [u8] {
    FILES[1].bytes
}

pub fn crf_model() -> &'static [u8] {
    FILES[2].bytes
}

pub fn domains_tsv() -> &'static str {
    std::str::from_utf8(FILES[3].bytes).expect("bundled domains.tsv is not UTF-8")
}

pub fn type_classifier_rf_json() -> &'static [u8] {
    FILES[4].bytes
}

/// Write embedded GECCO data into a stable cache directory and return that path.
pub fn materialize() -> io::Result<PathBuf> {
    let dir = std::env::temp_dir().join(format!(
        "gecco-bundled-data-{}-{}",
        env!("CARGO_PKG_VERSION"),
        bundle_size()
    ));
    std::fs::create_dir_all(&dir)?;

    for file in FILES {
        write_if_needed(&dir.join(file.name), file.bytes)?;
    }

    Ok(dir)
}

fn bundle_size() -> usize {
    FILES.iter().map(|file| file.bytes.len()).sum()
}

fn write_if_needed(path: &Path, bytes: &[u8]) -> io::Result<()> {
    if path
        .metadata()
        .map(|meta| meta.len() == bytes.len() as u64)
        .unwrap_or(false)
    {
        return Ok(());
    }

    let tmp = path.with_extension("tmp");
    std::fs::write(&tmp, bytes)?;
    std::fs::rename(tmp, path)?;
    Ok(())
}
