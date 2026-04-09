//! Resolve the GECCO data directory containing HMM, CRF model, and InterPro files.
//!
//! Resolution order:
//! 1. Explicit path (from `--data-dir` CLI flag or library API)
//! 2. `GECCO_DATA_DIR` environment variable
//! 3. `gecco_data/` next to the running binary
//! 4. `gecco_data/` in the current working directory

use std::path::PathBuf;

/// Default data directory name placed next to the binary.
const DEFAULT_DIR_NAME: &str = "gecco_data";

/// Resolve the data directory path.
///
/// If `explicit` is `Some`, returns that path directly.
/// Otherwise checks `GECCO_DATA_DIR` env var, then falls back to
/// `gecco_data/` next to the current executable.
pub fn resolve(explicit: Option<&PathBuf>) -> PathBuf {
    if let Some(dir) = explicit {
        return dir.clone();
    }

    if let Ok(dir) = std::env::var("GECCO_DATA_DIR") {
        return PathBuf::from(dir);
    }

    // Next to the binary
    if let Ok(exe) = std::env::current_exe() {
        if let Some(parent) = exe.parent() {
            let candidate = parent.join(DEFAULT_DIR_NAME);
            if candidate.is_dir() {
                return candidate;
            }
        }
    }

    // Current working directory
    PathBuf::from(DEFAULT_DIR_NAME)
}

/// Expected file paths within the data directory.
pub fn hmm_path(data_dir: &std::path::Path) -> PathBuf {
    data_dir.join("Pfam.h3m")
}

pub fn interpro_path(data_dir: &std::path::Path) -> PathBuf {
    data_dir.join("interpro.json")
}

pub fn crf_model_path(data_dir: &std::path::Path) -> PathBuf {
    data_dir.join("model.crfsuite")
}
