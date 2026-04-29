//! The `gecco build-data` subcommand — download HMM databases and data files.

use std::path::{Path, PathBuf};
use std::time::Duration;

use anyhow::{Context, Result};
use clap::Args;
use log::info;

const GECCO_VERSION: &str = "0.10.3";
const INTERPRO_URL: &str =
    "https://github.com/zellerlab/GECCO/raw/v0.10.3/gecco/interpro/interpro.json";
const MODEL_URL: &str = "https://github.com/henriksson-lab/gecco-rs/raw/main/model.crfsuite";

#[derive(Args)]
pub struct BuildDataArgs {
    /// Output directory for data files (default: ./gecco_data).
    #[arg(short, long)]
    pub output_dir: Option<PathBuf>,

    /// Force re-download even if files exist.
    #[arg(short, long)]
    pub force: bool,
}

impl BuildDataArgs {
    pub fn execute(&self) -> Result<()> {
        let output_dir = self
            .output_dir
            .clone()
            .unwrap_or_else(|| PathBuf::from("gecco_data"));
        std::fs::create_dir_all(&output_dir)?;

        // 1. Download Pfam HMM
        let h3m_path = output_dir.join("Pfam.h3m");
        if h3m_path.exists() && !self.force {
            info!(
                "{:?} already exists, skipping (use --force to re-download)",
                h3m_path
            );
        } else {
            // Try reading config from .ini if available, otherwise use defaults
            let ini_path = PathBuf::from("GECCO/gecco/hmmer/Pfam.ini");
            let (id, md5_expected, fallback_url) = if ini_path.exists() {
                let ini = std::fs::read_to_string(&ini_path)?;
                let cfg = parse_ini(&ini);
                (
                    cfg.get("id").cloned().unwrap_or_else(|| "Pfam".to_string()),
                    cfg.get("md5").cloned(),
                    cfg.get("url").cloned(),
                )
            } else {
                ("Pfam".to_string(), None, None)
            };
            download_hmm(
                &id,
                &h3m_path,
                md5_expected.as_deref(),
                fallback_url.as_deref(),
            )?;
        }

        // 2. Download InterPro metadata
        let interpro_path = output_dir.join("interpro.json");
        if interpro_path.exists() && !self.force {
            info!("{:?} already exists, skipping", interpro_path);
        } else {
            // Try copying from GECCO source tree first
            let source = PathBuf::from("GECCO/gecco/interpro/interpro.json");
            if source.exists() {
                info!("Copying InterPro metadata from {:?}", source);
                std::fs::copy(&source, &interpro_path)?;
            } else {
                info!("Downloading InterPro metadata");
                let data = ureq_get(INTERPRO_URL)?;
                std::fs::write(&interpro_path, &data)?;
            }
            info!("InterPro metadata ready at {:?}", interpro_path);
        }

        // 3. Copy CRF model if available
        let model_path = output_dir.join("model.crfsuite");
        if model_path.exists() && !self.force {
            info!("{:?} already exists, skipping", model_path);
        } else {
            // Check for existing .crfsuite model in common locations
            let candidates = [
                PathBuf::from("GECCO/gecco/crf/model.crfsuite"),
                PathBuf::from("model.crfsuite"),
                PathBuf::from("data/model.crfsuite"),
            ];
            let mut found = false;
            for src in &candidates {
                if src.exists() {
                    info!("Copying CRF model from {:?}", src);
                    std::fs::copy(src, &model_path)?;
                    found = true;
                    break;
                }
            }
            if !found {
                info!("Downloading CRF model from {}", MODEL_URL);
                let data = ureq_get(MODEL_URL)?;
                std::fs::write(&model_path, &data)?;
                info!("CRF model ready at {:?} ({} bytes)", model_path, data.len());
            }
        }

        // 4. Copy type-classifier data if available
        copy_local_data_file(
            &output_dir,
            "domains.tsv",
            &[
                PathBuf::from("GECCO/gecco/types/domains.tsv"),
                PathBuf::from("data/domains.tsv"),
            ],
            self.force,
        )?;
        copy_local_data_file(
            &output_dir,
            "type_classifier.rf.json",
            &[
                PathBuf::from("GECCO/gecco/types/type_classifier.rf.json"),
                PathBuf::from("data/type_classifier.rf.json"),
            ],
            self.force,
        )?;

        info!("Build data complete. Data directory: {:?}", output_dir);
        Ok(())
    }
}

fn copy_local_data_file(
    output_dir: &Path,
    file_name: &str,
    candidates: &[PathBuf],
    force: bool,
) -> Result<()> {
    let dst = output_dir.join(file_name);
    if dst.exists() && !force {
        info!("{:?} already exists, skipping", dst);
        return Ok(());
    }

    for src in candidates {
        if src.exists() {
            info!("Copying {} from {:?}", file_name, src);
            std::fs::copy(src, &dst)?;
            return Ok(());
        }
    }

    log::warn!(
        "{} not found in local GECCO data; type prediction will require this file",
        file_name
    );
    Ok(())
}

fn download_hmm(
    id: &str,
    output: &Path,
    expected_md5: Option<&str>,
    fallback_url: Option<&str>,
) -> Result<()> {
    // Download pre-built .h3m from GECCO GitHub releases
    let url = format!(
        "https://github.com/zellerlab/GECCO/releases/download/v{}/{}.h3m.gz",
        GECCO_VERSION, id
    );
    info!("Downloading {} from {}", id, url);

    match download_and_decompress(&url, output) {
        Ok(()) => {
            info!(
                "Downloaded {:?} ({} bytes)",
                output,
                std::fs::metadata(output)?.len()
            );
            if let Some(expected) = expected_md5 {
                let actual = md5_file(output)?;
                if actual != expected {
                    log::warn!(
                        "MD5 mismatch for {:?}: expected {}, got {} (this is normal for pre-filtered .h3m)",
                        output, expected, actual
                    );
                }
            }
            Ok(())
        }
        Err(e) => {
            log::warn!("GitHub download failed: {}. Trying fallback...", e);
            if let Some(url) = fallback_url {
                info!("Downloading from {}", url);
                download_and_decompress_hmm_filter(url, output)?;
                return Ok(());
            }
            Err(e)
        }
    }
}

fn download_and_decompress(url: &str, output: &Path) -> Result<()> {
    let response = ureq_get(url)?;

    // Decompress gzip
    let decoder = flate2::read::GzDecoder::new(response.as_slice());
    let mut reader = std::io::BufReader::new(decoder);
    let mut out = std::fs::File::create(output)?;

    std::io::copy(&mut reader, &mut out)?;
    Ok(())
}

fn ureq_get(url: &str) -> Result<Vec<u8>> {
    let agent = ureq::Agent::new_with_config(
        ureq::config::Config::builder()
            .timeout_global(Some(Duration::from_secs(300)))
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

fn download_and_decompress_hmm_filter(url: &str, output: &Path) -> Result<()> {
    info!("Downloading and filtering HMMs from source...");

    // Load domain whitelist
    let domains_path = "GECCO/gecco/types/domains.tsv";
    let whitelist: std::collections::HashSet<String> =
        if std::path::Path::new(domains_path).exists() {
            std::fs::read_to_string(domains_path)?
                .lines()
                .map(|l| l.trim().to_string())
                .filter(|l| !l.is_empty())
                .collect()
        } else {
            log::warn!("No domain whitelist found; downloading all HMMs");
            std::collections::HashSet::new()
        };

    let data = ureq_get(url)?;
    let decoder = flate2::read::GzDecoder::new(data.as_slice());
    let reader = std::io::BufReader::new(decoder);

    // Read all HMMs, filter by whitelist, write as text .hmm
    let all_hmms = hmmer::hmmfile::read_hmms(reader)
        .map_err(|e| anyhow::anyhow!("parsing HMM file: {}", e))?;

    info!("Read {} HMMs from source", all_hmms.len());

    let filtered: Vec<_> = if whitelist.is_empty() {
        all_hmms
    } else {
        all_hmms
            .into_iter()
            .filter(|h| {
                h.acc
                    .as_ref()
                    .map(|a| {
                        let base = a.split('.').next().unwrap_or(a);
                        whitelist.contains(base)
                    })
                    .unwrap_or(false)
            })
            .collect()
    };

    info!(
        "Filtered to {} HMMs matching domain whitelist",
        filtered.len()
    );

    // Write as text HMM format
    let mut out = std::fs::File::create(output)?;
    for hmm in &filtered {
        hmmer::hmmfile::write_hmm(&mut out, hmm)
            .map_err(|e| anyhow::anyhow!("writing HMM: {}", e))?;
    }

    Ok(())
}

fn md5_file(path: &Path) -> Result<String> {
    let data = std::fs::read(path)?;
    Ok(format!("{:x}", md5::compute(&data)))
}

fn parse_ini(content: &str) -> std::collections::HashMap<String, String> {
    let mut map = std::collections::HashMap::new();
    for line in content.lines() {
        let line = line.trim();
        if line.starts_with('[')
            || line.starts_with('#')
            || line.starts_with(';')
            || line.is_empty()
        {
            continue;
        }
        if let Some((key, val)) = line.split_once('=') {
            map.insert(key.trim().to_string(), val.trim().to_string());
        }
    }
    map
}
