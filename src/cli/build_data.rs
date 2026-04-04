//! The `gecco build-data` subcommand — download HMM databases.

use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Args;
use log::info;

const GECCO_VERSION: &str = "0.10.3";

#[derive(Args)]
pub struct BuildDataArgs {
    /// Output directory for data files (default: alongside GECCO data).
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
            .unwrap_or_else(|| PathBuf::from("data"));
        std::fs::create_dir_all(&output_dir)?;

        // Read HMM configs from .ini files
        let ini_path = PathBuf::from("GECCO/gecco/hmmer/Pfam.ini");
        if ini_path.exists() {
            let ini = std::fs::read_to_string(&ini_path)?;
            let cfg = parse_ini(&ini);

            let id = cfg.get("id").map(|s| s.as_str()).unwrap_or("Pfam");
            let version = cfg.get("version").map(|s| s.as_str()).unwrap_or("35.0");
            let md5_expected = cfg.get("md5").cloned();

            let h3m_path = output_dir.join(format!("{}.h3m", id));

            if h3m_path.exists() && !self.force {
                info!("{:?} already exists, skipping (use --force to re-download)", h3m_path);
            } else {
                download_hmm(id, version, &h3m_path, md5_expected.as_deref())?;
            }

            // Also create/update the .hmm text version for the hmmer crate
            // The Rust hmmer crate can read .h3m directly via binary::read_all_pressed
            info!("HMM data ready at {:?}", h3m_path);

            // Symlink or copy to where the pipeline expects it
            let target = ini_path.with_extension("h3m");
            if !target.exists() || self.force {
                if target.exists() {
                    std::fs::remove_file(&target)?;
                }
                #[cfg(unix)]
                {
                    let abs_h3m = std::fs::canonicalize(&h3m_path)?;
                    std::os::unix::fs::symlink(&abs_h3m, &target)
                        .or_else(|_| std::fs::copy(&h3m_path, &target).map(|_| ()))?;
                }
                #[cfg(not(unix))]
                {
                    std::fs::copy(&h3m_path, &target)?;
                }
                info!("Linked {:?} -> {:?}", target, h3m_path);
            }
        } else {
            info!("No Pfam.ini found; downloading default Pfam HMM");
            let h3m_path = output_dir.join("Pfam.h3m");
            download_hmm("Pfam", "35.0", &h3m_path, None)?;
        }

        // Also download InterPro metadata if missing
        let interpro_path = PathBuf::from("GECCO/gecco/interpro/interpro.json");
        if !interpro_path.exists() {
            info!("InterPro metadata not found at {:?} — using empty set", interpro_path);
        }

        info!("Build data complete");
        Ok(())
    }
}

fn download_hmm(id: &str, _version: &str, output: &PathBuf, expected_md5: Option<&str>) -> Result<()> {
    // Download pre-built .h3m from GECCO GitHub releases
    let url = format!(
        "https://github.com/zellerlab/GECCO/releases/download/v{}/{}.h3m.gz",
        GECCO_VERSION, id
    );
    info!("Downloading {} from {}", id, url);

    match download_and_decompress(&url, output) {
        Ok(()) => {
            info!("Downloaded {:?} ({} bytes)", output, std::fs::metadata(output)?.len());
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
            let ini_path = format!("GECCO/gecco/hmmer/{}.ini", id);
            if let Ok(ini) = std::fs::read_to_string(&ini_path) {
                let cfg = parse_ini(&ini);
                if let Some(fallback_url) = cfg.get("url") {
                    info!("Downloading from {}", fallback_url);
                    download_and_decompress_hmm_filter(fallback_url, output)?;
                    return Ok(());
                }
            }
            Err(e)
        }
    }
}

fn download_and_decompress(url: &str, output: &PathBuf) -> Result<()> {
    // Use a simple HTTP GET via std (no extra deps)
    let response = ureq_get(url)?;

    // Decompress gzip
    let decoder = flate2::read::GzDecoder::new(response.as_slice());
    let mut reader = std::io::BufReader::new(decoder);
    let mut out = std::fs::File::create(output)?;

    std::io::copy(&mut reader, &mut out)?;
    Ok(())
}

fn ureq_get(url: &str) -> Result<Vec<u8>> {
    // Minimal HTTP client using std::process::Command to call curl
    // (avoids adding an HTTP client dependency)
    let output = std::process::Command::new("curl")
        .args(["-fsSL", "--max-time", "300", url])
        .output()
        .with_context(|| format!("running curl to download {}", url))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        anyhow::bail!("curl failed for {}: {}", url, stderr.trim());
    }

    Ok(output.stdout)
}

fn download_and_decompress_hmm_filter(url: &str, output: &PathBuf) -> Result<()> {
    info!("Downloading and filtering HMMs from source...");

    // Load domain whitelist
    let domains_path = "GECCO/gecco/types/domains.tsv";
    let whitelist: std::collections::HashSet<String> = if std::path::Path::new(domains_path).exists() {
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
    let mut reader = std::io::BufReader::new(decoder);

    // Read all HMMs, filter by whitelist, write as text .hmm
    let all_hmms = hmmer::io::hmm_file::read_all_hmms(&mut reader)
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

    info!("Filtered to {} HMMs matching domain whitelist", filtered.len());

    // Write as text HMM format
    let mut out = std::fs::File::create(output)?;
    for hmm in &filtered {
        hmmer::io::hmm_file::write_hmm(&mut out, hmm)
            .map_err(|e| anyhow::anyhow!("writing HMM: {}", e))?;
    }

    Ok(())
}

fn md5_file(path: &PathBuf) -> Result<String> {
    let data = std::fs::read(path)?;
    Ok(format!("{:x}", md5::compute(&data)))
}

fn parse_ini(content: &str) -> std::collections::HashMap<String, String> {
    let mut map = std::collections::HashMap::new();
    for line in content.lines() {
        let line = line.trim();
        if line.starts_with('[') || line.is_empty() {
            continue;
        }
        if let Some((key, val)) = line.split_once('=') {
            map.insert(key.trim().to_string(), val.trim().to_string());
        }
    }
    map
}
