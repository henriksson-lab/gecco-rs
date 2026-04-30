use std::error::Error;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::time::Duration;

use sha2::{Digest, Sha256};

const DOWNLOAD_LIMIT: u64 = 200 * 1024 * 1024;

struct BundledDownload {
    url: &'static str,
    output: &'static str,
    sha256: &'static str,
    gzip: bool,
}

const DOWNLOADS: &[BundledDownload] = &[
    BundledDownload {
        url: "https://github.com/zellerlab/GECCO/releases/download/v0.10.3/Pfam.h3m.gz",
        output: "Pfam.h3m",
        sha256: "4383aa81f2ba9ce1f78334f6e249e0a8e69b9c0e8d4d9e69d144350597dac9c5",
        gzip: true,
    },
    BundledDownload {
        url:
            "https://raw.githubusercontent.com/zellerlab/GECCO/v0.10.3/gecco/interpro/interpro.json",
        output: "interpro.json",
        sha256: "f9ce8964b150918adc09eca80fc6aa0d641f183c322d023fc4ff051cd45716e9",
        gzip: false,
    },
    BundledDownload {
        url: "https://raw.githubusercontent.com/zellerlab/GECCO/v0.10.3/gecco/types/domains.tsv",
        output: "domains.tsv",
        sha256: "26fd7eb6c6818b3a6006408e73b0a9fba3c1879c1ec82568d42abdabac70fe06",
        gzip: false,
    },
    BundledDownload {
        url: "https://raw.githubusercontent.com/zellerlab/GECCO/v0.10.3/gecco/hmmer/Pfam.ini",
        output: "Pfam.ini",
        sha256: "542ee03dda0d9a7be012079b5a69f9271cf1d06accc0de28430082db14b229b0",
        gzip: false,
    },
];

fn main() -> Result<(), Box<dyn Error>> {
    println!("cargo:rerun-if-env-changed=CARGO_FEATURE_BUNDLED_DATA");

    if std::env::var_os("CARGO_FEATURE_BUNDLED_DATA").is_none() {
        return Ok(());
    }

    let out_dir = PathBuf::from(std::env::var_os("OUT_DIR").ok_or("OUT_DIR is not set")?);
    for download in DOWNLOADS {
        println!("cargo:rerun-if-changed=build.rs");
        ensure_downloaded(download, &out_dir)?;
    }

    Ok(())
}

fn ensure_downloaded(download: &BundledDownload, out_dir: &Path) -> Result<(), Box<dyn Error>> {
    let path = out_dir.join(download.output);
    if path.exists() && file_sha256(&path)? == download.sha256 {
        return Ok(());
    }

    println!(
        "cargo:warning=downloading bundled GECCO data from {}",
        download.url
    );
    let mut bytes = http_get(download.url)?;
    if download.gzip {
        let mut decoder = flate2::read::GzDecoder::new(bytes.as_slice());
        let mut decompressed = Vec::new();
        decoder.read_to_end(&mut decompressed)?;
        bytes = decompressed;
    }

    let actual = bytes_sha256(&bytes);
    if actual != download.sha256 {
        return Err(format!(
            "SHA256 mismatch for {}: expected {}, got {}",
            download.output, download.sha256, actual
        )
        .into());
    }

    let tmp = path.with_extension("tmp");
    fs::write(&tmp, bytes)?;
    fs::rename(tmp, path)?;
    Ok(())
}

fn http_get(url: &str) -> Result<Vec<u8>, Box<dyn Error>> {
    let agent = ureq::Agent::new_with_config(
        ureq::config::Config::builder()
            .timeout_global(Some(Duration::from_secs(300)))
            .build(),
    );
    let bytes = agent
        .get(url)
        .call()?
        .into_body()
        .with_config()
        .limit(DOWNLOAD_LIMIT)
        .read_to_vec()?;
    Ok(bytes)
}

fn file_sha256(path: &Path) -> Result<String, Box<dyn Error>> {
    Ok(bytes_sha256(&fs::read(path)?))
}

fn bytes_sha256(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    format!("{:x}", hasher.finalize())
}
