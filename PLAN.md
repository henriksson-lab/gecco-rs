# Rust Dependencies for GECCO-RS

## Context
Porting GECCO (Python BGC prediction tool) to pure Rust. User already has Rust solutions for CRFsuite and HMMER. Scope: full pipeline including training. Prodigal: using `orphos-core` crate.

## Dependency Map

### Already Covered by User
- **CRFsuite** — CRF model training and inference
- **HMMER** — Pure Rust HMMER at `/home/mahogny/github/claude/hmmer-pure-rs` (workspace with crates: `hmmer`, `hmmer-core`, `hmmer-io`, `hmmer-alphabet`, `hmmer-simd`, `hmmer-cli`)

### Core Algorithm Crates Needed

| Need | Python Dep | Rust Crate(s) | Notes |
|------|-----------|---------------|-------|
| Gene finding | pyrodigal | **`orphos-core`** ([github.com/FullHuman/orphos](https://github.com/FullHuman/orphos)) | Rust Prodigal equivalent. Eliminates the largest custom code effort |
| Random Forest | scikit-learn | **`smartcore`** or **`linfa-trees`** | Need `fit()` + `predict_proba()`. Multi-label: run N binary classifiers. Custom `TypeBinarizer` (~30 lines) |
| Dense arrays | numpy | **`ndarray`** | Probability arrays, domain composition matrices, element-wise ops |
| Sparse matrices | scipy.sparse | **`sprs`** + **`npyz`** + **`zip`** | Load `.npz` training data (CSR format). ~30 lines custom loader |
| Fisher exact test | scipy.stats | **`statrs`** (hypergeometric) | ~40 lines custom. Also need BH-FDR correction (~20 lines) |

### I/O Crates Needed

| Need | Python Dep | Rust Crate(s) | Notes |
|------|-----------|---------------|-------|
| FASTA reading | biopython | **`noodles-fasta`** | Standard, mature |
| GenBank read/write | biopython | **`gb-io`** | May need ~100 lines custom code for structured comment blocks (`##GECCO-Data-START##`) |
| DNA/protein types | biopython | **`bio`** (rust-bio) | Alphabets, codon tables. Or hand-roll for Prodigal |
| TSV tables | polars | **`csv`** + **`serde`** | Simple typed TSV read/write. No complex DataFrame ops needed |
| JSON metadata | json | **`serde_json`** | InterPro entries, GO terms. Trivial with derive macros |
| Compression | gzip/bz2/xz/lz4 | **`flate2`**, **`bzip2`**, **`lzma-rs`**, **`lz4_flex`** | Port `zopen()` magic-byte detection (~40 lines) |
| Serialization | pickle | **`serde`** + **`bincode`** | New format for CRF/RF models. Existing pickle models not compatible — need re-export or re-train |

### CLI & Utility Crates

| Need | Python Dep | Rust Crate(s) |
|------|-----------|---------------|
| CLI framework | argparse | **`clap`** (v4, derive) |
| Progress bars | rich | **`indicatif`** |
| Parallelism | multiprocessing | **`rayon`** |
| Error handling | — | **`anyhow`** or **`thiserror`** |
| Hashing | hashlib | **`md5`** |
| Logging | — | **`log`** + **`env_logger`** |
| Lazy statics | — | `std::sync::LazyLock` (Rust 1.80+) or **`once_cell`** |

## Key Custom Code Required (beyond crate usage)

1. **Prodigal integration** — Wire `orphos-core` into the pipeline (gene coords, strand, protein translation)
2. **NPZ sparse matrix loader** — unzip + parse .npy arrays → `sprs` CSR matrix
3. **Fisher exact test** — from hypergeometric distribution + Benjamini-Hochberg FDR
4. **GenBank structured comments** — `##GECCO-Data-START##` blocks
5. **`zopen` equivalent** — magic-byte-based transparent decompression
6. **Model migration** — Python pickle models must be re-exported to bincode or re-trained in Rust

## Complete `Cargo.toml` Dependencies (estimated)

```toml
[dependencies]
# Core algorithm
orphos-core = "*"          # Prodigal gene finding
hmmer = { path = "../path/to/hmmer-pure-rs/crates/hmmer" }  # Pure Rust HMMER from /home/mahogny/github/claude/hmmer-pure-rs
ndarray = "0.16"
smartcore = "0.3"          # or linfa-trees
sprs = "0.11"
statrs = "0.17"

# Bio I/O
noodles-fasta = "0.42"
gb-io = "0.7"
bio = "2.0"

# Data I/O
csv = "1.3"
serde = { version = "1", features = ["derive"] }
serde_json = "1"
bincode = "1.3"
npyz = "0.8"
zip = "2"

# Compression
flate2 = "1.0"
bzip2 = "0.5"
lzma-rs = "0.3"
lz4_flex = "0.11"

# CLI
clap = { version = "4", features = ["derive"] }
indicatif = "0.17"

# Utility
rayon = "1.10"
anyhow = "1"
md5 = "0.7"
log = "0.4"
env_logger = "0.11"
```

## Verification
- Compare output of Rust `gecco run` against Python `gecco run` on the same input genome
- Unit tests for each pipeline stage matching Python test suite in `GECCO/tests/`
- Verify model serialization round-trip (save + load produces identical predictions)
