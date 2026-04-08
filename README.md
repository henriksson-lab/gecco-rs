# GECCO-RS

## Overview

GECCO-RS is a Rust reimplementation of [GECCO](https://github.com/zellerlab/GECCO)
(Gene Cluster prediction with Conditional Random Fields), a fast and scalable
method for identifying putative novel Biosynthetic Gene Clusters (BGCs) in
genomic and metagenomic data using Conditional Random Fields (CRFs).

This port replaces all Python dependencies with native Rust crates, yielding
significant speedups (see [Benchmarks](#benchmarks) below) while preserving
compatibility with the original GECCO models and data formats.

## Building

GECCO-RS requires a Rust toolchain (1.70+). Build with:

```console
$ cargo build --release
```

Before first use, download the data files (Pfam HMM, InterPro metadata):

```console
$ gecco build-data
```

This creates a `gecco_data/` directory in the current working directory. At
runtime, GECCO looks for data files next to the binary (`gecco_data/` alongside
the `gecco` executable). You can override this with:

- `--data-dir /path/to/data` on any command
- The `GECCO_DATA_DIR` environment variable

## Usage

Run the full pipeline (gene finding, domain annotation, CRF prediction,
cluster refinement, and type classification) on a genome:

```console
$ gecco run --genome some_genome.fna -o output_dir --model model.crfsuite
```

To use data files from a custom location:

```console
$ gecco run --genome some_genome.fna --data-dir /opt/gecco_data
```

### Commands

| Command | Description |
|---------|-------------|
| `gecco run` | Full pipeline: gene finding -> annotation -> prediction -> clustering |
| `gecco annotate` | Domain annotation only (gene finding + HMMER) |
| `gecco predict` | Predict clusters from pre-annotated feature/gene tables |
| `gecco train` | Train a new CRF model from labeled data |
| `gecco cv` | K-fold or leave-one-type-out cross-validation |
| `gecco convert` | Format conversion (GenBank, FASTA, GFF3) |

### Global Options

| Flag | Description | Default |
|------|-------------|---------|
| `-v, --verbose` | Increase verbosity (repeat for more, e.g. `-vv`) | — |
| `-q, --quiet` | Reduce or disable console output | — |

### `gecco run`

Full pipeline: gene finding, domain annotation, CRF prediction, cluster
refinement, and type classification.

**Input:**

| Flag | Description | Default |
|------|-------------|---------|
| `-g, --genome` | Input genome file (FASTA or GenBank) | *required* |
| `--data-dir` | Data directory (HMM, CRF model, InterPro files) | `gecco_data/` next to binary |

**Gene calling:**

| Flag | Description | Default |
|------|-------------|---------|
| `-j, --jobs` | Number of threads (0 = auto-detect) | `0` |
| `-M, --mask` | Mask ambiguous nucleotides to prevent genes stretching across them | off |
| `--cds-feature` | Extract genes from annotated features instead of de-novo prediction | — |
| `--locus-tag` | Feature qualifier for naming extracted genes | `locus_tag` |

**Domain annotation:**

| Flag | Description | Default |
|------|-------------|---------|
| `--hmm` | Path to alternative HMM file(s) (can be repeated) | from data dir |
| `-e, --e-filter` | E-value cutoff for protein domains | — |
| `-p, --p-filter` | P-value cutoff for protein domains | `1e-9` |
| `--disentangle` | Disentangle overlapping domains | off |

**Cluster detection:**

| Flag | Description | Default |
|------|-------------|---------|
| `--model` | Path to alternative CRF model (`.crfsuite` format) | from data dir |
| `--no-pad` | Disable padding of short gene sequences | off |
| `-c, --cds` | Minimum coding sequences per cluster | `3` |
| `-m, --threshold` | Probability threshold for cluster membership | `0.8` |
| `-E, --edge-distance` | Minimum annotated genes separating cluster from sequence edge | `0` |
| `--no-trim` | Disable trimming genes without domain annotations at cluster edges | off |

**Output:**

| Flag | Description | Default |
|------|-------------|---------|
| `-o, --output-dir` | Output directory | `.` |
| `--force-tsv` | Write TSV files even when empty | off |
| `--merge-gbk` | Single GenBank file for all clusters instead of one per cluster | off |
| `--antismash-sideload` | Write AntiSMASH v6 sideload JSON file | off |

### `gecco annotate`

Run only gene finding and domain annotation (no cluster prediction).

| Flag | Description | Default |
|------|-------------|---------|
| `-g, --genome` | Input genome file | *required* |
| `-o, --output-dir` | Output directory | `.` |
| `--data-dir` | Data directory (HMM, InterPro files) | `gecco_data/` next to binary |
| `-j, --jobs` | Number of threads (0 = auto-detect) | `0` |
| `-M, --mask` | Mask ambiguous nucleotides | off |
| `--cds-feature` | Extract genes from annotated features | — |
| `--locus-tag` | Feature qualifier for naming genes | `locus_tag` |
| `--hmm` | Alternative HMM file(s) | from data dir |
| `-e, --e-filter` | E-value cutoff | — |
| `-p, --p-filter` | P-value cutoff | `1e-9` |
| `--disentangle` | Disentangle overlapping domains | off |
| `--force-tsv` | Write TSV files even when empty | off |

### `gecco predict`

Predict clusters from pre-annotated feature/gene tables.

| Flag | Description | Default |
|------|-------------|---------|
| `--genome` | Input genome file (for GenBank output) | *required* |
| `-g, --genes` | Gene coordinate table (TSV) | *required* |
| `-f, --features` | Domain annotation table(s) (can be repeated) | *required* |
| `-o, --output-dir` | Output directory | `.` |
| `--data-dir` | Data directory (CRF model) | `gecco_data/` next to binary |
| `-j, --jobs` | Number of threads (0 = auto-detect) | `0` |
| `-e, --e-filter` | E-value cutoff | — |
| `-p, --p-filter` | P-value cutoff | `1e-9` |
| `--model` | Alternative CRF model | from data dir |
| `--no-pad` | Disable feature padding | off |
| `-c, --cds` | Minimum coding sequences per cluster | `3` |
| `-m, --threshold` | Probability threshold | `0.8` |
| `-E, --edge-distance` | Minimum genes from sequence edge | `0` |
| `--no-trim` | Disable edge gene trimming | off |
| `--force-tsv` | Write TSV files even when empty | off |
| `--merge-gbk` | Single GenBank file for all clusters | off |
| `--antismash-sideload` | Write AntiSMASH sideload JSON | off |

### `gecco train`

Train a new CRF model from labeled annotation tables.

| Flag | Description | Default |
|------|-------------|---------|
| `-g, --genes` | Gene coordinate table (TSV) | *required* |
| `-f, --features` | Domain annotation table(s) (can be repeated) | *required* |
| `-c, --clusters` | Cluster annotation table (TSV) | *required* |
| `-o, --output-dir` | Output directory | `.` |
| `-e, --e-filter` | E-value cutoff | — |
| `-p, --p-filter` | P-value cutoff | `1e-9` |
| `--no-shuffle` | Disable data shuffling before fitting | off |
| `--seed` | Random number generator seed | `42` |
| `-W, --window-size` | CRF sliding window length | `5` |
| `--window-step` | CRF sliding window step | `1` |
| `--c1` | L1 regularization strength | `0.15` |
| `--c2` | L2 regularization strength | `0.15` |
| `--feature-type` | Feature extraction level (`protein` or `domain`) | `protein` |
| `--select` | Fraction of most significant features to select (0.0–1.0) | all |

### `gecco cv`

Cross-validation for model evaluation. Supports K-fold and Leave-One-Type-Out.

| Flag | Description | Default |
|------|-------------|---------|
| `-g, --genes` | Gene coordinate table (TSV) | *required* |
| `-f, --features` | Domain annotation table(s) (can be repeated) | *required* |
| `-c, --clusters` | Cluster annotation table (TSV) | *required* |
| `-o, --output` | Output file path | `cv.tsv` |
| `-e, --e-filter` | E-value cutoff | — |
| `-p, --p-filter` | P-value cutoff | `1e-9` |
| `--no-shuffle` | Disable data shuffling | off |
| `--seed` | Random number generator seed | `42` |
| `-W, --window-size` | CRF sliding window length | `5` |
| `--window-step` | CRF sliding window step | `1` |
| `--c1` | L1 regularization strength | `0.15` |
| `--c2` | L2 regularization strength | `0.15` |
| `--feature-type` | Feature extraction level (`protein` or `domain`) | `protein` |
| `--select` | Fraction of features to select | all |
| `--loto` | Use Leave-One-Type-Out instead of K-folds | off |
| `--splits` | Number of K-fold splits | `10` |

### `gecco convert`

Convert output files to other formats.

**`gecco convert gbk`** — Convert GenBank cluster files:

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input-dir` | Input directory containing `.gbk` files | *required* |
| `-o, --output-dir` | Output directory | same as input |
| `-f, --format` | Output format: `bigslice`, `fna`, or `faa` | *required* |

**`gecco convert clusters`** — Convert cluster tables:

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input-dir` | Input directory containing `.clusters.tsv` files | *required* |
| `-o, --output-dir` | Output directory | same as input |
| `-f, --format` | Output format: `gff` | *required* |

### `gecco build-data`

Download HMM databases and prepare data files for the pipeline.

| Flag | Description | Default |
|------|-------------|---------|
| `-o, --output-dir` | Output directory for data files | `gecco_data` |
| `-f, --force` | Force re-download even if files exist | off |

## Results

GECCO-RS produces the same output files as Python GECCO:

- `{genome}.genes.tsv` -- Predicted genes with per-gene BGC probabilities
- `{genome}.features.tsv` -- Identified protein domains in tabular format
- `{genome}.clusters.tsv` -- Predicted cluster coordinates and biosynthetic types
- `{genome}_cluster_{N}.gbk` -- GenBank file per cluster with annotated proteins and domains

## Architecture

The pipeline flows through 5 stages:

1. **Gene Finding** -- Predicts ORFs using a pure Rust port of Prodigal with rayon-parallel metagenomic model evaluation
2. **HMM Annotation** -- Annotates protein domains against Pfam/TIGRFAM using a pure Rust HMMER implementation (SIMD-accelerated SSE/AVX2)
3. **CRF Prediction** -- Sliding-window feature extraction + CRF marginal inference for per-gene biosynthetic probability
4. **Cluster Refinement** -- Groups contiguous high-probability genes into clusters, trims edges
5. **Type Classification** -- Random Forest predicts biosynthetic type (Polyketide, NRP, Terpene, RiPP, etc.)

### Key Differences from Python GECCO

| Component | Python GECCO | GECCO-RS |
|-----------|-------------|----------|
| Gene finding | pyrodigal (Cython/C) | Pure Rust Prodigal + rayon |
| HMMER | pyhmmer (C bindings) | Pure Rust HMMER (SSE/AVX2) |
| CRF | sklearn-crfsuite (pickle) | crfsuite-compliant-rs (CRFsuite binary format) |
| Random Forest | scikit-learn | smartcore |
| DataFrames | Polars | csv + serde |
| Parallelism | multiprocessing | rayon |

## Benchmarks

Benchmarked on a 5.3 Mbp bacterial genome (*Streptomyces* sp., GenBank CP157504.1, 5,401 predicted genes). Both tools run with `-j 4` on the same machine (Linux, x86_64).

### Performance

| Stage | Rust | Python | Speedup |
|-------|-----:|-------:|--------:|
| Gene finding | 5s | 9s | 1.8x |
| HMM annotation | 17s | 25s | 1.5x |
| CRF + clustering | 2s | 8s | 4.0x |
| **Total** | **25s** | **42s** | **1.7x** |

### Output Comparison

| Metric | Rust | Python | Match? |
|--------|-----:|-------:|--------|
| Predicted genes | 5,401 | 5,401 | Identical |
| Gene coordinates | -- | -- | Identical |
| Domain hits | 6,839 | 6,455 | +6% |
| Unique domain types | 1,655 | 1,623 | +2% |
| Clusters found | 9 | 3 | See below |

Gene finding produces identical results between the two implementations (same gene count, same coordinates). Rust finds slightly more domain hits due to minor numerical differences in the HMMER implementation.

Rust detects more clusters because its CRF marginal probabilities (from `crfsuite-compliant-rs`) run 5-15% higher than Python's `sklearn-crfsuite`, pushing additional regions above the default 0.8 threshold. All 3 Python clusters overlap with Rust clusters at matching genomic coordinates:

| Region | Rust | Python |
|--------|------|--------|
| 623 kbp | cluster_2 (avg 0.93) | cluster_1 (avg 0.91, Terpene) |
| 3,953 kbp | cluster_12 (avg 0.97) | cluster_3 (avg 0.97, Saccharide) |
| 4,138 kbp | cluster_14 (avg 0.98) | cluster_4 (avg 0.98, Unknown) |

### Running Benchmarks

```console
# Rust pipeline benchmark (per-stage timing)
$ cargo run --release --bin bench_pipeline

# Rust full pipeline benchmark (end-to-end)
$ cargo run --release --bin bench_full
```

## Dependencies

### Core Algorithm
- [crfsuite-compliant-rs](https://crates.io/crates/crfsuite-compliant-rs) -- CRF model loading, training, and marginal inference
- [prodigal](https://github.com/henriksson-lab/prodigal-rs) -- Gene prediction (pure Rust port with rayon parallelism)
- [hmmer-pure-rs](https://github.com/henriksson-lab/hmmer-pure-rs) -- HMMER3 domain search (SIMD-accelerated)
- [smartcore_proba](https://crates.io/crates/smartcore_proba) -- Random Forest type classifier

### Bio I/O
- [bio](https://crates.io/crates/bio) -- Bioinformatics utilities
- [noodles-fasta](https://crates.io/crates/noodles-fasta) -- FASTA parsing
- [gb-io](https://crates.io/crates/gb-io) -- GenBank I/O

## Reference

GECCO can be cited using the following publication:

> **Accurate de novo identification of biosynthetic gene clusters with GECCO**.
> Laura M Carroll, Martin Larralde, Jonas Simon Fleck, Ruby Ponnudurai, Alessio Milanese, Elisa Cappio Barazzone, Georg Zeller.
> bioRxiv 2021.05.03.442509; [doi:10.1101/2021.05.03.442509](https://doi.org/10.1101/2021.05.03.442509)

## License

This software is provided under the [GNU General Public License v3.0 or later](https://choosealicense.com/licenses/gpl-3.0/).
