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

Before first use, download the HMM databases:

```console
$ gecco build-data
```

## Usage

Run the full pipeline (gene finding, domain annotation, CRF prediction,
cluster refinement, and type classification) on a genome:

```console
$ gecco run --genome some_genome.fna -o output_dir --model model.crfsuite
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

### Key Options

| Flag | Description | Default |
|------|-------------|---------|
| `-g, --genome` | Input genome (FASTA or GenBank) | required |
| `-o, --output-dir` | Output directory | `.` |
| `-j, --jobs` | Number of threads (0 = auto) | `0` |
| `--model` | CRF model file (`.crfsuite`) | embedded |
| `-p, --p_filter` | P-value cutoff for domains | `1e-9` |
| `-m, --threshold` | Cluster membership probability threshold | `0.8` |
| `-c, --cds` | Minimum annotated CDS per cluster | `3` |
| `-M, --mask` | Mask ambiguous nucleotides | off |
| `--hmm` | Additional HMM file(s) | none |

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
