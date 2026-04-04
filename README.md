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
$ gecco run --genome some_genome.fna -o output_dir
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
| `gecco build-data` | Download and build required HMM databases |

### Key Options

- `-v` / `-q` -- Control verbosity (repeat for more: `-vv`, `-vvv`)
- `--jobs` -- Number of threads for parallelizable steps (default: auto-detect)
- `--cds` -- Minimum consecutive genes for a BGC region (default: 3)
- `--threshold` -- Minimum probability for a gene to be part of a BGC (default: 0.8)

## Results

GECCO-RS produces the same output files as Python GECCO:

- `{genome}.genes.tsv` -- Predicted genes with per-gene BGC probabilities
- `{genome}.features.tsv` -- Identified protein domains in tabular format
- `{genome}.clusters.tsv` -- Predicted cluster coordinates and biosynthetic types
- `{genome}_cluster_{N}.gbk` -- GenBank file per cluster with annotated proteins and domains

## Library Usage

GECCO-RS can be used as a Rust library. Add it to your `Cargo.toml`:

```toml
[dependencies]
gecco = { path = "path/to/gecco-rs" }
```

Then scan sequences for gene clusters:

```rust
use gecco::{Pipeline, SeqRecord};

fn main() -> anyhow::Result<()> {
    // Build a pipeline with model and HMM database paths
    let pipeline = Pipeline::builder()
        .crf_model("path/to/model.crfsuite")
        .hmm("path/to/Pfam.hmm")
        .threshold(0.8)
        .build()?;

    // Feed in DNA sequences
    let records = vec![SeqRecord {
        id: "contig_1".into(),
        seq: std::fs::read_to_string("genome.fna")?,
    }];

    // Run full pipeline: gene finding → HMMER → CRF → clustering
    let clusters = pipeline.scan(&records)?;

    for cluster in &clusters {
        println!("{}: {} genes, {:?}",
            cluster.id, cluster.genes.len(), cluster.cluster_type);
    }
    Ok(())
}
```

For finer control, individual stages are also available:

```rust
// Run stages separately
let mut genes = pipeline.find_genes(&records)?;
pipeline.annotate_domains(&mut genes)?;
let genes = pipeline.predict_probabilities(&genes)?;
let clusters = pipeline.extract_clusters(&genes);
```

Or skip gene finding and HMMER if you already have annotated genes:

```rust
use gecco::io::tables::FeatureTable;

let genes = FeatureTable::read_to_genes(std::fs::File::open("features.tsv")?)?;
let (genes, clusters) = pipeline.predict_from_genes(&genes)?;
```

## Architecture

The pipeline flows through 5 stages:

1. **Gene Finding** -- Predicts ORFs using `orphos-core` (Rust Prodigal)
2. **HMM Annotation** -- Annotates protein domains against Pfam/TIGRFAM using a pure Rust HMMER implementation
3. **CRF Prediction** -- Sliding-window feature extraction + CRF marginal inference for per-gene biosynthetic probability
4. **Cluster Refinement** -- Groups contiguous high-probability genes into clusters, trims edges
5. **Type Classification** -- Random Forest predicts biosynthetic type (Polyketide, NRP, Terpene, RiPP, etc.)

### Key Differences from Python GECCO

| Component | Python GECCO | GECCO-RS |
|-----------|-------------|----------|
| Gene finding | pyrodigal | orphos-core (pure Rust) |
| HMMER | pyhmmer (C bindings) | Pure Rust HMMER |
| CRF | sklearn-crfsuite (pickle) | crfs crate (CRFsuite binary format) |
| Random Forest | scikit-learn | smartcore |
| DataFrames | Polars | csv + serde |
| Parallelism | multiprocessing | rayon |

## Benchmarks

GECCO-RS includes benchmarks comparing Rust and Python implementations on
identical inputs. Both benchmarks use the same test genomes from the GECCO
test suite.

### Running Benchmarks

```console
# Rust pipeline benchmark (per-stage timing)
$ cargo run --release --bin bench_pipeline

# Rust full pipeline benchmark (end-to-end, excluding HMMER)
$ cargo run --release --bin bench_full

# Python comparison benchmark (requires GECCO installed)
$ python bench_python.py
```

### What Is Measured

**`bench_pipeline`** measures each stage individually:
- Gene finding on BGC0001737.fna (~34kb genome)
- HMMER annotation against a mini Pfam database
- CRF prediction (10-run average) on BGC0001866 features
- Feature extraction (1000-run average)

**`bench_full`** measures the full pipeline end-to-end (excluding HMMER annotation):
- Gene finding on BGC0001866.fna (~33kb genome)
- CRF prediction from pre-annotated features
- Cluster refinement
- Total wall time

**`bench_python.py`** runs the equivalent Python GECCO stages for direct comparison:
- Gene finding with PyrodigalFinder
- HMMER annotation with pyhmmer
- CRF prediction (10-run average)

## Reference

GECCO can be cited using the following publication:

> **Accurate de novo identification of biosynthetic gene clusters with GECCO**.
> Laura M Carroll, Martin Larralde, Jonas Simon Fleck, Ruby Ponnudurai, Alessio Milanese, Elisa Cappio Barazzone, Georg Zeller.
> bioRxiv 2021.05.03.442509; [doi:10.1101/2021.05.03.442509](https://doi.org/10.1101/2021.05.03.442509)

## License

This software is provided under the [GNU General Public License v3.0 or later](https://choosealicense.com/licenses/gpl-3.0/).
