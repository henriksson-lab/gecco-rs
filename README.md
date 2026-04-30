# GECCO-rs

## Overview

GECCO-rs is a Rust reimplementation of [GECCO](https://github.com/zellerlab/GECCO)
(Gene Cluster prediction with Conditional Random Fields), a fast and scalable
method for identifying putative novel Biosynthetic Gene Clusters (BGCs) in
genomic and metagenomic data using Conditional Random Fields (CRFs).

* 2026-04-29: Validated to give the same output on real life data. Still an early translation; compare with original gecco on your data before switching. On par or maybe 50% faster than original

## This is an LLM-mediated faithful (hopefully) translation, not the original code! 

Most users should probably first see if the existing original code works for them, unless they have reason otherwise. The original source
may have newer features and it has had more love in terms of fixing bugs. In fact, we aim to replicate bugs if they are present, for the
sake of reproducibility! (but then we might have added a few more in the process)

There are however cases when you might prefer this Rust version. We generally agree with [this manifesto](https://rewrites.bio/) but more specifically:
* We have had many issues with ensuring that our software works using existing containers (Docker, PodMan, Singularity). One size does not fit all and it eats our resources trying to keep up with every way of delivering software
* Common package managers do not work well. It was great when we had a few Linux distributions with stable procedures, but now there are just too many ecosystems (Homebrew, Conda). Conda has an NP-complete resolver which does not scale. Homebrew is only so-stable. And our dependencies in Python still break. These can no longer be considered professional serious options. Meanwhile, Cargo enables multiple versions of packages to be available, even within the same program(!)
* The future is the web. We deploy software in the web browser, and until now that has meant Javascript. This is a language where even the == operator is broken. Typescript is one step up, but a game changer is the ability to compile Rust code into webassembly, enabling performance and sharing of code with the backend. Translating code to Rust enables new ways of deployment and running code in the browser has especial benefits for science - researchers do not have deep pockets to run servers, so pushing compute to the user enables deployment that otherwise would be impossible
* Old CLI-based utilities are bad for the environment(!). A large amount of compute resources are spent creating and communicating via small files, which we can bypass by using code as libraries. Even better, we can avoid frequent reloading of databases by hoisting this stage, with up to 100x speedups in some cases. Less compute means faster compute and less electricity wasted
* LLM-mediated translations may actually be safer to use than the original code. This article shows that [running the same code on different operating systems can give somewhat different answers](https://doi.org/10.1038/nbt.3820). This is a gap that Rust+Cargo can reduce. Typesafe interfaces also reduce coding mistakes and error handling, as opposed to typical command-line scripting

But:

* **This approach should still be considered experimental**. The LLM technology is immature and has sharp corners. But there are opportunities to reap, and the genie is not going back into the bottle. This translation is as much aimed to learn how to improve the technology and get feedback on the results.
* Translations are not endorsed by the original authors unless otherwise noted. **Do not send bug reports to the original developers**. Use our Github issues page instead.
* **Do not trust the benchmarks on this page**. They are used to help evaluate the translation. If you want improved performance, you generally have to use this code as a library, and use the additional tricks it offers. We generally accept performance losses in order to reduce our dependency issues
* **Check the original Github pages for information about the package**. This README is kept sparse on purpose. It is not meant to be the primary source of information
* **If you are the author of the original code and wish to move to Rust, you can obtain ownership of this repository and crate**. Until then, our commitment is to offer an as-faithful-as-possible translation of a snapshot of your code. If we find serious bugs, we will report them to you. Otherwise we will just replicate them, to ensure comparability across studies that claim to use package XYZ v.666. Think of this like a fancy Ubuntu .deb-package of your software - that is how we treat it

This blurb might be out of date. Go to [this page](https://github.com/henriksson-lab/rustification) for the latest information and further information about how we approach translation


## Building

Build the command-line tool with:

```console
$ cargo build --release
```

Before first use, download the data files (Pfam HMM, InterPro metadata):

```console
$ gecco build-data
```

This creates a `gecco_data/` directory in the current working directory. At
runtime, GECCO looks for data files next to the binary (`gecco_data/` alongside
the `gecco` executable). The data directory also contains GECCO's exported type
classifier (`type_classifier.rf.json`) and its domain order (`domains.tsv`).
You can override this with:

- `--data-dir /path/to/data` on any command
- The `GECCO_DATA_DIR` environment variable

Alternatively, build with `--features bundled-data` to embed GECCO data in the
binary. The crate includes the Rust-specific converted models
(`model.crfsuite`, `type_classifier.rf.json`) and the build script downloads the
remaining original GECCO v0.10.3 assets from `zellerlab/GECCO` into Cargo's
`OUT_DIR`, verifies their SHA256 checksums, and embeds them at compile time.
This includes the large `Pfam.h3m.gz` release asset, so building this feature
requires network access.

## Command-Line Usage

Run the full pipeline on a genome:

```console
$ gecco run --genome genome.fna --output-dir output_dir
```

To use data files from a custom location:

```console
$ gecco run --genome genome.fna --data-dir /opt/gecco_data
```

Use `gecco <command> --help` for the full current option list.

| Command | Description |
|---------|-------------|
| `gecco run` | Full pipeline: gene finding, HMMER annotation, CRF prediction, clustering, and type classification |
| `gecco annotate` | Gene finding and HMMER annotation only |
| `gecco predict` | Predict clusters from pre-annotated feature/gene tables |
| `gecco train` | Train a new CRF model from labeled data |
| `gecco cv` | K-fold or leave-one-type-out cross-validation on labeled data |
| `gecco convert` | Format conversion (GenBank, FASTA, GFF3) |
| `gecco build-data` | Download and prepare the default data directory |
| `gecco update-interpro` | Rebuild InterPro metadata from upstream databases |

### Global Options

| Flag | Description | Default |
|------|-------------|---------|
| `-v, --verbose` | Increase verbosity (repeat for more, e.g. `-vv`) | — |
| `-q, --quiet` | Reduce or disable console output | — |

The `-f, --features` arguments below refer to GECCO domain feature tables, not
Cargo feature flags.

### `gecco run`

Run gene finding, domain annotation, CRF prediction, cluster refinement, and
type classification.

```console
$ gecco run --genome genome.fna --output-dir output_dir
```

| Flag | Description | Default |
|------|-------------|---------|
| `-g, --genome` | Input genome file (FASTA or GenBank) | *required* |
| `-o, --output-dir` | Output directory | `.` |
| `--data-dir` | Data directory (HMM, CRF model, InterPro files) | `gecco_data/` next to binary |
| `-j, --jobs` | Number of threads (0 = auto-detect) | `0` |
| `-M, --mask` | Mask ambiguous nucleotides | off |
| `--hmm` | Additional HMM file path; repeat for multiple databases | from data dir |
| `-e, --e-filter` | E-value cutoff for protein domains | — |
| `-p, --p-filter` | P-value cutoff for protein domains | `1e-9` |
| `--disentangle` | Disentangle overlapping domains | off |
| `--model` | Alternative CRF model file | from data dir |
| `--no-pad` | Disable padding of short gene sequences | off |
| `-c, --cds` | Minimum coding sequences per cluster | `3` |
| `-m, --threshold` | Probability threshold for cluster membership | `0.8` |
| `-E, --edge-distance` | Minimum genes separating a cluster from a sequence edge | `0` |
| `--no-trim` | Disable trimming genes without domain annotations | off |
| `--force-tsv` | Write TSV files even when empty | off |
| `--merge-gbk` | Write one merged GenBank file instead of one file per cluster | off |

### `gecco annotate`

Run only gene finding and domain annotation.

```console
$ gecco annotate --genome genome.fna --output-dir annotations
```

| Flag | Description | Default |
|------|-------------|---------|
| `-g, --genome` | Input genome file | *required* |
| `-o, --output-dir` | Output directory | `.` |
| `--data-dir` | Data directory (HMM, InterPro files) | `gecco_data/` next to binary |
| `-j, --jobs` | Number of threads (0 = auto-detect) | `0` |
| `-M, --mask` | Mask ambiguous nucleotides | off |
| `--hmm` | Additional HMM file path; repeat for multiple databases | from data dir |
| `-e, --e-filter` | E-value cutoff | — |
| `-p, --p-filter` | P-value cutoff | `1e-9` |
| `--disentangle` | Disentangle overlapping domains | off |
| `--force-tsv` | Write TSV files even when empty | off |

### `gecco predict`

Predict clusters from pre-annotated feature/gene tables.

```console
$ gecco predict --genome genome.fna --genes genome.genes.tsv --features genome.features.tsv --output-dir output_dir
```

| Flag | Description | Default |
|------|-------------|---------|
| `--genome` | Input genome file (for GenBank output) | *required* |
| `-g, --genes` | Gene coordinate table (TSV) | *required* |
| `-f, --features` | Domain annotation table(s); accepts multiple values | optional |
| `-o, --output-dir` | Output directory | `.` |
| `--data-dir` | Data directory (CRF model and type classifier data) | `gecco_data/` next to binary |
| `-j, --jobs` | Number of threads (0 = auto-detect) | `0` |
| `-e, --e-filter` | E-value cutoff | — |
| `-p, --p-filter` | P-value cutoff | `1e-9` |
| `--model` | Alternative CRF model | from data dir |
| `--no-pad` | Disable padding of short gene sequences | off |
| `-c, --cds` | Minimum coding sequences per cluster | `3` |
| `-m, --threshold` | Probability threshold | `0.8` |
| `-E, --edge-distance` | Minimum genes from sequence edge | `0` |
| `--no-trim` | Disable trimming genes without domain annotations | off |
| `--force-tsv` | Write TSV files even when empty | off |
| `--merge-gbk` | Single GenBank file for all clusters | off |

### `gecco train`

Train a new CRF model from labeled annotation tables.

```console
$ gecco train --genes training.genes.tsv --features training.features.tsv --clusters training.clusters.tsv --output-dir model_dir
```

| Flag | Description | Default |
|------|-------------|---------|
| `-g, --genes` | Gene coordinate table (TSV) | *required* |
| `-f, --features` | Domain annotation table(s); accepts multiple values | optional |
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

```console
$ gecco cv --genes training.genes.tsv --features training.features.tsv --clusters training.clusters.tsv --output cv.tsv
```

| Flag | Description | Default |
|------|-------------|---------|
| `-g, --genes` | Gene coordinate table (TSV) | *required* |
| `-f, --features` | Domain annotation table(s); accepts multiple values | optional |
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

```console
$ gecco build-data --output-dir gecco_data
```

| Flag | Description | Default |
|------|-------------|---------|
| `-o, --output-dir` | Output directory for data files | `gecco_data` |
| `-f, --force` | Force re-download even if files exist | off |

## Library Usage

GECCO-rs can be used as a Rust library. Add it to your `Cargo.toml`:

```toml
[dependencies]
gecco = { version = "0.5", default-features = false }
```

This pulls in only the core library without CLI dependencies (`clap`, `ureq`, etc.).

Then use the `Gecco` API to scan sequences for biosynthetic gene clusters:

```rust
use std::fs::File;
use std::path::Path;
use gecco::Gecco;
use gecco::io::genbank::read_sequences;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let pipeline = Gecco::builder()
        .threshold(0.8)
        .build()?;

    let records = read_sequences(Path::new("genome.fna"))?;

    // Runs gene finding, domain annotation, CRF prediction, and clustering.
    let results = pipeline.scan(&records)?;

    for cluster in &results.clusters {
        println!(
            "{}: {} genes, {}..{}",
            cluster.id,
            cluster.genes.len(),
            cluster.start(),
            cluster.end()
        );
    }

    results.write_gene_table(File::create("output.genes.tsv")?)?;
    results.write_feature_table(File::create("output.features.tsv")?)?;
    results.write_cluster_table(File::create("output.clusters.tsv")?)?;
    std::fs::create_dir_all("output_dir")?;
    results.write_cluster_gbks(Path::new("output_dir"))?;

    Ok(())
}
```

For more control, run individual pipeline stages separately:

```rust
let mut genes = pipeline.find_genes(&records)?;
pipeline.annotate_domains(&mut genes)?;
let genes = pipeline.predict_probabilities(&genes)?;
let clusters = pipeline.extract_clusters(&genes);
```

The builder supports many options — see the [GeccoBuilder](src/pipeline.rs) source for the full list:

```rust
let pipeline = Gecco::builder()
    .data_dir("/opt/gecco_data")
    .threshold(0.6)          // lower threshold → more clusters
    .jobs(4)                 // parallel threads
    .p_filter(1e-6)          // relaxed domain filtering
    .mask(true)              // mask ambiguous nucleotides
    .build()?;
```

## Results

GECCO-rs produces the same output files as Python GECCO:

- `{genome}.genes.tsv` -- Predicted genes with per-gene BGC probabilities
- `{genome}.features.tsv` -- Identified protein domains in tabular format
- `{genome}.clusters.tsv` -- Predicted cluster coordinates and biosynthetic types
- `{genome}_cluster_{N}.gbk` -- GenBank file per cluster with annotated proteins and domains

## Benchmarks

Benchmarked on a 5.3 Mbp bacterial genome (*Streptomyces* sp., GenBank CP157504.1, 5,401 predicted genes). Both tools run with `-j 4` on the same machine (Linux, x86_64).

### Performance

| Stage | Rust | Python | Speedup |
|-------|-----:|-------:|--------:|
| Gene finding | 5s | 9s | 1.8x |
| HMM annotation | 17s | 25s | 1.5x |
| CRF + clustering | 2s | 8s | 4.0x |
| **Total** | **25s** | **42s** | **1.7x** |

### Running Benchmarks

```console
# Rust pipeline benchmark (per-stage timing)
$ cargo run --release --features bench --bin bench_pipeline

# Rust full pipeline benchmark (end-to-end)
$ cargo run --release --features bench --bin bench_full
```

## Reference

GECCO can be cited using the following publication:

> **Accurate de novo identification of biosynthetic gene clusters with GECCO**.
> Laura M Carroll, Martin Larralde, Jonas Simon Fleck, Ruby Ponnudurai, Alessio Milanese, Elisa Cappio Barazzone, Georg Zeller.
> bioRxiv 2021.05.03.442509; [doi:10.1101/2021.05.03.442509](https://doi.org/10.1101/2021.05.03.442509)

## License

This software is provided under the [GNU General Public License v3.0 or later](https://choosealicense.com/licenses/gpl-3.0/).
