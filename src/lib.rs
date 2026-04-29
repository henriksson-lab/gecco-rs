//! GECCO: Gene Cluster prediction with Conditional Random Fields.
//!
//! This crate provides both a CLI tool and a library API for identifying
//! putative Biosynthetic Gene Clusters (BGCs) in genomic and metagenomic data.
//!
//! # Library usage
//!
//! The easiest way to use GECCO as a library is through the [`Gecco`] struct:
//!
//! ```no_run
//! use gecco::Gecco;
//! use gecco::orf::SeqRecord;
//!
//! let pipeline = Gecco::builder()
//!     .data_dir("gecco_data")
//!     .threshold(0.8)
//!     .build()
//!     .unwrap();
//!
//! let records = vec![SeqRecord {
//!     id: "contig_1".into(),
//!     seq: "ATGCCC...".into(),
//! }];
//!
//! let results = pipeline.scan(&records).unwrap();
//! for cluster in &results.clusters {
//!     println!("{}: {} genes", cluster.id, cluster.genes.len());
//! }
//! ```
//!
//! For finer control, individual pipeline stages are also available as public
//! methods on [`Gecco`], or through the lower-level modules directly.

#[cfg(feature = "bundled-data")]
pub mod bundled_data;
#[cfg(feature = "cli")]
pub mod cli;
pub mod crf;
pub mod data_dir;
pub mod hmmer;
pub mod interpro;
pub mod io;
pub mod model;
pub mod orf;
pub mod output;
pub mod pipeline;
pub mod refine;
pub mod sklearn_rf;
pub mod types;
pub mod util;

// Re-export key types at crate root for convenience.
pub use model::{Cluster, Domain, Gene, Protein};
pub use orf::SeqRecord;
pub use pipeline::{Gecco, GeccoBuilder, GeccoResults};
