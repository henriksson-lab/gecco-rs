//! GECCO: Gene Cluster prediction with Conditional Random Fields.
//!
//! This crate provides both a CLI tool and a library API for identifying
//! putative Biosynthetic Gene Clusters (BGCs) in genomic and metagenomic data.
//!
//! # Library usage
//!
//! The easiest way to use GECCO as a library is through the [`Pipeline`] struct:
//!
//! ```no_run
//! use gecco::Pipeline;
//! use gecco::orf::SeqRecord;
//!
//! let pipeline = Pipeline::builder()
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
//! let clusters = pipeline.scan(&records).unwrap();
//! for cluster in &clusters {
//!     println!("{}: {} genes", cluster.id, cluster.genes.len());
//! }
//! ```
//!
//! For finer control, individual pipeline stages are also available as public
//! methods on [`Pipeline`], or through the lower-level modules directly.

pub mod model;
pub mod interpro;
pub mod io;
pub mod orf;
pub mod hmmer;
pub mod crf;
pub mod refine;
pub mod types;
pub mod util;
pub mod pipeline;
pub mod cli;
pub mod data_dir;

// Re-export key types at crate root for convenience.
pub use model::{Cluster, Gene, Domain, Protein};
pub use orf::SeqRecord;
pub use pipeline::{Pipeline, PipelineBuilder};
