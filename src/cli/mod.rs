//! CLI commands for GECCO.

use clap::{Parser, Subcommand};

pub mod annotate;
pub mod build_data;
pub mod convert;
pub mod cv;
pub mod predict;
pub mod run;
pub mod train;
pub mod update_interpro;

#[derive(Parser)]
#[command(name = "gecco", version, about = "Gene Cluster prediction with Conditional Random Fields")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,

    /// Increase verbosity.
    #[arg(short, long, global = true, action = clap::ArgAction::Count)]
    pub verbose: u8,

    /// Reduce output.
    #[arg(short, long, global = true, action = clap::ArgAction::Count)]
    pub quiet: u8,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Run the full GECCO pipeline on a genome.
    Run(run::RunArgs),
    /// Annotate protein domains in a genome.
    Annotate(annotate::AnnotateArgs),
    /// Predict clusters from pre-annotated features.
    Predict(predict::PredictArgs),
    /// Train a new CRF model.
    Train(train::TrainArgs),
    /// Run cross-validation on training data.
    Cv(cv::CvArgs),
    /// Convert between output formats.
    Convert(convert::ConvertArgs),
    /// Download HMM databases and build data files.
    BuildData(build_data::BuildDataArgs),
    /// Rebuild interpro.json from upstream GO and InterPro databases.
    UpdateInterpro(update_interpro::UpdateInterProArgs),
}
