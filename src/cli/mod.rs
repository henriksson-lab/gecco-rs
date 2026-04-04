//! CLI commands for GECCO.

use clap::{Parser, Subcommand};

pub mod run;

#[derive(Parser)]
#[command(name = "gecco", version, about = "Gene Cluster prediction with Conditional Random Fields")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Run the full GECCO pipeline on a genome.
    Run(run::RunArgs),
    // TODO: Annotate, Predict, Train, Cv, Convert
}
