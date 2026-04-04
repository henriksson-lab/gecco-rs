use clap::Parser;
use gecco::cli::{Cli, Commands};

fn main() -> anyhow::Result<()> {
    env_logger::init();

    let cli = Cli::parse();
    match cli.command {
        Commands::Run(args) => args.execute(),
    }
}
