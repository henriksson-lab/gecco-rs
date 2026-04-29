use clap::Parser;
use gecco::cli::{Cli, Commands};

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    // Configure logging based on verbosity
    let log_level = match (cli.verbose, cli.quiet) {
        (0, 0) => "info",
        (1, _) => "debug",
        (v, _) if v >= 2 => "trace",
        (_, 1) => "warn",
        (_, _) => "error",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level)).init();

    match cli.command {
        Commands::Run(args) => args.execute(),
        Commands::Annotate(args) => args.execute(),
        Commands::Predict(args) => args.execute(),
        Commands::Train(args) => args.execute(),
        Commands::Cv(args) => args.execute(),
        Commands::Convert(args) => args.execute(),
        Commands::BuildData(args) => args.execute(),
        Commands::UpdateInterpro(args) => args.execute(),
    }
}
