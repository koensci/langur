mod cmd_delta;
mod cmd_index;
mod cmd_inspect;
mod cmd_quant;

use clap::{Parser, Subcommand};
use std::process::ExitCode;

#[derive(Parser)]
#[command(
    name = "langur",
    version,
    about = "Langur — fast, statistically calibrated splice-graph quantification and differential splicing engine"
)]
struct Cli {
    /// Number of worker threads (default: all logical cores)
    #[arg(short = 't', long = "threads", global = true)]
    threads: Option<usize>,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build a zero-copy mmap-able splice-graph + minimizer + selective-decoy index (.lidx)
    Index(cmd_index::IndexArgs),

    /// Inspect metadata, magic header, section layout, and load latency of a .lidx index
    Inspect(cmd_inspect::InspectArgs),

    /// Quantify single-end or paired-end FASTQ reads against a .lidx index
    Quant(cmd_quant::QuantArgs),

    /// Differential splicing analysis via covariate-aware Beta-Binomial GLM
    Delta(cmd_delta::DeltaArgs),
}

fn main() -> ExitCode {
    let cli = Cli::parse();

    if let Some(t) = cli.threads {
        if let Err(e) = rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
        {
            eprintln!("warning: failed to configure thread pool ({e})");
        }
    }

    let res = match cli.command {
        Commands::Index(args) => cmd_index::run(args),
        Commands::Inspect(args) => cmd_inspect::run(args),
        Commands::Quant(args) => cmd_quant::run(args),
        Commands::Delta(args) => cmd_delta::run(args),
    };

    match res {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
