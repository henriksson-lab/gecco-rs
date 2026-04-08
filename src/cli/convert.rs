//! The `gecco convert` subcommand — format conversion.

use std::io::Write;
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Args, Subcommand, ValueEnum};
use log::info;

#[derive(Args)]
pub struct ConvertArgs {
    #[command(subcommand)]
    pub command: ConvertCommand,
}

#[derive(Subcommand)]
pub enum ConvertCommand {
    /// Convert GenBank cluster files to other formats.
    Gbk(ConvertGbkArgs),
    /// Convert cluster TSV tables to other formats.
    Clusters(ConvertClustersArgs),
}

/// Output format for GenBank conversion.
#[derive(Debug, Clone, ValueEnum)]
pub enum GbkFormat {
    /// BiG-SLiCE feature table format.
    Bigslice,
    /// Nucleotide FASTA.
    Fna,
    /// Protein FASTA.
    Faa,
}

/// Output format for cluster table conversion.
#[derive(Debug, Clone, ValueEnum)]
pub enum ClusterFormat {
    /// GFF3 annotation format.
    Gff,
}

#[derive(Args)]
pub struct ConvertGbkArgs {
    /// Input directory containing GenBank files.
    #[arg(short, long)]
    pub input_dir: PathBuf,

    /// Output directory (defaults to input directory).
    #[arg(short, long)]
    pub output_dir: Option<PathBuf>,

    /// Output format.
    #[arg(short, long)]
    pub format: GbkFormat,
}

#[derive(Args)]
pub struct ConvertClustersArgs {
    /// Input directory containing cluster TSV files.
    #[arg(short, long)]
    pub input_dir: PathBuf,

    /// Output directory (defaults to input directory).
    #[arg(short, long)]
    pub output_dir: Option<PathBuf>,

    /// Output format.
    #[arg(short, long)]
    pub format: ClusterFormat,
}

impl ConvertArgs {
    pub fn execute(&self) -> Result<()> {
        match &self.command {
            ConvertCommand::Gbk(args) => convert_gbk(args),
            ConvertCommand::Clusters(args) => convert_clusters(args),
        }
    }
}

fn convert_gbk(args: &ConvertGbkArgs) -> Result<()> {
    let out_dir = args.output_dir.as_ref().unwrap_or(&args.input_dir);
    std::fs::create_dir_all(out_dir)?;

    info!(
        "Converting GenBank files from {:?} to {:?} ({:?})",
        args.input_dir, out_dir, args.format
    );

    // Find all .gbk files in input directory
    let entries: Vec<_> = std::fs::read_dir(&args.input_dir)?
        .filter_map(|e| e.ok())
        .filter(|e| {
            e.path()
                .extension()
                .map(|ext| ext == "gbk")
                .unwrap_or(false)
        })
        .collect();

    for entry in &entries {
        let path = entry.path();
        let stem = path.file_stem().unwrap_or_default().to_string_lossy();

        let gbk_data = std::fs::read(&path)?;
        let seqs: Vec<_> = gb_io::reader::SeqReader::new(&gbk_data[..])
            .filter_map(|r| r.ok())
            .collect();

        match args.format {
            GbkFormat::Fna => {
                let out_path = out_dir.join(format!("{}.fna", stem));
                let mut out = std::fs::File::create(&out_path)?;
                for seq in &seqs {
                    let name = seq.name.as_deref().unwrap_or("unknown");
                    writeln!(out, ">{}", name)?;
                    let dna = String::from_utf8_lossy(&seq.seq);
                    for chunk in dna.as_bytes().chunks(80) {
                        out.write_all(chunk)?;
                        writeln!(out)?;
                    }
                }
                info!("Wrote {:?}", out_path);
            }
            GbkFormat::Faa => {
                let out_path = out_dir.join(format!("{}.faa", stem));
                let mut out = std::fs::File::create(&out_path)?;
                for seq in &seqs {
                    for feat in &seq.features {
                        if feat.kind != std::borrow::Cow::from("CDS") {
                            continue;
                        }
                        let locus_tag = feat
                            .qualifier_values("locus_tag")
                            .next()
                            .unwrap_or("unknown");
                        let translation = feat
                            .qualifier_values("translation")
                            .next()
                            .unwrap_or("");
                        if !translation.is_empty() {
                            writeln!(out, ">{}", locus_tag)?;
                            for chunk in translation.as_bytes().chunks(80) {
                                out.write_all(chunk)?;
                                writeln!(out)?;
                            }
                        }
                    }
                }
                info!("Wrote {:?}", out_path);
            }
            GbkFormat::Bigslice => {
                // BiG-SLiCE expects region-named files
                let out_path = out_dir.join(format!("{}.region001.gbk", stem));
                std::fs::copy(&path, &out_path)?;
                info!("Wrote {:?}", out_path);
            }
        }
    }

    info!("Converted {} file(s)", entries.len());
    Ok(())
}

fn convert_clusters(args: &ConvertClustersArgs) -> Result<()> {
    let out_dir = args.output_dir.as_ref().unwrap_or(&args.input_dir);
    std::fs::create_dir_all(out_dir)?;

    info!(
        "Converting cluster tables from {:?} to {:?} ({:?})",
        args.input_dir, out_dir, args.format
    );

    // Find all .clusters.tsv files
    let entries: Vec<_> = std::fs::read_dir(&args.input_dir)?
        .filter_map(|e| e.ok())
        .filter(|e| {
            e.path()
                .to_string_lossy()
                .ends_with(".clusters.tsv")
        })
        .collect();

    for entry in &entries {
        let path = entry.path();
        let full_name = path.to_string_lossy();
        let base = full_name
            .strip_suffix(".clusters.tsv")
            .unwrap_or(&full_name);
        let base_name = std::path::Path::new(base)
            .file_name()
            .unwrap_or_default()
            .to_string_lossy();

        match args.format {
            ClusterFormat::Gff => {
                let out_path = out_dir.join(format!("{}.clusters.gff", base_name));
                let rdr = std::fs::File::open(&path)
                    .with_context(|| format!("opening {:?}", path))?;
                let mut tsv_rdr = csv::ReaderBuilder::new()
                    .delimiter(b'\t')
                    .from_reader(rdr);

                let mut out = std::fs::File::create(&out_path)?;
                writeln!(out, "##gff-version 3")?;

                for result in tsv_rdr.deserialize::<crate::io::tables::ClusterRow>() {
                    let row = result?;
                    writeln!(
                        out,
                        "{}\tGECCO\tgene_cluster\t{}\t{}\t.\t.\t.\tID={};type={}",
                        row.sequence_id,
                        row.start,
                        row.end,
                        row.cluster_id,
                        row.r#type,
                    )?;
                }
                info!("Wrote {:?}", out_path);
            }
        }
    }

    info!("Converted {} file(s)", entries.len());
    Ok(())
}
