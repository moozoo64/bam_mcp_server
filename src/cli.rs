use std::path::PathBuf;

use anyhow::bail;
use clap::Parser;

/// BAM Pileup MCP Server — serves a genomic pileup tool over the Model Context Protocol.
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    /// Path to sorted, indexed BAM file (.bai index must exist alongside)
    #[arg(short = 'b', long)]
    pub bam: PathBuf,

    /// Path to FASTA reference file (.fai index must exist alongside)
    #[arg(short = 'r', long)]
    pub reference: PathBuf,

    /// Default window half-width in bp around query position
    #[arg(short = 'w', long, default_value_t = 75)]
    pub window: u32,

    /// Maximum reads to display in pileup
    #[arg(short = 'd', long, default_value_t = 50)]
    pub max_depth: u32,

    /// Minimum mapping quality filter
    #[arg(short = 'q', long, default_value_t = 0)]
    pub min_mapq: u8,

    /// Minimum base quality to show as uppercase
    #[arg(short = 'Q', long, default_value_t = 20)]
    pub min_baseq: u8,

    /// Run as HTTP server on the given address instead of stdio MCP transport
    /// (e.g. 127.0.0.1:8090)
    #[arg(long, value_name = "ADDR:PORT")]
    pub sse: Option<String>,

    /// Enable debug logging to stderr
    #[arg(long)]
    pub debug: bool,

    /// Write debug log to file instead of stderr
    #[arg(long)]
    pub log_file: Option<PathBuf>,
}

/// Validated, ready-to-use configuration derived from CLI args.
#[derive(Debug, Clone)]
pub struct AppConfig {
    pub bam_path: PathBuf,
    pub reference_path: PathBuf,
    pub window: u32,
    pub max_depth: u32,
    pub min_mapq: u8,
    pub min_baseq: u8,
    pub sse: Option<String>,
    pub debug: bool,
    pub log_file: Option<PathBuf>,
}

impl TryFrom<Args> for AppConfig {
    type Error = anyhow::Error;

    fn try_from(args: Args) -> anyhow::Result<Self> {
        // Validate BAM file
        if !args.bam.exists() {
            bail!("BAM file not found: {}", args.bam.display());
        }
        let bai = args.bam.with_extension(format!(
            "{}.bai",
            args.bam
                .extension()
                .and_then(|e| e.to_str())
                .unwrap_or("bam")
        ));
        // Support both <file>.bam.bai and <file>.bai
        let bai_alt = {
            let mut p = args.bam.clone();
            let stem = p
                .file_name()
                .and_then(|n| n.to_str())
                .map(|n| format!("{}.bai", n))
                .unwrap_or_default();
            p.set_file_name(stem);
            p
        };
        if !bai.exists() && !bai_alt.exists() {
            bail!(
                "BAM index not found: expected {} or {}",
                bai.display(),
                bai_alt.display()
            );
        }

        // Validate FASTA reference
        if !args.reference.exists() {
            bail!("Reference FASTA not found: {}", args.reference.display());
        }
        let fai = args.reference.with_extension(format!(
            "{}.fai",
            args.reference
                .extension()
                .and_then(|e| e.to_str())
                .unwrap_or("fa")
        ));
        let fai_alt = {
            let mut p = args.reference.clone();
            let stem = p
                .file_name()
                .and_then(|n| n.to_str())
                .map(|n| format!("{}.fai", n))
                .unwrap_or_default();
            p.set_file_name(stem);
            p
        };
        if !fai.exists() && !fai_alt.exists() {
            bail!(
                "FASTA index not found: expected {} or {}",
                fai.display(),
                fai_alt.display()
            );
        }

        if args.window == 0 {
            bail!("--window must be greater than 0");
        }
        if args.max_depth == 0 {
            bail!("--max-depth must be greater than 0");
        }

        Ok(AppConfig {
            bam_path: args.bam,
            reference_path: args.reference,
            window: args.window,
            max_depth: args.max_depth,
            min_mapq: args.min_mapq,
            min_baseq: args.min_baseq,
            sse: args.sse,
            debug: args.debug,
            log_file: args.log_file,
        })
    }
}
