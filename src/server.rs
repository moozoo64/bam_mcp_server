use std::sync::Arc;
use std::time::Instant;

use rmcp::handler::server::wrapper::Parameters;
use rmcp::{tool, tool_router};
use schemars::JsonSchema;
use serde::Deserialize;
use tracing::debug;

use noodles::bam;

use crate::cli::AppConfig;
use crate::pileup::expand_reads;
use crate::query::{contig_length, query_region};
use crate::reference::ReferenceReader;
use crate::render::{RenderOpts, render_pileup};

/// Input parameters for the `query_pileup` MCP tool.
#[derive(Debug, Deserialize, JsonSchema)]
pub struct QueryPileupParams {
    /// Chromosome name exactly as it appears in the BAM header (e.g. chr1, 1, chrM)
    pub chrom: String,

    /// 1-based genomic position to query
    pub position: i64,

    /// Half-width of region in bp (10–500): this many bases are shown on each side of
    /// the queried base, 2*window+1 in total. Overrides --window CLI default.
    pub window: Option<u32>,

    /// Minimum mapping quality (0–60). Overrides --min-mapq CLI default.
    pub min_mapq: Option<u8>,

    /// Show strand (+/-) column per read
    pub show_strand: Option<bool>,

    /// Show mapping quality column per read
    pub show_mapq: Option<bool>,
}

/// Allowed range for the per-call `window` parameter (bp).
const MIN_WINDOW: u32 = 10;
const MAX_WINDOW: u32 = 500;

/// Run a synchronous tool body on Tokio's blocking pool and turn its result
/// into the tool's text output.
///
/// The work does blocking file I/O, so it must not run on an async worker.
/// Running it in its own task also contains panics: a panic becomes an error
/// reply instead of a request that never gets answered.
async fn run_blocking_tool<F>(name: &str, f: F) -> String
where
    F: FnOnce() -> anyhow::Result<String> + Send + 'static,
{
    match tokio::task::spawn_blocking(f).await {
        Ok(Ok(output)) => output,
        Ok(Err(e)) => {
            tracing::error!("{name}: {e:#}");
            format!("Error: {e:#}")
        }
        Err(join_err) => {
            let msg = if join_err.is_panic() {
                let payload = join_err.into_panic();
                payload
                    .downcast_ref::<&str>()
                    .map(|s| s.to_string())
                    .or_else(|| payload.downcast_ref::<String>().cloned())
                    .unwrap_or_else(|| "unknown panic".to_string())
            } else {
                join_err.to_string()
            };
            tracing::error!("{name}: internal error: {msg}");
            format!("Error: internal error in {name}: {msg}")
        }
    }
}

#[derive(Clone)]
pub struct PileupServer {
    pub config: Arc<AppConfig>,
}

impl PileupServer {
    pub fn new(config: AppConfig) -> Self {
        Self {
            config: Arc::new(config),
        }
    }

    pub fn from_arc(config: Arc<AppConfig>) -> Self {
        Self { config }
    }
}

#[tool_router(server_handler)]
impl PileupServer {
    /// Query a BAM file at a genomic position and return a text pileup showing aligned
    /// short reads with dots for reference matches and explicit bases for mismatches,
    /// insertions, and deletions.
    #[tool(description = "Query a BAM file at a genomic position and return a text pileup")]
    async fn query_pileup(&self, Parameters(params): Parameters<QueryPileupParams>) -> String {
        let this = self.clone();
        run_blocking_tool("query_pileup", move || this.do_query_pileup(params)).await
    }

    /// Return the README documentation for this server, including usage instructions,
    /// tool descriptions, pileup symbol reference, and client configuration examples.
    #[tool(description = "Return the README documentation for this MCP server")]
    async fn get_documentation(&self) -> String {
        include_str!("../README.md").to_string()
    }

    /// Return the SAM-format header of the BAM file, including @HD, @SQ, @RG, @PG lines.
    #[tool(description = "Return the SAM-format header of the BAM file")]
    async fn get_header(&self) -> String {
        let this = self.clone();
        run_blocking_tool("get_header", move || this.do_get_header()).await
    }
}

impl PileupServer {
    fn do_get_header(&self) -> anyhow::Result<String> {
        use anyhow::Context as _;
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&self.config.bam_path)
            .with_context(|| format!("Failed to open BAM: {}", self.config.bam_path.display()))?;
        let header = reader.read_header().context("Failed to read BAM header")?;
        let mut buf = Vec::new();
        let mut writer = noodles::sam::io::Writer::new(&mut buf);
        writer
            .write_header(&header)
            .context("Failed to write SAM header")?;
        String::from_utf8(buf).context("SAM header is not valid UTF-8")
    }

    fn do_query_pileup(&self, params: QueryPileupParams) -> anyhow::Result<String> {
        let config = &self.config;

        // Resolve per-call overrides against CLI defaults.
        if let Some(w) = params.window
            && !(MIN_WINDOW..=MAX_WINDOW).contains(&w)
        {
            anyhow::bail!(
                "window {w} is out of range: must be between {MIN_WINDOW} and {MAX_WINDOW} bp"
            );
        }
        let window = params.window.unwrap_or(config.window);
        let min_mapq = params.min_mapq.unwrap_or(config.min_mapq);
        let show_strand = params.show_strand.unwrap_or(true);
        let show_mapq = params.show_mapq.unwrap_or(true);

        // Validate the position before touching the reference or BAM. Out-of-range
        // positions would otherwise index outside the fetched region while rendering.
        let contig_len = contig_length(&config.bam_path, &params.chrom)? as i64;
        if params.position < 1 || params.position > contig_len {
            anyhow::bail!(
                "position {} is out of range for {}: must be between 1 and {} (1-based)",
                params.position,
                params.chrom,
                contig_len
            );
        }

        // Convert 1-based input position to 0-based internal coords.
        let query_pos_0 = params.position - 1;
        let region_start = (query_pos_0 - window as i64).max(0);
        // `window` is a half-width: `window` bases each side of the queried base,
        // 2*window+1 in total (fewer only where clipped at a contig end).
        let region_end = (query_pos_0 + window as i64 + 1).min(contig_len);

        debug!(
            chrom = %params.chrom,
            position = params.position,
            window,
            min_mapq,
            "query_pileup called"
        );
        debug!("Region: {}:{}-{}", params.chrom, region_start, region_end);

        // Fetch reference sequence.
        let ref_seq = ReferenceReader::open(&config.reference_path)?.fetch(
            &params.chrom,
            region_start,
            region_end,
        )?;
        debug!(
            "Reference fetched: {}bp starting {:?}",
            ref_seq.len(),
            &ref_seq[..ref_seq.len().min(8)]
        );

        // Query BAM.
        let query_result = query_region(
            &config.bam_path,
            &params.chrom,
            region_start,
            region_end,
            query_pos_0,
            min_mapq,
            config.max_depth as usize,
        )?;
        debug!(
            "BAM query: {} passing reads ({} kept for display); filtered MAPQ={} flags={}",
            query_result.total_passed,
            query_result.reads.len(),
            query_result.filtered_mapq,
            query_result.filtered_flags,
        );

        // Expand CIGAR into reference-coordinate-aligned arrays.
        let t0 = Instant::now();
        let expanded = expand_reads(&query_result.reads, &ref_seq, region_start, region_end);
        debug!(
            "CIGAR expansion: {} reads in {:.1}ms",
            expanded.len(),
            t0.elapsed().as_secs_f64() * 1000.0
        );

        // Build render options.
        let bam_name = config
            .bam_path
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .into_owned();
        let opts = RenderOpts {
            chrom: params.chrom.clone(),
            bam_name,
            show_strand,
            show_mapq,
            min_baseq: config.min_baseq,
        };

        // Render to string.
        let output = render_pileup(
            &expanded,
            &ref_seq,
            region_start,
            query_pos_0,
            &query_result,
            &opts,
        );
        debug!(
            "Render complete: {} lines, {} chars",
            output.lines().count(),
            output.len()
        );

        Ok(output)
    }
}

#[cfg(test)]
mod tests {
    use super::run_blocking_tool;

    #[tokio::test]
    async fn run_blocking_tool_returns_output() {
        let out = run_blocking_tool("t", || Ok("hello".to_string())).await;
        assert_eq!(out, "hello");
    }

    #[tokio::test]
    async fn run_blocking_tool_formats_errors() {
        let out = run_blocking_tool("t", || anyhow::bail!("boom")).await;
        assert_eq!(out, "Error: boom");
    }

    #[tokio::test]
    async fn run_blocking_tool_turns_panic_into_error() {
        let out = run_blocking_tool("t", || -> anyhow::Result<String> {
            panic!("deliberate test panic")
        })
        .await;
        assert!(out.starts_with("Error: internal error in t:"), "{out}");
        assert!(out.contains("deliberate test panic"), "{out}");
    }
}
