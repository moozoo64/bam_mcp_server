use std::sync::Arc;
use std::time::Instant;

use rmcp::handler::server::wrapper::Parameters;
use rmcp::{tool, tool_router};
use schemars::JsonSchema;
use serde::Deserialize;
use tracing::debug;

use crate::cli::AppConfig;
use crate::pileup::expand_reads;
use crate::query::query_region;
use crate::reference::ReferenceReader;
use crate::render::{RenderOpts, render_pileup};

/// Input parameters for the `query_pileup` MCP tool.
#[derive(Debug, Deserialize, JsonSchema)]
pub struct QueryPileupParams {
    /// Chromosome name exactly as it appears in the BAM header (e.g. chr1, 1, chrM)
    pub chrom: String,

    /// 1-based genomic position to query
    pub position: i64,

    /// Half-width of region around position in bp (10–500). Overrides --window CLI default.
    pub window: Option<u32>,

    /// Minimum mapping quality (0–60). Overrides --min-mapq CLI default.
    pub min_mapq: Option<u8>,

    /// Show strand (+/-) column per read
    pub show_strand: Option<bool>,

    /// Show mapping quality column per read
    pub show_mapq: Option<bool>,
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
        match self.do_query_pileup(params).await {
            Ok(output) => output,
            Err(e) => {
                tracing::error!("{:#}", e);
                format!("Error: {:#}", e)
            }
        }
    }
}

impl PileupServer {
    async fn do_query_pileup(&self, params: QueryPileupParams) -> anyhow::Result<String> {
        let config = &self.config;

        // Resolve per-call overrides against CLI defaults.
        let window = params.window.unwrap_or(config.window);
        let min_mapq = params.min_mapq.unwrap_or(config.min_mapq);
        let show_strand = params.show_strand.unwrap_or(true);
        let show_mapq = params.show_mapq.unwrap_or(true);

        // Convert 1-based input position to 0-based internal coords.
        let query_pos_0 = params.position - 1;
        let region_start = (query_pos_0 - window as i64).max(0);
        let region_end = query_pos_0 + window as i64;

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
            min_mapq,
            config.max_depth as usize,
        )?;
        debug!(
            "BAM query returned {} records; filtered MAPQ={} flags={}",
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
