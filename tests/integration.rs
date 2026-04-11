//! End-to-end integration tests for the full BAM pileup pipeline.
//!
//! These tests exercise the complete chain:
//!   ReferenceReader::fetch → query_region → expand_reads → render_pileup
//!
//! They require the sample files in `samples/` and are therefore gated on
//! those files existing (they are not `#[ignore]`-d — the CI environment is
//! expected to provide them).

use std::path::PathBuf;

use bam_mcp_server::{
    pileup::expand_reads,
    query::query_region,
    reference::ReferenceReader,
    render::{RenderOpts, render_pileup},
};

fn sample_bam() -> PathBuf {
    PathBuf::from("samples/SampleHuman-30x-WGS.bam")
}

fn sample_fasta() -> PathBuf {
    PathBuf::from("samples/GRCh38.primary_assembly.genome.fa.gz")
}

fn default_opts(chrom: &str) -> RenderOpts {
    RenderOpts {
        chrom: chrom.to_string(),
        bam_name: "SampleHuman-30x-WGS.bam".to_string(),
        show_strand: true,
        show_mapq: true,
        min_baseq: 20,
    }
}

/// Full pipeline against the first chromosome in the BAM at a position with
/// expected 30× WGS coverage.
#[test]
fn test_pipeline_chr1() {
    // Discover the first chromosome name dynamically (handles chr-prefix variants).
    use noodles::bam;
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(sample_bam())
        .expect("open BAM");
    let header = reader.read_header().expect("read header");
    let first_chrom = header
        .reference_sequences()
        .keys()
        .next()
        .expect("at least one reference sequence")
        .to_string();
    let chrom_len = header
        .reference_sequences()
        .get(first_chrom.as_bytes())
        .map(|rs| rs.length().get() as i64)
        .unwrap_or(1_000_000);
    drop(reader);

    // Use a stable well-covered region near position 10,000.
    let query_pos_0: i64 = 10_075_i64.min(chrom_len - 200);
    let window: i64 = 75;
    let region_start = (query_pos_0 - window).max(0);
    let region_end = query_pos_0 + window;

    // 1. Reference.
    let ref_seq = ReferenceReader::open(&sample_fasta())
        .expect("open reference")
        .fetch(&first_chrom, region_start, region_end)
        .expect("fetch reference");
    assert_eq!(
        ref_seq.len(),
        (region_end - region_start) as usize,
        "ref_seq length should match window"
    );

    // 2. BAM query.
    let query_result = query_region(&sample_bam(), &first_chrom, region_start, region_end, 0, 50)
        .expect("query_region should succeed");
    assert!(
        !query_result.reads.is_empty(),
        "expected reads in a 30× WGS region at {}:{}-{}",
        first_chrom,
        region_start,
        region_end
    );

    // 3. CIGAR expansion.
    let expanded = expand_reads(&query_result.reads, &ref_seq, region_start, region_end);
    assert_eq!(
        expanded.len(),
        query_result.reads.len(),
        "one AlignedRead per ReadRecord"
    );
    for ar in &expanded {
        assert_eq!(
            ar.bases.len(),
            ref_seq.len(),
            "bases length must equal region width"
        );
        assert_eq!(
            ar.ins_after.len(),
            ref_seq.len(),
            "ins_after length must equal region width"
        );
    }

    // 4. Render.
    let opts = default_opts(&first_chrom);
    let output = render_pileup(
        &expanded,
        &ref_seq,
        region_start,
        query_pos_0,
        &query_result,
        &opts,
    );

    // Structural assertions on the rendered text.
    assert!(
        output.contains("Region:"),
        "output should contain Region header"
    );
    assert!(output.contains("REF   "), "output should contain REF line");
    assert!(
        output.contains("R01   "),
        "output should contain at least one read line"
    );
    assert!(
        output.contains("COV"),
        "output should contain coverage footer"
    );
    assert!(
        output.contains("Reads shown:"),
        "output should contain filter summary"
    );
    assert!(
        output.contains(&first_chrom),
        "output should mention the queried chromosome"
    );
    assert!(
        output.contains("["),
        "query position should be bracketed in output"
    );
}

/// The pipeline must return a clear error when given a chromosome name that
/// does not exist in the BAM header.
#[test]
fn test_pipeline_invalid_chrom() {
    let err = query_region(&sample_bam(), "nonexistent_chrXYZ99", 0, 100, 0, 50)
        .expect_err("should fail for unknown chromosome");
    let msg = err.to_string();
    assert!(
        msg.contains("Unknown chromosome"),
        "error should mention 'Unknown chromosome', got: {msg}"
    );
}

/// A region with no overlapping reads produces a valid pileup showing 0× coverage.
#[test]
fn test_pipeline_empty_region() {
    // Use a small window at the very start of the reference (position 0–10).
    // At 30× WGS the first few bases are typically uncovered.
    use noodles::bam;
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(sample_bam())
        .expect("open BAM");
    let header = reader.read_header().expect("read header");
    let first_chrom = header
        .reference_sequences()
        .keys()
        .next()
        .expect("at least one reference sequence")
        .to_string();
    drop(reader);

    let region_start: i64 = 0;
    let region_end: i64 = 10;
    let query_pos_0: i64 = 5;

    let ref_seq = ReferenceReader::open(&sample_fasta())
        .expect("open reference")
        .fetch(&first_chrom, region_start, region_end)
        .expect("fetch reference");

    let query_result = query_region(&sample_bam(), &first_chrom, region_start, region_end, 0, 50)
        .expect("query should succeed even with zero reads");

    // Whether reads are present or not, the pipeline must not panic.
    let expanded = expand_reads(&query_result.reads, &ref_seq, region_start, region_end);
    let opts = default_opts(&first_chrom);
    let output = render_pileup(
        &expanded,
        &ref_seq,
        region_start,
        query_pos_0,
        &query_result,
        &opts,
    );

    // Output must always contain the structural sections.
    assert!(output.contains("Region:"), "must always render a header");
    assert!(output.contains("REF   "), "must always render a REF line");
    assert!(output.contains("COV"), "must always render a footer");
    // If no reads, coverage should be 0.
    if query_result.reads.is_empty() {
        assert!(
            output.contains("COV   0x"),
            "zero reads should yield 0× coverage, got:\n{output}"
        );
    }
}
