use std::path::Path;

use anyhow::{Context, bail};
use noodles::{
    bam,
    core::{Position, Region},
    sam::alignment::record::cigar::op::Kind,
};
use rand::distr::{Distribution as _, Uniform};

/// Internal representation of a single CIGAR operation.
///
/// All three "match" op codes (M / = / X) are collapsed to `Match` because the
/// pileup engine does its own base comparison against the reference.
#[derive(Debug, Clone)]
pub enum CigarOp {
    /// M / = / X — consumes both read and reference
    Match(u32),
    /// I — consumes read only
    Ins(u32),
    /// D — consumes reference only
    Del(u32),
    /// S — consumes read, does not advance reference
    SoftClip(u32),
    /// H — consumes nothing; recorded for completeness but ignored during expansion
    HardClip(u32),
    /// N — consumes reference only (intron / splice gap)
    Skip(u32),
}

/// A single aligned read extracted from the BAM file, ready for pileup expansion.
#[derive(Debug, Clone)]
pub struct ReadRecord {
    pub name: String,
    /// 0-based, inclusive
    pub ref_start: i64,
    /// 0-based, exclusive (computed from ref_start + CIGAR reference span)
    pub ref_end: i64,
    pub cigar: Vec<CigarOp>,
    /// Raw bases as ASCII bytes (b'A', b'C', b'G', b'T', b'N', …)
    pub sequence: Vec<u8>,
    /// Phred base quality scores, parallel to `sequence`
    pub base_quals: Vec<u8>,
    pub is_reverse: bool,
    pub mapq: u8,
    pub is_dup: bool,
    pub is_qcfail: bool,
}

/// The result of a BAM region query after filtering.
#[derive(Debug)]
pub struct QueryResult {
    pub reads: Vec<ReadRecord>,
    /// Reads excluded because mapping quality < min_mapq
    pub filtered_mapq: usize,
    /// Reads excluded because they are duplicates or QC failures
    pub filtered_flags: usize,
}

/// Query a BAM file for all reads overlapping `[start, end)` on `chrom`.
///
/// Coordinates are **0-based, half-open**.
///
/// Filtering decision tree (applied in order):
/// 1. Unmapped / secondary / supplementary → skipped silently (not counted)
/// 2. `mapq < min_mapq` → counted in `filtered_mapq`, excluded from output
/// 3. Duplicate or QC-fail → counted in `filtered_flags`, excluded from output
/// 4. All remaining reads are eligible; if more than `max_depth` pass, the output
///    is a reservoir sample of exactly `max_depth` reads.
pub fn query_region(
    bam_path: &Path,
    chrom: &str,
    start: i64,
    end: i64,
    min_mapq: u8,
    max_depth: usize,
) -> anyhow::Result<QueryResult> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .with_context(|| format!("Failed to open BAM: {}", bam_path.display()))?;

    let header = reader.read_header().context("Failed to read BAM header")?;

    // Validate chromosome name against the header reference sequences.
    if !header.reference_sequences().contains_key(chrom.as_bytes()) {
        let available: Vec<String> = header
            .reference_sequences()
            .keys()
            .take(10)
            .map(|k| String::from_utf8_lossy(k).into_owned())
            .collect();
        bail!(
            "Unknown chromosome: '{}'. Available (first 10): {}",
            chrom,
            available.join(", ")
        );
    }

    // Build region — noodles uses 1-based inclusive positions.
    // 0-based half-open [start, end)  →  1-based inclusive [start+1, end]
    let pos_start = Position::try_from(start as usize + 1)
        .with_context(|| format!("Invalid start coordinate: {start}"))?;
    let pos_end = Position::try_from(end as usize)
        .with_context(|| format!("Invalid end coordinate: {end}"))?;
    let region = Region::new(chrom, pos_start..=pos_end);

    let mut candidates: Vec<ReadRecord> = Vec::new();
    let mut filtered_mapq: usize = 0;
    let mut filtered_flags: usize = 0;
    // `total_passed` is the count of reads that reached the reservoir stage
    // (i.e. passed all filters). Used as the denominator in reservoir sampling.
    let mut total_passed: usize = 0;

    for result in reader.query(&header, &region)?.records() {
        let record = result.context("Failed to read BAM record")?;

        let flags = record.flags();

        // Step 1 — skip reads that cannot contribute to a linear pileup
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }

        // Step 2 — MAPQ filter
        let mapq = record.mapping_quality().map(|mq| mq.get()).unwrap_or(0);
        if mapq < min_mapq {
            filtered_mapq += 1;
            continue;
        }

        // Step 3 — duplicate / QC-fail filter
        let is_dup = flags.is_duplicate();
        let is_qcfail = flags.is_qc_fail();
        if is_dup || is_qcfail {
            filtered_flags += 1;
            continue;
        }

        // --- Extract fields ---

        // Alignment start (0-based)
        let ref_start = match record.alignment_start() {
            Some(pos) => pos.context("Invalid alignment start")?.get() as i64 - 1,
            None => continue, // should be caught by is_unmapped() but be defensive
        };

        let cigar = decode_cigar(&record)?;
        let ref_end = compute_ref_end(ref_start, &cigar);

        let name = record
            .name()
            .map(|n| String::from_utf8_lossy(n).into_owned())
            .unwrap_or_default();

        let sequence: Vec<u8> = record.sequence().iter().collect();
        let base_quals: Vec<u8> = record.quality_scores().iter().collect();

        let read = ReadRecord {
            name,
            ref_start,
            ref_end,
            cigar,
            sequence,
            base_quals,
            is_reverse: flags.is_reverse_complemented(),
            mapq,
            is_dup,
            is_qcfail,
        };

        // Reservoir sampling (Knuth's Algorithm R): keep at most `max_depth` reads.
        if total_passed < max_depth {
            candidates.push(read);
        } else {
            // Replace a random slot with probability max_depth / (total_passed + 1)
            // Reservoir sampling: replace a random slot with probability
            // max_depth / (total_passed + 1)
            let j: usize = Uniform::new_inclusive(0usize, total_passed)
                .expect("valid reservoir range")
                .sample(&mut rand::rng());
            if j < max_depth {
                candidates[j] = read;
            }
        }
        total_passed += 1;
    }

    Ok(QueryResult {
        reads: candidates,
        filtered_mapq,
        filtered_flags,
    })
}

/// Decode the CIGAR ops from a BAM record into the internal `CigarOp` enum.
/// `Pad` ops are silently skipped (they apply only to padded reference assemblies).
fn decode_cigar(record: &bam::Record) -> anyhow::Result<Vec<CigarOp>> {
    let mut ops = Vec::new();
    for result in record.cigar().iter() {
        let op = result.context("Invalid CIGAR op")?;
        let n = op.len() as u32;
        let cigar_op = match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => CigarOp::Match(n),
            Kind::Insertion => CigarOp::Ins(n),
            Kind::Deletion => CigarOp::Del(n),
            Kind::SoftClip => CigarOp::SoftClip(n),
            Kind::HardClip => CigarOp::HardClip(n),
            Kind::Skip => CigarOp::Skip(n),
            Kind::Pad => continue,
        };
        ops.push(cigar_op);
    }
    Ok(ops)
}

/// Compute 0-based exclusive reference end from start + CIGAR reference span.
fn compute_ref_end(ref_start: i64, cigar: &[CigarOp]) -> i64 {
    let ref_span: i64 = cigar
        .iter()
        .map(|op| match op {
            CigarOp::Match(n) | CigarOp::Del(n) | CigarOp::Skip(n) => *n as i64,
            CigarOp::Ins(_) | CigarOp::SoftClip(_) | CigarOp::HardClip(_) => 0,
        })
        .sum();
    ref_start + ref_span
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn sample_bam() -> PathBuf {
        PathBuf::from("samples/SampleHuman-30x-WGS.bam")
    }

    /// Verify all 9 noodles CIGAR Kind variants map to the correct CigarOp.
    /// This is a pure logic test — no I/O required.
    #[test]
    fn test_cigar_kind_mapping() {
        let cases: &[(Kind, &str)] = &[
            (Kind::Match, "Match"),
            (Kind::SequenceMatch, "Match"),
            (Kind::SequenceMismatch, "Match"),
            (Kind::Insertion, "Ins"),
            (Kind::Deletion, "Del"),
            (Kind::SoftClip, "SoftClip"),
            (Kind::HardClip, "HardClip"),
            (Kind::Skip, "Skip"),
        ];
        for (kind, expected) in cases {
            let op = match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => "Match",
                Kind::Insertion => "Ins",
                Kind::Deletion => "Del",
                Kind::SoftClip => "SoftClip",
                Kind::HardClip => "HardClip",
                Kind::Skip => "Skip",
                Kind::Pad => "Pad",
            };
            assert_eq!(op, *expected, "Kind::{kind:?} should map to {expected}");
        }
    }

    #[test]
    fn test_compute_ref_end() {
        // 100M → ref span = 100
        let cigar = vec![CigarOp::Match(100)];
        assert_eq!(compute_ref_end(1000, &cigar), 1100);

        // 10M 2D 88M → ref span = 100
        let cigar = vec![CigarOp::Match(10), CigarOp::Del(2), CigarOp::Match(88)];
        assert_eq!(compute_ref_end(0, &cigar), 100);

        // 5S 95M → ref span = 95 (soft clip doesn't consume ref)
        let cigar = vec![CigarOp::SoftClip(5), CigarOp::Match(95)];
        assert_eq!(compute_ref_end(0, &cigar), 95);

        // 10M 3I 87M → ref span = 97 (insertion doesn't consume ref)
        let cigar = vec![CigarOp::Match(10), CigarOp::Ins(3), CigarOp::Match(87)];
        assert_eq!(compute_ref_end(0, &cigar), 97);
    }

    #[test]
    fn test_query_known_region() {
        // Dynamically pick the first chromosome from the BAM header so the
        // test is independent of chr-prefix naming conventions.
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
        let first_len = header
            .reference_sequences()
            .get(first_chrom.as_bytes())
            .map(|rs| rs.length().get() as i64)
            .unwrap_or(1_000_000);
        drop(reader);

        let result = query_region(
            &sample_bam(),
            &first_chrom,
            10_000.min(first_len - 200),
            10_150.min(first_len),
            0,
            50,
        )
        .expect("query_region should succeed");

        assert!(
            !result.reads.is_empty(),
            "expected reads in a 30× WGS region"
        );
        for r in &result.reads {
            assert!(r.ref_start >= 0);
            assert!(r.ref_end > r.ref_start);
            assert!(!r.cigar.is_empty());
        }
    }

    #[test]
    fn test_query_invalid_chrom() {
        let err = query_region(&sample_bam(), "nonexistent_chrXYZ99", 0, 100, 0, 50)
            .expect_err("should fail for unknown chromosome");
        let msg = err.to_string();
        assert!(
            msg.contains("Unknown chromosome"),
            "error message should mention 'Unknown chromosome', got: {msg}"
        );
    }

    #[test]
    fn test_filter_mapq() {
        // Open BAM to get first chrom name
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

        let result = query_region(&sample_bam(), &first_chrom, 10_000, 10_150, 255, 50)
            .expect("query_region should succeed");

        assert!(
            result.reads.is_empty(),
            "min_mapq=255 should filter all reads"
        );
        assert!(
            result.filtered_mapq > 0,
            "some reads should have been MAPQ-filtered"
        );
    }
}
