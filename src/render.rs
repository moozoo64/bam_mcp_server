use crate::pileup::{AlignedRead, PileupBase};
use crate::query::QueryResult;

/// Options that control how the pileup is rendered.
#[derive(Debug, Clone)]
pub struct RenderOpts {
    /// Chromosome / contig name (for the header line).
    pub chrom: String,
    /// Basename of the BAM file (for the header line).
    pub bam_name: String,
    /// Whether to show strand orientation in the per-read suffix.
    pub show_strand: bool,
    /// Whether to show mapping quality in the per-read suffix.
    pub show_mapq: bool,
    /// Minimum base quality; mismatches below this are shown in lowercase.
    pub min_baseq: u8,
}

// ── Section renderers ────────────────────────────────────────────────────────

fn render_header(opts: &RenderOpts, region_start: i64, region_end: i64) -> String {
    let len = region_end - region_start;
    format!(
        "Region: {}:{}-{}  ({}bp window)  BAM: {}",
        opts.chrom,
        format_with_commas(region_start),
        format_with_commas(region_end),
        len,
        opts.bam_name,
    )
}

fn render_pos_ruler(region_start: i64, region_end: i64, query_offset: usize) -> String {
    let region_len = (region_end - region_start) as usize;
    // Content width = region_len + 2 (room for the `[X]` bracket at query_offset)
    let content_width = region_len + 2;

    let start_str = format!("{}", region_start);
    let end_str = format!("{}", region_end);

    // Mark the query position with `[*]` on the POS ruler line.
    // Build a character array for the content.
    let mut chars: Vec<u8> = vec![b' '; content_width];

    // Place start label at position 0 (accounting for bracket shift).
    for (i, b) in start_str.bytes().enumerate() {
        if i < content_width {
            chars[i] = b;
        }
    }

    // Place `[*]` at the query offset (shifted by 1 if query_offset > 0 due to bracket).
    // The content column for ref index i is i when i < query_offset, or i+2 when i > query_offset,
    // and i..i+2 when i == query_offset.
    let q_col = query_offset; // bracket occupies q_col, q_col+1, q_col+2
    if q_col < content_width.saturating_sub(2) {
        chars[q_col] = b'[';
        chars[q_col + 1] = b'*';
        chars[q_col + 2] = b']';
    }

    // Right-align end label at the far right of the content.
    let end_bytes = end_str.as_bytes();
    let end_start = content_width.saturating_sub(end_bytes.len());
    for (i, b) in end_bytes.iter().enumerate() {
        let col = end_start + i;
        if col < content_width {
            chars[col] = *b;
        }
    }

    format!("POS   {}", String::from_utf8_lossy(&chars))
}

fn render_ref_line(ref_seq: &[u8], query_offset: usize) -> String {
    let mut out = String::with_capacity(ref_seq.len() + 8);
    out.push_str("REF   ");
    for (i, &b) in ref_seq.iter().enumerate() {
        if i == query_offset {
            out.push('[');
            out.push(b as char);
            out.push(']');
        } else {
            out.push(b as char);
        }
    }
    out
}

fn render_read_line(
    idx: usize,
    ar: &AlignedRead,
    query_offset: usize,
    opts: &RenderOpts,
) -> String {
    let region_len = ar.bases.len();
    let mut out = String::with_capacity(region_len + 20);

    // Name prefix: R01, R02, … padded to 6 chars
    out.push_str(&format!("R{:02}   ", idx + 1));

    for i in 0..region_len {
        // Build the symbol(s) for this slot.
        let sym = match &ar.bases[i] {
            None => " ".to_string(),
            Some(PileupBase::Match) => ".".to_string(),
            Some(PileupBase::Deletion) => "-".to_string(),
            Some(PileupBase::Mismatch { base, qual }) => {
                let ch = if *qual >= opts.min_baseq {
                    (*base as char).to_ascii_uppercase()
                } else {
                    (*base as char).to_ascii_lowercase()
                };
                ch.to_string()
            }
        };

        // Append insertion annotation if any.
        let ins = ar.ins_after[i];
        let slot = if ins > 0 {
            format!("{}^{}", sym, ins)
        } else {
            sym
        };

        // Wrap in brackets at the query position.
        if i == query_offset {
            out.push('[');
            out.push_str(&slot);
            out.push(']');
        } else {
            out.push_str(&slot);
        }
    }

    // Suffix.
    if opts.show_strand || opts.show_mapq {
        let strand = if ar.record.is_reverse { '-' } else { '+' };
        if opts.show_strand && opts.show_mapq {
            out.push_str(&format!("  {} MQ{}", strand, ar.record.mapq));
        } else if opts.show_strand {
            out.push_str(&format!("  {}", strand));
        } else {
            out.push_str(&format!("  MQ{}", ar.record.mapq));
        }
    }

    out
}

fn render_footer(
    reads: &[AlignedRead],
    ref_seq: &[u8],
    query_offset: usize,
    query_pos: i64,
    query_result: &QueryResult,
    opts: &RenderOpts,
) -> String {
    // Coverage = reads that have a base at query_offset.
    let covered: Vec<&AlignedRead> = reads
        .iter()
        .filter(|ar| ar.bases[query_offset].is_some())
        .collect();
    let coverage = covered.len();

    // Allele counts at query position.
    let mut allele_a = 0usize;
    let mut allele_c = 0usize;
    let mut allele_g = 0usize;
    let mut allele_t = 0usize;
    let mut allele_del = 0usize;

    let ref_base = ref_seq[query_offset].to_ascii_uppercase();

    for ar in &covered {
        match &ar.bases[query_offset] {
            Some(PileupBase::Match) => {
                // Tally the reference base.
                match ref_base {
                    b'A' => allele_a += 1,
                    b'C' => allele_c += 1,
                    b'G' => allele_g += 1,
                    b'T' => allele_t += 1,
                    _ => {}
                }
            }
            Some(PileupBase::Mismatch { base, .. }) => match base.to_ascii_uppercase() {
                b'A' => allele_a += 1,
                b'C' => allele_c += 1,
                b'G' => allele_g += 1,
                b'T' => allele_t += 1,
                _ => {}
            },
            Some(PileupBase::Deletion) => allele_del += 1,
            None => {}
        }
    }

    let mut allele_lines: Vec<String> = Vec::new();
    for (label, count) in [
        ("A", allele_a),
        ("C", allele_c),
        ("G", allele_g),
        ("T", allele_t),
        ("del", allele_del),
    ] {
        if count > 0 {
            let pct = 100.0 * count as f64 / coverage.max(1) as f64;
            allele_lines.push(format!("{}: {} ({:.1}%)", label, count, pct));
        }
    }

    // Strand balance.
    let fwd = covered.iter().filter(|ar| !ar.record.is_reverse).count();
    let rev = covered.len() - fwd;

    // MQ stats across all reads returned by the query.
    let (mq_mean, mq_min) = if query_result.reads.is_empty() {
        (0.0f64, 0u8)
    } else {
        let sum: usize = query_result.reads.iter().map(|r| r.mapq as usize).sum();
        let mean = sum as f64 / query_result.reads.len() as f64;
        let min = query_result.reads.iter().map(|r| r.mapq).min().unwrap_or(0);
        (mean, min)
    };

    let mut lines: Vec<String> = Vec::new();
    lines.push(format!(
        "COV   {}x at {}:{}",
        coverage,
        opts.chrom,
        format_with_commas(query_pos),
    ));
    if !allele_lines.is_empty() {
        lines.push(format!("      {}", allele_lines.join("  ")));
    }
    lines.push(format!("      +strand: {}  -strand: {}", fwd, rev));
    lines.push(format!("      MQ_mean: {:.1}  MQ_min: {}", mq_mean, mq_min));
    lines.push(format!(
        "      Reads shown: {}  Reads filtered (MAPQ): {}  Reads filtered (dup/qcfail): {}",
        query_result.reads.len(),
        query_result.filtered_mapq,
        query_result.filtered_flags,
    ));

    lines.join("\n")
}

// ── Public entry point ───────────────────────────────────────────────────────

/// Render a full pileup text block.
///
/// * `reads`        — CIGAR-expanded reads from `expand_reads()`
/// * `ref_seq`      — uppercase ASCII bytes; length == `region_end - region_start`
/// * `region_start` — 0-based inclusive left bound of the window
/// * `query_pos`    — 0-based position the user queried (centre of window)
/// * `query_result` — raw query result for filter counts and MQ stats
/// * `opts`         — display options
pub fn render_pileup(
    reads: &[AlignedRead],
    ref_seq: &[u8],
    region_start: i64,
    query_pos: i64,
    query_result: &QueryResult,
    opts: &RenderOpts,
) -> String {
    let region_end = region_start + ref_seq.len() as i64;
    let query_offset = (query_pos - region_start) as usize;

    let header = render_header(opts, region_start, region_end);
    let ruler = render_pos_ruler(region_start, region_end, query_offset);
    let ref_line = render_ref_line(ref_seq, query_offset);
    let read_lines: Vec<String> = reads
        .iter()
        .enumerate()
        .map(|(i, ar)| render_read_line(i, ar, query_offset, opts))
        .collect();
    let footer = render_footer(reads, ref_seq, query_offset, query_pos, query_result, opts);

    let mut parts: Vec<String> = Vec::new();
    parts.push(header);
    parts.push(String::new()); // blank line
    parts.push(ruler);
    parts.push(ref_line);
    parts.push(String::new()); // blank line
    parts.extend(read_lines);
    parts.push(String::new()); // blank line
    parts.push(footer);

    parts.join("\n")
}

// ── Utilities ────────────────────────────────────────────────────────────────

fn format_with_commas(n: i64) -> String {
    let s = n.abs().to_string();
    let mut out = String::new();
    for (i, ch) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }
    let mut result: String = out.chars().rev().collect();
    if n < 0 {
        result.insert(0, '-');
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::{AlignedRead, PileupBase};
    use crate::query::{QueryResult, ReadRecord};

    fn dummy_record(is_reverse: bool, mapq: u8) -> ReadRecord {
        ReadRecord {
            name: "r".to_string(),
            ref_start: 0,
            ref_end: 0,
            cigar: vec![],
            sequence: vec![],
            base_quals: vec![],
            is_reverse,
            mapq,
            is_dup: false,
            is_qcfail: false,
        }
    }

    fn make_ar(bases: Vec<Option<PileupBase>>, is_reverse: bool, mapq: u8) -> AlignedRead {
        let len = bases.len();
        AlignedRead {
            record: dummy_record(is_reverse, mapq),
            bases,
            ins_after: vec![0u32; len],
        }
    }

    fn default_opts() -> RenderOpts {
        RenderOpts {
            chrom: "chr1".to_string(),
            bam_name: "test.bam".to_string(),
            show_strand: true,
            show_mapq: true,
            min_baseq: 20,
        }
    }

    fn dummy_query_result(reads: Vec<ReadRecord>) -> QueryResult {
        QueryResult {
            reads,
            filtered_mapq: 0,
            filtered_flags: 0,
        }
    }

    // ── tests ────────────────────────────────────────────────────────────────

    #[test]
    fn test_basic_render() {
        // 5-base window, query at index 2. One matching read.
        let ref_seq = b"ACGTA";
        let region_start = 10i64;
        let query_pos = 12i64; // offset 2
        let ar = make_ar(
            vec![
                Some(PileupBase::Match),
                Some(PileupBase::Match),
                Some(PileupBase::Match),
                Some(PileupBase::Match),
                Some(PileupBase::Match),
            ],
            false,
            60,
        );
        let qr = dummy_query_result(vec![dummy_record(false, 60)]);
        let output = render_pileup(
            &[ar],
            ref_seq,
            region_start,
            query_pos,
            &qr,
            &default_opts(),
        );

        // REF line must contain [G] at the query position (ref_seq[2] = 'G')
        assert!(
            output.contains("[G]"),
            "REF line should contain [G]: {}",
            output
        );
        // Read line must contain [.] at the query position
        assert!(
            output.contains("[.]"),
            "Read line should contain [.]: {}",
            output
        );
    }

    #[test]
    fn test_mismatch_high_qual() {
        let ref_seq = b"ACGTA";
        let ar = make_ar(
            vec![
                Some(PileupBase::Mismatch {
                    base: b'T',
                    qual: 30,
                }),
                Some(PileupBase::Match),
                Some(PileupBase::Match),
                Some(PileupBase::Match),
                Some(PileupBase::Match),
            ],
            false,
            60,
        );
        let qr = dummy_query_result(vec![dummy_record(false, 60)]);
        let output = render_pileup(&[ar], ref_seq, 10, 10, &qr, &default_opts());
        // query_offset=0, so position 0 is bracketed: [T]
        assert!(
            output.contains("[T]"),
            "High-qual mismatch should be uppercase: {}",
            output
        );
    }

    #[test]
    fn test_mismatch_low_qual() {
        let ref_seq = b"ACGTA";
        let ar = make_ar(
            vec![
                Some(PileupBase::Mismatch {
                    base: b'T',
                    qual: 10,
                }), // below min_baseq=20
                Some(PileupBase::Match),
                Some(PileupBase::Match),
                Some(PileupBase::Match),
                Some(PileupBase::Match),
            ],
            false,
            60,
        );
        let qr = dummy_query_result(vec![dummy_record(false, 60)]);
        let output = render_pileup(&[ar], ref_seq, 10, 10, &qr, &default_opts());
        // Should contain lowercase [t]
        assert!(
            output.contains("[t]"),
            "Low-qual mismatch should be lowercase: {}",
            output
        );
    }

    #[test]
    fn test_deletion_render() {
        let ref_seq = b"ACGTA";
        let ar = make_ar(
            vec![
                Some(PileupBase::Match),
                Some(PileupBase::Match),
                Some(PileupBase::Deletion),
                Some(PileupBase::Match),
                Some(PileupBase::Match),
            ],
            false,
            60,
        );
        let qr = dummy_query_result(vec![dummy_record(false, 60)]);
        // query at offset 2 (the deletion)
        let output = render_pileup(&[ar], ref_seq, 10, 12, &qr, &default_opts());
        assert!(
            output.contains("[-]"),
            "Deletion at query pos should be [-]: {}",
            output
        );
    }

    #[test]
    fn test_insertion_render() {
        let ref_seq = b"ACGTA";
        let len = ref_seq.len();
        // Insertion of 3 after position 1
        let mut ar = make_ar(vec![Some(PileupBase::Match); len], false, 60);
        ar.ins_after[1] = 3;
        let qr = dummy_query_result(vec![dummy_record(false, 60)]);
        let output = render_pileup(&[ar], ref_seq, 10, 10, &qr, &default_opts());
        // Position 1 is not the query (query=offset 0), so plain .^3
        assert!(
            output.contains(".^3"),
            "Insertion annotation should appear as .^3: {}",
            output
        );
    }

    #[test]
    fn test_partial_overlap_spaces() {
        let ref_seq = b"ACGTA";
        // Read only covers positions 2..5 (None at 0 and 1)
        let ar = make_ar(
            vec![
                None,
                None,
                Some(PileupBase::Match),
                Some(PileupBase::Match),
                Some(PileupBase::Match),
            ],
            false,
            60,
        );
        let qr = dummy_query_result(vec![dummy_record(false, 60)]);
        let output = render_pileup(&[ar], ref_seq, 10, 12, &qr, &default_opts());
        // read line should start with "R01   " followed by two spaces then [.]...
        let read_line = output.lines().find(|l| l.starts_with("R01")).unwrap();
        // The first two content chars after "R01   " should be spaces
        let content = &read_line["R01   ".len()..];
        assert_eq!(
            &content[..2],
            "  ",
            "Un-covered positions should be spaces: {:?}",
            content
        );
    }

    #[test]
    fn test_footer_allele_counts() {
        let ref_seq = b"ACGTA"; // ref at offset 2 = G
        let region_start = 10i64;
        let query_pos = 12i64; // offset 2

        // 2 reads with G (match), 1 read with T (mismatch high qual)
        let ar1 = make_ar(vec![Some(PileupBase::Match); 5], false, 60);
        let ar2 = make_ar(vec![Some(PileupBase::Match); 5], false, 60);
        let mut bases3 = vec![Some(PileupBase::Match); 5];
        bases3[2] = Some(PileupBase::Mismatch {
            base: b'T',
            qual: 30,
        });
        let ar3 = make_ar(bases3, true, 55);

        let qr = dummy_query_result(vec![
            dummy_record(false, 60),
            dummy_record(false, 60),
            dummy_record(true, 55),
        ]);
        let output = render_pileup(
            &[ar1, ar2, ar3],
            ref_seq,
            region_start,
            query_pos,
            &qr,
            &default_opts(),
        );

        // Footer should say 3x coverage
        assert!(
            output.contains("COV   3x"),
            "Expected 3x coverage: {}",
            output
        );
        // G: 2 reads (66.7%), T: 1 read (33.3%)
        assert!(
            output.contains("G: 2"),
            "Expected G: 2 in footer: {}",
            output
        );
        assert!(
            output.contains("T: 1"),
            "Expected T: 1 in footer: {}",
            output
        );
    }
}
