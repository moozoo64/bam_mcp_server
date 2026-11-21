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

// Coordinates are passed around 0-based half-open, but every label in the
// output is 1-based (inclusive), matching the `position` input, VCF, IGV and
// samtools. For a 0-based half-open `[start, end)` that is `start+1 ..= end`.

fn render_header(opts: &RenderOpts, region_start: i64, region_end: i64) -> String {
    let len = region_end - region_start;
    format!(
        "Region: {}:{}-{}  ({}bp window)  BAM: {}",
        opts.chrom,
        format_with_commas(region_start + 1),
        format_with_commas(region_end),
        len,
        opts.bam_name,
    )
}

fn render_pos_ruler(region_start: i64, region_end: i64, query_offset: usize) -> String {
    let region_len = (region_end - region_start) as usize;
    // Content width = region_len + 2 (room for the `[X]` bracket at query_offset)
    let content_width = region_len + 2;

    let start_str = format!("{}", region_start + 1);
    let end_str = format!("{}", region_end);

    let mut chars: Vec<u8> = vec![b' '; content_width];

    // `[*]` marks the query column and always wins. Ref index i sits at column
    // i before the query and i+2 after it, so the bracket spans
    // query_offset..=query_offset+2.
    let q_col = query_offset;
    let (q_first, q_last) = (q_col, q_col + 2);
    if q_last < content_width {
        chars[q_first..=q_last].copy_from_slice(b"[*]");
    }

    // Start label, left-aligned, only if it leaves a space before `[*]`.
    // Otherwise it is omitted rather than overwritten; the Region line always
    // shows the full range.
    let start_fits = start_str.len() < q_first;
    if start_fits {
        chars[..start_str.len()].copy_from_slice(start_str.as_bytes());
    }

    // End label, right-aligned, only if it leaves a space after `[*]` and
    // after the start label.
    let occupied_until = if start_fits {
        q_last.max(start_str.len() - 1)
    } else {
        q_last
    };
    if let Some(end_col) = content_width.checked_sub(end_str.len())
        && end_col > occupied_until + 1
    {
        chars[end_col..].copy_from_slice(end_str.as_bytes());
    }

    format!("POS   {}", String::from_utf8_lossy(&chars).trim_end())
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
            Some(PileupBase::Match { .. }) => ".".to_string(),
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

fn render_bq_section(
    reads: &[AlignedRead],
    query_offset: usize,
    query_pos: i64,
    opts: &RenderOpts,
) -> String {
    // Collect (read_idx, is_reverse, bq) for reads with a base call at query_offset.
    let mut entries: Vec<(usize, bool, u8)> = reads
        .iter()
        .enumerate()
        .filter_map(|(i, ar)| match &ar.bases[query_offset] {
            Some(PileupBase::Match { qual }) => Some((i, ar.record.is_reverse, *qual)),
            Some(PileupBase::Mismatch { qual, .. }) => Some((i, ar.record.is_reverse, *qual)),
            _ => None,
        })
        .collect();

    // Sort: reverse-strand (-) first, then forward (+); stable preserves read-index order within group.
    entries.sort_by_key(|&(_, is_rev, _)| if is_rev { 0u8 } else { 1u8 });

    let header = format!(
        "BQ    {}:{}  (shown reads spanning query position, sorted by strand)",
        opts.chrom,
        format_with_commas(query_pos + 1),
    );

    if entries.is_empty() {
        return format!("{}\n\n      (no reads)", header);
    }

    let bq_values: Vec<u8> = entries.iter().map(|&(_, _, bq)| bq).collect();
    let mean = bq_values.iter().map(|&b| b as f64).sum::<f64>() / bq_values.len() as f64;
    let min = *bq_values.iter().min().unwrap();
    let max = *bq_values.iter().max().unwrap();

    let col_header = "      Read    Strand   BQ-at-query";
    let rows: Vec<String> = entries
        .iter()
        .map(|&(idx, is_rev, bq)| {
            let name = format!("R{:02}", idx + 1);
            let strand = if is_rev { '-' } else { '+' };
            format!("      {:<8}{:<9}{}", name, strand, bq)
        })
        .collect();

    let stats = format!("      mean: {:.1}   min: {}   max: {}", mean, min, max);

    let mut lines = vec![header, String::new(), col_header.to_string()];
    lines.extend(rows);
    lines.push(stats);
    lines.join("\n")
}

fn render_footer(
    reads: &[AlignedRead],
    query_offset: usize,
    query_pos: i64,
    query_result: &QueryResult,
    opts: &RenderOpts,
) -> String {
    // Depth, alleles and strand counts come from every read that passed the
    // filters, not just the (possibly sampled) reads drawn above.
    let site = &query_result.site;

    let mut allele_lines: Vec<String> = Vec::new();
    for (label, count) in [
        ("A", site.a),
        ("C", site.c),
        ("G", site.g),
        ("T", site.t),
        ("del", site.del),
    ] {
        if count > 0 {
            let pct = 100.0 * count as f64 / site.depth.max(1) as f64;
            allele_lines.push(format!("{}: {} ({:.1}%)", label, count, pct));
        }
    }

    let shown = query_result.reads.len();
    let reads_shown = if shown < query_result.total_passed {
        format!(
            "Reads shown: {} of {} passing filters (random sample, capped by --max-depth)",
            shown, query_result.total_passed
        )
    } else {
        format!("Reads shown: {}", shown)
    };

    let mut lines: Vec<String> = Vec::new();
    lines.push(format!(
        "COV   {}x at {}:{}",
        site.depth,
        opts.chrom,
        format_with_commas(query_pos + 1),
    ));
    if !allele_lines.is_empty() {
        lines.push(format!("      {}", allele_lines.join("  ")));
    }
    for (kind, prefix, count) in [
        ("ins after", '+', site.ins_after),
        ("del after", '-', site.del_after),
    ] {
        if count > 0 {
            lines.push(format!(
                "      {}: {} ({:.1}%)  {}",
                kind,
                count,
                100.0 * count as f64 / site.depth.max(1) as f64,
                indel_breakdown(&site.indels, prefix),
            ));
        }
    }
    lines.push(format!(
        "      +strand: {}  -strand: {}",
        site.fwd, site.rev
    ));
    lines.push(format!(
        "      MQ_mean: {:.1}  MQ_min: {}",
        query_result.mapq_mean(),
        query_result.mapq_min.unwrap_or(0)
    ));
    lines.push(format!(
        "      {}  Reads filtered (MAPQ): {}  Reads filtered (dup/qcfail): {}",
        reads_shown, query_result.filtered_mapq, query_result.filtered_flags,
    ));

    lines.push(String::new());
    lines.push(render_bq_section(reads, query_offset, query_pos, opts));

    lines.join("\n")
}

/// Format the anchored indels whose label starts with `prefix` (`+` or `-`),
/// most common first, e.g. `+A: 24  +AA: 1`. Shows at most 4 and sums the rest.
fn indel_breakdown(indels: &std::collections::BTreeMap<String, usize>, prefix: char) -> String {
    let mut entries: Vec<(&String, &usize)> = indels
        .iter()
        .filter(|(label, _)| label.starts_with(prefix))
        .collect();
    // Most common first; BTreeMap order breaks ties deterministically.
    entries.sort_by(|a, b| b.1.cmp(a.1));
    const MAX_SHOWN: usize = 4;
    let mut parts: Vec<String> = entries
        .iter()
        .take(MAX_SHOWN)
        .map(|(label, n)| format!("{label}: {n}"))
        .collect();
    let rest: usize = entries.iter().skip(MAX_SHOWN).map(|(_, n)| **n).sum();
    if rest > 0 {
        parts.push(format!("other: {rest}"));
    }
    parts.join("  ")
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
    let footer = render_footer(reads, query_offset, query_pos, query_result, opts);

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
    use crate::query::{QueryResult, ReadRecord, SiteStats};

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
            total_passed: reads.len(),
            mapq_sum: reads.iter().map(|r| r.mapq as u64).sum(),
            mapq_min: reads.iter().map(|r| r.mapq).min(),
            reads,
            ..Default::default()
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
                Some(PileupBase::Match { qual: 30 }),
                Some(PileupBase::Match { qual: 30 }),
                Some(PileupBase::Match { qual: 30 }),
                Some(PileupBase::Match { qual: 30 }),
                Some(PileupBase::Match { qual: 30 }),
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
                Some(PileupBase::Match { qual: 30 }),
                Some(PileupBase::Match { qual: 30 }),
                Some(PileupBase::Match { qual: 30 }),
                Some(PileupBase::Match { qual: 30 }),
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
                Some(PileupBase::Match { qual: 30 }),
                Some(PileupBase::Match { qual: 30 }),
                Some(PileupBase::Match { qual: 30 }),
                Some(PileupBase::Match { qual: 30 }),
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
                Some(PileupBase::Match { qual: 30 }),
                Some(PileupBase::Match { qual: 30 }),
                Some(PileupBase::Deletion),
                Some(PileupBase::Match { qual: 30 }),
                Some(PileupBase::Match { qual: 30 }),
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
        let mut ar = make_ar(vec![Some(PileupBase::Match { qual: 30 }); len], false, 60);
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
                Some(PileupBase::Match { qual: 30 }),
                Some(PileupBase::Match { qual: 30 }),
                Some(PileupBase::Match { qual: 30 }),
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
        let ar1 = make_ar(vec![Some(PileupBase::Match { qual: 30 }); 5], false, 60);
        let ar2 = make_ar(vec![Some(PileupBase::Match { qual: 30 }); 5], false, 60);
        let mut bases3 = vec![Some(PileupBase::Match { qual: 30 }); 5];
        bases3[2] = Some(PileupBase::Mismatch {
            base: b'T',
            qual: 30,
        });
        let ar3 = make_ar(bases3, true, 55);

        let mut qr = dummy_query_result(vec![
            dummy_record(false, 60),
            dummy_record(false, 60),
            dummy_record(true, 55),
        ]);
        qr.site = SiteStats {
            depth: 3,
            g: 2,
            t: 1,
            fwd: 2,
            rev: 1,
            ..Default::default()
        };
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

    /// When the displayed reads are a capped sample, COV, allele fractions,
    /// strand counts and MQ must come from all passing reads, not the sample.
    #[test]
    fn test_footer_uses_full_depth_when_sampled() {
        let ref_seq = b"ACGTA";
        // Only 2 reads drawn...
        let ar1 = make_ar(vec![Some(PileupBase::Match { qual: 30 }); 5], false, 60);
        let ar2 = make_ar(vec![Some(PileupBase::Match { qual: 30 }); 5], true, 60);
        let mut qr = dummy_query_result(vec![dummy_record(false, 60), dummy_record(true, 60)]);
        // ...out of 120 passing reads, 40 of which cover the query position.
        qr.total_passed = 120;
        qr.mapq_sum = 120 * 50;
        qr.mapq_min = Some(7);
        qr.site = SiteStats {
            depth: 40,
            g: 30,
            t: 10,
            fwd: 22,
            rev: 18,
            ..Default::default()
        };

        let output = render_pileup(&[ar1, ar2], ref_seq, 10, 12, &qr, &default_opts());

        assert!(output.contains("COV   40x"), "{output}");
        assert!(output.contains("G: 30 (75.0%)"), "{output}");
        assert!(output.contains("T: 10 (25.0%)"), "{output}");
        assert!(output.contains("+strand: 22  -strand: 18"), "{output}");
        assert!(output.contains("MQ_mean: 50.0  MQ_min: 7"), "{output}");
        assert!(
            output.contains("Reads shown: 2 of 120 passing filters"),
            "{output}"
        );
    }

    /// All labels are 1-based inclusive: 0-based region [10, 15) with the query
    /// at 0-based 12 is displayed as 11-15 with the query at 13.
    #[test]
    fn test_labels_are_one_based() {
        let ref_seq = b"ACGTA";
        let ar = make_ar(vec![Some(PileupBase::Match { qual: 30 }); 5], false, 60);
        let mut qr = dummy_query_result(vec![dummy_record(false, 60)]);
        qr.site = SiteStats {
            depth: 1,
            g: 1,
            fwd: 1,
            ..Default::default()
        };
        let output = render_pileup(&[ar], ref_seq, 10, 12, &qr, &default_opts());

        assert!(
            output.contains("Region: chr1:11-15  (5bp window)"),
            "{output}"
        );
        // POS ruler numbers are covered by the test_pos_ruler_* tests; this
        // 5bp window is too narrow to show them.
        assert!(output.contains("COV   1x at chr1:13"), "{output}");
        assert!(output.contains("BQ    chr1:13"), "{output}");
    }

    #[test]
    fn test_pos_ruler_normal() {
        // 0-based [100, 121), query at offset 10 → labels 101 .. 121.
        let line = render_pos_ruler(100, 121, 10);
        assert_eq!(
            line,
            format!("POS   101{}[*]{}121", " ".repeat(7), " ".repeat(7))
        );
        // Same width as the REF line content (region_len + 2).
        assert_eq!(line.len(), "POS   ".len() + 21 + 2);
    }

    /// Clipped window at the end of chr1: labels must not overwrite each other
    /// or the `[*]` marker. Previously this rendered "POS   248956248956422".
    #[test]
    fn test_pos_ruler_contig_end_no_overlap() {
        // 0-based [248956409, 248956422), query at offset 10 (1-based 248,956,420).
        let line = render_pos_ruler(248_956_409, 248_956_422, 10);
        let content = &line["POS   ".len()..];
        assert_eq!(&content[10..13], "[*]", "{line:?}");
        // The 9-digit start label fits before `[*]` with a space; the end
        // label has no room after `[*]`, so it is omitted.
        assert_eq!(line, "POS   248956410 [*]");
    }

    #[test]
    fn test_pos_ruler_query_at_first_base() {
        // Query at offset 0: no room for the start label; end label still shown.
        let line = render_pos_ruler(0, 13, 0);
        assert_eq!(line, format!("POS   [*]{}13", " ".repeat(10)));
    }

    #[test]
    fn test_pos_ruler_keeps_one_label_when_only_one_fits() {
        // Room before `[*]` for the start label but not after it for the end label.
        let line = render_pos_ruler(999, 1013, 12);
        let content = &line["POS   ".len()..];
        assert!(content.starts_with("1000 "), "{line:?}");
        assert_eq!(&content[12..15], "[*]", "{line:?}");
        assert!(!content.contains("1013"), "{line:?}");
    }

    #[test]
    fn test_footer_anchored_indels() {
        let ref_seq = b"ACGTA";
        let ar = make_ar(vec![Some(PileupBase::Match { qual: 30 }); 5], false, 60);
        let mut qr = dummy_query_result(vec![dummy_record(false, 60)]);
        qr.site = SiteStats {
            depth: 40,
            g: 40,
            fwd: 20,
            rev: 20,
            ins_after: 25,
            del_after: 2,
            indels: [
                ("+A", 20),
                ("+AA", 2),
                ("+C", 1),
                ("+G", 1),
                ("+T", 1),
                ("-1bp", 2),
            ]
            .into_iter()
            .map(|(l, n)| (l.to_string(), n))
            .collect(),
            ..Default::default()
        };
        let output = render_pileup(&[ar], ref_seq, 10, 12, &qr, &default_opts());
        assert!(
            output.contains("ins after: 25 (62.5%)  +A: 20  +AA: 2  +C: 1  +G: 1  other: 1"),
            "{output}"
        );
        assert!(output.contains("del after: 2 (5.0%)  -1bp: 2"), "{output}");
    }

    #[test]
    fn test_footer_no_indel_lines_without_indels() {
        let ref_seq = b"ACGTA";
        let ar = make_ar(vec![Some(PileupBase::Match { qual: 30 }); 5], false, 60);
        let qr = dummy_query_result(vec![dummy_record(false, 60)]);
        let output = render_pileup(&[ar], ref_seq, 10, 12, &qr, &default_opts());
        assert!(!output.contains("ins after"), "{output}");
        assert!(!output.contains("del after"), "{output}");
    }

    #[test]
    fn test_bq_section() {
        let ref_seq = b"ACGTA";
        let region_start = 10i64;
        let query_pos = 12i64; // offset 2

        // One reverse-strand read (bq=38), one forward-strand read (bq=40), both Match at offset 2.
        let mut bases_rev = vec![Some(PileupBase::Match { qual: 30 }); 5];
        bases_rev[2] = Some(PileupBase::Match { qual: 38 });
        let mut bases_fwd = vec![Some(PileupBase::Match { qual: 30 }); 5];
        bases_fwd[2] = Some(PileupBase::Match { qual: 40 });
        let ar_rev = make_ar(bases_rev, true, 60); // R01, reverse
        let ar_fwd = make_ar(bases_fwd, false, 60); // R02, forward

        let qr = dummy_query_result(vec![dummy_record(true, 60), dummy_record(false, 60)]);
        let output = render_pileup(
            &[ar_rev, ar_fwd],
            ref_seq,
            region_start,
            query_pos,
            &qr,
            &default_opts(),
        );

        // BQ header must be present.
        assert!(
            output.contains("BQ    chr1:13"),
            "BQ header missing: {}",
            output
        );
        // Reverse-strand read (R01, bq=38) should appear before forward-strand (R02, bq=40).
        let bq_header_pos = output.find("BQ    ").unwrap();
        let rev_pos = output.find("R01     -").unwrap();
        let fwd_pos = output.find("R02     +").unwrap();
        assert!(rev_pos > bq_header_pos, "Rev row should be after BQ header");
        assert!(rev_pos < fwd_pos, "Rev row should come before fwd row");
        // Stats line: mean=39.0, min=38, max=40.
        assert!(
            output.contains("mean: 39.0"),
            "Expected mean: 39.0: {}",
            output
        );
        assert!(output.contains("min: 38"), "Expected min: 38: {}", output);
        assert!(output.contains("max: 40"), "Expected max: 40: {}", output);
    }
}
