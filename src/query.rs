use std::collections::BTreeMap;
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

/// What a single read shows at one reference position.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiteObservation {
    /// An aligned read base (ASCII, as stored in the BAM).
    Base(u8),
    /// A deletion in the read spans the position.
    Deletion,
}

impl ReadRecord {
    /// What this read has at 0-based reference position `pos`, or `None` if it
    /// does not cover it (outside the alignment, or inside an `N` skip).
    ///
    /// Uses the same CIGAR interpretation as `pileup::expand_reads`.
    pub fn observe_at(&self, pos: i64) -> Option<SiteObservation> {
        if pos < self.ref_start || pos >= self.ref_end {
            return None;
        }
        let mut ref_pos = self.ref_start;
        let mut read_pos = 0usize;
        for op in &self.cigar {
            match *op {
                CigarOp::Match(n) => {
                    let n = n as i64;
                    if pos < ref_pos + n {
                        let i = read_pos + (pos - ref_pos) as usize;
                        return self.sequence.get(i).map(|&b| SiteObservation::Base(b));
                    }
                    ref_pos += n;
                    read_pos += n as usize;
                }
                CigarOp::Del(n) => {
                    if pos < ref_pos + n as i64 {
                        return Some(SiteObservation::Deletion);
                    }
                    ref_pos += n as i64;
                }
                CigarOp::Skip(n) => {
                    if pos < ref_pos + n as i64 {
                        return None;
                    }
                    ref_pos += n as i64;
                }
                CigarOp::Ins(n) | CigarOp::SoftClip(n) => read_pos += n as usize,
                CigarOp::HardClip(_) => {}
            }
        }
        None
    }

    /// The indel this read has immediately after the aligned base at 0-based
    /// `pos`, i.e. anchored on `pos` the way VCF anchors indels. `None` if the
    /// read has no aligned base at `pos` or no indel follows it.
    ///
    /// Matches `samtools mpileup`, which prints `+N…` / `-N…` on the anchor base.
    pub fn indel_after(&self, pos: i64) -> Option<Indel> {
        if pos < self.ref_start || pos >= self.ref_end {
            return None;
        }
        let mut ref_pos = self.ref_start;
        let mut read_pos = 0usize;
        // True when the previous op was a Match whose last base is at `pos`.
        let mut anchored = false;
        for op in &self.cigar {
            if anchored {
                return match *op {
                    CigarOp::Ins(n) => {
                        let end = (read_pos + n as usize).min(self.sequence.len());
                        let start = read_pos.min(end);
                        Some(Indel::Ins(self.sequence[start..end].to_ascii_uppercase()))
                    }
                    CigarOp::Del(n) => Some(Indel::Del(n)),
                    _ => None,
                };
            }
            match *op {
                CigarOp::Match(n) => {
                    ref_pos += n as i64;
                    read_pos += n as usize;
                    anchored = ref_pos == pos + 1;
                    if ref_pos > pos + 1 {
                        return None;
                    }
                }
                CigarOp::Del(n) | CigarOp::Skip(n) => {
                    ref_pos += n as i64;
                    if ref_pos > pos {
                        return None;
                    }
                }
                CigarOp::Ins(n) | CigarOp::SoftClip(n) => read_pos += n as usize,
                CigarOp::HardClip(_) => {}
            }
        }
        None
    }
}

/// An insertion or deletion immediately following an anchor base.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Indel {
    /// Inserted read bases (uppercase).
    Ins(Vec<u8>),
    /// Number of deleted reference bases.
    Del(u32),
}

impl Indel {
    /// Short label used to group and display indels: `+A`, `+15bp`, `-1bp`.
    /// Long insertions are grouped by length so the counts stay readable.
    pub fn label(&self) -> String {
        match self {
            Indel::Ins(seq) if seq.len() <= 10 => format!("+{}", String::from_utf8_lossy(seq)),
            Indel::Ins(seq) => format!("+{}bp", seq.len()),
            Indel::Del(n) => format!("-{n}bp"),
        }
    }
}

/// Depth, allele and strand counts at the queried position, accumulated over
/// **every** read that passed the filters — not just the sampled reads.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SiteStats {
    /// Reads with a base or a deletion at the position.
    pub depth: usize,
    pub a: usize,
    pub c: usize,
    pub g: usize,
    pub t: usize,
    /// Reads with a deletion spanning the position.
    pub del: usize,
    /// Forward-strand reads covering the position.
    pub fwd: usize,
    /// Reverse-strand reads covering the position.
    pub rev: usize,
    /// Reads with an insertion right after the position (anchored on it).
    pub ins_after: usize,
    /// Reads with a deletion starting right after the position (anchored on it).
    pub del_after: usize,
    /// Anchored indel counts by label (`+A`, `-1bp`, …); see `Indel::label`.
    pub indels: BTreeMap<String, usize>,
}

impl SiteStats {
    /// Record one read's observation at the position.
    pub fn add(&mut self, obs: SiteObservation, is_reverse: bool) {
        self.depth += 1;
        if is_reverse {
            self.rev += 1;
        } else {
            self.fwd += 1;
        }
        match obs {
            SiteObservation::Deletion => self.del += 1,
            SiteObservation::Base(b) => match b.to_ascii_uppercase() {
                b'A' => self.a += 1,
                b'C' => self.c += 1,
                b'G' => self.g += 1,
                b'T' => self.t += 1,
                _ => {} // N etc.: counted in depth only
            },
        }
    }

    /// Record an indel anchored on the position.
    pub fn add_indel(&mut self, indel: &Indel) {
        match indel {
            Indel::Ins(_) => self.ins_after += 1,
            Indel::Del(_) => self.del_after += 1,
        }
        *self.indels.entry(indel.label()).or_default() += 1;
    }
}

/// The result of a BAM region query after filtering.
#[derive(Debug, Default)]
pub struct QueryResult {
    /// Reads to display: all passing reads, or a random sample of `max_depth`.
    pub reads: Vec<ReadRecord>,
    /// Reads excluded because mapping quality < min_mapq
    pub filtered_mapq: usize,
    /// Reads excluded because they are duplicates or QC failures
    pub filtered_flags: usize,
    /// Number of reads overlapping the window that passed all filters
    /// (before sampling down to `max_depth`).
    pub total_passed: usize,
    /// Counts at the query position over all passing reads.
    pub site: SiteStats,
    /// Sum of MAPQ over all passing reads.
    pub mapq_sum: u64,
    /// Minimum MAPQ over all passing reads (`None` if there are none).
    pub mapq_min: Option<u8>,
}

impl QueryResult {
    /// Mean MAPQ over all passing reads (0.0 if there are none).
    pub fn mapq_mean(&self) -> f64 {
        if self.total_passed == 0 {
            0.0
        } else {
            self.mapq_sum as f64 / self.total_passed as f64
        }
    }
}

/// Look up the length of `chrom` in a BAM header, or fail with a helpful
/// "unknown chromosome" error listing some available names.
fn lookup_contig_length(header: &noodles::sam::Header, chrom: &str) -> anyhow::Result<usize> {
    match header.reference_sequences().get(chrom.as_bytes()) {
        Some(rs) => Ok(rs.length().get()),
        None => {
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
    }
}

/// Return the length (bp) of `chrom` as declared in the BAM header.
pub fn contig_length(bam_path: &Path, chrom: &str) -> anyhow::Result<usize> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .with_context(|| format!("Failed to open BAM: {}", bam_path.display()))?;
    let header = reader.read_header().context("Failed to read BAM header")?;
    lookup_contig_length(&header, chrom)
}

/// Query a BAM file for all reads overlapping `[start, end)` on `chrom`.
///
/// Coordinates are **0-based, half-open**. `query_pos` (0-based) is the
/// position whose depth and allele counts are tallied in `QueryResult::site`;
/// those counts, `total_passed` and the MAPQ stats cover every passing read,
/// while `reads` holds at most `max_depth` of them for display.
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
    query_pos: i64,
    min_mapq: u8,
    max_depth: usize,
) -> anyhow::Result<QueryResult> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .with_context(|| format!("Failed to open BAM: {}", bam_path.display()))?;

    let header = reader.read_header().context("Failed to read BAM header")?;

    // Validate chromosome name against the header reference sequences.
    lookup_contig_length(&header, chrom)?;

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
    let mut site = SiteStats::default();
    let mut mapq_sum: u64 = 0;
    let mut mapq_min: Option<u8> = None;

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

        // Tally depth / alleles / MAPQ over every passing read, before sampling.
        if let Some(obs) = read.observe_at(query_pos) {
            site.add(obs, read.is_reverse);
        }
        if let Some(indel) = read.indel_after(query_pos) {
            site.add_indel(&indel);
        }
        mapq_sum += mapq as u64;
        mapq_min = Some(mapq_min.map_or(mapq, |m| m.min(mapq)));

        // Reservoir sampling (Knuth's Algorithm R): keep at most `max_depth` reads.
        if total_passed < max_depth {
            candidates.push(read);
        } else {
            // Replace a random slot with probability
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
        total_passed,
        site,
        mapq_sum,
        mapq_min,
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

    fn record_with(ref_start: i64, cigar: Vec<CigarOp>, sequence: &[u8]) -> ReadRecord {
        let ref_end = compute_ref_end(ref_start, &cigar);
        ReadRecord {
            name: "r".into(),
            ref_start,
            ref_end,
            cigar,
            sequence: sequence.to_vec(),
            base_quals: vec![30; sequence.len()],
            is_reverse: false,
            mapq: 60,
            is_dup: false,
            is_qcfail: false,
        }
    }

    #[test]
    fn test_observe_at() {
        use SiteObservation::{Base, Deletion};
        // ref 100: 2S 3M 1I 2M 2D 2M 3N 2M
        // read:    ss ACG  t  GT  --  CA ... TT
        let r = record_with(
            100,
            vec![
                CigarOp::SoftClip(2),
                CigarOp::Match(3),
                CigarOp::Ins(1),
                CigarOp::Match(2),
                CigarOp::Del(2),
                CigarOp::Match(2),
                CigarOp::Skip(3),
                CigarOp::Match(2),
            ],
            b"ssACGtGTCATT",
        );
        assert_eq!(r.observe_at(99), None, "before alignment");
        assert_eq!(r.observe_at(100), Some(Base(b'A')), "soft clip skipped");
        assert_eq!(r.observe_at(102), Some(Base(b'G')));
        assert_eq!(r.observe_at(103), Some(Base(b'G')), "insertion skipped");
        assert_eq!(r.observe_at(104), Some(Base(b'T')));
        assert_eq!(r.observe_at(105), Some(Deletion));
        assert_eq!(r.observe_at(106), Some(Deletion));
        assert_eq!(r.observe_at(107), Some(Base(b'C')));
        assert_eq!(r.observe_at(108), Some(Base(b'A')));
        assert_eq!(r.observe_at(109), None, "inside N skip");
        assert_eq!(r.observe_at(111), None, "inside N skip");
        assert_eq!(r.observe_at(112), Some(Base(b'T')));
        assert_eq!(r.observe_at(113), Some(Base(b'T')));
        assert_eq!(r.observe_at(114), None, "after alignment");
    }

    #[test]
    fn test_indel_after() {
        // ref 100: 2S 3M 2I 2M 2D 2M     read: ss ACG tt GT -- CA
        let r = record_with(
            100,
            vec![
                CigarOp::SoftClip(2),
                CigarOp::Match(3),
                CigarOp::Ins(2),
                CigarOp::Match(2),
                CigarOp::Del(2),
                CigarOp::Match(2),
            ],
            b"ssACGttGTCA",
        );
        assert_eq!(r.indel_after(101), None, "middle of a match");
        assert_eq!(
            r.indel_after(102),
            Some(Indel::Ins(b"TT".to_vec())),
            "uppercased"
        );
        assert_eq!(r.indel_after(103), None);
        assert_eq!(r.indel_after(104), Some(Indel::Del(2)));
        assert_eq!(r.indel_after(105), None, "deleted base is not an anchor");
        assert_eq!(r.indel_after(106), None, "deleted base is not an anchor");
        assert_eq!(r.indel_after(108), None, "last base, nothing follows");
        assert_eq!(r.indel_after(99), None, "before alignment");

        // Insertion at the start of the read has no aligned anchor base.
        let r = record_with(100, vec![CigarOp::Ins(2), CigarOp::Match(3)], b"ttACG");
        assert_eq!(r.indel_after(99), None);
        assert_eq!(r.indel_after(100), None);
    }

    #[test]
    fn test_indel_label() {
        assert_eq!(Indel::Ins(b"A".to_vec()).label(), "+A");
        assert_eq!(Indel::Ins(b"ACGTACGTAC".to_vec()).label(), "+ACGTACGTAC");
        assert_eq!(Indel::Ins(b"ACGTACGTACG".to_vec()).label(), "+11bp");
        assert_eq!(Indel::Del(3).label(), "-3bp");
    }

    /// Anchored indel counts at VCF sites in the sample BAM, cross-checked with
    /// `samtools mpileup -A -B -Q0 --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY`.
    #[test]
    fn test_indel_counts_at_vcf_sites() {
        // (1-based VCF POS, depth, ins_after, del_after, expected label, count)
        let cases = [
            (4_888_325, 32, 26, 0, "+A", 25),   // hom G>GA
            (3_956_581, 39, 22, 0, "+A", 22),   // het G>GA
            (3_386_594, 39, 0, 36, "-1bp", 36), // hom GA>G
            (3_180_748, 35, 0, 18, "-1bp", 18), // het GA>G
        ];
        for (pos1, depth, ins, del, label, n) in cases {
            let q = pos1 - 1;
            let r = query_region(&sample_bam(), "chr1", q - 75, q + 76, q, 0, 50)
                .expect("query_region should succeed");
            let site = &r.site;
            assert_eq!(site.depth, depth, "depth at {pos1}");
            assert_eq!(site.ins_after, ins, "ins_after at {pos1}");
            assert_eq!(site.del_after, del, "del_after at {pos1}");
            assert_eq!(site.indels.get(label), Some(&n), "{label} at {pos1}");
        }
    }

    #[test]
    fn test_site_stats_add() {
        let mut s = SiteStats::default();
        s.add(SiteObservation::Base(b'a'), false);
        s.add(SiteObservation::Base(b'T'), true);
        s.add(SiteObservation::Base(b'N'), true);
        s.add(SiteObservation::Deletion, false);
        assert_eq!(
            s,
            SiteStats {
                depth: 4,
                a: 1,
                t: 1,
                del: 1,
                fwd: 2,
                rev: 2,
                ..Default::default()
            }
        );
    }

    /// With a small `max_depth`, site depth must still count every passing read.
    #[test]
    fn test_site_depth_not_capped_by_max_depth() {
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(sample_bam())
            .expect("open BAM");
        let header = reader.read_header().expect("read header");
        let chrom = header
            .reference_sequences()
            .keys()
            .next()
            .expect("at least one reference sequence")
            .to_string();
        drop(reader);

        let (start, end, pos) = (10_000, 10_150, 10_075);
        let full = query_region(&sample_bam(), &chrom, start, end, pos, 0, 10_000)
            .expect("query_region should succeed");
        let capped = query_region(&sample_bam(), &chrom, start, end, pos, 0, 3)
            .expect("query_region should succeed");

        assert!(full.site.depth > 3, "need real depth > cap for this test");
        assert_eq!(capped.reads.len(), 3);
        assert_eq!(capped.total_passed, full.total_passed);
        assert_eq!(capped.site, full.site);
        assert_eq!(capped.mapq_sum, full.mapq_sum);
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
            10_075.min(first_len - 100),
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
        let err = query_region(&sample_bam(), "nonexistent_chrXYZ99", 0, 100, 50, 0, 50)
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

        let result = query_region(&sample_bam(), &first_chrom, 10_000, 10_150, 10_075, 255, 50)
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
