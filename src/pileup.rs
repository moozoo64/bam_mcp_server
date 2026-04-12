use crate::query::{CigarOp, ReadRecord};

/// The alignment state of a single read at one reference position.
#[derive(Debug, Clone, PartialEq)]
pub enum PileupBase {
    /// Read base matches the reference.
    Match { qual: u8 },
    /// Read base differs from the reference.
    Mismatch { base: u8, qual: u8 },
    /// A deletion in the read spans this reference position.
    Deletion,
}

/// A single read expanded into reference-coordinate space for `[region_start, region_end)`.
#[derive(Debug)]
pub struct AlignedRead {
    pub record: ReadRecord,
    /// One slot per reference position in `[region_start, region_end)`.
    /// `None` = the read does not cover this position.
    pub bases: Vec<Option<PileupBase>>,
    /// Number of inserted bases immediately *after* position `i`.
    /// `0` = no insertion follows. Parallel to `bases`.
    pub ins_after: Vec<u32>,
}

/// Expand a slice of reads into reference-coordinate-aligned arrays.
///
/// * `ref_seq`      — uppercase ASCII bytes; length must equal `region_end - region_start`.
/// * `region_start` — 0-based inclusive left bound of the window.
/// * `region_end`   — 0-based exclusive right bound of the window.
pub fn expand_reads(
    reads: &[ReadRecord],
    ref_seq: &[u8],
    region_start: i64,
    region_end: i64,
) -> Vec<AlignedRead> {
    let region_len = (region_end - region_start) as usize;
    reads
        .iter()
        .map(|record| expand_one(record, ref_seq, region_start, region_len))
        .collect()
}

fn expand_one(
    record: &ReadRecord,
    ref_seq: &[u8],
    region_start: i64,
    region_len: usize,
) -> AlignedRead {
    let mut bases: Vec<Option<PileupBase>> = vec![None; region_len];
    let mut ins_after: Vec<u32> = vec![0u32; region_len];

    let mut ref_pos: i64 = record.ref_start;
    let mut read_pos: usize = 0;

    for op in &record.cigar {
        match op {
            CigarOp::Match(n) => {
                for _ in 0..*n {
                    let ri = ref_pos - region_start;
                    if ri >= 0 && ri < region_len as i64 {
                        let ri = ri as usize;
                        let base = record.sequence[read_pos];
                        let qual = record.base_quals[read_pos];
                        let ref_base = ref_seq[ri];
                        bases[ri] = Some(
                            if base.to_ascii_uppercase() == ref_base.to_ascii_uppercase() {
                                PileupBase::Match { qual }
                            } else {
                                PileupBase::Mismatch { base, qual }
                            },
                        );
                    }
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            CigarOp::Del(n) => {
                for _ in 0..*n {
                    let ri = ref_pos - region_start;
                    if ri >= 0 && ri < region_len as i64 {
                        bases[ri as usize] = Some(PileupBase::Deletion);
                    }
                    ref_pos += 1;
                }
            }
            CigarOp::Ins(n) => {
                // Insertion is attributed to the preceding reference position.
                let prev_ri = ref_pos - region_start - 1;
                if prev_ri >= 0 && prev_ri < region_len as i64 {
                    ins_after[prev_ri as usize] += *n;
                }
                read_pos += *n as usize;
                // ref_pos is unchanged
            }
            CigarOp::SoftClip(n) => {
                read_pos += *n as usize;
                // ref_pos is unchanged
            }
            CigarOp::HardClip(_) => {
                // consumes nothing
            }
            CigarOp::Skip(n) => {
                ref_pos += *n as i64;
                // covered positions remain None
            }
        }
    }

    AlignedRead {
        record: record.clone(),
        bases,
        ins_after,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::query::{CigarOp, ReadRecord};

    // Shared region: reference positions [10, 15), ref bases = b"ACGTA"
    const REGION_START: i64 = 10;
    const REGION_END: i64 = 15;
    const REF: &[u8] = b"ACGTA";

    fn make_record(ref_start: i64, cigar: Vec<CigarOp>, sequence: &[u8]) -> ReadRecord {
        ReadRecord {
            name: "read".to_string(),
            ref_start,
            ref_end: ref_start, // not used by expand_one
            cigar,
            sequence: sequence.to_vec(),
            base_quals: vec![30u8; sequence.len()],
            is_reverse: false,
            mapq: 60,
            is_dup: false,
            is_qcfail: false,
        }
    }

    #[test]
    fn test_pure_match() {
        // 5M starting at region_start: every position should match ACGTA
        let rec = make_record(10, vec![CigarOp::Match(5)], b"ACGTA");
        let result = expand_reads(&[rec], REF, REGION_START, REGION_END);
        let ar = &result[0];
        assert_eq!(ar.bases.len(), 5);
        for b in &ar.bases {
            assert_eq!(b, &Some(PileupBase::Match { qual: 30 }));
        }
        assert!(ar.ins_after.iter().all(|&x| x == 0));
    }

    #[test]
    fn test_snp() {
        // Position 12 (region index 2): read has T where ref has G
        let mut rec = make_record(10, vec![CigarOp::Match(5)], b"ACTTA");
        rec.base_quals[2] = 25;
        let result = expand_reads(&[rec], REF, REGION_START, REGION_END);
        let ar = &result[0];
        assert_eq!(ar.bases[0], Some(PileupBase::Match { qual: 30 })); // A==A
        assert_eq!(ar.bases[1], Some(PileupBase::Match { qual: 30 })); // C==C
        assert_eq!(
            ar.bases[2],
            Some(PileupBase::Mismatch {
                base: b'T',
                qual: 25
            })
        );
        assert_eq!(ar.bases[3], Some(PileupBase::Match { qual: 30 })); // T==T
        assert_eq!(ar.bases[4], Some(PileupBase::Match { qual: 30 })); // A==A
    }

    #[test]
    fn test_deletion() {
        // 2M 1D 2M: position 12 (index 2) is a deletion
        // read bases consumed: A,C then T,A (4 bytes)
        let rec = make_record(
            10,
            vec![CigarOp::Match(2), CigarOp::Del(1), CigarOp::Match(2)],
            b"ACTA",
        );
        let result = expand_reads(&[rec], REF, REGION_START, REGION_END);
        let ar = &result[0];
        assert_eq!(ar.bases[0], Some(PileupBase::Match { qual: 30 })); // pos10: A==A
        assert_eq!(ar.bases[1], Some(PileupBase::Match { qual: 30 })); // pos11: C==C
        assert_eq!(ar.bases[2], Some(PileupBase::Deletion)); // pos12
        assert_eq!(ar.bases[3], Some(PileupBase::Match { qual: 30 })); // pos13: T==T
        assert_eq!(ar.bases[4], Some(PileupBase::Match { qual: 30 })); // pos14: A==A
    }

    #[test]
    fn test_insertion() {
        // 2M 3I 3M: insertion after region index 1 (ref_pos 11)
        // read bytes: A,C (2M) + N,N,N (3I) + G,T,A (3M) = 8 bytes
        let rec = make_record(
            10,
            vec![CigarOp::Match(2), CigarOp::Ins(3), CigarOp::Match(3)],
            b"ACNNNGTA",
        );
        let result = expand_reads(&[rec], REF, REGION_START, REGION_END);
        let ar = &result[0];
        assert_eq!(ar.bases[0], Some(PileupBase::Match { qual: 30 })); // pos10: A==A
        assert_eq!(ar.bases[1], Some(PileupBase::Match { qual: 30 })); // pos11: C==C
        assert_eq!(ar.ins_after[1], 3); // 3 inserted bases after region index 1
        assert_eq!(ar.bases[2], Some(PileupBase::Match { qual: 30 })); // pos12: G==G
        assert_eq!(ar.bases[3], Some(PileupBase::Match { qual: 30 })); // pos13: T==T
        assert_eq!(ar.bases[4], Some(PileupBase::Match { qual: 30 })); // pos14: A==A
    }

    #[test]
    fn test_soft_clip() {
        // 2S 5M: soft clip eats 2 read bases without advancing ref_pos
        // seq: NN + ACGTA (7 bytes total)
        let rec = make_record(
            10,
            vec![CigarOp::SoftClip(2), CigarOp::Match(5)],
            b"NNACGTA",
        );
        let result = expand_reads(&[rec], REF, REGION_START, REGION_END);
        let ar = &result[0];
        for b in &ar.bases {
            assert_eq!(b, &Some(PileupBase::Match { qual: 30 }));
        }
    }

    #[test]
    fn test_partial_overlap_left() {
        // Read starts at ref_pos 8, 10M: covers 8..18; only [10,15) is in region
        // seq positions 2..7 land in region (indices 0..5 of seq = NN + ACGTA + YYY)
        let rec = make_record(8, vec![CigarOp::Match(10)], b"NNACGTAYYY");
        let result = expand_reads(&[rec], REF, REGION_START, REGION_END);
        let ar = &result[0];
        for b in &ar.bases {
            assert_eq!(b, &Some(PileupBase::Match { qual: 30 }));
        }
    }

    #[test]
    fn test_partial_overlap_right() {
        // Read starts at ref_pos 13, 5M: covers 13..18; only [13,15) is in region
        // seq[0]=T matches ref[3]=T; seq[1]=A matches ref[4]=A
        let rec = make_record(13, vec![CigarOp::Match(5)], b"TAABC");
        let result = expand_reads(&[rec], REF, REGION_START, REGION_END);
        let ar = &result[0];
        assert_eq!(ar.bases[0], None); // pos10
        assert_eq!(ar.bases[1], None); // pos11
        assert_eq!(ar.bases[2], None); // pos12
        assert_eq!(ar.bases[3], Some(PileupBase::Match { qual: 30 })); // pos13: T==T
        assert_eq!(ar.bases[4], Some(PileupBase::Match { qual: 30 })); // pos14: A==A
    }
}
