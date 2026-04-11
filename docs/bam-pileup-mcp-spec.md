# BAM Pileup MCP Server — Specification & Implementation Plan

## Overview

A Model Context Protocol (MCP) server written in Rust (2024 edition) that accepts a genomic
coordinate, queries a BAM file, and returns a text-based pileup representation suitable for
reasoning by an LLM. The server runs as a stdio-based MCP process, is stateless per call, and
requires a sorted, indexed BAM file and an accompanying FASTA reference.

---

## 1. Command-Line Interface

```
bam-pileup-mcp [OPTIONS] --bam <PATH> --reference <PATH>

Options:
  -b, --bam <PATH>           Path to sorted, indexed BAM file (.bam + .bai must exist)
  -r, --reference <PATH>     Path to FASTA reference file (.fa/.fasta + .fai must exist)
  -w, --window <INT>         Default window half-width in bp around query position [default: 75]
  -d, --max-depth <INT>      Maximum reads to display in pileup [default: 50]
  -q, --min-mapq <INT>       Minimum mapping quality filter [default: 0]
  -Q, --min-baseq <INT>      Minimum base quality to show as uppercase [default: 20]
      --debug                Write debug output to stderr (MCP uses stdout for protocol)
      --log-file <PATH>      Write debug log to file instead of stderr
  -h, --help                 Print help
  -V, --version              Print version
```

### Notes
- `--debug` must never write to stdout — MCP protocol uses stdout exclusively for JSON-RPC
- `.bai` index is assumed to be at `<bam>.bai`; `.fai` index at `<reference>.fai`
- `--window` can be overridden per-call by the MCP tool arguments
- All paths are validated at startup and the process exits with a clear error if any are missing

---

## 2. MCP Tool Specification

### Tool: `query_pileup`

**Description:** Query a BAM file at a genomic position and return a text pileup showing
aligned short reads with dots for reference matches and explicit bases for mismatches,
insertions, and deletions.

#### Input Schema

```json
{
  "type": "object",
  "properties": {
    "chrom": {
      "type": "string",
      "description": "Chromosome name exactly as it appears in the BAM header (e.g. chr1, 1, chrM)"
    },
    "position": {
      "type": "integer",
      "description": "1-based genomic position to query"
    },
    "window": {
      "type": "integer",edition download
      "description": "Half-width of region around position in bp. Overrides --window CLI default.",
      "minimum": 10,
      "maximum": 500
    },
    "min_mapq": {
      "type": "integer",
      "description": "Minimum mapping quality. Overrides --min-mapq CLI default.",
      "minimum": 0,
      "maximum": 60
    },
    "show_strand": {
      "type": "boolean",
      "description": "Show strand (+/-) column per read",
      "default": true
    },
    "show_mapq": {
      "type": "boolean",
      "description": "Show mapping quality column per read",
      "default": true
    }
  },
  "required": ["chrom", "position"]
}
```

#### Output

Returns a single `text` content block containing the pileup as a UTF-8 string.

---

## 3. Output Format

### 3.1 Header

```
Region: chr7:117,548,920-117,549,070  (150bp window)  Reference: GRCh38  BAM: sample.bam
```

### 3.2 Reference Line

```
POS   117548920                   [*]                              117549070
REF   TGCAATCCGAATCGGCATGCCTACG...[A]...CATGCATCGAATCGGCATGCCT
```

The queried position is bracketed with `[X]` in the reference line.

### 3.3 Read Lines

```
R01   ..............[G]..................................  + MQ60
R02   ...................................................  - MQ60
R03        ..........................[A]...........         + MQ59
```

#### Per-base encoding

| Situation                        | Symbol         |
|----------------------------------|----------------|
| Match to reference               | `.`            |
| Mismatch (high base quality)     | Base letter (ACGT) uppercase |
| Mismatch (low base quality)      | Base letter lowercase |
| Deletion in read                 | `-`            |
| Insertion after this position    | `^N` where N = inserted bases count |
| Soft clip (not shown by default) | omitted, read starts/ends at alignment start |
| No coverage (read doesn't span)  | ` ` (space)    |
| Queried position                 | `[X]` bracketed regardless of match/mismatch |

#### Per-read suffix (if enabled)

```
  + MQ60    →   forward strand, mapping quality 60
  - MQ23    →   reverse strand, mapping quality 23
```

### 3.4 Footer / Summary

```
COV   47x at chr7:117,548,975
      A: 24 (51.1%)  G: 23 (48.9%)
      +strand: 25  -strand: 22
      MQ_mean: 59.1  MQ_min: 23
      Reads shown: 47  Reads filtered (MAPQ): 2  Reads filtered (dup/qcfail): 1
```

### 3.5 Full Example Output

```
Region: chr7:117,548,920-117,549,070  (150bp window)  BAM: SQ3J6L38.bam

POS   117548920                                            117549070
REF   TGCAATCCGAATCGGCATGCCTACGATTTACGGCATGCATCGAATCGGCAT[A]GCATCGGAATCGCATGCCTACGATTTACGGCATGCATCGAATCGGCATGCCT

R01   .............................................................................................................  + MQ60
R02   ...............................................................[G]............................................ + MQ60
R03          ..........................................................................................             + MQ59
R04   ...............................................................[G]............................................ - MQ60
R05                 ......................................................................                         - MQ58
R06   ...............................................................[G]............................................ + MQ60
R07   ...............................................................[A]............................................ - MQ60
R08   ...............................................................[G]............................................ + MQ60

COV   8x at chr7:117,548,975
      A: 4 (50.0%)  G: 4 (50.0%)
      +strand: 5  -strand: 3
      MQ_mean: 59.4  MQ_min: 58
      Reads shown: 8  Reads filtered (MAPQ): 0  Reads filtered (dup/qcfail): 0
```

---

## 4. Crate Dependencies

```toml
[package]
name    = "bam-pileup-mcp"
version = "0.1.0"
edition = "2024"

[dependencies]
# MCP server framework
rmcp = { version = "0.1", features = ["server", "transport-io"] }

# Genomics I/O
noodles-bam  = "0.72"
noodles-bai  = "0.52"
noodles-fasta = "0.40"
noodles-fai  = "0.25"
noodles-core = "0.16"
noodles-sam  = "0.66"    # for record/cigar types shared with BAM

# Async runtime (rmcp is async)
tokio = { version = "1", features = ["full"] }

# CLI
clap = { version = "4", features = ["derive"] }

# Serialization (MCP JSON-RPC)
serde       = { version = "1", features = ["derive"] }
serde_json  = "1"

# Base64 (if image content ever added later)
base64 = "0.22"

# Error handling
anyhow  = "1"
thiserror = "1"

# Logging (to stderr only)
tracing            = "0.1"
tracing-subscriber = { version = "0.3", features = ["env-filter"] }
```

---

## 5. Module Structure

```
bam-pileup-mcp/
├── Cargo.toml
├── Cargo.lock
├── README.md
└── src/
    ├── main.rs          # Entry point: parse CLI, init logging, start MCP server
    ├── cli.rs           # Clap CLI struct and validation
    ├── server.rs        # rmcp Tool impl, request dispatch
    ├── query.rs         # BAM region query, read collection, MAPQ/flag filtering
    ├── pileup.rs        # Core pileup logic: CIGAR expansion, row alignment
    ├── reference.rs     # FASTA/FAI reference sequence fetching
    ├── render.rs        # Text rendering: dots-for-matches, brackets, summary
    └── error.rs         # Error types
```

---

## 6. Implementation Plan

### Phase 1 — Scaffold & CLI (Day 1 morning)

**`cli.rs`**
- Define `Args` struct with clap `#[derive(Parser)]`
- Validate at startup:
  - BAM file exists and `.bai` sidecar exists
  - FASTA file exists and `.fai` sidecar exists
  - window > 0, max_depth > 0
- Store validated paths in an `AppConfig` struct passed to the server

**`main.rs`**
- Init `tracing_subscriber` directed to stderr (never stdout)
- Parse `Args`, build `AppConfig`
- Construct `PileupServer` and hand to rmcp's stdio transport

```rust
#[tokio::main]
async fn main() -> anyhow::Result<()> {
    let args = Args::parse();
    init_logging(&args);
    let config = AppConfig::try_from(args)?;
    let server = PileupServer::new(config);
    let transport = rmcp::transport::io::stdio();
    rmcp::serve(server, transport).await?;
    Ok(())
}
```

---

### Phase 2 — Reference Fetching (Day 1 afternoon)

**`reference.rs`**

```rust
pub struct ReferenceReader {
    path: PathBuf,
    index: fai::Index,
}

impl ReferenceReader {
    pub fn open(fasta: &Path) -> anyhow::Result<Self>
    pub fn fetch(&self, chrom: &str, start: usize, end: usize) -> anyhow::Result<Vec<u8>>
}
```

- Use `noodles_fasta` with `fai` index for random access
- Return raw bytes (uppercase); caller converts to `Vec<char>` for rendering
- Validate that requested region is within sequence bounds

---

### Phase 3 — BAM Query & Filtering (Day 1 afternoon)

**`query.rs`**

```rust
pub struct ReadRecord {
    pub name:       String,
    pub ref_start:  i64,       // 0-based
    pub ref_end:    i64,       // 0-based exclusive
    pub cigar:      Vec<CigarOp>,
    pub sequence:   Vec<u8>,   // raw bases
    pub base_quals: Vec<u8>,
    pub is_reverse: bool,
    pub mapq:       u8,
    pub is_dup:     bool,
    pub is_qcfail:  bool,
}

pub struct QueryResult {
    pub reads:           Vec<ReadRecord>,
    pub filtered_mapq:   usize,
    pub filtered_flags:  usize,
}

pub fn query_region(
    bam_path: &Path,
    chrom: &str,
    start: i64,   // 0-based
    end: i64,     // 0-based exclusive
    min_mapq: u8,
) -> anyhow::Result<QueryResult>
```

- Open BAM with `noodles_bam::io::indexed_reader::Builder`
- Parse header to validate chromosome name exists; return clear error if not
- Build `noodles_core::Region` from chrom/start/end
- Iterate query results, extract fields, apply filters:
  - Skip reads below `min_mapq`
  - Skip secondary, supplementary alignments
  - Track (but don't display) duplicates and QC-fail reads separately
- Cap at `max_depth` reads (reservoir sample if over cap)

**CIGAR op enum (internal)**

```rust
pub enum CigarOp {
    Match(u32),      // M, = (ref consumed, read consumed)
    Ins(u32),        // I (read consumed only)
    Del(u32),        // D (ref consumed only)
    SoftClip(u32),   // S (read consumed, not shown)
    HardClip(u32),   // H (ignored)
    Skip(u32),       // N (intron/large gap)
}
```

---

### Phase 4 — Pileup Construction (Day 2 morning)

**`pileup.rs`**

The core of the work. Expands each read's CIGAR into a reference-coordinate-aligned sequence.

```rust
pub struct AlignedRead {
    pub record:      ReadRecord,
    pub ref_bases:   Vec<Option<PileupBase>>,  // indexed by ref position offset
}

pub enum PileupBase {
    Match,
    Mismatch { base: char, qual: u8 },
    Deletion,
    Insertion { count: u32 },  // insertion after this ref position
}

pub fn expand_reads(
    reads: &[ReadRecord],
    ref_seq: &[u8],
    region_start: i64,
    region_end: i64,
) -> Vec<AlignedRead>
```

**CIGAR expansion algorithm:**

```
read_pos  = 0
ref_pos   = record.ref_start

for each cigar op:
  match op:
    Match(n) | SeqMatch(n) | SeqMismatch(n):
      for i in 0..n:
        ref_base  = ref_seq[ref_pos - region_start]
        read_base = record.sequence[read_pos]
        emit Match or Mismatch at ref_pos
        ref_pos++; read_pos++

    Del(n):
      for i in 0..n:
        emit Deletion at ref_pos
        ref_pos++

    Ins(n):
      emit Insertion{count: n} at current ref_pos (after)
      read_pos += n

    SoftClip(n):
      read_pos += n   // skip, do not advance ref_pos

    Skip(n):
      ref_pos += n    // splice gap, emit spaces
```

---

### Phase 5 — Text Rendering (Day 2 morning)

**`render.rs`**

```rust
pub fn render_pileup(
    reads:        &[AlignedRead],
    ref_seq:      &[u8],
    region_start: i64,
    query_pos:    i64,
    query_result: &QueryResult,
    opts:         &RenderOpts,
) -> String
```

**Rendering steps:**

1. **Header line** — region, window size, BAM filename
2. **POS ruler** — start and end positions, query position marker
3. **REF line** — reference bases with `[X]` at query position
4. **Separator** — `|||...|||`
5. **Read lines** — for each `AlignedRead`:
   - Left-pad with spaces to `ref_start - region_start`
   - For each ref position in window:
     - space if read doesn't span
     - `.` if match
     - uppercase base if mismatch + high qual
     - lowercase base if mismatch + low qual
     - `-` if deletion
     - `^N` if insertion follows (N = count)
     - `[X]` bracket at query position regardless
   - Append strand and MAPQ suffix
6. **Separator**
7. **Summary footer** — coverage, allele counts, strand balance, filter counts

**Allele counting** at query position:
- Walk all reads, extract base at `query_pos`
- Count each distinct base (A/C/G/T/del)
- Compute percentages

---

### Phase 6 — MCP Server (Day 2 afternoon)

**`server.rs`**

```rust
use rmcp::{ServerHandler, Tool, ToolResult, Content};

pub struct PileupServer {
    config: Arc<AppConfig>,
}

#[rmcp::tool(name = "query_pileup")]
impl PileupServer {
    async fn query_pileup(&self, params: QueryPileupParams) -> ToolResult {
        // 1. Resolve window (param override or config default)
        // 2. Compute region start/end from position ± window
        // 3. Fetch reference sequence
        // 4. Query BAM
        // 5. Expand CIGAR into AlignedReads
        // 6. Render to string
        // 7. Return as Content::text(rendered)
    }
}
```

rmcp handles:
- JSON-RPC framing over stdio
- `initialize` / `tools/list` / `tools/call` protocol messages
- Serialization of input params
- Error responses for invalid params

---

## 7. Error Handling Strategy

| Error situation                        | Response                                      |
|----------------------------------------|-----------------------------------------------|
| Chromosome not in BAM header           | MCP error: `"Unknown chromosome: chrX"`       |
| Position out of reference bounds       | MCP error with valid range                    |
| BAM index missing at runtime           | MCP error (caught at startup ideally)         |
| Zero reads in region                   | Valid response: pileup with "No reads in region" message |
| Region too large (> 2× max window)     | MCP error: suggest smaller window             |
| CIGAR parse failure on a read          | Log warning to stderr, skip that read, continue |

All user-facing errors are `anyhow::Result` propagated to rmcp as tool errors.
Internal unexpected errors use `tracing::error!` to stderr before returning.

---

## 8. Debug Logging

With `--debug` or `RUST_LOG=debug`:

```
[DEBUG] BAM: /data/SQ3J6L38.bam  Reference: /data/GRCh38.fa
[DEBUG] query_pileup called: chr7:117548975 window=75 min_mapq=0
[DEBUG] Region: chr7:117548900-117549050
[DEBUG] Reference fetched: 150bp starting TGCAATCC...
[DEBUG] BAM query returned 49 records
[DEBUG] Filtered: 2 MAPQ, 0 dup, 0 qcfail → 47 displayed
[DEBUG] CIGAR expansion: 47 reads expanded in 1.2ms
[DEBUG] Render complete: 54 lines, 8240 chars
```

All to stderr via `tracing`. Never touches stdout.

---

## 9. Testing Plan

### Unit tests

- `reference.rs`: fetch known sequence from test FASTA, assert bytes match
- `pileup.rs`: hand-crafted CIGAR strings → assert correct expansion
  - Pure match read
  - Read with SNP at known position
  - Read with 2bp deletion
  - Read with 3bp insertion
  - Soft-clipped read
- `render.rs`: known `AlignedRead` set → assert specific lines in output string

### Integration tests

- Use a small synthetic BAM (generate with `samtools view` from NA12878 or similar)
- Known het SNP position → assert allele counts in footer
- Known homozygous position → assert 100% reference dots
- Region with no reads → assert graceful empty output
- Invalid chromosome → assert MCP error response

### Manual smoke test

```bash
# Start server
bam-pileup-mcp --bam sample.bam --reference GRCh38.fa --debug

# In another terminal via MCP client or raw JSON-RPC:
{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"query_pileup","arguments":{"chrom":"chr7","position":117548975}}}
```

---

## 10. Known Limitations & Future Work

| Item | Notes |
|------|-------|
| Multi-allelic indels | Complex indel rendering is hard; v1 shows `-` for all deletions |
| Long reads (ONT/PacBio) | Window size may need to be much larger; error rate rendering gets noisy |
| CRAM support | noodles supports CRAM; add `--cram` flag later |
| Phasing | Reads could be coloured by HP tag in a future version |
| Multiple BAM files | Tumour/normal pair comparison would be a natural extension |
| Strand bias test | Could compute Fisher exact p-value and include in footer |
| Base quality histogram | Per-position BQ summary useful for diagnosing systematic errors |

---

## 11. Estimated Effort

| Phase | Estimated time |
|-------|----------------|
| Phase 1: Scaffold & CLI | 2–3 hours |
| Phase 2: Reference fetching | 2–3 hours |
| Phase 3: BAM query & filtering | 3–4 hours |
| Phase 4: CIGAR expansion (pileup core) | 4–6 hours |
| Phase 5: Text rendering | 3–4 hours |
| Phase 6: MCP server wiring | 2–3 hours |
| Testing & debugging | 4–6 hours |
| **Total** | **~20–29 hours** |

The CIGAR expansion (Phase 4) is the highest-risk phase — edge cases around indels at region
boundaries and soft-clip handling tend to surface bugs. Budget extra time there.
Given your existing `dna-vcf-mcp` / `manta-vcf-mcp` codebase, the MCP wiring (Phase 6) should
be very fast — mostly pattern-matching what you've already built.
