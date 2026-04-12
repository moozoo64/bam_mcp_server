# BAM Pileup MCP Server — Specification

## Overview

A Model Context Protocol (MCP) server written in Rust (2024 edition) that accepts a genomic
coordinate, queries a BAM file, and returns a text-based pileup representation suitable for
reasoning by an LLM. The server is stateless per call and requires a sorted, indexed BAM file
and an accompanying FASTA reference.

Two transports are supported:
- **stdio** (default) — standard MCP process wired over stdin/stdout
- **HTTP** (`--sse`) — streamable-HTTP MCP transport served by axum; responses are SSE
  (`text/event-stream`). Clients must include `Accept: application/json, text/event-stream`.

---

## 1. Command-Line Interface

```
bam_mcp_server [OPTIONS] --bam <PATH> --reference <PATH>

Options:
  -b, --bam <PATH>           Path to sorted, indexed BAM file (.bam + .bai must exist)
  -r, --reference <PATH>     Path to FASTA reference file (.fai index must exist alongside)
  -w, --window <INT>         Default window half-width in bp around query position [default: 75]
  -d, --max-depth <INT>      Maximum reads to display in pileup [default: 50]
  -q, --min-mapq <INT>       Minimum mapping quality filter [default: 0]
  -Q, --min-baseq <INT>      Minimum base quality to show as uppercase [default: 20]
      --sse <ADDR:PORT>       Run as HTTP MCP server instead of stdio (e.g. 127.0.0.1:8090)
      --allowed-host <HOST>  Allow this Host header value (repeatable; default: localhost/127.0.0.1/::1)
      --allow-all-hosts      Disable Host header checking entirely (HTTP mode only)
      --debug                Write debug output to stderr
      --log-file <PATH>      Write debug log to file instead of stderr
  -h, --help                 Print help
  -V, --version              Print version
```

### Notes
- `--debug` never writes to stdout — stdio MCP protocol uses stdout exclusively for JSON-RPC
- BAM index: both `<file>.bam.bai` and `<file>.bai` sidecars are accepted
- FASTA index: both `<file>.fa.gz.fai` and `<file>.fai` sidecars are accepted; BGZF `.fa.gz`
  requires an additional `.fa.gz.gzi` block index
- `--window` can be overridden per-call by the MCP tool arguments
- All paths are validated at startup; the process exits with a clear error if any are missing
- **Host filtering (HTTP mode):** by default only `localhost`, `127.0.0.1`, and `::1` are
  accepted via the `Host` header (DNS rebinding protection). Use `--allowed-host <HOST>`
  (repeatable) to permit additional values, or `--allow-all-hosts` to disable host checking
  entirely. The effective list is logged at DEBUG level (`--debug`).

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
      "type": ["integer", "null"],
      "description": "Half-width of region around position in bp. Overrides --window CLI default.",
      "minimum": 0
    },
    "min_mapq": {
      "type": ["integer", "null"],
      "description": "Minimum mapping quality. Overrides --min-mapq CLI default.",
      "minimum": 0,
      "maximum": 255
    },
    "show_strand": {
      "type": ["boolean", "null"],
      "description": "Show strand (+/-) column per read",
      "default": true
    },
    "show_mapq": {
      "type": ["boolean", "null"],
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
Region: chr7:117,548,920-117,549,070  (150bp window)  BAM: sample.bam
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
name    = "bam_mcp_server"
version = "0.1.0"
edition = "2024"
license = "MIT"

[dependencies]
noodles = { version = "0.109.0", features = ["bgzf", "core", "fasta", "bam", "sam"] }
rand = "0.10.1"
rmcp = { version = "1.4.0", features = [
  "server",
  "macros",
  "transport-io",
  "transport-streamable-http-server",
] }
tokio = { version = "1.51.1", features = ["full"] }
tokio-util = "0.7.18"
serde = { version = "1.0.228", features = ["derive"] }
serde_json = "1.0.149"
clap = { version = "4.6.0", features = ["derive"] }
hyper = { version = "1.9.0", features = ["server", "http1"] }
hyper-util = { version = "0.1.20", features = ["tokio", "server", "server-auto"] }
axum = "0.8.8"
anyhow = "1.0.102"
uuid = { version = "1.23.0", features = ["v4"] }
tracing = "0.1"
tracing-subscriber = { version = "0.3", features = ["env-filter"] }
schemars = "1"

[dev-dependencies]
criterion = "0.8.2"
tempfile = "3.27.0"
reqwest = { version = "0.12", default-features = false, features = ["json"] }
```

---

## 5. Module Structure

```
bam_mcp_server/
├── Cargo.toml
├── Cargo.lock
├── README.md
├── .github/
│   └── copilot-instructions.md
├── docs/
│   └── bam-pileup-mcp-spec.md   (this file)
├── benches/
│   └── bam_queries.rs
├── tests/
│   ├── integration.rs            # end-to-end pipeline tests
│   └── http.rs                   # HTTP transport tests
└── src/
    ├── main.rs       — Entry point: parse CLI, init logging, start stdio or HTTP server
    ├── lib.rs        — Library root; all modules are pub for integration test access
    ├── cli.rs        — Clap Args struct + AppConfig validation
    ├── server.rs     — PileupServer: rmcp #[tool_router], wires the pipeline
    ├── query.rs      — BAM region query, MAPQ/flag filtering, reservoir sampling
    ├── pileup.rs     — CIGAR expansion into reference-coordinate-aligned arrays
    ├── reference.rs  — BGZF/plain FASTA random-access via noodles IndexedReader
    ├── render.rs     — Text rendering: dots, brackets, footer, coverage stats
    └── error.rs      — anyhow re-exports
```

---

## 6. Implementation Details

### 6.1 CLI & Configuration (`cli.rs`)

`Args` is a `clap::Parser` derive struct. `AppConfig::try_from(Args)` validates at startup:
- BAM file exists; both `<file>.bam.bai` and `<file>.bai` sidecar paths are checked
- FASTA file exists; both `<file>.fai` and `<file>.<ext>.fai` checked
- All missing files produce a clear error before any MCP traffic

### 6.2 Entry point (`main.rs`)

`tracing_subscriber` is directed to stderr exclusively. Transport selection:

```rust
if let Some(ref addr) = config.sse {
    // HTTP: StreamableHttpService over axum, mount at /mcp
    let http_config = if config.allow_all_hosts {
        StreamableHttpServerConfig::default().disable_allowed_hosts()
    } else if config.allowed_hosts.is_empty() {
        StreamableHttpServerConfig::default()
    } else {
        StreamableHttpServerConfig::default()
            .with_allowed_hosts(config.allowed_hosts.iter().cloned())
    };
    let service = StreamableHttpService::new(
        move || Ok(PileupServer::from_arc(Arc::clone(&config_arc))),
        LocalSessionManager::default().into(),
        http_config,
    );
    let router = Router::new().nest_service("/mcp", service);
    axum::serve(listener, router).await?;
} else {
    // stdio: rmcp's built-in IO transport
    let service = server.serve(stdio()).await?;
    service.waiting().await?;
}
```

### 6.3 Reference Fetching (`reference.rs`)

```rust
pub struct ReferenceReader {
    inner: noodles::fasta::io::IndexedReader<File>,
}

impl ReferenceReader {
    pub fn open(fasta: &Path) -> anyhow::Result<Self>
    pub fn fetch(&mut self, chrom: &str, start: i64, end: i64) -> anyhow::Result<Vec<u8>>
}
```

- Uses `noodles::fasta::io::indexed_reader::Builder` for transparent BGZF / plain handling
- Coordinates are **0-based half-open** `[start, end)`; converted to noodles 1-based inclusive
  positions internally
- Returns raw uppercase ASCII bytes

### 6.4 BAM Query & Filtering (`query.rs`)

```rust
pub enum CigarOp { Match(u32), Ins(u32), Del(u32), SoftClip(u32), HardClip(u32), Skip(u32) }

pub struct ReadRecord {
    pub name: String,
    pub ref_start: i64,   // 0-based inclusive
    pub ref_end: i64,     // 0-based exclusive
    pub cigar: Vec<CigarOp>,
    pub sequence: Vec<u8>,
    pub base_quals: Vec<u8>,
    pub is_reverse: bool,
    pub mapq: u8,
    pub is_dup: bool,
    pub is_qcfail: bool,
}

pub struct QueryResult {
    pub reads: Vec<ReadRecord>,
    pub filtered_mapq: usize,
    pub filtered_flags: usize,
}

pub fn query_region(
    bam_path: &Path, chrom: &str,
    start: i64, end: i64,   // 0-based half-open
    min_mapq: u8, max_depth: usize,
) -> anyhow::Result<QueryResult>
```

Filtering decision tree (applied in order):
1. Unmapped / secondary / supplementary → skipped silently (not counted)
2. `mapq < min_mapq` → counted in `filtered_mapq`, excluded from output
3. Duplicate or QC-fail → counted in `filtered_flags`, excluded from output
4. If eligible reads exceed `max_depth`, a **reservoir sample** of exactly `max_depth` is kept

All three "match" CIGAR op codes (M / = / X) are collapsed to `CigarOp::Match` — the pileup
engine does its own base comparison against the reference.

### 6.5 CIGAR Expansion (`pileup.rs`)

```rust
pub enum PileupBase {
    Match,
    Mismatch { base: u8, qual: u8 },
    Deletion,
}

pub struct AlignedRead {
    pub record: ReadRecord,
    pub bases: Vec<Option<PileupBase>>,  // one slot per ref position in [region_start, region_end)
    pub ins_after: Vec<u32>,             // insertion count after position i; parallel to bases
}

pub fn expand_reads(
    reads: &[ReadRecord], ref_seq: &[u8],
    region_start: i64, region_end: i64,
) -> Vec<AlignedRead>
```

CIGAR expansion algorithm:

```
ref_pos  = record.ref_start
read_pos = 0

for each cigar op:
  Match(n):
    compare read base to ref base; emit Match or Mismatch at ref_pos
    ref_pos += 1; read_pos += 1  (×n)

  Del(n):
    emit Deletion at ref_pos; ref_pos += 1  (×n)

  Ins(n):
    set ins_after[ref_pos - 1] = n (clamped to window); read_pos += n

  SoftClip(n):
    read_pos += n  (ref_pos unchanged)

  Skip(n) | HardClip(n):
    ref_pos += n  (or ignored)
```

Slots outside `[region_start, region_end)` are ignored. `None` means the read does not span
that reference position.

### 6.6 Text Rendering (`render.rs`)

```rust
pub struct RenderOpts {
    pub chrom: String,
    pub bam_name: String,
    pub show_strand: bool,
    pub show_mapq: bool,
    pub min_baseq: u8,
}

pub fn render_pileup(
    reads: &[AlignedRead], ref_seq: &[u8],
    region_start: i64, query_pos: i64,
    query_result: &QueryResult, opts: &RenderOpts,
) -> String
```

Rendering steps:
1. **Header** — `Region: chrom:start-end  (Nbp window)  BAM: filename`
2. **POS ruler** — start/end coords, `[*]` marker at query position
3. **REF line** — reference bases; query position bracketed as `[X]`
4. *(blank line)*
5. **Read lines** — `R01 … RNN` rows, one per `AlignedRead`:
   - `None` slot → space
   - `Match` → `.`
   - `Mismatch`, qual ≥ `min_baseq` → uppercase base
   - `Mismatch`, qual < `min_baseq` → lowercase base
   - `Deletion` → `-`
   - `ins_after[i] > 0` → `^N` appended to that slot's symbol
   - Query position → `[X]` brackets wrap the slot
   - Suffix: `  + MQ60` (strand and/or MAPQ if enabled)
6. *(blank line)*
7. **Footer** — coverage, allele counts at query position, strand balance, MQ stats, filter counts

### 6.7 MCP Server (`server.rs`)

`PileupServer` is a `Clone + Send + Sync` struct wrapping `Arc<AppConfig>`. The `#[tool_router]`
macro wires `query_pileup` into the MCP `tools/list` and `tools/call` handlers. Tool-level
errors are returned as plain-text content strings, not JSON-RPC protocol errors.

Per-call pipeline:
1. Resolve `window` / `min_mapq` overrides (tool params take precedence over CLI defaults)
2. Convert 1-based input position to 0-based `query_pos_0`; compute `[region_start, region_end)`
3. `ReferenceReader::fetch` → `ref_seq: Vec<u8>`
4. `query_region` → `QueryResult`
5. `expand_reads` → `Vec<AlignedRead>`
6. `render_pileup` → `String` returned as MCP text content

---

## 7. Error Handling Strategy

| Error situation                        | Response                                                    |
|----------------------------------------|-------------------------------------------------------------|
| Chromosome not in BAM header           | Tool content: `"Error: Unknown chromosome: chrX"`           |
| BAM / FASTA file not found at startup  | Process exits with clear message before MCP traffic begins  |
| Zero reads in region                   | Valid pileup: REF line + empty read section + `0x` footer   |
| CIGAR parse failure on a read          | `tracing::error!` to stderr; read is skipped, others shown  |
| Invalid 0-based coordinate             | `anyhow::bail!` propagated as tool-level error text         |

All errors inside `do_query_pileup` are `anyhow::Result` propagated to the `query_pileup`
handler which returns them as plain-text content rather than JSON-RPC protocol errors.

---

## 8. Debug Logging

With `--debug` or `RUST_LOG=debug` (all output to stderr via `tracing`):

```
[DEBUG] query_pileup called: chr7:117548975 window=75 min_mapq=0
[DEBUG] Region: chr7:117548900-117549050
[DEBUG] Reference fetched: 150bp starting TGCAATCC...
[DEBUG] BAM query returned 49 records; filtered MAPQ=2 flags=0
[DEBUG] CIGAR expansion: 47 reads in 1.2ms
[DEBUG] Render complete: 54 lines, 8240 chars
```

---

## 9. Testing

### Test suites

| Suite | File | Count | Requires samples/ |
|-------|------|-------|-------------------|
| Unit tests | `src/**` | 17 (7 `#[ignore]`) | 7 ignored tests only |
| HTTP transport | `tests/http.rs` | 4 | Yes |
| End-to-end pipeline | `tests/integration.rs` | 3 | Yes |

### Running tests

```bash
cargo test                             # unit + HTTP + end-to-end
cargo test -- --include-ignored        # also runs #[ignore]-gated sample-file unit tests
cargo test --test integration          # end-to-end pipeline only
cargo test --test http                 # HTTP transport only
```

### HTTP test approach

`tests/http.rs` spins up a real axum server on a random port (`TcpListener::bind("127.0.0.1:0")`)
for each test — no mocking. The server requires `Accept: application/json, text/event-stream`;
all responses are SSE. The `parse_sse_json()` helper in that file extracts the JSON-RPC payload
from the `data:` lines.

Tests cover:
- `test_initialize_response` — MCP handshake, session ID header, protocol version
- `test_tools_list_contains_query_pileup` — tool listing with correct schema
- `test_tools_call_query_pileup_valid` — real pileup output structure
- `test_tools_call_query_pileup_invalid_chrom` — tool-level error message format

### Manual smoke test (stdio)

```bash
bam_mcp_server --bam sample.bam --reference GRCh38.fa.gz --debug
# In another terminal, via raw JSON-RPC on stdin:
{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"query_pileup","arguments":{"chrom":"chr7","position":117548975}}}
```

### Manual smoke test (HTTP)

```bash
bam_mcp_server --bam sample.bam --reference GRCh38.fa.gz --sse 127.0.0.1:8090
curl -X POST http://127.0.0.1:8090/mcp \
  -H 'Accept: application/json, text/event-stream' \
  -H 'Content-Type: application/json' \
  -d '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{},"clientInfo":{"name":"test","version":"0.1"}}}'
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


