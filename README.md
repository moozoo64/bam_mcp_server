# bam-pileup-mcp: Exposing DNA BAM Format files to LLMs for Analysis

> **Early Release** — this project is in early development and may have significant bugs,
> incomplete features, or breaking changes without notice. Use with appropriate caution.

A [Model Context Protocol](https://modelcontextprotocol.io/) (MCP) server written in Rust that accepts a genomic coordinate, queries a BAM file, and returns a text-based pileup representation suitable for reasoning by an LLM.

The server supports two transports:
- **stdio** (default) — for direct MCP client integration
- **HTTP** (`--sse`) — streamable-HTTP MCP transport served by axum, suitable for networked or multi-client use

## Example Output

```
Region: chr7:117,548,920-117,549,070  (150bp window)  BAM: sample.bam

POS   117548920                                            117549070
REF   TGCAATCCGAATCGGCATGCCTACGATTTACGGCATGCATCGAATCGGCAT[A]GCATCGGAATCGCATGCCTACGATTTACGGCATGCATCGAATCGGCATGCCT

R01   .............................................................................................................  + MQ60
R02   ...............................................................[G]............................................ + MQ60
R03          ..........................................................................................             + MQ59
R04   ...............................................................[G]............................................ - MQ60
R05                 ......................................................................                         - MQ58

COV   5x at chr7:117,548,975
      A: 2 (40.0%)  G: 3 (60.0%)
      +strand: 3  -strand: 2
      MQ_mean: 59.4  MQ_min: 58
      Reads shown: 5  Reads filtered (MAPQ): 0  Reads filtered (dup/qcfail): 0
```

## Requirements

- Sorted, indexed BAM file (`.bam` + `.bam.bai`)
- FASTA reference file with index (`.fa.gz` + `.fa.gz.fai` + `.fa.gz.gzi`, or plain `.fa` + `.fa.fai`)
- Rust 1.85+ (2024 edition)

## Building

```bash
cargo build --release
```

The binary is at `target/release/bam_mcp_server`.

## Usage

```
bam_mcp_server --bam <PATH> --reference <PATH> [OPTIONS]

Options:
  -b, --bam <PATH>           Path to sorted, indexed BAM file
  -r, --reference <PATH>     Path to FASTA reference (.fai index must exist alongside)
  -w, --window <INT>         Default window half-width in bp around query position [default: 75]
  -d, --max-depth <INT>      Maximum reads to display in pileup [default: 50]
  -q, --min-mapq <INT>       Minimum mapping quality filter [default: 0]
  -Q, --min-baseq <INT>      Minimum base quality to show as uppercase [default: 20]
      --sse <ADDR:PORT>       Run as HTTP MCP server instead of stdio (e.g. 127.0.0.1:8090)
      --allowed-host <HOST>  Allow this Host header value for remote access (repeatable; default: localhost/127.0.0.1/::1)
      --allow-all-hosts      Disable Host header checking entirely — allow any host (HTTP mode only)
      --debug                Write debug output to stderr
      --log-file <PATH>      Write debug log to file instead of stderr
  -h, --help                 Print help
  -V, --version              Print version
```

All debug/log output goes to **stderr**. The MCP JSON-RPC protocol uses **stdout** (stdio mode) or the HTTP response body (HTTP mode).

## MCP Tools

### `query_pileup`

**Input parameters:**

| Parameter     | Type    | Required | Description |
|---------------|---------|----------|-------------|
| `chrom`       | string  | ✓        | Chromosome name exactly as in the BAM header (e.g. `chr1`, `1`, `chrM`) |
| `position`    | integer | ✓        | 1-based genomic position to query |
| `window`      | integer |          | Half-width of region in bp (overrides `--window`) |
| `min_mapq`    | integer |          | Minimum mapping quality (overrides `--min-mapq`) |
| `show_strand` | boolean |          | Show strand orientation per read [default: true] |
| `show_mapq`   | boolean |          | Show mapping quality per read [default: true] |--allow-all-hosts

**Output:** A single `text` content block containing the pileup as a UTF-8 string.

---

### `get_header`

Returns the SAM-format header of the BAM file (`@HD`, `@SQ`, `@RG`, `@PG` lines), useful for discovering available chromosomes, read groups, and the programs used to generate the file.

**Input parameters:** none

**Output:** A `text` content block containing the SAM header as a UTF-8 string.

---

### `get_documentation`

Returns the full README for this server, including usage instructions, tool descriptions, and client configuration examples.

**Input parameters:** none

**Output:** A `text` content block containing the README as a UTF-8 string.

---

### Pileup symbols

| Symbol | Meaning |
|--------|---------|
| `.`    | Base matches reference |
| `A/C/G/T` (uppercase) | Mismatch, base quality ≥ min-baseq |
| `a/c/g/t` (lowercase) | Mismatch, base quality < min-baseq |
| `-`    | Deletion in read |
| `^N`   | N inserted bases follow this position |
| ` `    | Read does not span this position |
| `[X]`  | Queried position (regardless of match/mismatch) |

## MCP Client Configuration

### Claude Desktop — stdio (`claude_desktop_config.json`)

```json
{
  "mcpServers": {
    "bam-pileup": {
      "command": "/path/to/bam_mcp_server",
      "args": [
        "--bam", "/path/to/sample.bam",
        "--reference", "/path/to/reference.fa.gz"
      ]
    }
  }
}
```

### Claude Desktop — remote via `mcp-remote` (`claude_desktop_config.json`)

[`mcp-remote`](https://github.com/geelen/mcp-remote) acts as a local stdio proxy that forwards to a remote HTTP MCP server, bridging Claude Desktop's stdio-only support to an HTTP endpoint.

Start the server on the remote host:

```bash
# Allow a specific hostname (recommended)
bam_mcp_server --bam /path/to/sample.bam --reference /path/to/reference.fa.gz --sse 0.0.0.0:8090 --allowed-host <remote-host>

# Or disable host checking entirely (e.g. when the client sends an unexpected Host value)
bam_mcp_server --bam /path/to/sample.bam --reference /path/to/reference.fa.gz --sse 0.0.0.0:8090 --allow-all-hosts
```

Then configure Claude Desktop on your local machine:

```json
{
  "mcpServers": {
    "bam-pileup": {
      "command": "npx",
      "args": [
        "mcp-remote@latest",
        "http://<remote-host>:8090/mcp"
      ]
    }
  }
}
```

> **Use `mcp-remote@latest`:** This server uses the [streamable-HTTP MCP transport](https://modelcontextprotocol.io/docs/concepts/transports#streamable-http), which supersedes the older SSE transport. Older versions of `mcp-remote` only supported the legacy SSE transport and will fail to connect. Using `mcp-remote@latest` ensures you have support for the streamable-HTTP transport.

> **Security note:** The server has no built-in authentication. If exposed beyond localhost, place it behind a reverse proxy (e.g. nginx or Caddy) with TLS and authentication.

### VS Code — stdio (`.vscode/mcp.json`)

```json
{
  "servers": {
    "bam-pileup": {
      "type": "stdio",
      "command": "/path/to/bam_mcp_server",
      "args": [
        "--bam", "/path/to/sample.bam",
        "--reference", "/path/to/reference.fa.gz"
      ]
    }
  }
}
```

### VS Code — HTTP (`.vscode/mcp.json`)

Start the server with `--sse 127.0.0.1:8090`, then point VS Code at it:

```json
{
  "servers": {
    "bam-pileup": {
      "type": "sse",
      "url": "http://127.0.0.1:8090/mcp"
    }
  }
}
```

### Testing the HTTP transport with curl

The MCP protocol requires a 3-step handshake before sending requests. The server returns SSE-formatted responses (`text/event-stream`) — look for the `data: {…}` line for the JSON payload.

```bash
# Step 1: Initialize — capture the session ID from the response header
SESSION_ID=$(curl -s -D - -X POST http://localhost:8090/mcp \
  -H "Content-Type: application/json" \
  -H "Accept: application/json, text/event-stream" \
  -d '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{},"clientInfo":{"name":"curl","version":"0.1"}}}' \
  | grep -i "mcp-session-id" | awk '{print $2}' | tr -d '\r')

# Step 2: Send the initialized notification (server responds 202, no body)
curl -s -X POST http://localhost:8090/mcp \
  -H "Content-Type: application/json" \
  -H "Accept: application/json, text/event-stream" \
  -H "Mcp-Session-Id: $SESSION_ID" \
  -d '{"jsonrpc":"2.0","method":"notifications/initialized","params":{}}'

# Step 3: List available tools
curl -s -X POST http://localhost:8090/mcp \
  -H "Content-Type: application/json" \
  -H "Accept: application/json, text/event-stream" \
  -H "Mcp-Session-Id: $SESSION_ID" \
  -d '{"jsonrpc":"2.0","id":2,"method":"tools/list","params":{}}'
```

## Architecture

```
src/
├── main.rs       — Entry point: parse CLI, init logging, start MCP server
├── lib.rs        — Library root (exposes modules for integration tests)
├── cli.rs        — Clap CLI args and AppConfig validation
├── server.rs     — rmcp tool handler; wires the pipeline together
├── reference.rs  — BGZF/plain FASTA random-access reader (noodles)
├── query.rs      — BAM region query, filtering, reservoir sampling (noodles)
├── pileup.rs     — CIGAR expansion into reference-coordinate-aligned arrays
├── render.rs     — Text rendering: dots, brackets, footer
└── error.rs      — Error type re-exports (anyhow)
```

The pipeline per call:

1. Resolve `window` / `min_mapq` overrides from tool params vs CLI defaults
2. Compute 0-based `[region_start, region_end)` from the 1-based query position
3. Fetch reference sub-sequence (`ReferenceReader::fetch`)
4. Query BAM for overlapping reads with filtering (`query_region`)
5. Expand each read's CIGAR into per-reference-position slots (`expand_reads`)
6. Render to text (`render_pileup`)

## Testing

```bash
# Unit + HTTP + end-to-end integration tests (requires sample files)
cargo test

# Including sample-file unit tests that are marked #[ignore]
cargo test -- --include-ignored

# Individual suites
cargo test --test integration   # end-to-end pipeline tests
cargo test --test http          # HTTP transport tests
```

31 tests total: 17 unit tests (10 pure logic + 7 sample-file, the latter marked `#[ignore]`), 4 HTTP transport tests, and 3 end-to-end pipeline integration tests.

## License

MIT — see [LICENSE](LICENSE) for the full text.

## Credits

Vibe coded by [Michael Simmons](https://github.com/michael-simmons) using
[GitHub Copilot](https://github.com/features/copilot) powered by Claude Sonnet 4.6.
