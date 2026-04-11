# bam-pileup-mcp — Copilot Instructions

## Project Overview

A Rust MCP server that queries BAM/CRAM files at genomic coordinates and returns
text pileups for LLM reasoning. Supports stdio and HTTP (streamable-HTTP via axum)
transports.

## Build & Test

```bash
cargo build --release                  # binary → target/release/bam_mcp_server
cargo test                             # unit + HTTP + end-to-end (requires samples/)
cargo test -- --include-ignored        # also runs #[ignore]-gated sample-file unit tests
cargo test --test integration          # end-to-end pipeline only
cargo test --test http                 # HTTP transport tests only
```

Sample files live in `samples/` (git-ignored). Tests that require them are in
`tests/integration.rs` and `tests/http.rs`.

## Architecture

```
src/
├── main.rs       — CLI entry point; starts stdio or HTTP transport
├── lib.rs        — Library root; all modules are pub so tests can reach them
├── cli.rs        — Clap args + AppConfig validation
├── server.rs     — rmcp tool handler (PileupServer); wires the pipeline
├── reference.rs  — BGZF/plain FASTA random-access via noodles
├── query.rs      — BAM region query, MAPQ/flag filtering, reservoir sampling
├── pileup.rs     — CIGAR expansion into reference-coordinate-aligned arrays
├── render.rs     — Text renderer (dots, brackets, coverage footer)
└── error.rs      — anyhow re-exports
```

Per-call pipeline: `ReferenceReader::fetch` → `query_region` → `expand_reads` → `render_pileup`

## Key Conventions

- **Coordinates**: all internal coords are 0-based half-open `[start, end)`;
  the MCP tool input is 1-based and converted immediately in `server.rs`.
- **Error handling**: `anyhow::Result` throughout; tool-level errors are
  returned as plain-text `String` content (not JSON-RPC errors) by `PileupServer`.
- **HTTP transport**: uses rmcp's `StreamableHttpService` over axum. Responses
  are always `text/event-stream` SSE — clients must include
  `Accept: application/json, text/event-stream`. Notifications return `202 Accepted` with no body.
- **No `mod server;` in main.rs**: `server.rs` lives in the library crate
  (`lib.rs` exposes it as `pub mod server`), so integration tests can construct
  `PileupServer` directly.
- **Edition 2024**: use `use` import conventions and `async fn` in traits as
  stabilised in Rust 2024.

## Dependencies to Know

| Crate           | Role                                                       |
| --------------- | ---------------------------------------------------------- |
| `noodles`       | BAM/FASTA/BGZF I/O (features: bgzf, bam, fasta, sam, core) |
| `rmcp`          | MCP protocol + `#[tool]` / `#[tool_router]` macros         |
| `axum`          | HTTP server for `--sse` transport                          |
| `schemars`      | JSON Schema generation for tool input structs              |
| `anyhow`        | Error propagation                                          |
| `clap`          | CLI arg parsing (derive feature)                           |
| `reqwest` (dev) | HTTP client used in `tests/http.rs`                        |

## Testing Notes

- `tests/http.rs` spins up a real axum server on a random port per test using
  `TcpListener::bind("127.0.0.1:0")` — no mocking needed.
- SSE response parsing is handled by `parse_sse_json()` in that file; reuse it
  for new HTTP tests rather than duplicating the logic.
- Integration and HTTP tests all require `samples/` to be present.
