//! HTTP transport integration tests for the MCP server.
//!
//! These tests spin up the axum-based streamable-HTTP MCP server on a random
//! port, perform the MCP JSON-RPC handshake, and verify tool listing and
//! invocation.  They require the sample files in `samples/` — the same
//! prerequisite as the pipeline integration tests.

use std::net::SocketAddr;
use std::path::PathBuf;
use std::sync::Arc;

use axum::Router;
use bam_mcp_server::cli::AppConfig;
use bam_mcp_server::server::PileupServer;
use rmcp::transport::streamable_http_server::{
    StreamableHttpServerConfig, StreamableHttpService, session::local::LocalSessionManager,
};
use serde_json::{Value, json};

// ── helpers ─────────────────────────────────────────────────────────────────

fn make_config() -> AppConfig {
    AppConfig {
        bam_path: PathBuf::from("samples/SampleHuman-30x-WGS.bam"),
        reference_path: PathBuf::from("samples/GRCh38.primary_assembly.genome.fa.gz"),
        window: 75,
        max_depth: 50,
        min_mapq: 0,
        min_baseq: 20,
        sse: None,
        allowed_hosts: vec![],
        debug: false,
        log_file: None,
    }
}

/// Bind to an OS-assigned port and serve the MCP HTTP service in a background
/// task.  Returns the bound socket address.
async fn start_server() -> SocketAddr {
    let config = Arc::new(make_config());

    let service = StreamableHttpService::new(
        move || Ok(PileupServer::from_arc(Arc::clone(&config))),
        LocalSessionManager::default().into(),
        StreamableHttpServerConfig::default(),
    );

    let router = Router::new().nest_service("/mcp", service);
    let listener = tokio::net::TcpListener::bind("127.0.0.1:0")
        .await
        .expect("bind to random port");
    let addr = listener.local_addr().expect("local_addr");

    tokio::spawn(async move {
        axum::serve(listener, router).await.expect("server error");
    });

    addr
}

/// Extract the JSON-RPC payload from an SSE response body.
///
/// The streamable-HTTP transport always returns `text/event-stream`.  Each
/// response event looks like:
///
/// ```text
/// data:
/// id: 0
/// retry: 3000
///
/// data: {"jsonrpc":"2.0","id":1,"result":{…}}
/// id: 0/1
/// ```
///
/// We scan lines for `data: {…}` and return the first parseable JSON object.
fn parse_sse_json(body: &str) -> Value {
    for line in body.lines() {
        if let Some(payload) = line.strip_prefix("data:") {
            let payload = payload.trim();
            if payload.starts_with('{') {
                return serde_json::from_str(payload).expect("SSE data line must be valid JSON");
            }
        }
    }
    panic!("no JSON object found in SSE body:\n{body}");
}

/// Manages a single MCP session: tracks the session ID and auto-increments
/// JSON-RPC message IDs.
struct McpSession {
    client: reqwest::Client,
    url: String,
    session_id: Option<String>,
    next_id: u32,
}

impl McpSession {
    /// Connect to the server and complete the MCP initialize handshake.
    async fn connect(addr: SocketAddr) -> Self {
        let client = reqwest::Client::new();
        let url = format!("http://{}/mcp", addr);
        let mut s = Self {
            client,
            url,
            session_id: None,
            next_id: 1,
        };
        s.handshake().await;
        s
    }

    fn take_id(&mut self) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        id
    }

    /// Send a JSON-RPC request and return the parsed response body.
    async fn request(&mut self, body: Value) -> Value {
        let resp = self.send_raw(&body).await;
        let sid: Option<String> = resp
            .headers()
            .get("Mcp-Session-Id")
            .and_then(|v| v.to_str().ok())
            .map(|s| s.to_owned());
        if let Some(s) = sid {
            self.session_id = Some(s);
        }
        let text = resp.text().await.expect("read response body");
        parse_sse_json(&text)
    }

    /// Send a JSON-RPC notification (no response body expected, server returns
    /// 202 Accepted with an empty body).
    async fn notify(&mut self, body: Value) {
        let resp = self.send_raw(&body).await;
        assert_eq!(
            resp.status().as_u16(),
            202,
            "notification must be acknowledged with 202"
        );
    }

    async fn send_raw(&self, body: &Value) -> reqwest::Response {
        let mut b = self
            .client
            .post(&self.url)
            // The streamable-HTTP transport requires SSE to be acceptable.
            .header("Accept", "application/json, text/event-stream")
            .header("Content-Type", "application/json");
        if let Some(ref sid) = self.session_id {
            b = b.header("Mcp-Session-Id", sid);
        }
        b.json(body).send().await.expect("HTTP send failed")
    }

    /// Perform the MCP initialize + initialized handshake.
    async fn handshake(&mut self) {
        let id = self.take_id();
        let resp = self
            .request(json!({
                "jsonrpc": "2.0",
                "id": id,
                "method": "initialize",
                "params": {
                    "protocolVersion": "2024-11-05",
                    "capabilities": {},
                    "clientInfo": {"name": "test-client", "version": "0.1.0"}
                }
            }))
            .await;
        assert_eq!(
            resp["jsonrpc"], "2.0",
            "initialize response must be JSON-RPC 2.0"
        );
        assert!(
            resp["result"]["protocolVersion"].is_string(),
            "initialize result must contain protocolVersion"
        );

        // The initialized notification has no id and no response body.
        self.notify(json!({
            "jsonrpc": "2.0",
            "method": "notifications/initialized",
            "params": {}
        }))
        .await;
    }
}

// ── tests ────────────────────────────────────────────────────────────────────

/// The initialize response must carry the protocol version and server info,
/// and include an `Mcp-Session-Id` response header.
#[tokio::test]
async fn test_initialize_response() {
    let addr = start_server().await;
    let client = reqwest::Client::new();
    let url = format!("http://{}/mcp", addr);

    let resp = client
        .post(&url)
        .header("Accept", "application/json, text/event-stream")
        .header("Content-Type", "application/json")
        .json(&json!({
            "jsonrpc": "2.0",
            "id": 1,
            "method": "initialize",
            "params": {
                "protocolVersion": "2024-11-05",
                "capabilities": {},
                "clientInfo": {"name": "test-client", "version": "0.1.0"}
            }
        }))
        .send()
        .await
        .expect("send initialize");

    assert!(resp.status().is_success(), "initialize should return 2xx");
    assert!(
        resp.headers().contains_key("mcp-session-id"),
        "response must include Mcp-Session-Id header"
    );

    let body = parse_sse_json(&resp.text().await.expect("read body"));
    assert_eq!(body["jsonrpc"], "2.0");
    assert_eq!(body["id"], 1);
    assert!(
        body["result"]["protocolVersion"].is_string(),
        "result.protocolVersion must be a string"
    );
    assert!(
        body["result"]["serverInfo"]["name"].is_string(),
        "result.serverInfo.name must be a string"
    );
}

/// After initialization, `tools/list` must return an array that includes the
/// `query_pileup` tool with the expected input schema.
#[tokio::test]
async fn test_tools_list_contains_query_pileup() {
    let addr = start_server().await;
    let mut session = McpSession::connect(addr).await;

    let id = session.take_id();
    let resp = session
        .request(json!({
            "jsonrpc": "2.0",
            "id": id,
            "method": "tools/list",
            "params": {}
        }))
        .await;

    assert_eq!(resp["jsonrpc"], "2.0");
    let tools = resp["result"]["tools"]
        .as_array()
        .expect("result.tools must be an array");

    let tool = tools
        .iter()
        .find(|t| t["name"] == "query_pileup")
        .expect("query_pileup tool must be listed");

    assert!(
        tool["description"].is_string(),
        "query_pileup must have a description"
    );
    assert!(
        tool["inputSchema"].is_object(),
        "query_pileup must have an inputSchema"
    );
}

/// Calling `query_pileup` with valid coordinates returns a pileup that
/// contains the standard section headers produced by the render layer.
#[tokio::test]
async fn test_tools_call_query_pileup_valid() {
    let addr = start_server().await;
    let mut session = McpSession::connect(addr).await;

    let id = session.take_id();
    let resp = session
        .request(json!({
            "jsonrpc": "2.0",
            "id": id,
            "method": "tools/call",
            "params": {
                "name": "query_pileup",
                "arguments": {
                    "chrom": "chr1",
                    "position": 10076
                }
            }
        }))
        .await;

    assert_eq!(resp["jsonrpc"], "2.0");
    let content = resp["result"]["content"]
        .as_array()
        .expect("result.content must be an array");
    assert!(!content.is_empty(), "content must not be empty");

    let text = content[0]["text"]
        .as_str()
        .expect("first content item must have a text field");

    assert!(
        text.contains("Region:"),
        "pileup output must contain 'Region:'"
    );
    assert!(
        text.contains("REF"),
        "pileup output must contain the REF line"
    );
    assert!(
        text.contains("COV"),
        "pileup output must contain the COV footer"
    );
}

/// Calling `query_pileup` with an unknown chromosome must return a tool-level
/// error message (not a JSON-RPC protocol error).
#[tokio::test]
async fn test_tools_call_query_pileup_invalid_chrom() {
    let addr = start_server().await;
    let mut session = McpSession::connect(addr).await;

    let id = session.take_id();
    let resp = session
        .request(json!({
            "jsonrpc": "2.0",
            "id": id,
            "method": "tools/call",
            "params": {
                "name": "query_pileup",
                "arguments": {
                    "chrom": "nonexistent_chrXYZ99",
                    "position": 1000
                }
            }
        }))
        .await;

    assert_eq!(
        resp["jsonrpc"], "2.0",
        "must still be a valid JSON-RPC response"
    );
    // The server wraps tool errors as content text rather than a protocol error.
    let content = resp["result"]["content"]
        .as_array()
        .expect("result.content must be an array");
    let text = content[0]["text"]
        .as_str()
        .expect("first content item must have text");
    assert!(
        text.contains("Error"),
        "error text must mention 'Error', got: {text}"
    );
}
