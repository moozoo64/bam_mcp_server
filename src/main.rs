use std::sync::Arc;

use axum::extract::Request;
use axum::response::Response;
use axum::{Router, middleware};
use bam_mcp_server::cli::{AppConfig, Args};
use bam_mcp_server::server::PileupServer;
use clap::Parser;
use rmcp::transport::streamable_http_server::{
    StreamableHttpServerConfig, StreamableHttpService, session::local::LocalSessionManager,
};
use rmcp::{ServiceExt, transport::stdio};

/// Axum middleware that logs the incoming Host header and request line at DEBUG level.
/// This makes it easy to diagnose `403 Forbidden: Host header is not allowed` errors
/// by showing exactly what Host value the client is sending.
async fn log_host_header(req: Request, next: middleware::Next) -> Response {
    let host = req
        .headers()
        .get("host")
        .and_then(|v| v.to_str().ok())
        .unwrap_or("<missing>");
    tracing::debug!(
        method = %req.method(),
        path = %req.uri().path(),
        host = %host,
        "incoming HTTP request"
    );
    next.run(req).await
}

fn init_logging(args: &Args) {
    use tracing_subscriber::{EnvFilter, fmt};

    let builder = fmt::Subscriber::builder()
        .with_env_filter(EnvFilter::try_from_default_env().unwrap_or_else(|_| {
            if args.debug {
                EnvFilter::new("debug")
            } else {
                EnvFilter::new("warn")
            }
        }))
        // MCP protocol uses stdout exclusively — all logs must go to stderr
        .with_writer(std::io::stderr);

    if let Some(log_path) = &args.log_file {
        match std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(log_path)
        {
            Ok(file) => {
                builder
                    .with_writer(move || -> Box<dyn std::io::Write> {
                        Box::new(file.try_clone().expect("clone log file handle"))
                    })
                    .init();
            }
            Err(e) => {
                // Fall back to stderr if log file can't be opened
                builder.init();
                eprintln!(
                    "Warning: could not open log file {}: {e}",
                    log_path.display()
                );
            }
        }
    } else {
        builder.init();
    }
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    let args = Args::parse();
    init_logging(&args);

    let config = AppConfig::try_from(args)?;

    tracing::info!(
        bam = %config.bam_path.display(),
        reference = %config.reference_path.display(),
        window = config.window,
        max_depth = config.max_depth,
        "Starting bam-pileup-mcp server"
    );

    if let Some(ref addr) = config.sse {
        // ── HTTP / Streamable-HTTP transport ─────────────────────────────────
        let addr = addr.clone();

        let http_config = if config.allow_all_hosts {
            StreamableHttpServerConfig::default().disable_allowed_hosts()
        } else if config.allowed_hosts.is_empty() {
            StreamableHttpServerConfig::default()
        } else {
            StreamableHttpServerConfig::default()
                .with_allowed_hosts(config.allowed_hosts.iter().cloned())
        };

        let config_arc = Arc::new(config);

        let service = StreamableHttpService::new(
            {
                let config_arc = Arc::clone(&config_arc);
                move || Ok(PileupServer::from_arc(Arc::clone(&config_arc)))
            },
            LocalSessionManager::default().into(),
            http_config,
        );

        let router = Router::new()
            .nest_service("/mcp", service)
            .layer(middleware::from_fn(log_host_header));
        let listener = tokio::net::TcpListener::bind(&addr).await?;
        let bound = listener.local_addr()?;
        // Always print the effective allowed-hosts list to stderr so the user can
        // diagnose "403 Forbidden: Host header is not allowed" without --debug.
        let effective_hosts = if config_arc.allow_all_hosts {
            "ALL (host checking disabled)".to_string()
        } else if config_arc.allowed_hosts.is_empty() {
            "localhost, 127.0.0.1, ::1 (default)".to_string()
        } else {
            config_arc.allowed_hosts.join(", ")
        };
        tracing::info!("HTTP MCP server listening on http://{}/mcp", bound);
        eprintln!("HTTP MCP server listening on http://{}/mcp", bound);
        tracing::debug!("Allowed Host headers: {effective_hosts}");
        tracing::debug!(
            "Use --allowed-host <HOST> to allow additional hosts (repeatable), or --allow-all-hosts to disable host checking."
        );

        axum::serve(listener, router)
            .with_graceful_shutdown(async {
                tokio::signal::ctrl_c()
                    .await
                    .expect("failed to install Ctrl+C handler");
                tracing::info!("Shutting down");
            })
            .await?;
    } else {
        // ── stdio transport (default MCP) ─────────────────────────────────────
        let server = PileupServer::new(config);
        let service = server.serve(stdio()).await?;
        service.waiting().await?;
    }

    Ok(())
}
