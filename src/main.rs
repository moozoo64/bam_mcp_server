use std::sync::Arc;

use axum::Router;
use bam_mcp_server::cli::{AppConfig, Args};
use bam_mcp_server::server::PileupServer;
use clap::Parser;
use rmcp::transport::streamable_http_server::{
    StreamableHttpServerConfig, StreamableHttpService, session::local::LocalSessionManager,
};
use rmcp::{ServiceExt, transport::stdio};

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
        let config_arc = Arc::new(config);

        let service = StreamableHttpService::new(
            move || Ok(PileupServer::from_arc(Arc::clone(&config_arc))),
            LocalSessionManager::default().into(),
            StreamableHttpServerConfig::default(),
        );

        let router = Router::new().nest_service("/mcp", service);
        let listener = tokio::net::TcpListener::bind(&addr).await?;
        let bound = listener.local_addr()?;
        tracing::info!("HTTP MCP server listening on http://{}/mcp", bound);
        eprintln!("HTTP MCP server listening on http://{}/mcp", bound);

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
