use std::sync::Arc;

use axum::extract::Request;
use axum::response::Response;
use axum::{Router, middleware};
use bam_mcp_server::cli::{AppConfig, Args};
use bam_mcp_server::server::PileupServer;
use clap::Parser;
use rmcp::transport::streamable_http_server::{
    StreamableHttpServerConfig, StreamableHttpService,
    session::local::{LocalSessionManager, SessionConfig},
};
use rmcp::{ServiceExt, transport::stdio};
use tokio_util::sync::CancellationToken;

/// Resolves when the process receives Ctrl+C (SIGINT) or, on Unix, SIGTERM.
async fn shutdown_signal() {
    let ctrl_c = async {
        if let Err(e) = tokio::signal::ctrl_c().await {
            tracing::error!("failed to install Ctrl+C handler: {e}");
            std::future::pending::<()>().await;
        }
    };

    #[cfg(unix)]
    let terminate = async {
        match tokio::signal::unix::signal(tokio::signal::unix::SignalKind::terminate()) {
            Ok(mut sig) => {
                sig.recv().await;
            }
            Err(e) => {
                tracing::error!("failed to install SIGTERM handler: {e}");
                std::future::pending::<()>().await;
            }
        }
    };
    #[cfg(not(unix))]
    let terminate = std::future::pending::<()>();

    tokio::select! {
        _ = ctrl_c => tracing::info!("Received Ctrl+C, shutting down"),
        _ = terminate => tracing::info!("Received SIGTERM, shutting down"),
    }
}

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

fn main() -> anyhow::Result<()> {
    let runtime = tokio::runtime::Builder::new_multi_thread()
        .enable_all()
        .build()?;
    let result = runtime.block_on(run());
    // Tokio's stdin is read on a blocking thread that cannot be interrupted, so a
    // normal runtime drop would hang until the client sends more input. All async
    // cleanup has already finished by this point, so don't wait for it.
    runtime.shutdown_background();
    result
}

async fn run() -> anyhow::Result<()> {
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

        // Cancelling this token closes all open MCP sessions / SSE streams, which
        // lets axum's graceful shutdown finish instead of waiting on them forever.
        let shutdown_token = CancellationToken::new();

        let http_config = if config.allow_all_hosts {
            StreamableHttpServerConfig::default().disable_allowed_hosts()
        } else if config.allowed_hosts.is_empty() {
            StreamableHttpServerConfig::default()
        } else {
            StreamableHttpServerConfig::default()
                .with_allowed_hosts(config.allowed_hosts.iter().cloned())
        }
        .with_cancellation_token(shutdown_token.clone());

        let config_arc = Arc::new(config);

        let mut session_config = SessionConfig::default();
        session_config.keep_alive = Some(std::time::Duration::from_secs(36000)); // 10 hours
        let mut session_manager = LocalSessionManager::default();
        session_manager.session_config = session_config;

        let service = StreamableHttpService::new(
            {
                let config_arc = Arc::clone(&config_arc);
                move || Ok(PileupServer::from_arc(Arc::clone(&config_arc)))
            },
            Arc::new(session_manager),
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

        tokio::spawn({
            let shutdown_token = shutdown_token.clone();
            async move {
                shutdown_signal().await;
                eprintln!("Shutting down HTTP MCP server...");
                shutdown_token.cancel();
                // A second Ctrl+C forces an immediate exit if draining hangs.
                if tokio::signal::ctrl_c().await.is_ok() {
                    eprintln!("Forced shutdown");
                    std::process::exit(130);
                }
            }
        });

        axum::serve(listener, router)
            .with_graceful_shutdown(shutdown_token.cancelled_owned())
            .await?;
        tracing::info!("HTTP MCP server stopped");
    } else {
        // ── stdio transport (default MCP) ─────────────────────────────────────
        let server = PileupServer::new(config);
        // Install the signal handler before the handshake: `serve` blocks until
        // the client sends `initialize`, and Ctrl+C must work during that wait too.
        let shutdown = shutdown_signal();
        tokio::pin!(shutdown);
        let service = tokio::select! {
            res = server.serve(stdio()) => res?,
            _ = &mut shutdown => {
                tracing::info!("stdio MCP server stopped before initialization");
                return Ok(());
            }
        };
        let token = service.cancellation_token();
        let waiting = service.waiting();
        tokio::pin!(waiting);
        tokio::select! {
            res = &mut waiting => {
                res?;
            }
            _ = &mut shutdown => {
                token.cancel();
                // Let the service finish its cleanup after cancellation.
                waiting.await?;
            }
        }
        tracing::info!("stdio MCP server stopped");
    }

    Ok(())
}
