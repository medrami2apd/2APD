#!/usr/bin/env Rscript
source("scripts/00_setup/setup.R")
start_logging("11_draws_inference")

stopifnot(
  exists("propagate_draws"), exists("mcmc_diag_derived"), exists("export_table"))

if (!file.exists(CONFIG$PATHS$posterior_draws))
  stop("posterior_draws.rds not found. Run 02_estimation.R first (it stores per-draw paths).")

prep  <- readRDS(CONFIG$PATHS$processed)
conn  <- readRDS(CONFIG$PATHS$gfevd)
draws <- readRDS(CONFIG$PATHS$posterior_draws)
fit   <- if (file.exists(CONFIG$PATHS$fit)) readRDS(CONFIG$PATHS$fit) else NULL

k <- conn$k; p <- conn$p; sectors <- conn$sectors
nchains <- as.integer((fit$diagnostics$n_chains %||% CONFIG$N_CHAINS) %||% 1L)

prop <- propagate_draws(prep = prep, draws = draws, conn_dates = conn$dates,
                        sectors = sectors, k = k, p = p, cfg = CONFIG)

export_table(prop$h1,                "table_36_h1_tci_trend_draws")
export_table(prop$tci_bands,         "table_37_tci_credible_band")
export_table(prop$h3,                "table_38_h3_draws_ci")
export_table(prop$h4_episode,        "table_39_h4_break_posterior")
export_table(prop$h4_count,          "table_40_h4_break_count_posterior")
if (nrow(prop$h4_year_inclusion) > 0)
  export_table(prop$h4_year_inclusion, "table_41_h4_break_year_inclusion")

conv_derived <- tryCatch(
  mcmc_diag_derived(draws = draws, cfg = CONFIG, sectors = sectors,
                    k = k, p = p, nchains = nchains),
  error = function(e) tibble::tibble(Note = paste("derived diagnostics failed:", e$message)))
export_table(conv_derived, "table_42_convergence_derived")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  FIG_DIR <- CONFIG$PATHS$figures
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  pl <- ggplot2::ggplot(prop$tci_bands, ggplot2::aes(x = Date)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = TCI_lwr, ymax = TCI_upr),
                         alpha = 0.25, fill = "#2c7fb8") +
    ggplot2::geom_line(ggplot2::aes(y = TCI_med), color = "#08306b", linewidth = 0.7) +
    ggplot2::labs(title = "Total Connectedness Index with 95% posterior credible band",
                  subtitle = sprintf("Draws-mode propagation over %d posterior draws", prop$n_prop),
                  x = NULL, y = "TCI (%)") + ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(FIG_DIR, "figure_30_tci_credible_band.png"), pl,
                  width = 9, height = 5, dpi = 150, bg = "white")

  yi <- prop$h4_year_inclusion
  if (!is.null(yi) && nrow(yi) > 0) {
    ycol <- if ("Year" %in% names(yi)) "Year" else names(yi)[1]
    pcol <- if ("Post_prob" %in% names(yi)) "Post_prob" else
            if ("PostProb" %in% names(yi)) "PostProb" else
            names(yi)[vapply(yi, is.numeric, logical(1))][which(names(yi)[vapply(yi, is.numeric, logical(1))] != ycol)[1]]
    bd <- data.frame(Year = as.integer(yi[[ycol]]), Post_prob = as.numeric(yi[[pcol]]))
    episodes <- as.integer(CONFIG$PREREGISTERED_EPISODES)
    p41 <- ggplot2::ggplot(bd, ggplot2::aes(x = Year, y = Post_prob)) +
      ggplot2::geom_col(fill = "#2c7fb8", width = 0.8) +
      ggplot2::geom_vline(xintercept = episodes, linetype = "dashed",
                          color = "#c0392b", linewidth = 0.6) +
      ggplot2::labs(
        title = "Posterior density of structural break dates (H4)",
        subtitle = "Per-year posterior inclusion probability; dashed lines = pre-registered historical episodes",
        x = NULL, y = "Posterior inclusion probability") +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    ggplot2::ggsave(file.path(FIG_DIR, "figure_41_break_posterior_density.png"), p41,
                    width = 11, height = 5, dpi = 300, bg = "white")
    tryCatch(ggplot2::ggsave(file.path(FIG_DIR, "figure_41_break_posterior_density.pdf"), p41,
                             width = 11, height = 5, device = "pdf"), error = function(e) NULL)
    cat("  saved figure_41_break_posterior_density.png\n")
  }
}

dir.create(dirname(CONFIG$PATHS$draws_inference), recursive = TRUE, showWarnings = FALSE)
saveRDS(prop, CONFIG$PATHS$draws_inference)

cat("\n=== DRAWS-MODE INFERENCE (posterior uncertainty propagated) ===\n")
cat("\n--- H1 (TCI trend) ---\n"); print(prop$h1)
cat("\n--- H3 (centrality -> resilience, 95% credible intervals) ---\n"); print(prop$h3)
cat("\n--- H4 (pre-registered episode posterior inclusion) ---\n"); print(prop$h4_episode)
cat("\n--- H4 (posterior over number of breaks) ---\n"); print(prop$h4_count)
cat("\n--- Convergence on TCI_t / beta_t ---\n"); print(conv_derived)

end_logging()
