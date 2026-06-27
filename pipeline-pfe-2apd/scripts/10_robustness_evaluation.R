#!/usr/bin/env Rscript
source("scripts/00_setup/setup.R")
start_logging("10_robustness_evaluation")

suppressWarnings(suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble)
}))

stopifnot(
  exists("robustness_checklist"), exists("mcmc_diag"), exists("tv_sv_relevance"),
  exists("model_comparison"), exists("build_xy"), exists("ols_var"),
  exists("predict_metrics"), exists("waic_dic"), exists("rolling_vs_tvp_tci"),
  exists("export_table"))

cfg   <- CONFIG$ROBUSTNESS %||% list()
prep  <- readRDS(CONFIG$PATHS$processed)
post  <- readRDS(CONFIG$PATHS$posterior)
conn  <- readRDS(CONFIG$PATHS$gfevd)
fit   <- tryCatch(readRDS(CONFIG$PATHS$fit),             error = function(e) NULL)
draws <- tryCatch(readRDS(CONFIG$PATHS$posterior_draws), error = function(e) NULL)

FIG_DIR <- CONFIG$PATHS$figures
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
LINE_COL <- "#2c3e50"; FILL_COL <- "#3498db"; ALT_COL <- "#c0392b"

save_eval_plot <- function(p, name, w = 11, h = 6, dpi = 300) {
  ggplot2::ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = dpi, bg = "white")
  tryCatch(ggplot2::ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf"),
           error = function(e) NULL)
  cat(sprintf("  saved %s.png\n", name))
}

results <- list()

chk <- robustness_checklist()
export_table(chk, "table_30_robustness_checklist"); results$checklist <- chk
cat(sprintf("Checklist: %d items triaged for this model.\n", nrow(chk)))

if (!is.null(fit)) {
  diag_tbl <- mcmc_diag(fit, post$sectors, cfg)
  export_table(diag_tbl, "table_31_mcmc_diagnostics"); results$mcmc <- diag_tbl
  if ("Rhat_pass" %in% names(diag_tbl))
    cat(sprintf("MCMC diagnostics: %d/%d params pass Rhat < %.3f.\n",
                sum(diag_tbl$Rhat_pass, na.rm = TRUE), nrow(diag_tbl), cfg$RHAT_THRESHOLD %||% 1.01))
} else cat("No fit object; skipped MCMC diagnostics.\n")

rel <- tv_sv_relevance(post, post$sectors)
export_table(rel, "table_32_tvp_sv_relevance"); results$relevance <- rel

mc <- model_comparison(prep, post, cfg)
base_rmse <- mc$table$RMSE[mc$table$Model == "Constant BVAR (std)"]
base_lpds <- mc$table$LPDS[mc$table$Model == "Constant BVAR (std)"]
mc$table$RMSE_gain_vs_std_pct <- round((base_rmse - mc$table$RMSE) / base_rmse * 100, 2)
mc$table$LPDS_gain_vs_std     <- round(mc$table$LPDS - base_lpds, 3)
export_table(mc$table, "table_33_model_comparison"); results$model_comparison <- mc$table
print(mc$table)

wd <- tryCatch(waic_dic(prep, post, draws, cfg), error = function(e) { warning(conditionMessage(e)); NULL })
if (!is.null(wd)) {
  wd_tbl <- tibble::tibble(
    Criterion = c("WAIC", "p_WAIC", "lppd", "DIC", "pD", "LOOIC", "n_draws_used"),
    Value     = c(wd$waic, wd$p_waic, wd$lppd, wd$dic, wd$pD, wd$looic, wd$n_draws))
  export_table(wd_tbl, "table_34_information_criteria"); results$information_criteria <- wd_tbl
  cat(sprintf("WAIC = %.1f | DIC = %.1f (full TVP-BVAR-SV)\n", wd$waic, wd$dic))
} else cat("No posterior draws; skipped WAIC/DIC/LOO.\n")

rt <- tryCatch(rolling_vs_tvp_tci(prep, conn, post, cfg),
               error = function(e) { warning(conditionMessage(e)); NULL })
if (!is.null(rt)) {
  export_table(rt, "table_35_rolling_vs_tvp_tci"); results$rolling <- rt
  cc <- suppressWarnings(stats::cor(rt$TVP_TCI, rt$Rolling_TCI, use = "complete.obs"))
  cat(sprintf("Rolling-window vs TVP TCI correlation: %.3f\n", cc))
}

if (!is.null(mc$pit)) {
  pit_df <- tibble::tibble(PIT = mc$pit[is.finite(mc$pit)])
  nb <- as.integer(cfg$PIT_BINS %||% 10L)
  p_pit <- ggplot2::ggplot(pit_df, ggplot2::aes(PIT)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = nb, fill = FILL_COL, color = "white", boundary = 0) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = ALT_COL, linewidth = 0.9) +
    ggplot2::labs(title = "PIT - TVP-BVAR-SV one-step predictive calibration",
                  subtitle = "Flat at density 1 => well-calibrated predictive density",
                  x = "Probability integral transform", y = "Density") +
    ggplot2::theme_minimal(base_size = 13)
  save_eval_plot(p_pit, "figure_eval_01_pit_full", 9, 6)
}
if (!is.null(rt) && any(is.finite(rt$Rolling_TCI))) {
  rl <- tidyr::pivot_longer(rt, c(TVP_TCI, Rolling_TCI), names_to = "Estimator", values_to = "TCI")
  p_rt <- ggplot2::ggplot(rl, ggplot2::aes(Date, TCI, color = Estimator)) +
    ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(TVP_TCI = LINE_COL, Rolling_TCI = ALT_COL),
                                labels = c(TVP_TCI = "TVP-VAR", Rolling_TCI = "Rolling-window VAR")) +
    ggplot2::labs(title = "Robustness: rolling-window VAR vs TVP-VAR connectedness",
                  subtitle = sprintf("Rolling window = %d periods", as.integer(cfg$ROLLING_WINDOW %||% 40L)),
                  x = NULL, y = "Total connectedness (%)", color = NULL) +
    ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(legend.position = "bottom")
  save_eval_plot(p_rt, "figure_eval_02_rolling_vs_tvp", 12, 6)
}

if (!is.null(fit) && !is.null(fit$sv_params)) {
  nchains <- as.integer((fit$diagnostics$n_chains %||% CONFIG$N_CHAINS) %||% 1L)
  if (nchains < 1L) nchains <- 1L
  sectors <- post$sectors; k <- post$k

  avg_sigma <- colMeans(post$sigma_array, na.rm = TRUE)
  ord       <- order(avg_sigma, decreasing = TRUE)
  hi_idx    <- head(ord, 3L)
  lo_idx    <- tail(ord, 3L)
  sel_idx   <- unique(c(hi_idx, lo_idx))

  sv <- fit$sv_params
  nd_sv <- dim(sv)[1]; kpc_sv <- nd_sv %/% nchains

  tidy_trace <- function(vec, kpc, nch, label) {
    use <- kpc * nch
    if (use < 1L) return(NULL)
    v <- vec[seq_len(use)]
    tibble::tibble(
      Iteration = rep(seq_len(kpc), times = nch),
      Chain     = factor(rep(seq_len(nch), each = kpc)),
      Value     = v, Parameter = label)
  }

  trace_rows <- list()
  grp_of <- function(i) if (i %in% hi_idx) "high-sigma" else "low-sigma"
  for (i in sel_idx) {
    trace_rows[[length(trace_rows) + 1L]] <-
      tidy_trace(sv[, i, 1], kpc_sv, nchains, sprintf("mu[%s] (%s)", sectors[i], grp_of(i)))
    trace_rows[[length(trace_rows) + 1L]] <-
      tidy_trace(sv[, i, 2], kpc_sv, nchains, sprintf("phi[%s] (%s)", sectors[i], grp_of(i)))
  }

  if (!is.null(draws) && !is.null(draws$beta_draws)) {
    mfg <- which(sectors == "MFG"); if (length(mfg) == 0L) mfg <- min(3L, k)
    bd <- draws$beta_draws
    D  <- dim(bd)[1]; n_eff <- dim(bd)[2]
    kpc_b <- D %/% nchains
    beta_vec <- bd[, n_eff, mfg, 1]
    trace_rows[[length(trace_rows) + 1L]] <-
      tidy_trace(beta_vec, kpc_b, nchains,
                 sprintf("beta_t[%s<-lag1, t=T]", sectors[mfg]))
  }

  trace_df <- dplyr::bind_rows(trace_rows)
  if (!is.null(trace_df) && nrow(trace_df) > 0) {
    trace_df$Parameter <- factor(trace_df$Parameter, levels = unique(trace_df$Parameter))
    p_tr <- ggplot2::ggplot(trace_df, ggplot2::aes(Iteration, Value, color = Chain)) +
      ggplot2::geom_line(linewidth = 0.3, alpha = 0.8, na.rm = TRUE) +
      ggplot2::facet_wrap(~ Parameter, scales = "free_y", ncol = 2) +
      ggplot2::scale_color_viridis_d(option = "D", end = 0.9) +
      ggplot2::labs(
        title = "MCMC trace plots for key parameters",
        subtitle = "SV mu/phi for the 3 highest- and 3 lowest-volatility sectors + a representative time-varying beta element",
        x = "Iteration (within chain)", y = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")
    save_eval_plot(p_tr, "figure_eval_03_trace_plots", 12, max(7, 1.6 * length(unique(trace_df$Parameter)) / 2))
  }
} else cat("No SV draws in fit; skipped trace plots (figure_eval_03).\n")

{
  xy   <- build_xy(prep, post$p); Xmat <- xy$Xmat; Ydep <- xy$Ydep; n <- xy$n
  kk   <- post$k
  ols  <- ols_var(Xmat, Ydep)
  per_rows <- list()
  for (i in 1:kk) {
    tvp_mu   <- rowSums(Xmat * post$beta_array[, i, ])
    tvp_sd   <- pmax(post$sigma_array[, i], 1e-8)
    const_mu <- as.numeric(Xmat %*% ols$B[i, ])
    const_sd <- pmax(ols$sigma[i], 1e-8)
    yi <- Ydep[i, ]
    rmse_tvp   <- sqrt(mean((yi - tvp_mu)^2))
    rmse_const <- sqrt(mean((yi - const_mu)^2))
    lpd_tvp    <- mean(dnorm(yi, tvp_mu,   tvp_sd,   log = TRUE))
    lpd_const  <- mean(dnorm(yi, const_mu, const_sd, log = TRUE))
    per_rows[[i]] <- tibble::tibble(
      Sector = post$sectors[i],
      RMSE_TVP_SV = round(rmse_tvp, 4),  RMSE_Const = round(rmse_const, 4),
      RMSE_reduction_pct = round((rmse_const - rmse_tvp) / rmse_const * 100, 2),
      LPD_TVP_SV = round(lpd_tvp, 4),    LPD_Const = round(lpd_const, 4),
      LPD_gain = round(lpd_tvp - lpd_const, 4))
  }
  per_sector <- dplyr::bind_rows(per_rows) %>% dplyr::arrange(dplyr::desc(RMSE_reduction_pct))
  export_table(per_sector, "table_33b_model_comparison_by_sector")
  results$model_comparison_by_sector <- per_sector

  ps <- per_sector
  ps$Sector <- factor(ps$Sector, levels = ps$Sector[order(ps$RMSE_reduction_pct)])
  p_ps <- ggplot2::ggplot(ps, ggplot2::aes(Sector, RMSE_reduction_pct,
                                           fill = RMSE_reduction_pct >= 0)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_hline(yintercept = 0, color = LINE_COL, linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = FILL_COL, `FALSE` = ALT_COL), guide = "none") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Per-sector one-step RMSE reduction: TVP-BVAR-SV vs constant BVAR",
      subtitle = "Positive = TVP-BVAR-SV lowers in-sample one-step RMSE relative to the constant BVAR",
      x = NULL, y = "RMSE reduction (%)") +
    ggplot2::theme_minimal(base_size = 12)
  save_eval_plot(p_ps, "figure_eval_04_predictive_by_sector", 10, 7)
}

dir.create(dirname(CONFIG$PATHS$robustness), recursive = TRUE, showWarnings = FALSE)
saveRDS(results, CONFIG$PATHS$robustness)
cat(sprintf("\nRobustness & evaluation complete -> %s\n", CONFIG$PATHS$robustness))
end_logging()
