#!/usr/bin/env Rscript
source("scripts/00_setup/setup.R")
start_logging("11b_credible_bands")

suppressWarnings(suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) }))

stopifnot(
  exists("compute_gfevd_array"), exists("compute_network_tv"),
  exists("compute_resilience_tv"), exists("export_table"))

if (!file.exists(CONFIG$PATHS$posterior_draws))
  stop("posterior_draws.rds not found. Run 02_estimation.R first.")

draws <- readRDS(CONFIG$PATHS$posterior_draws)
post  <- readRDS(CONFIG$PATHS$posterior)
conn  <- readRDS(CONFIG$PATHS$gfevd)

stopifnot(!is.null(draws$beta_draws), !is.null(draws$sigma_draws))

k       <- conn$k; p <- conn$p
sectors <- conn$sectors
dates   <- conn$dates
n_eff   <- length(dates)
dates_v <- as.Date(dates)
labels  <- CONFIG$SECTOR_LABELS
lab_vec <- ifelse(sectors %in% names(labels), unname(labels[sectors]), sectors)

FIG_DIR <- CONFIG$PATHS$figures
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

probs <- CONFIG$CREDIBLE_BANDS$QUANTILE_PROBS
stopifnot(length(probs) == 3L)

D     <- dim(draws$beta_draws)[1]
ndraw <- min(as.integer(CONFIG$CREDIBLE_BANDS$N_DRAWS), D)
sel   <- unique(round(seq(1, D, length.out = ndraw)))
ndraw <- length(sel)
cat(sprintf("Credible bands: propagating %d of %d posterior draws through GFEVD->network->resilience.\n",
            ndraw, D))

cent_metrics <- c("OutDegree", "InDegree", "Betweenness", "Closeness")
res_metrics  <- c("Resistance", "RecoverySpeed", "Reconfiguration")

mk_acc <- function() array(NA_real_, dim = c(k, n_eff, ndraw))
cent_acc <- setNames(lapply(cent_metrics, function(...) mk_acc()), cent_metrics)
res_acc  <- setNames(lapply(res_metrics,  function(...) mk_acc()), res_metrics)

as_mat <- function(df, metric) {
  d <- df %>% dplyr::mutate(.s = factor(Sector, levels = sectors)) %>%
    dplyr::arrange(TimeIndex, .s)
  matrix(d[[metric]], nrow = k, ncol = n_eff)
}

for (di in seq_along(sel)) {
  d <- sel[di]
  beta_d  <- draws$beta_draws[d, , , ]
  sigma_d <- draws$sigma_draws[d, , ]

  gfevd_d <- compute_gfevd_array(beta_d, sigma_d, k, p, CONFIG$HORIZON_GFEVD, CONFIG$EPS)
  net_d   <- compute_network_tv(gfevd_d, sectors, dates,
                                threshold_type = CONFIG$NETWORK_THRESHOLD_TYPE,
                                threshold_value = CONFIG$NETWORK_THRESHOLD)
  res_d   <- compute_resilience_tv(
    beta_array = beta_d, sigma_array = sigma_d, gfevd_array = gfevd_d,
    sectors = sectors, dates = dates, k = k, p = p,
    horizon = CONFIG$HORIZON_GIRF, resist_h = CONFIG$RESIST_HORIZON,
    eps = CONFIG$RESIST_EPS)

  for (m in cent_metrics) cent_acc[[m]][, , di] <- as_mat(net_d$centrality_tv, m)
  for (m in res_metrics)  res_acc[[m]][, , di]  <- as_mat(res_d, m)

  if (di %% 25L == 0L || di == ndraw)
    cat(sprintf("  ... %d/%d draws done\n", di, ndraw))
}

band_long <- function(acc_list, metric_names) {
  out <- list()
  for (m in metric_names) {
    arr <- acc_list[[m]]
    q <- apply(arr, c(1, 2), function(v) stats::quantile(v, probs = probs, na.rm = TRUE))
    lower <- q[1, , ]; med <- q[2, , ]; upper <- q[3, , ]
    out[[m]] <- tibble::tibble(
      Sector    = rep(sectors, times = n_eff),
      SectorLabel = rep(lab_vec, times = n_eff),
      TimeIndex = rep(seq_len(n_eff), each = k),
      Date      = rep(dates_v, each = k),
      Metric    = m,
      Lower     = as.numeric(lower),
      Median    = as.numeric(med),
      Upper     = as.numeric(upper))
  }
  dplyr::bind_rows(out)
}

centrality_bands <- band_long(cent_acc, cent_metrics)
resilience_bands <- band_long(res_acc,  res_metrics)

export_table(centrality_bands, "table_43_centrality_credible_bands")
export_table(resilience_bands, "table_44_resilience_credible_bands")

ribbon_fig <- function(bands, metric, fig_name, ylab) {
  df <- dplyr::filter(bands, Metric == metric)
  df$SectorLabel <- factor(df$SectorLabel, levels = lab_vec)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Date)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper),
                         alpha = 0.25, fill = "#2c7fb8") +
    ggplot2::geom_line(ggplot2::aes(y = Median), color = "#08306b", linewidth = 0.5) +
    ggplot2::facet_wrap(~ SectorLabel, scales = "free_y", ncol = 3) +
    ggplot2::labs(
      title = sprintf("%s: posterior median with 95%% credible band", metric),
      subtitle = sprintf("Time-varying band over %d posterior draws (GFEVD->network->resilience per draw)", ndraw),
      x = NULL, y = ylab) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(FIG_DIR, paste0(fig_name, ".png")), p,
                  width = 12, height = 9, dpi = 300, bg = "white")
  tryCatch(ggplot2::ggsave(file.path(FIG_DIR, paste0(fig_name, ".pdf")), p,
                           width = 12, height = 9, device = "pdf"), error = function(e) NULL)
  cat(sprintf("  saved %s.png\n", fig_name))
}

ribbon_fig(centrality_bands, "OutDegree",       "figure_31_outdegree_credible",      "Out-degree (to-others)")
ribbon_fig(centrality_bands, "InDegree",        "figure_32_indegree_credible",       "In-degree (from-others)")
ribbon_fig(centrality_bands, "Betweenness",     "figure_33_betweenness_credible",    "Betweenness")
ribbon_fig(centrality_bands, "Closeness",       "figure_34_closeness_credible",      "Closeness")
ribbon_fig(resilience_bands, "Resistance",      "figure_35_resistance_credible",     "Resistance")
ribbon_fig(resilience_bands, "RecoverySpeed",   "figure_36_recovery_credible",       "Recovery speed")
ribbon_fig(resilience_bands, "Reconfiguration", "figure_37_reconfiguration_credible","Reconfiguration")

sig <- post$sigma_array
stopifnot(!is.null(sig))
vol_long <- tibble::tibble(
  Sector      = rep(sectors, each = n_eff),
  SectorLabel = rep(lab_vec, each = n_eff),
  TimeIndex   = rep(seq_len(n_eff), times = k),
  Date        = rep(dates_v, times = k),
  Sigma       = as.numeric(sig))
vol_long$SectorLabel <- factor(vol_long$SectorLabel, levels = lab_vec)

p_vol <- ggplot2::ggplot(vol_long, ggplot2::aes(Date, Sigma)) +
  ggplot2::geom_line(color = "#b2182b", linewidth = 0.5, na.rm = TRUE) +
  ggplot2::facet_wrap(~ SectorLabel, scales = "free_y", ncol = 3) +
  ggplot2::labs(
    title = expression("Stochastic volatility by sector: " * sigma[it] == exp(h[it]/2)),
    subtitle = "Posterior-mean time-varying standard deviation of the reduced-form shocks",
    x = NULL, y = expression(sigma[it])) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
ggplot2::ggsave(file.path(FIG_DIR, "figure_38_volatility_by_sector.png"), p_vol,
                width = 12, height = 9, dpi = 300, bg = "white")
tryCatch(ggplot2::ggsave(file.path(FIG_DIR, "figure_38_volatility_by_sector.pdf"), p_vol,
                         width = 12, height = 9, device = "pdf"), error = function(e) NULL)
cat("  saved figure_38_volatility_by_sector.png\n")

vol_summary <- vol_long %>% dplyr::group_by(Sector, SectorLabel) %>%
  dplyr::summarise(
    Mean_sigma = mean(Sigma, na.rm = TRUE),
    Min_sigma  = min(Sigma, na.rm = TRUE),
    Max_sigma  = max(Sigma, na.rm = TRUE),
    SD_sigma   = stats::sd(Sigma, na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(Mean_sigma))
export_table(vol_summary, "table_43b_volatility_summary")

dir.create(dirname(CONFIG$PATHS$credible_bands), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(centrality_bands = centrality_bands,
             resilience_bands = resilience_bands,
             volatility_long = vol_long, volatility_summary = vol_summary,
             n_draws_used = ndraw, quantile_probs = probs,
             sectors = sectors, dates = dates_v),
        CONFIG$PATHS$credible_bands)
cat(sprintf("Credible bands done. Centrality rows: %d | Resilience rows: %d | draws used: %d.\n",
            nrow(centrality_bands), nrow(resilience_bands), ndraw))

end_logging()
