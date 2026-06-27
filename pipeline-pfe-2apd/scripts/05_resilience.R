#!/usr/bin/env Rscript
source("scripts/00_setup/setup.R")
start_logging("05_resilience")

suppressWarnings(suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) }))

stopifnot(
  exists("compute_resilience_tv"), exists("export_table"))

FIG_DIR <- CONFIG$PATHS$figures
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

post <- readRDS(CONFIG$PATHS$posterior)
conn <- readRDS(CONFIG$PATHS$gfevd)

if (!is.null(conn$n_unstable) && conn$n_unstable > 0)
  warning(sprintf("Resilience computed over %d period(s) with spectral radius >= 1 (see gfevd_results.rds$unstable).",
                  conn$n_unstable))

resilience_tv <- compute_resilience_tv(
  beta_array = post$beta_array, sigma_array = post$sigma_array,
  gfevd_array = conn$gfevd, sectors = post$sectors, dates = post$dates,
  k = post$k, p = post$p, horizon = CONFIG$HORIZON_GIRF,
  resist_h = CONFIG$RESIST_HORIZON, eps = CONFIG$RESIST_EPS)

resilience_avg <- resilience_tv %>% dplyr::group_by(Sector) %>%
  dplyr::summarise(Resistance = mean(Resistance, na.rm = TRUE),
                   RecoverySpeed = mean(RecoverySpeed, na.rm = TRUE),
                   Reconfiguration = mean(Reconfiguration, na.rm = TRUE),
                   .groups = "drop")

export_table(resilience_tv, "table_12_resilience_time_varying")
export_table(resilience_avg, "table_13_resilience_average")

to_wide <- function(panel, value_col) {
  panel %>% dplyr::arrange(TimeIndex) %>%
    tidyr::pivot_wider(id_cols = c(TimeIndex, Date),
                       names_from = Sector, values_from = dplyr::all_of(value_col))
}
recovery_wide       <- to_wide(resilience_tv, "RecoverySpeed")
resistance_wide     <- to_wide(resilience_tv, "Resistance")
reconfiguration_wide<- to_wide(resilience_tv, "Reconfiguration")
export_table(recovery_wide,        "table_14_recovery_speed_by_sector_tv")
export_table(resistance_wide,      "table_15_resistance_by_sector_tv")
export_table(reconfiguration_wide, "table_16_reconfiguration_by_sector_tv")
cat(sprintf("Recovery-speed index: %d periods x %d sectors.\n",
            nrow(recovery_wide), ncol(recovery_wide) - 2L))

eps_grid    <- CONFIG$RESIST_EPS_GRID
eps_default <- CONFIG$RESIST_EPS
res_metrics <- c("Resistance", "RecoverySpeed", "Reconfiguration")
sectors     <- post$sectors
n_sec       <- length(sectors)
dates_v     <- as.Date(post$dates)
n_eff       <- length(dates_v)

res_by_eps <- list()
for (e in eps_grid) {
  res_by_eps[[as.character(e)]] <-
    if (isTRUE(all.equal(e, eps_default))) resilience_tv
    else compute_resilience_tv(
      beta_array = post$beta_array, sigma_array = post$sigma_array,
      gfevd_array = conn$gfevd, sectors = sectors, dates = post$dates,
      k = post$k, p = post$p, horizon = CONFIG$HORIZON_GIRF,
      resist_h = CONFIG$RESIST_HORIZON, eps = e)
}
res_default <- res_by_eps[[as.character(eps_default)]]
if (is.null(res_default)) res_default <- resilience_tv

spearman_at_t <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4L) return(NA_real_)
  if (stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}
as_mat <- function(panel, metric) {
  d <- panel %>% dplyr::arrange(TimeIndex, Sector)
  matrix(d[[metric]], nrow = n_sec)
}

eps_rows <- list()
for (metric in res_metrics) {
  def_mat <- as_mat(res_default, metric)
  for (e in eps_grid) {
    alt_mat <- as_mat(res_by_eps[[as.character(e)]], metric)
    rho_t <- vapply(seq_len(n_eff), function(t) spearman_at_t(def_mat[, t], alt_mat[, t]), numeric(1))
    eps_rows[[length(eps_rows) + 1L]] <- tibble::tibble(
      Metric = metric, Eps = e, TimeIndex = seq_len(n_eff),
      Date = dates_v, Spearman = rho_t)
  }
}
eps_stability_tv <- dplyr::bind_rows(eps_rows)

resist_eps_sensitivity <- eps_stability_tv %>%
  dplyr::group_by(Metric, Eps) %>%
  dplyr::summarise(
    Median_Spearman = stats::median(Spearman, na.rm = TRUE),
    Mean_Spearman   = mean(Spearman, na.rm = TRUE),
    Min_Spearman    = suppressWarnings(min(Spearman, na.rm = TRUE)),
    SD_Spearman     = stats::sd(Spearman, na.rm = TRUE),
    N_periods       = sum(is.finite(Spearman)),
    .groups = "drop") %>%
  dplyr::mutate(Is_default = abs(Eps - eps_default) < 1e-12) %>%
  dplyr::arrange(Metric, Eps)
export_table(resist_eps_sensitivity, "table_45_resist_eps_sensitivity")

plot_df <- eps_stability_tv %>%
  dplyr::filter(abs(Eps - eps_default) > 1e-12) %>%
  dplyr::mutate(Eps = factor(format(Eps, scientific = TRUE)))
if (nrow(plot_df) > 0) {
  p39 <- ggplot2::ggplot(plot_df, ggplot2::aes(Date, Spearman, color = Eps)) +
    ggplot2::geom_line(linewidth = 0.6, na.rm = TRUE) +
    ggplot2::facet_wrap(~ Metric, ncol = 1) +
    ggplot2::scale_color_viridis_d(option = "C", end = 0.85, name = expression(epsilon)) +
    ggplot2::coord_cartesian(ylim = c(-0.2, 1)) +
    ggplot2::labs(
      title = "RESIST_EPS sensitivity: time-varying ranking stability (H3)",
      subtitle = sprintf("Spearman rank correlation of sector resilience vs default epsilon (%s), per quarter",
                         format(eps_default, scientific = TRUE)),
      x = NULL, y = "Spearman rank correlation") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
  ggplot2::ggsave(file.path(FIG_DIR, "figure_39_resist_eps_stability.png"),
                  p39, width = 11, height = 8, dpi = 300, bg = "white")
  tryCatch(ggplot2::ggsave(file.path(FIG_DIR, "figure_39_resist_eps_stability.pdf"),
                           p39, width = 11, height = 8, device = "pdf"), error = function(e) NULL)
  cat("  saved figure_39_resist_eps_stability.png\n")
}
cat(sprintf("RESIST_EPS sensitivity done over grid {%s}; median Spearman (default excl.) summarised.\n",
            paste(format(eps_grid, scientific = TRUE), collapse = ", ")))

dir.create(dirname(CONFIG$PATHS$resilience), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(resilience_tv = resilience_tv, resilience_avg = resilience_avg,
             recovery_speed_by_sector_tv = recovery_wide,
             resistance_by_sector_tv = resistance_wide,
             reconfiguration_by_sector_tv = reconfiguration_wide,
             resist_eps_sensitivity = resist_eps_sensitivity,
             resist_eps_stability_tv = eps_stability_tv,
             horizon = CONFIG$HORIZON_GIRF, sectors = post$sectors),
        CONFIG$PATHS$resilience)
cat(sprintf("Resilience done. %d panel obs.\n", nrow(resilience_tv)))

end_logging()
