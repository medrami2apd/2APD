#!/usr/bin/env Rscript
source("scripts/00_setup/setup.R")
start_logging("04_network")

suppressWarnings(suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) }))

stopifnot(
  exists("compute_network_tv"), exists("export_table"))

FIG_DIR <- CONFIG$PATHS$figures
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
labels <- CONFIG$SECTOR_LABELS
lab_of <- function(s) ifelse(s %in% names(labels), unname(labels[s]), s)

save_plot <- function(p, name, w = 11, h = 6, dpi = 300) {
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = dpi, bg = "white")
  tryCatch(ggplot2::ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf"),
           error = function(e) NULL)
  cat(sprintf("  saved %s.png\n", name))
}

conn <- readRDS(CONFIG$PATHS$gfevd)

net <- compute_network_tv(conn$gfevd, conn$sectors, conn$dates,
                          threshold_type = CONFIG$NETWORK_THRESHOLD_TYPE,
                          threshold_value = CONFIG$NETWORK_THRESHOLD)

thr_grid    <- CONFIG$NETWORK_THRESHOLD_GRID
thr_default <- CONFIG$NETWORK_THRESHOLD
path_metrics <- c("Betweenness", "Closeness")

cent_by_thr <- list()
for (thr in thr_grid) {
  cent_by_thr[[as.character(thr)]] <-
    if (isTRUE(all.equal(thr, thr_default))) net$centrality_tv
    else compute_network_tv(conn$gfevd, conn$sectors, conn$dates,
                            threshold_type = CONFIG$NETWORK_THRESHOLD_TYPE,
                            threshold_value = thr)$centrality_tv
}
cent_default <- if (!is.null(cent_by_thr[[as.character(thr_default)]]))
  cent_by_thr[[as.character(thr_default)]] else net$centrality_tv

thr_rows <- list()
for (thr in thr_grid) {
  avg <- cent_by_thr[[as.character(thr)]] %>% dplyr::group_by(Sector) %>%
    dplyr::summarise(Betweenness = mean(Betweenness, na.rm = TRUE),
                     Closeness   = mean(Closeness,   na.rm = TRUE), .groups = "drop")
  avg$Threshold <- thr
  thr_rows[[length(thr_rows) + 1L]] <- avg
}
export_table(dplyr::bind_rows(thr_rows), "table_46_centrality_threshold_sensitivity")

n_eff   <- length(conn$dates)
dates_v <- as.Date(conn$dates)

spearman_at_t <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4L) return(NA_real_)
  if (stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}

stab_rows <- list()
for (metric in path_metrics) {
  d_def <- cent_default %>% dplyr::arrange(TimeIndex, Sector)
  def_mat <- matrix(d_def[[metric]], nrow = length(conn$sectors))
  for (thr in thr_grid) {
    d_alt <- cent_by_thr[[as.character(thr)]] %>% dplyr::arrange(TimeIndex, Sector)
    alt_mat <- matrix(d_alt[[metric]], nrow = length(conn$sectors))
    rho_t <- vapply(seq_len(n_eff), function(t) spearman_at_t(def_mat[, t], alt_mat[, t]), numeric(1))
    stab_rows[[length(stab_rows) + 1L]] <- tibble::tibble(
      Metric = metric, Threshold = thr, TimeIndex = seq_len(n_eff),
      Date = dates_v, Spearman = rho_t)
  }
}
stability_tv <- dplyr::bind_rows(stab_rows)

threshold_sensitivity_revised <- stability_tv %>%
  dplyr::group_by(Metric, Threshold) %>%
  dplyr::summarise(
    Mean_Spearman   = mean(Spearman, na.rm = TRUE),
    Median_Spearman = stats::median(Spearman, na.rm = TRUE),
    Min_Spearman    = suppressWarnings(min(Spearman, na.rm = TRUE)),
    SD_Spearman     = stats::sd(Spearman, na.rm = TRUE),
    N_periods       = sum(is.finite(Spearman)),
    .groups = "drop") %>%
  dplyr::mutate(Is_default = abs(Threshold - thr_default) < 1e-12) %>%
  dplyr::arrange(Metric, Threshold)
export_table(threshold_sensitivity_revised, "table_46_centrality_threshold_sensitivity_revised")

plot_df <- stability_tv %>%
  dplyr::filter(abs(Threshold - thr_default) > 1e-12) %>%
  dplyr::mutate(Threshold = factor(sprintf("%.2f", Threshold)))
if (nrow(plot_df) > 0) {
  p40 <- ggplot2::ggplot(plot_df, ggplot2::aes(Date, Spearman, color = Threshold)) +
    ggplot2::geom_line(linewidth = 0.6, na.rm = TRUE) +
    ggplot2::facet_wrap(~ Metric, ncol = 1) +
    ggplot2::scale_color_viridis_d(option = "D", end = 0.9, name = "Alt. threshold") +
    ggplot2::coord_cartesian(ylim = c(-0.2, 1)) +
    ggplot2::labs(
      title = "Network Threshold Sensitivity: time-varying ranking stability (H3)",
      subtitle = sprintf("Spearman rank correlation of sector %s vs default threshold (%.2f), per quarter",
                         paste(path_metrics, collapse = " / "), thr_default),
      x = NULL, y = "Spearman rank correlation") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
  save_plot(p40, "figure_40_threshold_stability", 12, 8)
}

centrality_avg <- net$centrality_tv %>% dplyr::group_by(Sector) %>%
  dplyr::summarise(dplyr::across(c(OutDegree, InDegree, NetDegree, Betweenness,
                                   Closeness, Eigenvector, PageRank),
                                 ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(OutDegree))

export_table(net$centrality_tv, "table_08_centrality_time_varying")
export_table(centrality_avg, "table_09_centrality_average")
export_table(net$density_tv, "table_10_network_density")

dir.create(dirname(CONFIG$PATHS$network), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(centrality_tv = net$centrality_tv, centrality_avg = centrality_avg,
             density_tv = net$density_tv, sectors = conn$sectors,
             threshold_stability_tv = stability_tv,
             threshold_sensitivity_revised = threshold_sensitivity_revised),
        CONFIG$PATHS$network)
cat("Network analysis done. Centrality measures: OutDegree, InDegree, NetDegree, Betweenness, Closeness, Eigenvector, PageRank\n")
cat(sprintf("Time-resolved threshold sensitivity computed for %s across %d thresholds.\n",
            paste(path_metrics, collapse = "/"), length(thr_grid)))

end_logging()
