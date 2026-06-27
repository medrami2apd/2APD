.theme_thesis <- function() ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold"))

.save_plot <- function(p, name, w = 9, h = 5) {
  ggplot2::ggsave(file.path("outputs/figures", paste0(name, ".png")),
                  p, width = w, height = h, dpi = 200)
  invisible(p)
}

break_dates_from <- function(bp) {
  if (is.null(bp) || is.null(bp$regimes) || nrow(bp$regimes) < 2) return(NULL)
  as.Date(bp$regimes$EndDate[-nrow(bp$regimes)])
}
regime_segments <- function(bp) {
  if (is.null(bp) || is.null(bp$regimes)) return(NULL)
  data.frame(x = as.Date(bp$regimes$StartDate), xend = as.Date(bp$regimes$EndDate),
             y = bp$regimes$Mean)
}

plot_tci_trend <- function(total_conn, dates) {
  d <- total_conn; d$Date <- as.Date(dates)
  fit <- stats::lm(TotalConnectedness ~ as.numeric(Date), data = d)
  d$Fitted <- stats::predict(fit)
  p <- ggplot2::ggplot(d, ggplot2::aes(Date, TotalConnectedness)) +
    ggplot2::geom_line(color = "#2c3e50", linewidth = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = Fitted), color = "#c0392b",
                       linetype = "dashed", linewidth = 0.8) +
    ggplot2::labs(title = "H1: Total Connectedness Index with linear trend",
                  x = NULL, y = "TCI (%)") + .theme_thesis()
  .save_plot(p, "fig_h1_tci_trend")
}

plot_tci <- function(total_conn, dates, breaks = NULL) {
  d <- total_conn; d$Date <- as.Date(dates)
  p <- ggplot2::ggplot(d, ggplot2::aes(Date, TotalConnectedness)) +
    ggplot2::geom_line(color = "#2c3e50", linewidth = 0.5)
  if (!is.null(breaks) && length(breaks) > 0L)
    p <- p + ggplot2::geom_vline(xintercept = as.numeric(as.Date(breaks)),
                                 linetype = "dashed", color = "#c0392b")
  p <- p + ggplot2::labs(title = "Total Connectedness Index (TCI)",
                         x = NULL, y = "TCI (%)") + .theme_thesis()
  .save_plot(p, "fig_tci")
}

plot_series_breaks <- function(series_df, value_col, dates, bp, title, name,
                               ylab = NULL) {
  d <- series_df; d$Date <- as.Date(dates); d$.val <- d[[value_col]]
  p <- ggplot2::ggplot(d, ggplot2::aes(Date, .val)) +
    ggplot2::geom_line(color = "#2c3e50", linewidth = 0.5)
  segs <- regime_segments(bp)
  if (!is.null(segs)) p <- p + ggplot2::geom_segment(
    data = segs, ggplot2::aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE, color = "#2980b9", linewidth = 1)
  bd <- break_dates_from(bp)
  if (!is.null(bd)) p <- p + ggplot2::geom_vline(xintercept = as.numeric(bd),
                                                 linetype = "dashed", color = "#c0392b")
  p <- p + ggplot2::labs(title = title, x = NULL, y = ylab %||% value_col) + .theme_thesis()
  .save_plot(p, name)
}

plot_net_directional <- function(directional_df) {
  agg <- directional_df %>% dplyr::group_by(Sector) %>%
    dplyr::summarise(Net = mean(NetSpillover, na.rm = TRUE), .groups = "drop")
  p <- ggplot2::ggplot(agg, ggplot2::aes(stats::reorder(Sector, Net), Net,
                                          fill = Net > 0)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#27ae60", "FALSE" = "#c0392b"),
                               labels = c("TRUE" = "Transmitter", "FALSE" = "Absorber"),
                               name = NULL) +
    ggplot2::labs(title = "H2: Mean net directional connectedness",
                  x = NULL, y = "Net (To - From)") + .theme_thesis()
  .save_plot(p, "fig_h2_net_directional")
}

plot_spillover_heatmap <- function(spillover, sectors) {
  df <- as.data.frame(as.table(as.matrix(spillover)))
  names(df) <- c("Receiver", "Source", "Share")
  p <- ggplot2::ggplot(df, ggplot2::aes(Source, Receiver, fill = Share)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08306b", name = "%") +
    ggplot2::labs(title = "Average pairwise spillovers (row = receiver, col = source)",
                  x = "Source (transmitter)", y = "Receiver (absorber)") +
    .theme_thesis() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  .save_plot(p, "fig_h2_spillover_heatmap", w = 8, h = 7)
}

plot_net_heatmap <- function(directional_df, dates) {
  d <- directional_df
  d$Date <- as.Date(dates)[d$TimeIndex]
  ord <- d %>% dplyr::group_by(Sector) %>%
    dplyr::summarise(m = mean(NetSpillover, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(m) %>% dplyr::pull(Sector)
  d$Sector <- factor(d$Sector, levels = ord)
  p <- ggplot2::ggplot(d, ggplot2::aes(Date, Sector, fill = NetSpillover)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#c0392b", mid = "white", high = "#27ae60",
                                  midpoint = 0, name = "Net") +
    ggplot2::labs(title = "H2: Net spillover by sector over time", x = NULL, y = NULL) +
    .theme_thesis()
  .save_plot(p, "fig_h2_net_heatmap", w = 10, h = 5)
}

plot_role_timeline <- function(directional_df, dates, band = 0) {
  d <- directional_df
  d$Date <- as.Date(dates)[d$TimeIndex]
  d$Role <- dplyr::case_when(d$NetSpillover >  band ~ "Transmitter",
                             d$NetSpillover < -band ~ "Absorber",
                             TRUE ~ "Neutral")
  ord <- d %>% dplyr::group_by(Sector) %>%
    dplyr::summarise(m = mean(NetSpillover, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(m) %>% dplyr::pull(Sector)
  d$Sector <- factor(d$Sector, levels = ord)
  p <- ggplot2::ggplot(d, ggplot2::aes(Date, Sector, fill = Role)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = c(Transmitter = "#27ae60",
                                          Absorber = "#c0392b", Neutral = "#ecf0f1")) +
    ggplot2::labs(title = "H2: Sector role timeline", x = NULL, y = NULL) +
    .theme_thesis()
  .save_plot(p, "fig_h2_role_timeline", w = 10, h = 5)
}

plot_centrality_series <- function(centrality_tv) {
  long <- centrality_tv %>%
    dplyr::select(Date, Sector, OutDegree, InDegree, Betweenness, Closeness) %>%
    tidyr::pivot_longer(-c(Date, Sector), names_to = "Measure", values_to = "Value")
  p <- ggplot2::ggplot(long, ggplot2::aes(Date, Value, group = Sector)) +
    ggplot2::geom_line(alpha = 0.35, color = "#34495e", linewidth = 0.3) +
    ggplot2::facet_wrap(~Measure, scales = "free_y") +
    ggplot2::labs(title = "H3: Time-varying centrality (all sectors)", x = NULL, y = NULL) +
    .theme_thesis()
  .save_plot(p, "fig_h3_centrality_series", w = 10, h = 6)
}

plot_centrality_resilience_scatter <- function(centrality_tv, resilience_tv) {
  panel <- dplyr::inner_join(centrality_tv, resilience_tv,
                             by = c("Sector", "TimeIndex"))
  cl <- panel %>% tidyr::pivot_longer(c(OutDegree, InDegree, Betweenness, Closeness),
                                      names_to = "Centrality", values_to = "CentVal")
  both <- cl %>% tidyr::pivot_longer(c(Resistance, RecoverySpeed, Reconfiguration),
                                     names_to = "Resilience", values_to = "ResVal")
  p <- ggplot2::ggplot(both, ggplot2::aes(CentVal, ResVal)) +
    ggplot2::geom_point(alpha = 0.15, size = 0.5, color = "#34495e") +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "#c0392b",
                         linewidth = 0.7, formula = y ~ x) +
    ggplot2::facet_grid(Resilience ~ Centrality, scales = "free") +
    ggplot2::labs(title = "H3: Centrality vs resilience (pooled, with OLS fit)",
                  x = "Centrality", y = "Resilience") + .theme_thesis()
  .save_plot(p, "fig_h3_centrality_resilience_scatter", w = 11, h = 8)
}

plot_density <- function(density_df) {
  p <- ggplot2::ggplot(density_df, ggplot2::aes(Date, Density)) +
    ggplot2::geom_line(color = "#8e44ad", linewidth = 0.5) +
    ggplot2::labs(title = "Network density over time", x = NULL, y = "Density") +
    .theme_thesis()
  .save_plot(p, "fig_density")
}

plot_network_graph <- function(gfevd_array, sectors, dates, target_date,
                               name, q = 0.85) {
  dates <- as.Date(dates)
  idx <- which.min(abs(dates - as.Date(target_date)))
  adj <- build_transmission_adjacency(gfevd_array[, , idx], sectors)
  pos <- adj[adj > 0]
  thr <- if (length(pos)) stats::quantile(pos, q, na.rm = TRUE) else 0
  a <- adj; a[a < thr] <- 0
  g <- igraph::graph_from_adjacency_matrix(a, mode = "directed", weighted = TRUE)
  out_strength <- rowSums(adj)
  rng <- range(out_strength, finite = TRUE)
  rsc <- if (diff(rng) > 0) (out_strength - rng[1]) / diff(rng) else rep(0.5, length(out_strength))
  vsize <- 10 + 20 * rsc
  ew <- if (igraph::ecount(g) > 0) scales::rescale(igraph::E(g)$weight, to = c(0.4, 4)) else 1
  qlab <- sprintf("%s-Q%d", format(dates[idx], "%Y"),
                  as.integer(format(dates[idx], "%m")) %/% 3L + 1L)
  grDevices::png(file.path("outputs/figures", paste0(name, ".png")),
                 width = 1500, height = 1300, res = 200)
  on.exit(grDevices::dev.off())
  set.seed(1)
  plot(g,
       layout = igraph::layout_in_circle,
       vertex.size = vsize,
       vertex.color = "#aed6f1",
       vertex.label = sectors,
       vertex.label.cex = 0.8,
       vertex.label.color = "#1b2631",
       edge.width = ew,
       edge.color = grDevices::adjustcolor("#566573", 0.6),
       edge.arrow.size = 0.25,
       main = sprintf("Spillover network - %s", qlab))
  invisible(name)
}
