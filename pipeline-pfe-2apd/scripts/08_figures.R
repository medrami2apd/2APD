#!/usr/bin/env Rscript
source("scripts/00_setup/setup.R")
start_logging("08_figures")

suppressWarnings(suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble); library(scales)
}))
has_zoo <- requireNamespace("zoo", quietly = TRUE)

stopifnot(
  exists("compute_spillover_matrix"), exists("build_transmission_adjacency"),
  exists("export_matrix"))

FIG_DIR <- CONFIG$PATHS$figures
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

BR_COL <- "#c0392b"; LINE_COL <- "#2c3e50"; FILL_COL <- "#3498db"; EP_COL <- "#16a085"

theme_pub <- function(base = 13) {
  ggplot2::theme_minimal(base_size = base) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = base + 2),
      plot.subtitle = ggplot2::element_text(color = "grey30", size = base - 2),
      axis.title    = ggplot2::element_text(face = "bold"),
      strip.text    = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank())
}

save_plot <- function(p, name, w = 11, h = 7, dpi = 300) {
  ggplot2::ggsave(file.path(FIG_DIR, paste0(name, ".png")), p,
                  width = w, height = h, dpi = dpi, bg = "white")
  tryCatch(ggplot2::ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p,
                           width = w, height = h, device = "pdf"),
           error = function(e) NULL)
  cat(sprintf("  saved %s.png\n", name))
}

fig_by_sector <- function(panel, value_col, title, subtitle, ylab, name,
                          color = FILL_COL, ncol = 4, w = 14, h = 9) {
  if (is.null(panel) || !(value_col %in% names(panel)) ||
      !all(c("Sector", "Date") %in% names(panel))) return(invisible(NULL))
  df <- panel
  df$Date <- as.Date(df$Date)
  df$SectorLab <- lab_of(df$Sector)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$Date, .data[[value_col]])) +
    ggplot2::geom_line(color = color, linewidth = 0.6, na.rm = TRUE) +
    ggplot2::facet_wrap(~ SectorLab, ncol = ncol, scales = "free_y") +
    theme_pub(11) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = ylab)
  save_plot(p, name, w, h)
}

conn   <- readRDS(CONFIG$PATHS$gfevd)
net    <- readRDS(CONFIG$PATHS$network)
res    <- tryCatch(readRDS(CONFIG$PATHS$resilience), error = function(e) NULL)
breaks <- tryCatch(readRDS(CONFIG$PATHS$breaks),     error = function(e) NULL)

sectors <- conn$sectors
labels  <- CONFIG$SECTOR_LABELS
n_eff   <- length(conn$dates)
date_lk <- tibble::tibble(TimeIndex = 1:n_eff, Date = as.Date(conn$dates))
lab_of  <- function(s) ifelse(s %in% names(labels), unname(labels[s]), s)
tci     <- dplyr::left_join(conn$total_conn, date_lk, by = "TimeIndex")

break_dates_of <- function(b) {
  if (is.null(b) || is.null(b$regimes) || nrow(b$regimes) < 2) return(as.Date(character(0)))
  as.Date(b$regimes$EndDate[-nrow(b$regimes)])
}
add_breaks <- function(p, bd) {
  if (length(bd) > 0)
    p <- p + ggplot2::geom_vline(xintercept = bd, linetype = "dashed",
                                 color = BR_COL, linewidth = 0.8)
  p
}
episode_df <- function() {
  ep <- CONFIG$HISTORICAL_EPISODES
  if (is.null(ep) || length(ep) == 0) return(NULL)
  d <- tibble::tibble(Date = as.Date(paste0(names(ep), "-01-01")), Label = unname(ep))
  rng <- range(as.Date(conn$dates))
  dplyr::filter(d, Date >= rng[1], Date <= rng[2])
}

cat("\n==== Generating H1-H4 figure suite ====\n")

cat("[H1] figure_01_tci ...\n")
bd_tci  <- break_dates_of(breaks$TCI)
tci_pl  <- tci
if (has_zoo) {
  w   <- min(12, max(4, floor(n_eff / 8)))
  sdv <- zoo::rollapply(tci_pl$TotalConnectedness, width = w, FUN = sd, fill = NA, align = "right")
  tci_pl$lo <- tci_pl$TotalConnectedness - 1.96 * sdv
  tci_pl$hi <- tci_pl$TotalConnectedness + 1.96 * sdv
  sub1 <- sprintf("Shaded band: +/-1.96 rolling SD (%d-period). Red dashed: structural breaks (H4).", w)
} else {
  tci_pl$lo <- NA; tci_pl$hi <- NA
  sub1 <- "Red dashed: structural breaks (H4)."
}
p1 <- ggplot(tci_pl, aes(Date, TotalConnectedness))
if (has_zoo) p1 <- p1 + geom_ribbon(aes(ymin = lo, ymax = hi), fill = FILL_COL, alpha = 0.20, na.rm = TRUE)
p1 <- p1 + geom_line(color = LINE_COL, linewidth = 1.1) +
  labs(title = "Total Intersectoral Connectedness Index (H1)", subtitle = sub1,
       x = NULL, y = "Total connectedness (%)") + theme_pub()
p1 <- add_breaks(p1, bd_tci)
save_plot(p1, "figure_01_tci", 12, 6)

cat("[H1] figure_02_tci_regimes ...\n")
if (!is.null(breaks$TCI$regimes) && nrow(breaks$TCI$regimes) >= 1) {
  reg <- breaks$TCI$regimes
  reg$StartDate <- as.Date(reg$StartDate); reg$EndDate <- as.Date(reg$EndDate)
  p2 <- ggplot(tci, aes(Date, TotalConnectedness)) +
    geom_line(color = LINE_COL, linewidth = 0.7, alpha = 0.6) +
    geom_segment(data = reg, inherit.aes = FALSE,
                 aes(x = StartDate, xend = EndDate, y = Mean, yend = Mean),
                 color = BR_COL, linewidth = 1.3) +
    labs(title = "TCI with Regime Means (H1 segmented by H4 breaks)",
         subtitle = "Red segments: mean TCI within each Bai-Perron regime",
         x = NULL, y = "Total connectedness (%)") + theme_pub()
  p2 <- add_breaks(p2, bd_tci)
  save_plot(p2, "figure_02_tci_regimes", 12, 6)
}

cat("[H1] figure_10_network_density ...\n")
dens <- net$density_tv %>% dplyr::mutate(Date = as.Date(Date))
p10 <- ggplot(dens, aes(Date, Density)) +
  geom_line(color = "#27ae60", linewidth = 1) +
  labs(title = "Directed Network Density Over Time (H1)",
       subtitle = "Density of the time-varying transmission network",
       x = NULL, y = "Density") + theme_pub()
p10 <- add_breaks(p10, bd_tci)
save_plot(p10, "figure_10_network_density", 12, 6)

cat("[H2] figure_03_net_spillover_ts ...\n")
dir_ts <- dplyr::left_join(conn$directional, date_lk, by = "TimeIndex")
dir_ts$SectorLab <- lab_of(dir_ts$Sector)
p3 <- ggplot(dir_ts, aes(Date, NetSpillover)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_line(color = "#8e44ad", linewidth = 0.6) +
  facet_wrap(~ SectorLab, ncol = 4, scales = "free_y") +
  labs(title = "Net Directional Spillovers by Sector (H2)",
       subtitle = "Positive = net transmitter; negative = net absorber",
       x = NULL, y = "Net spillover (%)") + theme_pub(11)
save_plot(p3, "figure_03_net_spillover_ts", 14, 9)

cat("[H2] figure_04_net_directional_bar ...\n")
dir_avg <- conn$directional %>% dplyr::group_by(Sector) %>%
  dplyr::summarise(To = mean(ToOthers, na.rm = TRUE), From = mean(FromOthers, na.rm = TRUE),
                   Net = mean(NetSpillover, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(SectorLab = lab_of(Sector),
                Role = ifelse(Net > 0, "Transmitter", "Absorber")) %>%
  dplyr::arrange(Net)
dir_avg$SectorLab <- factor(dir_avg$SectorLab, levels = dir_avg$SectorLab)
p4 <- ggplot(dir_avg, aes(SectorLab, Net, fill = Role)) +
  geom_col(color = "black", alpha = 0.85) + coord_flip() +
  scale_fill_manual(values = c(Absorber = "#2980b9", Transmitter = "#e67e22")) +
  labs(title = "Average Net Directional Connectedness (H2)",
       x = NULL, y = "Mean net spillover (%)", fill = NULL) + theme_pub()
save_plot(p4, "figure_04_net_directional_bar", 10, 7)

cat("[H2] figure_05_spillover_matrix ...\n")
sm <- conn$spillover
if (is.null(sm)) sm <- compute_spillover_matrix(conn$gfevd, sectors)
sm <- as.matrix(sm)
if (is.null(dimnames(sm))) dimnames(sm) <- list(sectors, sectors)
sm_df <- as.data.frame(as.table(sm)); names(sm_df) <- c("From", "To", "Value")
p5 <- ggplot(sm_df, aes(To, From, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "mako", name = "Spillover (%)") +
  labs(title = "Average Pairwise Spillover Matrix (H2)",
       subtitle = "Row = receiver (From), column = transmitter (To); diagonal removed",
       x = "Transmitter", y = "Receiver") + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p5, "figure_05_spillover_matrix", 9, 8)

cat("[H3] figure_06_centrality_evolution ...\n")
cent_flow <- net$centrality_tv %>%
  dplyr::select(Sector, Date, OutDegree, InDegree) %>%
  dplyr::mutate(Date = as.Date(Date)) %>%
  tidyr::pivot_longer(c(OutDegree, InDegree), names_to = "Flow", values_to = "Value") %>%
  dplyr::mutate(Flow = dplyr::recode(Flow,
                  OutDegree = "To others (transmission)",
                  InDegree  = "From others (absorption)"),
                SectorLab = lab_of(Sector))
p6 <- ggplot(cent_flow, aes(Date, Value, color = Flow)) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  facet_wrap(~ SectorLab, ncol = 4, scales = "free_y") +
  scale_color_manual(values = c("To others (transmission)" = "#e67e22",
                                "From others (absorption)" = "#2980b9")) +
  labs(title = "Directional Centrality by Sector: Transmission vs Absorption (H3)",
       subtitle = "Per-sector small multiples - orange = to others (out-degree); blue = from others (in-degree)",
       x = NULL, y = "Connectedness (%)", color = NULL) +
  theme_pub(11)
save_plot(p6, "figure_06_centrality_evolution", 14, 9)

cat("[H3] figure_07_centrality_heatmap ...\n")
cm_long <- net$centrality_avg %>%
  dplyr::select(Sector, OutDegree, InDegree, Betweenness, Closeness, Eigenvector, PageRank) %>%
  tidyr::pivot_longer(-Sector, names_to = "Metric", values_to = "Value") %>%
  dplyr::group_by(Metric) %>%
  dplyr::mutate(Z = as.numeric(scale(Value))) %>% dplyr::ungroup()
cm_long$SectorLab <- lab_of(cm_long$Sector)
p7 <- ggplot(cm_long, aes(Metric, SectorLab, fill = Z)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2980b9", mid = "white", high = BR_COL, midpoint = 0, name = "z-score") +
  labs(title = "Average Centrality Profile by Sector (H3)",
       subtitle = "Column-standardised (per-measure z-scores)", x = NULL, y = NULL) +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p7, "figure_07_centrality_heatmap", 9, 8)

cat("[H2/H3] figure_19_role_map ...\n")
role_df <- conn$directional %>% dplyr::group_by(Sector) %>%
  dplyr::summarise(To = mean(ToOthers, na.rm = TRUE), From = mean(FromOthers, na.rm = TRUE),
                   Net = mean(NetSpillover, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(SectorLab = lab_of(Sector), Total = To + From,
                Role = ifelse(Net > 0, "Net transmitter", "Net absorber"))
p19 <- ggplot(role_df, aes(From, To)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(aes(color = Role, size = Total), alpha = 0.85) +
  geom_text(aes(label = SectorLab), size = 3, vjust = -0.9, color = LINE_COL) +
  scale_color_manual(values = c("Net transmitter" = "#e67e22", "Net absorber" = "#2980b9")) +
  scale_size_continuous(range = c(3, 11), name = "Total degree (%)") +
  labs(title = "Sectoral Roles in the Spillover Network (H2/H3)",
       subtitle = "Above the 45-degree line = net transmitter; below = net absorber; point size = total connectedness",
       x = "From others (absorption, %)", y = "To others (transmission, %)", color = NULL) +
  theme_pub()
save_plot(p19, "figure_19_role_map", 10, 9)

cent_metrics <- c(
  OutDegree   = "Out-degree centrality - transmission to others",
  InDegree    = "In-degree centrality - absorption from others",
  NetDegree   = "Net directional centrality (transmitter minus absorber)",
  Betweenness = "Weighted betweenness centrality",
  Closeness   = "Out-closeness centrality",
  Eigenvector = "Eigenvector centrality",
  PageRank    = "PageRank centrality")
cent_idx <- c(OutDegree = 20L, InDegree = 21L, NetDegree = 22L, Betweenness = 23L,
              Closeness = 24L, Eigenvector = 25L, PageRank = 26L)
for (m in names(cent_metrics)) {
  cat(sprintf("[H3] figure_%02d_%s_by_sector ...\n", cent_idx[[m]], tolower(m)))
  fig_by_sector(net$centrality_tv, m,
    sprintf("%s by Sector Over Time (H3)", m), cent_metrics[[m]], "Centrality",
    sprintf("figure_%02d_%s_by_sector", cent_idx[[m]], tolower(m)),
    color = "#2c3e50")
}

if (!is.null(res) && !is.null(res$resilience_tv)) {
  cat("[H3] figure_27_resistance_by_sector ...\n")
  fig_by_sector(res$resilience_tv, "Resistance",
    "Resistance Resilience Index by Sector (H3)",
    "GIRF inverse worst-case impact (higher = more resistant)", "Resistance",
    "figure_27_resistance_by_sector", color = "#16a085")
  cat("[H3] figure_28_reconfiguration_by_sector ...\n")
  fig_by_sector(res$resilience_tv, "Reconfiguration",
    "Reconfiguration Index by Sector (H3)",
    "Long-run cumulative multiplier at horizon H", "Reconfiguration",
    "figure_28_reconfiguration_by_sector", color = "#8e44ad")
}

if (!is.null(res)) {
  cat("[H3] figure_08_resilience_dimensions ...\n")
  res_long <- res$resilience_tv %>%
    dplyr::select(Date, Resistance, RecoverySpeed, Reconfiguration) %>%
    dplyr::mutate(Date = as.Date(Date)) %>%
    tidyr::pivot_longer(-Date, names_to = "Dimension", values_to = "Value") %>%
    dplyr::group_by(Date, Dimension) %>%
    dplyr::summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")
  p8 <- ggplot(res_long, aes(Date, Value, color = Dimension)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ Dimension, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c(Resistance = "#16a085", RecoverySpeed = "#2980b9", Reconfiguration = "#e67e22")) +
    labs(title = "Resilience Dimensions Over Time (H3)",
         subtitle = "Cross-sector mean of GIRF-based resilience indicators",
         x = NULL, y = NULL) + theme_pub(12) + theme(legend.position = "none")
  save_plot(p8, "figure_08_resilience_dimensions", 11, 9)

  cat("[H3] figure_09_centrality_resilience_scatter ...\n")
  cr <- dplyr::inner_join(net$centrality_avg, res$resilience_avg, by = "Sector")
  if (nrow(cr) > 0) {
    cr$SectorLab <- lab_of(cr$Sector)
    p9 <- ggplot(cr, aes(OutDegree, Resistance)) +
      geom_smooth(method = "lm", se = TRUE, color = BR_COL, fill = "#f5b7b1", linewidth = 0.8) +
      geom_point(size = 3, color = LINE_COL) +
      geom_text(aes(label = SectorLab), vjust = -0.7, size = 3.1) +
      labs(title = "Centrality vs Resilience (H3)",
           subtitle = "Sector averages: out-degree centrality vs resistance",
           x = "Mean out-degree centrality", y = "Mean resistance") + theme_pub()
    save_plot(p9, "figure_09_centrality_resilience_scatter", 10, 7)

    cat("[H3] figure_29_h3_sector_profile_map ...\n")
    qd <- cr
    cx <- stats::median(qd$OutDegree,  na.rm = TRUE)
    cy <- stats::median(qd$Resistance, na.rm = TRUE)
    qd$Quadrant <- ifelse(qd$OutDegree >= cx & qd$Resistance >= cy, "Central & resilient",
                   ifelse(qd$OutDegree >= cx & qd$Resistance <  cy, "Central, less resilient",
                   ifelse(qd$OutDegree <  cx & qd$Resistance >= cy, "Peripheral, resilient",
                                                                    "Peripheral, less resilient")))
    rho <- suppressWarnings(stats::cor(qd$OutDegree, qd$Resistance,
                                       method = "spearman", use = "complete.obs"))
    p29 <- ggplot(qd, aes(OutDegree, Resistance)) +
      geom_vline(xintercept = cx, linetype = "dashed", color = "grey60") +
      geom_hline(yintercept = cy, linetype = "dashed", color = "grey60") +
      geom_point(aes(color = Quadrant), size = 3.5) +
      geom_text(aes(label = SectorLab), vjust = -0.8, size = 3.1, color = LINE_COL) +
      labs(title = "H3 Sector Profile Map: Centrality vs Resilience",
           subtitle = sprintf("Sector averages, split at medians; Spearman rho = %.2f (source: table_23_h3_sector_profile)", rho),
           x = "Mean out-degree centrality", y = "Mean resistance", color = NULL) +
      theme_pub() + theme(legend.position = "bottom")
    save_plot(p29, "figure_29_h3_sector_profile_map", 10, 7.5)
  }
}

cat("[H4] figure_11_tci_breaks_episodes ...\n")
ep <- episode_df()
p11 <- ggplot(tci, aes(Date, TotalConnectedness))
if (!is.null(breaks$TCI$regimes) && nrow(breaks$TCI$regimes) >= 1) {
  reg <- breaks$TCI$regimes
  reg$StartDate <- as.Date(reg$StartDate); reg$EndDate <- as.Date(reg$EndDate)
  reg$Shade <- ifelse(reg$Regime %% 2 == 0, "a", "b")
  p11 <- p11 + geom_rect(data = reg, inherit.aes = FALSE,
                         aes(xmin = StartDate, xmax = EndDate, ymin = -Inf, ymax = Inf, fill = Shade),
                         alpha = 0.12) +
    scale_fill_manual(values = c(a = "#7f8c8d", b = "#ffffff"), guide = "none")
}
p11 <- p11 + geom_line(color = LINE_COL, linewidth = 1.1)
p11 <- add_breaks(p11, bd_tci)
if (!is.null(ep) && nrow(ep) > 0) {
  ymax <- max(tci$TotalConnectedness, na.rm = TRUE)
  p11 <- p11 +
    geom_vline(data = ep, inherit.aes = FALSE, aes(xintercept = Date),
               linetype = "dotted", color = EP_COL, linewidth = 0.7) +
    geom_text(data = ep, inherit.aes = FALSE, aes(x = Date, y = ymax, label = Label),
              angle = 90, hjust = 1, vjust = -0.3, size = 2.8, color = EP_COL)
}
p11 <- p11 + labs(title = "Structural Reconfiguration of Connectedness (H4)",
                  subtitle = "Shaded: Bai-Perron regimes | red dashed: breaks | green dotted: full historical-episode timeline (context; H4 tests the pre-registered subset - see table_26/table_39)",
                  x = NULL, y = "Total connectedness (%)") + theme_pub()
save_plot(p11, "figure_11_tci_breaks_episodes", 13, 7)

make_break_panel <- function(comp, series_df, name, ylab) {
  if (is.null(comp) || is.null(series_df) || nrow(series_df) == 0) return(invisible(NULL))
  bd <- break_dates_of(comp)
  p <- ggplot(series_df, aes(Date, Value)) + geom_line(color = LINE_COL, linewidth = 0.8)
  if (!is.null(comp$regimes) && nrow(comp$regimes) >= 1) {
    reg <- comp$regimes; reg$StartDate <- as.Date(reg$StartDate); reg$EndDate <- as.Date(reg$EndDate)
    p <- p + geom_segment(data = reg, inherit.aes = FALSE,
                          aes(x = StartDate, xend = EndDate, y = Mean, yend = Mean),
                          color = BR_COL, linewidth = 1.1)
  }
  p <- add_breaks(p, bd)
  p <- p + labs(title = sprintf("Structural Breaks: %s (H4)", comp$label),
                subtitle = sprintf("%d break(s) detected; red segments = regime means", comp$n_breaks %||% 0L),
                x = NULL, y = ylab) + theme_pub()
  save_plot(p, name, 12, 6)
}

cat("[H4] figure_12_breaks_centrality ...\n")
cent_series <- net$centrality_tv %>% dplyr::group_by(TimeIndex) %>%
  dplyr::summarise(Date = as.Date(dplyr::first(Date)), Value = mean(OutDegree, na.rm = TRUE), .groups = "drop")
make_break_panel(breaks$Centrality, cent_series, "figure_12_breaks_centrality", "Mean out-degree centrality")

if (!is.null(res)) {
  cat("[H4] figure_13_breaks_resilience ...\n")
  res_series <- res$resilience_tv %>% dplyr::group_by(TimeIndex) %>%
    dplyr::summarise(Date = as.Date(dplyr::first(Date)), Value = mean(Resistance, na.rm = TRUE), .groups = "drop")
  make_break_panel(breaks$Resilience, res_series, "figure_13_breaks_resilience", "Mean resistance")
}

cat("[H4] figure_14_break_alignment ...\n")
if (!is.null(breaks)) {
  align_all <- dplyr::bind_rows(lapply(names(breaks), function(nm) {
    a <- breaks[[nm]]$alignment
    if (!is.null(a) && nrow(a) > 0) { a$Series <- breaks[[nm]]$label; a } else NULL
  }))
  if (!is.null(align_all) && nrow(align_all) > 0)
    align_all <- align_all[align_all$Episode != "(no match)", , drop = FALSE]
  if (!is.null(align_all) && nrow(align_all) > 0) {
    align_all$BreakDate <- as.Date(align_all$BreakDate)
    p14 <- ggplot(align_all, aes(BreakDate, Series, color = Episode)) +
      geom_point(size = 4) +
      geom_text(aes(label = Episode), vjust = -1, size = 3, show.legend = FALSE) +
      scale_color_viridis_d(option = "turbo") +
      labs(title = "Detected Breaks vs Historical Episodes (H4)",
           subtitle = sprintf("Each point = a detected break aligned to a pre-registered episode (within +/-%d year(s)); unmatched breaks omitted", CONFIG$BREAK_ALIGN_TOL),
           x = NULL, y = NULL) + theme_pub() + theme(legend.position = "none")
    save_plot(p14, "figure_14_break_alignment", 12, 5)
  }
}

if (!is.null(res) && !is.null(res$resilience_tv)) {
  reb_df <- res$resilience_tv
  reb_df$Date <- as.Date(reb_df$Date)
  reb_df$SectorLab <- vapply(reb_df$Sector, lab_of, character(1))
  p_reb <- ggplot2::ggplot(reb_df, ggplot2::aes(Date, RecoverySpeed)) +
    ggplot2::geom_line(color = FILL_COL, linewidth = 0.6, na.rm = TRUE) +
    ggplot2::facet_wrap(~ SectorLab, scales = "free_y") +
    theme_pub(11) +
    ggplot2::labs(title = "Recovery-speed resilience index by sector over time",
                  subtitle = "GIRF post-peak half-life of each sector's response (lower = faster recovery)",
                  x = NULL, y = "Recovery speed (half-life, periods)")
  save_plot(p_reb, "figure_18_recovery_by_sector", 13, 9)
}

has_igraph <- requireNamespace("igraph", quietly = TRUE)

if (has_igraph) {
  draw_network <- function(adj, sectors, title, q = 0.80) {
    diag(adj) <- 0
    pos <- adj[adj > 0]
    thr <- if (length(pos)) stats::quantile(pos, q, na.rm = TRUE) else 0
    a <- adj; a[a < thr] <- 0
    g <- igraph::graph_from_adjacency_matrix(a, mode = "directed", weighted = TRUE)
    out_s <- rowSums(adj); net_s <- rowSums(adj) - colSums(adj)
    rng <- range(out_s, finite = TRUE)
    rsc <- if (diff(rng) > 0) (out_s - rng[1]) / diff(rng) else rep(0.5, length(out_s))
    vsize <- 12 + 22 * rsc
    vcol  <- ifelse(net_s > 0, "#e67e22", "#2980b9")
    ew <- if (igraph::ecount(g) > 0) scales::rescale(igraph::E(g)$weight, to = c(0.4, 4.5)) else 1
    set.seed(1)
    igraph::plot.igraph(g, layout = igraph::layout_in_circle,
      vertex.size = vsize, vertex.color = vcol, vertex.frame.color = "white",
      vertex.label = sectors, vertex.label.cex = 0.85, vertex.label.color = "#1b2631",
      edge.width = ew, edge.color = grDevices::adjustcolor("#566573", 0.55),
      edge.arrow.size = 0.25, edge.curved = 0.15, main = title)
    graphics::legend("bottomleft", legend = c("Net transmitter", "Net absorber"),
                     pch = 21, pt.bg = c("#e67e22", "#2980b9"), pt.cex = 1.6,
                     bty = "n", cex = 0.9)
  }
  save_base_fig <- function(drawer, name, w = 1700, h = 1500, res = 200) {
    grDevices::png(file.path(FIG_DIR, paste0(name, ".png")), width = w, height = h, res = res)
    drawer(); grDevices::dev.off()
    tryCatch({
      grDevices::pdf(file.path(FIG_DIR, paste0(name, ".pdf")), width = w / res, height = h / res)
      drawer(); grDevices::dev.off()
    }, error = function(e) NULL)
    cat(sprintf("  saved %s.png\n", name))
  }

  avg_gfevd <- apply(conn$gfevd, c(1, 2), mean)
  adj_avg <- build_transmission_adjacency(avg_gfevd, sectors)
  cat("[NET] figure_15_network_average ...\n")
  save_base_fig(function() draw_network(adj_avg, sectors,
    "Average sectoral spillover network\n(orange = net transmitter, blue = net absorber)"),
    "figure_15_network_average")

  dts <- as.Date(conn$dates); yrs <- as.integer(format(dts, "%Y"))
  ep_years <- as.integer(names(CONFIG$HISTORICAL_EPISODES))
  ep_in <- ep_years[ep_years >= min(yrs) & ep_years <= max(yrs)]
  if (length(ep_in) >= 2) {
    cat("[NET] figure_16_network_evolution ...\n")
    save_base_fig(function() {
      old <- graphics::par(mfrow = c(ceiling(length(ep_in) / 2), 2), mar = c(1, 1, 3, 1))
      on.exit(graphics::par(old))
      for (ey in ep_in) {
        idx <- which.min(abs(yrs - ey))
        adj <- build_transmission_adjacency(conn$gfevd[, , idx], sectors)
        draw_network(adj, sectors,
          sprintf("%d - %s", ey, unname(CONFIG$HISTORICAL_EPISODES[as.character(ey)])))
      }
    }, "figure_16_network_evolution", w = 1700, h = 850 * ceiling(length(ep_in) / 2))
  }

  if (!is.null(breaks$TCI$regimes) && nrow(breaks$TCI$regimes) >= 2) {
    cat("[NET] figure_17_network_regimes ...\n")
    reg <- breaks$TCI$regimes
    save_base_fig(function() {
      nr <- nrow(reg)
      old <- graphics::par(mfrow = c(ceiling(nr / 2), 2), mar = c(1, 1, 3, 1))
      on.exit(graphics::par(old))
      for (r in seq_len(nr)) {
        w <- reg$StartIdx[r]:reg$EndIdx[r]
        adj <- build_transmission_adjacency(
          apply(conn$gfevd[, , w, drop = FALSE], c(1, 2), mean), sectors)
        draw_network(adj, sectors,
          sprintf("Regime %d (%s to %s)", reg$Regime[r],
                  format(as.Date(reg$StartDate[r]), "%Y"), format(as.Date(reg$EndDate[r]), "%Y")))
      }
    }, "figure_17_network_regimes", w = 1700, h = 850 * ceiling(nrow(breaks$TCI$regimes) / 2))
  }
} else {
  cat("[NET] igraph not available; skipped network diagrams.\n")
}

article_map <- list(
  "article_H1_densification"         = c("figure_01_tci", "figure_02_tci_regimes", "figure_10_network_density", "figure_15_network_average"),
  "article_H2_directional_roles"     = c("figure_03_net_spillover_ts", "figure_04_net_directional_bar", "figure_05_spillover_matrix", "figure_15_network_average", "figure_19_role_map"),
  "article_H3_centrality_resilience" = c("figure_06_centrality_evolution", "figure_07_centrality_heatmap", "figure_08_resilience_dimensions", "figure_09_centrality_resilience_scatter", "figure_29_h3_sector_profile_map", "figure_18_recovery_by_sector", "figure_19_role_map", "figure_20_outdegree_by_sector", "figure_21_indegree_by_sector", "figure_22_netdegree_by_sector", "figure_23_betweenness_by_sector", "figure_24_closeness_by_sector", "figure_25_eigenvector_by_sector", "figure_26_pagerank_by_sector", "figure_27_resistance_by_sector", "figure_28_reconfiguration_by_sector"),
  "article_H4_reconfiguration"       = c("figure_11_tci_breaks_episodes", "figure_12_breaks_centrality", "figure_13_breaks_resilience", "figure_14_break_alignment", "figure_16_network_evolution", "figure_17_network_regimes")
)
n_png <- length(list.files(FIG_DIR, pattern = "\\.png$"))
thesis_dir <- file.path(FIG_DIR, "thesis")
if (!dir.exists(thesis_dir)) dir.create(thesis_dir, recursive = TRUE)
for (f in list.files(FIG_DIR, pattern = "\\.(png|pdf)$", full.names = TRUE))
  file.copy(f, file.path(thesis_dir, basename(f)), overwrite = TRUE)
for (art in names(article_map)) {
  adir <- file.path(FIG_DIR, art)
  if (!dir.exists(adir)) dir.create(adir, recursive = TRUE)
  for (b in article_map[[art]]) for (ext in c("png", "pdf")) {
    src <- file.path(FIG_DIR, paste0(b, ".", ext))
    if (file.exists(src)) file.copy(src, file.path(adir, paste0(b, ".", ext)), overwrite = TRUE)
  }
}
cat(sprintf("\nFigures complete. %d top-level PNG files in %s (PDF companions alongside).\n", n_png, FIG_DIR))
cat("Organized into per-article folders (article_H1..H4) + thesis/.\n")
end_logging()
