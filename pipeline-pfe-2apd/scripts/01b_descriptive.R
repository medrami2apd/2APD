#!/usr/bin/env Rscript
source("scripts/00_setup/setup.R")
start_logging("01b_descriptive")

suppressWarnings(suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble); library(scales)
}))

FIG_DIR <- CONFIG$PATHS$figures
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
sectors <- CONFIG$SECTORS; labels <- CONFIG$SECTOR_LABELS
lab_of  <- function(s) ifelse(s %in% names(labels), unname(labels[s]), s)

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
melt_mat <- function(mat, dates) {
  df <- as.data.frame(mat); df$Date <- dates
  tidyr::pivot_longer(df, -Date, names_to = "Sector", values_to = "Value") %>%
    dplyr::mutate(SectorLab = lab_of(Sector))
}

raw  <- load_raw(CONFIG$PATHS$raw_data, sectors)
prep <- readRDS(CONFIG$PATHS$processed)

levels_mat <- as.matrix(raw$Y); colnames(levels_mat) <- sectors
lvl_dates  <- as.Date(raw$dates)
growth_mat <- sweep(sweep(prep$Y, 2, prep$scale, "*"), 2, prep$center, "+")
colnames(growth_mat) <- sectors
g_dates <- as.Date(prep$dates)

cat("figure_d1_sectoral_levels ...\n")
lvl_long <- melt_mat(levels_mat, lvl_dates)
pD1 <- ggplot(lvl_long, aes(Date, Value, color = Sector)) +
  geom_line(linewidth = 0.5) + scale_color_viridis_d(option = "turbo") +
  labs(title = "Sectoral Value Added in Morocco (1980-2025)",
       x = NULL, y = "Value added (levels)", color = NULL) +
  theme_pub() + theme(legend.position = "right")
save_plot(pD1, "figure_d1_sectoral_levels", 12, 7)

cat("figure_d2_growth_rates ...\n")
g_long <- melt_mat(growth_mat, g_dates)
pD2 <- ggplot(g_long, aes(Date, Value, color = Sector)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  facet_wrap(~ SectorLab, scales = "free_y", ncol = 4) +
  labs(title = "Sectoral Growth Rates (%)",
       subtitle = "Morocco; log-difference of value added",
       x = NULL, y = "Growth (%)") + theme_pub(11) + theme(legend.position = "none")
save_plot(pD2, "figure_d2_growth_rates", 14, 9)

cat("figure_d3_sector_shares ...\n")
tot    <- rowSums(levels_mat)
shares <- sweep(levels_mat, 1, tot, "/") * 100
sh_long <- melt_mat(shares, lvl_dates)
pD3 <- ggplot(sh_long, aes(Date, Value, fill = SectorLab)) +
  geom_area(position = "stack", alpha = 0.85) + scale_fill_viridis_d(option = "D") +
  labs(title = "Sectoral Composition of Value Added",
       x = NULL, y = "Share of total (%)", fill = NULL) +
  theme_pub() + theme(legend.position = "right")
save_plot(pD3, "figure_d3_sector_shares", 12, 7)

cat("figure_d4_correlation_heatmap ...\n")
cor_mat <- cor(growth_mat, use = "complete.obs")
cor_df  <- as.data.frame(as.table(cor_mat)); names(cor_df) <- c("S1", "S2", "R")
pD4 <- ggplot(cor_df, aes(S2, S1, fill = R)) + geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", R)), size = 2.6, color = "grey15") +
  scale_fill_gradient2(low = "#2980b9", mid = "white", high = "#c0392b",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  labs(title = "Cross-Sector Growth Correlations", x = NULL, y = NULL) +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(pD4, "figure_d4_correlation_heatmap", 10, 9)
export_matrix(cor_mat, "table_02_correlation_matrix")

cat("table_01_summary_statistics ...\n")
summ <- data.frame(
  Sector = sectors, Label = unname(lab_of(sectors)),
  Obs  = nrow(levels_mat),
  Mean = round(colMeans(levels_mat, na.rm = TRUE), 2),
  SD   = round(apply(levels_mat, 2, sd,  na.rm = TRUE), 2),
  Min  = round(apply(levels_mat, 2, min, na.rm = TRUE), 2),
  Max  = round(apply(levels_mat, 2, max, na.rm = TRUE), 2))
export_table(summ, "table_01_summary_statistics")

cat("table_03_sector_shares_by_decade ...\n")
sh_df <- as.data.frame(shares)
sh_df$Decade <- paste0(floor(as.integer(format(lvl_dates, "%Y")) / 10) * 10, "s")
shares_decade <- sh_df %>% dplyr::group_by(Decade) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(sectors), ~ round(mean(.x, na.rm = TRUE), 2)),
                   .groups = "drop")
export_table(shares_decade, "table_03_sector_shares_by_decade")

cat("table_04_growth_volatility ...\n")
vol <- g_long %>% dplyr::group_by(Sector) %>%
  dplyr::summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE),
                   CV = SD / abs(Mean) * 100,
                   Min = min(Value, na.rm = TRUE), Max = max(Value, na.rm = TRUE),
                   .groups = "drop") %>% dplyr::arrange(dplyr::desc(SD))
export_table(vol, "table_04_growth_volatility")

cat("Descriptive analysis complete.\n")
end_logging()
