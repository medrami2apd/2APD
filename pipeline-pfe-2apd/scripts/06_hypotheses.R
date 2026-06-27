#!/usr/bin/env Rscript
source("scripts/00_setup/setup.R")
start_logging("06_hypotheses")

stopifnot(
  exists("test_h1"), exists("test_h2"), exists("test_h3"),
  exists("h3_sector_cross_section"), exists("export_table"))

prep <- readRDS(CONFIG$PATHS$processed)
conn <- readRDS(CONFIG$PATHS$gfevd)
net  <- readRDS(CONFIG$PATHS$network)
res  <- readRDS(CONFIG$PATHS$resilience)

h1 <- test_h1(conn$total_conn$TotalConnectedness, Y = prep$Y, p = CONFIG$P_LAGS,
              horizon = CONFIG$HORIZON_GFEVD, B = CONFIG$BOOT_B,
              block_len = CONFIG$BLOCK_LEN, eps = CONFIG$EPS)
export_table(h1$table, "table_20_h1_densification")

h2 <- test_h2(net$centrality_tv, neutral_band = CONFIG$ROLE_NEUTRAL_BAND, B = CONFIG$BOOT_B)
export_table(h2$table, "table_21_h2_persistence")

h3 <- test_h3(net$centrality_tv, res$resilience_tv,
              centrality_vars = c("OutDegree","InDegree","Betweenness","Closeness"),
              resilience_vars = c("Resistance","RecoverySpeed","Reconfiguration"),
              cluster = CONFIG$PANEL_CLUSTER, padj = CONFIG$MULTIPLE_TEST)
export_table(h3, "table_22_h3_panel")

h3_xs <- h3_sector_cross_section(
  net$centrality_avg, res$resilience_avg,
  centrality_vars = c("OutDegree","InDegree","Betweenness","Closeness"),
  resilience_vars = c("Resistance","RecoverySpeed","Reconfiguration"))
export_table(h3_xs$profile,     "table_23_h3_sector_profile")
export_table(h3_xs$correlation, "table_23b_h3_cross_sectional_corr")

hyp_path <- file.path(CONFIG$PATHS$hypotheses, "hypothesis_results.rds")
dir.create(dirname(hyp_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(h1 = h1, h2 = h2, h3 = h3,
             h3_sector_profile = h3_xs$profile,
             h3_cross_sectional_corr = h3_xs$correlation),
        hyp_path)

cat("\n--- H1 ---\n"); print(h1$table); cat(sprintf("H1 verdict: %s\n", h1$verdict))
cat("\n--- H2 ---\n"); print(h2$table); cat(sprintf("H2 verdict: %s\n", h2$verdict))
cat("\n--- H3 ---\n"); print(h3)
cat("\n--- H3 per-sector profile ---\n"); print(h3_xs$profile)
cat("\n--- H3 cross-sectional correlation ---\n"); print(h3_xs$correlation)

end_logging()
