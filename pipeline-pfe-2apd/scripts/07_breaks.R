#!/usr/bin/env Rscript
source("scripts/00_setup/setup.R")
start_logging("07_breaks")

conn <- readRDS(CONFIG$PATHS$gfevd)
net  <- readRDS(CONFIG$PATHS$network)
res  <- readRDS(CONFIG$PATHS$resilience)

prereg_episodes <- CONFIG$HISTORICAL_EPISODES[as.character(CONFIG$PREREGISTERED_EPISODES)]
run_break <- function(series, dates, label) {
  bp <- detect_breaks(series, dates, CONFIG$BREAK_MAX, CONFIG$BREAK_MIN_SIZE)
  hits <- align_to_history(bp$break_dates, prereg_episodes, tol = CONFIG$BREAK_ALIGN_TOL)
  list(label = label, n_breaks = bp$n_breaks, regimes = bp$regimes, alignment = hits)
}

results <- list()
results$TCI <- run_break(conn$total_conn$TotalConnectedness, conn$dates, "TCI")

cent_series <- net$centrality_tv %>% dplyr::group_by(TimeIndex) %>%
  dplyr::summarise(Date = dplyr::first(Date), v = mean(OutDegree, na.rm = TRUE), .groups = "drop")
results$Centrality <- run_break(cent_series$v, cent_series$Date, "Mean OutDegree")

res_series <- res$resilience_tv %>% dplyr::group_by(TimeIndex) %>%
  dplyr::summarise(Date = dplyr::first(Date), v = mean(Resistance, na.rm = TRUE), .groups = "drop")
results$Resilience <- run_break(res_series$v, res_series$Date, "Mean Resistance")

align_all <- dplyr::bind_rows(lapply(names(results), function(nm) {
  a <- results[[nm]]$alignment; if (nrow(a) > 0) a$Series <- results[[nm]]$label; a
}))
if (nrow(align_all) > 0) export_table(align_all, "table_25_break_alignment")

summary_tbl <- tibble::tibble(
  Series = vapply(results, function(x) x$label, character(1)),
  N_Breaks = vapply(results, function(x) x$n_breaks, integer(1)))
export_table(summary_tbl, "table_24_break_summary")

n_series           <- length(results)
total_breaks       <- sum(vapply(results, function(x) x$n_breaks, integer(1)))
series_with_breaks <- sum(vapply(results, function(x) as.integer(x$n_breaks > 0), integer(1)))
aligned <- sum(vapply(names(results), function(nm) {
  a <- results[[nm]]$alignment
  if (is.null(a) || nrow(a) == 0) 0L else sum(a$Episode != "(no match)")
}, integer(1)))
episodes_hit <- unique(unlist(lapply(results, function(x) {
  a <- x$alignment
  if (is.null(a) || nrow(a) == 0) character(0) else a$Episode[a$Episode != "(no match)"]
})))
verdict <- if (total_breaks == 0) "Not supported (no breaks detected)" else
           if (aligned > 0) "Supported (breaks align to historical episodes)" else
           "Partial (breaks present, none aligned)"
h4_tbl <- tibble::tibble(
  Hypothesis = "H4: Structural Reconfiguration",
  Metric = c("Total structural breaks", "Series exhibiting breaks",
             "Breaks aligned to historical episodes", "H4 verdict"),
  Value  = c(total_breaks, series_with_breaks, aligned, aligned),
  Result = c(sprintf("across %d monitored series", n_series),
             sprintf("of %d series", n_series),
             if (length(episodes_hit)) paste(episodes_hit, collapse = "; ") else "none within tolerance",
             verdict))
export_table(h4_tbl, "table_26_h4_reconfiguration")
cat("\n--- H4 ---\n"); print(h4_tbl)

saveRDS(results, CONFIG$PATHS$breaks)
cat("Structural breaks done.\n"); print(summary_tbl)

end_logging()
