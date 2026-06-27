#!/usr/bin/env Rscript
source("scripts/00_setup/setup.R")
start_logging("12_kappa_sensitivity")

stopifnot(
  exists("compute_gfevd_array"), exists("compute_tci"),
  exists("compute_network_tv"), exists(".hac_trend"), exists("test_h2"),
  exists("detect_breaks"), exists("tv_sv_relevance"), exists("export_table"))

grid <- CONFIG$KAPPA_Q_GRID
sens <- CONFIG$SENSITIVITY
base_paths <- CONFIG$PATHS
rows <- list()

SENS_ROOT <- CONFIG$PATHS$sensitivity %||% "outputs/sensitivity"
dir.create(SENS_ROOT, recursive = TRUE, showWarnings = FALSE)

for (kap in grid) {
  tag   <- gsub("\\.", "p", format(kap, trim = TRUE))
  sroot <- file.path(SENS_ROOT, paste0("kappa_", tag))
  dir.create(sroot, recursive = TRUE, showWarnings = FALSE)
  paths <- base_paths
  paths$posterior       <- file.path(sroot, "posterior_mean_arrays.rds")
  paths$posterior_draws <- file.path(sroot, "posterior_draws.rds")
  paths$fit             <- file.path(sroot, "tvpvar_fit.rds")

  cat(sprintf("\n================ KAPPA_Q = %g ================\n", kap))
  .CONFIG_OVERRIDE <<- list(KAPPA_Q = kap, N_ITER = sens$N_ITER,
                            N_BURN = sens$N_BURN, THIN = sens$THIN, PATHS = paths)
  ok <- tryCatch({ source("scripts/02_estimation.R", local = new.env()); TRUE },
                 error = function(e) { message("estimation failed: ", e$message); FALSE })
  rm(.CONFIG_OVERRIDE, envir = .GlobalEnv)
  if (!ok) next

  post <- readRDS(paths$posterior)
  g    <- compute_gfevd_array(post$beta_array, post$sigma_array, post$k, post$p,
                              CONFIG$HORIZON_GFEVD, CONFIG$EPS)
  tci  <- compute_tci(g)$TotalConnectedness
  net  <- compute_network_tv(g, post$sectors, post$dates,
                             threshold_type = CONFIG$NETWORK_THRESHOLD_TYPE,
                             threshold_value = CONFIG$NETWORK_THRESHOLD)
  h1 <- .hac_trend(tci)
  h2     <- test_h2(net$centrality_tv, neutral_band = CONFIG$ROLE_NEUTRAL_BAND, B = CONFIG$BOOT_B)
  h2_tbl <- if (is.list(h2) && !is.data.frame(h2) && !is.null(h2$table)) h2$table else h2
  p_stay <- h2_tbl$Value[h2_tbl$Metric == "Mean P(stay)"]
  perm_p <- h2_tbl$Value[h2_tbl$Metric == "Permutation p-value"]
  bp <- detect_breaks(tci, post$dates, CONFIG$BREAK_MAX, CONFIG$BREAK_MIN_SIZE)
  tvr <- mean(tv_sv_relevance(post, post$sectors)$TV_ratio, na.rm = TRUE)

  rows[[length(rows) + 1L]] <- tibble::tibble(
    KAPPA_Q = kap, Mean_TCI = mean(tci),
    H1_trend_slope = unname(h1["slope"]), H1_trend_p = unname(h1["p"]),
    H2_P_stay = p_stay, H2_perm_p = perm_p,
    H4_n_breaks = bp$n_breaks, TV_ratio = round(tvr, 4))
}

sens_tbl <- dplyr::bind_rows(rows)
export_table(sens_tbl, "table_45_kappa_sensitivity")
kappa_rds <- file.path(dirname(CONFIG$PATHS$breaks), "kappa_sensitivity.rds")
dir.create(dirname(kappa_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(sens_tbl, kappa_rds)
cat("\n=== KAPPA_Q SENSITIVITY (H1/H2/H4 across the drift-prior grid) ===\n")
print(sens_tbl)
cat("\nIf H1 slope sign, H2 persistence, or H4 break count flip across the grid,\n")
cat("the corresponding finding is prior-driven, not data-driven.\n")

end_logging()
