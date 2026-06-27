#!/usr/bin/env Rscript
source("scripts/00_setup/setup.R")
start_logging("03_connectedness")

stopifnot(
  exists("compute_gfevd_array"), exists("companion_spectral_radius"),
  exists("compute_tci"), exists("compute_directional"),
  exists("compute_spillover_matrix"), exists("export_table"), exists("export_matrix"))

post <- readRDS(CONFIG$PATHS$posterior)
k <- post$k; p <- post$p; sectors <- post$sectors
beta_array <- post$beta_array; sigma_array <- post$sigma_array
dates <- post$dates; n_eff <- post$n_eff

gfevd_array <- compute_gfevd_array(beta_array, sigma_array, k, p,
                                   CONFIG$HORIZON_GFEVD, CONFIG$EPS)

rad <- vapply(1:n_eff, function(t) companion_spectral_radius(beta_array[t, , ], k, p), numeric(1))
unstable <- rad >= 1
n_unstable <- sum(unstable)
if (n_unstable > 0)
  warning(sprintf("%d/%d periods have spectral radius >= 1 (IRF measures unreliable there).",
                  n_unstable, n_eff))

diag_dir <- dirname(CONFIG$PATHS$breaks)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  data.frame(TimeIndex = 1:n_eff, Date = dates, SpectralRadius = rad, Unstable = unstable),
  file.path(diag_dir, "companion_stability.csv"), row.names = FALSE)

total_conn  <- compute_tci(gfevd_array)
directional <- compute_directional(gfevd_array, sectors)
spillover   <- compute_spillover_matrix(gfevd_array, sectors)

export_table(total_conn, "table_05_total_connectedness")
export_matrix(spillover, "table_06_spillover_matrix")

conn_dir <- dirname(CONFIG$PATHS$gfevd)
dir.create(conn_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(directional, file.path(conn_dir, "directional.csv"), row.names = FALSE)

dir.create(dirname(CONFIG$PATHS$gfevd), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(gfevd = gfevd_array, total_conn = total_conn,
             directional = directional, spillover = spillover,
             dates = dates, sectors = sectors, k = k, p = p,
             spectral_radius = rad, unstable = unstable, n_unstable = n_unstable),
        CONFIG$PATHS$gfevd)
cat(sprintf("Connectedness done. Mean TCI = %.2f%% | unstable periods: %d/%d\n",
            mean(total_conn$TotalConnectedness), n_unstable, n_eff))

end_logging()
