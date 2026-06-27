#!/usr/bin/env Rscript
source("scripts/00_setup/setup.R")
start_logging("01c_lag_selection")

prep <- readRDS(CONFIG$PATHS$processed)
lag_max <- 8L

if (!requireNamespace("vars", quietly = TRUE)) {
  warning("Package 'vars' not available; skipping lag selection.")
} else {
  vs <- vars::VARselect(prep$Y, lag.max = lag_max, type = "const")
  crit <- as.data.frame(t(vs$criteria))
  names(crit) <- c("AIC", "HQ", "SC_BIC", "FPE")
  crit$Lag <- seq_len(nrow(crit))
  crit <- crit[, c("Lag", "AIC", "HQ", "SC_BIC", "FPE")]
  sel <- vs$selection
  crit$Selected_by <- vapply(crit$Lag, function(L) {
    paste(names(sel)[sel == L], collapse = "+")
  }, character(1))
  export_table(crit, "table_43_lag_selection")

  cat("\n=== VAR LAG SELECTION (information criteria) ===\n")
  print(crit)
  cat(sprintf("\nAIC->p=%d  HQ->p=%d  SC/BIC->p=%d  FPE->p=%d   (CONFIG$P_LAGS = %d)\n",
              sel["AIC(n)"], sel["HQ(n)"], sel["SC(n)"], sel["FPE(n)"], CONFIG$P_LAGS))
  if (!identical(as.integer(sel["SC(n)"]), as.integer(CONFIG$P_LAGS)))
    cat("NOTE: BIC-preferred lag differs from CONFIG$P_LAGS; report p=2 as robustness.\n")
}

end_logging()
