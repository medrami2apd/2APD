#!/usr/bin/env Rscript
if (!exists(".PIPELINE_ROOT")) .PIPELINE_ROOT <- getwd()

.run_kappa <- !identical(Sys.getenv("RUN_KAPPA_SENSITIVITY", "1"), "0")

stages <- c(
  "scripts/01_data.R",
  "scripts/01b_descriptive.R",
  "scripts/01c_lag_selection.R",
  "scripts/02_estimation.R",
  "scripts/03_connectedness.R",
  "scripts/04_network.R",
  "scripts/05_resilience.R",
  "scripts/06_hypotheses.R",
  "scripts/07_breaks.R",
  "scripts/08_figures.R",
  "scripts/09_tables.R",
  "scripts/10_robustness_evaluation.R",
  "scripts/11_draws_inference.R",
  "scripts/11b_credible_bands.R",
  if (.run_kappa) "scripts/12_kappa_sensitivity.R" else NULL
)

for (s in stages) {
  cat(sprintf("\n================ RUNNING %s ================\n", s))
  ok <- tryCatch({ source(s, local = new.env()); TRUE },
                 error = function(e) { message(sprintf("STAGE FAILED (%s): %s", s, e$message)); FALSE })
  if (!ok) { message("Halting pipeline."); quit(status = 1L) }
}
cat("\nPipeline complete. See outputs/.\n")
