options(stringsAsFactors = FALSE, scipen = 999)

if (!exists(".PIPELINE_ROOT")) {
  .PIPELINE_ROOT <- tryCatch(normalizePath(file.path(dirname(sys.frame(1)$ofile), "..", "..")),
                             error = function(e) getwd())
}
setwd(.PIPELINE_ROOT)

source("scripts/00_setup/config.R")
set.seed(CONFIG$SEED)

.required <- c(
  "tibble","dplyr","tidyr","readr","purrr","stringr",
  "vars","strucchange","zoo",
  "MASS","posterior",
  "igraph",
  "plm","lmtest","sandwich",
  "Matrix","expm","matrixStats",
  "ggplot2","scales","openxlsx","xtable"
)
.load_pkg <- function(p) suppressWarnings(suppressMessages(
  requireNamespace(p, quietly = TRUE) && library(p, character.only = TRUE, logical.return = TRUE)))
for (p in .required) if (!.load_pkg(p)) warning(sprintf("Package not available: %s", p))

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

.func_files <- c(
  "functions/data_functions.R",
  "functions/connectedness_functions.R",
  "functions/network_functions.R",
  "functions/resilience_functions.R",
  "functions/hypothesis_functions.R",
  "functions/breaks_functions.R",
  "functions/draws_functions.R",
  "functions/robustness_functions.R",
  "functions/export_functions.R",
  "functions/plotting_functions.R"
)
for (f in .func_files) if (file.exists(f)) source(f) else warning(sprintf("Missing function file: %s", f))

for (d in c("data/raw","data/processed",
            "outputs/processed","outputs/posterior","outputs/connectedness",
            "outputs/network","outputs/resilience","outputs/diagnostics",
            "outputs/hypotheses","outputs/tables","outputs/tables/latex",
            "outputs/figures","outputs/logs")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat(sprintf("[setup] root=%s  k=%d  p=%d  mode=%s\n",
            .PIPELINE_ROOT, CONFIG$K, CONFIG$P_LAGS, CONFIG$UNCERTAINTY_MODE))
