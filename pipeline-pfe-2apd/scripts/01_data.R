#!/usr/bin/env Rscript
source("scripts/00_setup/setup.R")
start_logging("01_data")

if (!file.exists(CONFIG$PATHS$raw_data)) {
  stop(sprintf(paste0("Raw data not found at %s.\n",
    "Provide a CSV with a Date/Year column and one column per sector: %s"),
    CONFIG$PATHS$raw_data, paste(CONFIG$SECTORS, collapse = ", ")))
}

prep <- prepare_data(CONFIG, standardize = TRUE)

export_table(prep$stationarity, "table_01_stationarity")
cat(sprintf("Saved processed data: %d obs x %d sectors -> %s\n",
            nrow(prep$Y), ncol(prep$Y), CONFIG$PATHS$processed))

end_logging()
