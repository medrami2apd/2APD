parse_period <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) {
    if (all(x >= 1700 & x <= 2300, na.rm = TRUE)) return(as.Date(paste0(x, "-12-31")))
    return(as.Date(x, origin = "1899-12-30"))
  }
  s <- trimws(as.character(x))
  if (any(grepl("[Qq]", s))) {
    yr <- as.integer(sub(".*?(\\d{4}).*", "\\1", s))
    qq <- as.integer(sub(".*[Qq]\\s*([1-4]).*", "\\1", s))
    mon <- (qq - 1L) * 3L + 1L
    return(as.Date(sprintf("%04d-%02d-01", yr, mon)))
  }
  if (all(grepl("^\\d{4}$", s)))
    return(as.Date(paste0(s, "-12-31")))
  if (all(grepl("^\\d{4}[-/]\\d{1,2}$", s)))
    return(as.Date(paste0(sub("/", "-", s), "-01")))
  d <- suppressWarnings(as.Date(s))
  if (all(is.na(d))) stop("Could not parse period/date column into dates.")
  d
}

load_raw <- function(path, sectors) {
  if (!file.exists(path)) stop(sprintf("Raw data not found: %s", path))
  df <- readr::read_csv(path, show_col_types = FALSE)
  date_col <- intersect(c("Date","date","DATE","Year","year","Period","period",
                          "Quarter","quarter","Month","month","Time","time"),
                        names(df))[1]
  if (is.na(date_col)) stop(sprintf("No Date/Year/Period column found. Columns: %s",
                                    paste(names(df), collapse = ", ")))
  missing <- setdiff(sectors, names(df))
  if (length(missing) > 0) stop(sprintf("Missing sector columns: %s", paste(missing, collapse=", ")))
  dates <- parse_period(df[[date_col]])
  list(dates = as.Date(dates), Y = as.matrix(df[, sectors, drop = FALSE]))
}

to_growth <- function(Y, log_diff = TRUE) {
  if (log_diff && all(Y > 0, na.rm = TRUE)) {
    g <- diff(log(Y)) * 100
  } else {
    g <- diff(Y)
  }
  g
}

stationarity_report <- function(Ymat, sectors) {
  has_urca <- requireNamespace("urca", quietly = TRUE)
  out <- lapply(seq_along(sectors), function(j) {
    x <- Ymat[, j]; x <- x[is.finite(x)]
    if (!has_urca || length(x) < 10) return(data.frame(Sector=sectors[j], ADF_stat=NA, KPSS_stat=NA))
    adf <- tryCatch(urca::ur.df(x, type="drift", selectlags="AIC")@teststat[1], error=function(e) NA)
    kpss <- tryCatch(urca::ur.kpss(x, type="mu")@teststat[1], error=function(e) NA)
    data.frame(Sector=sectors[j], ADF_stat=adf, KPSS_stat=kpss)
  })
  dplyr::bind_rows(out)
}

prepare_data <- function(cfg, standardize = TRUE) {
  raw <- load_raw(cfg$PATHS$raw_data, cfg$SECTORS)
  g <- to_growth(raw$Y)
  g_dates <- raw$dates[-1]
  keep <- stats::complete.cases(g)
  Y <- g[keep, , drop = FALSE]
  dates <- g_dates[keep]
  if (!is.null(cfg$SUBSET_START) || !is.null(cfg$SUBSET_END)) {
    lo <- if (is.null(cfg$SUBSET_START)) min(dates) else as.Date(cfg$SUBSET_START)
    hi <- if (is.null(cfg$SUBSET_END))   max(dates) else as.Date(cfg$SUBSET_END)
    sel <- dates >= lo & dates <= hi
    Y <- Y[sel, , drop = FALSE]; dates <- dates[sel]
    message(sprintf("Sub-sample filter applied: %s to %s (%d obs)", lo, hi, nrow(Y)))
  }
  center <- colMeans(Y); scale <- apply(Y, 2, stats::sd)
  scale[scale == 0 | !is.finite(scale)] <- 1
  if (standardize) Y <- sweep(sweep(Y, 2, center, "-"), 2, scale, "/")
  out <- list(Y = Y, dates = dates, sectors = cfg$SECTORS,
              center = center, scale = scale,
              stationarity = stationarity_report(Y, cfg$SECTORS))
  saveRDS(out, cfg$PATHS$processed)
  message(sprintf("Prepared data: %d obs x %d sectors", nrow(Y), ncol(Y)))
  out
}
