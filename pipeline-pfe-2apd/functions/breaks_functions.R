detect_breaks <- function(series, dates = NULL, max_breaks = 5L,
                          min_size = 0.10) {
  n <- length(series)
  ts_df <- data.frame(y = series, t = seq_len(n))
  bp <- tryCatch(strucchange::breakpoints(y ~ t, data = ts_df,
                                          h = min_size, breaks = max_breaks),
                 error = function(e) NULL)
  if (is.null(bp)) return(list(n_breaks = 0L, break_indices = integer(0),
                               break_dates = NULL, regimes = NULL, bic_path = NULL))

  bic_path <- tryCatch(summary(bp)$RSS["BIC", ], error = function(e) NULL)
  if (is.null(bic_path) || all(is.na(bic_path))) {
    bic_path <- tryCatch(stats::BIC(bp), error = function(e) NULL)
  }
  if (is.null(bic_path) || all(is.na(bic_path))) {
    opt <- bp
    bi <- stats::na.omit(bp$breakpoints)
  } else {
    n_opt <- as.integer(names(which.min(bic_path))) %||% 0L
    if (is.na(n_opt)) n_opt <- which.min(bic_path) - 1L
    bi <- if (n_opt >= 1) strucchange::breakpoints(bp, breaks = n_opt)$breakpoints else integer(0)
    bi <- stats::na.omit(bi)
  }

  bi <- as.integer(bi[is.finite(bi)])
  n_breaks <- length(bi)
  break_dates <- if (!is.null(dates) && n_breaks > 0) dates[bi] else NULL

  bounds <- c(0L, bi, n)
  regimes <- do.call(rbind, lapply(seq_len(length(bounds) - 1), function(r) {
    start_i <- bounds[r] + 1L
    end_i   <- bounds[r + 1L]
    data.frame(Regime = r, StartIdx = start_i, EndIdx = end_i,
               StartDate = if (!is.null(dates)) dates[start_i] else NA,
               EndDate   = if (!is.null(dates)) dates[end_i] else NA,
               Mean = mean(series[start_i:end_i], na.rm = TRUE))
  }))

  list(n_breaks = n_breaks, break_indices = bi, break_dates = break_dates,
       regimes = regimes, bic_path = bic_path)
}

align_to_history <- function(break_dates, episodes, tol = 1L) {
  if (is.null(break_dates) || length(break_dates) == 0)
    return(tibble::tibble(BreakDate = as.Date(character()), Year = integer(),
                          Episode = character(), Distance = integer()))
  ep_years <- as.integer(names(episodes))
  out <- lapply(break_dates, function(d) {
    yr <- as.integer(format(as.Date(d), "%Y"))
    dist <- abs(ep_years - yr); j <- which.min(dist)
    tibble::tibble(BreakDate = as.Date(d), Year = yr,
                   Episode = if (dist[j] <= tol) unname(episodes[j]) else "(no match)",
                   Distance = dist[j])
  })
  dplyr::bind_rows(out)
}
