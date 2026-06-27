.hac_trend <- function(y) {
  n <- length(y); d <- data.frame(y = as.numeric(y), t = seq_len(n))
  d <- d[is.finite(d$y), , drop = FALSE]
  if (nrow(d) < 10L) return(c(slope = NA_real_, p = NA_real_))
  m <- stats::lm(y ~ t, data = d)
  V <- tryCatch(sandwich::NeweyWest(m, prewhite = FALSE, adjust = TRUE),
                error = function(e) sandwich::vcovHC(m, type = "HC1"))
  se <- sqrt(diag(V))["t"]; sl <- stats::coef(m)["t"]
  if (!is.finite(se) || se <= 0) return(c(slope = unname(sl), p = NA_real_))
  c(slope = unname(sl), p = unname(2 * stats::pnorm(-abs(sl / se))))
}

.panel_coef <- function(df, dv, iv) {
  keep <- is.finite(df[[dv]]) & is.finite(df[[iv]])
  d <- df[keep, c("Sector", "TimeIndex", dv, iv), drop = FALSE]
  if (nrow(d) < 12L || length(unique(d$Sector)) < 3L ||
      length(unique(d$TimeIndex)) < 3L) return(NA_real_)
  sdx <- stats::sd(d[[iv]]); sdy <- stats::sd(d[[dv]])
  if (!is.finite(sdx) || !is.finite(sdy) || sdx == 0 || sdy == 0) return(NA_real_)
  d$.y <- (d[[dv]] - mean(d[[dv]])) / sdy
  d$.x <- (d[[iv]] - mean(d[[iv]])) / sdx
  pd <- tryCatch(plm::pdata.frame(d, index = c("Sector", "TimeIndex")),
                 error = function(e) NULL)
  if (is.null(pd)) return(NA_real_)
  m <- tryCatch(plm::plm(.y ~ .x, data = pd, model = "within", effect = "twoways"),
                error = function(e) NULL)
  if (is.null(m)) return(NA_real_)
  tryCatch(unname(stats::coef(m)[".x"]), error = function(e) NA_real_)
}

propagate_draws <- function(prep, draws, conn_dates, sectors, k, p, cfg,
                            centrality_vars = c("OutDegree","InDegree","Betweenness","Closeness"),
                            resilience_vars = c("Resistance","RecoverySpeed","Reconfiguration")) {
  stopifnot(!is.null(draws$beta_draws), !is.null(draws$sigma_draws))
  D     <- dim(draws$beta_draws)[1]
  n_eff <- dim(draws$beta_draws)[2]
  nprop <- min(D, as.integer(cfg$N_DRAWS_PROP %||% 500L))
  sel   <- unique(round(seq(1, D, length.out = nprop)))
  nprop <- length(sel)
  dates <- as.Date(conn_dates)

  thr_type   <- cfg$NETWORK_THRESHOLD_TYPE %||% "none"
  thr_val    <- cfg$NETWORK_THRESHOLD %||% 0
  resist_h   <- cfg$RESIST_HORIZON %||% 4L
  resist_eps <- cfg$RESIST_EPS %||% 1e-2
  hg         <- cfg$HORIZON_GFEVD; hr <- cfg$HORIZON_GIRF
  eps        <- cfg$EPS %||% 1e-8

  pairs   <- expand.grid(Resilience = resilience_vars, Centrality = centrality_vars,
                         stringsAsFactors = FALSE)
  tci_mat  <- matrix(NA_real_, nprop, n_eff)
  h1_slope <- rep(NA_real_, nprop)
  h3_coef  <- matrix(NA_real_, nprop, nrow(pairs))
  nbreaks  <- rep(NA_integer_, nprop)
  byears   <- vector("list", nprop)

  cat(sprintf("[draws] propagating %d posterior draws (of %d) ...\n", nprop, D))
  for (di in seq_along(sel)) {
    d  <- sel[di]
    ba <- draws$beta_draws[d, , , ]
    sa <- draws$sigma_draws[d, , ]
    rd <- tryCatch(suppressWarnings(suppressMessages({
      g   <- compute_gfevd_array(ba, sa, k, p, hg, eps)
      tci <- compute_tci(g)$TotalConnectedness
      net <- compute_network_tv(g, sectors, dates,
                                threshold_type = thr_type, threshold_value = thr_val)
      rez <- compute_resilience_tv(beta_array = ba, sigma_array = sa, gfevd_array = g,
                                   sectors = sectors, dates = dates, k = k, p = p,
                                   horizon = hr, resist_h = resist_h, eps = resist_eps)
      list(tci = tci, cent = net$centrality_tv, rez = rez)
    })), error = function(e) NULL)
    if (is.null(rd)) next
    tci_mat[di, ] <- rd$tci
    h1_slope[di]  <- .hac_trend(rd$tci)["slope"]
    panel <- dplyr::inner_join(rd$cent, rd$rez, by = c("Sector", "TimeIndex"))
    for (pp in seq_len(nrow(pairs)))
      h3_coef[di, pp] <- .panel_coef(panel, pairs$Resilience[pp], pairs$Centrality[pp])
    bp <- tryCatch(detect_breaks(rd$tci, dates, cfg$BREAK_MAX, cfg$BREAK_MIN_SIZE),
                   error = function(e) NULL)
    if (!is.null(bp)) {
      nbreaks[di] <- bp$n_breaks
      byears[[di]] <- if (!is.null(bp$break_dates))
        as.integer(format(as.Date(bp$break_dates), "%Y")) else integer(0)
    }
    if (di %% 50L == 0L) cat(sprintf("[draws]   %d/%d done\n", di, nprop))
  }

  q3 <- function(x) stats::quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE)

  tci_bands <- tibble::tibble(
    TimeIndex = seq_len(n_eff), Date = dates,
    TCI_mean = colMeans(tci_mat, na.rm = TRUE),
    TCI_lwr  = apply(tci_mat, 2, function(c) stats::quantile(c, 0.025, na.rm = TRUE)),
    TCI_med  = apply(tci_mat, 2, function(c) stats::quantile(c, 0.5,   na.rm = TRUE)),
    TCI_upr  = apply(tci_mat, 2, function(c) stats::quantile(c, 0.975, na.rm = TRUE)))
  s_ci <- q3(h1_slope)
  h1 <- tibble::tibble(
    Hypothesis = "H1: Densification (TCI trend)",
    Slope_mean = mean(h1_slope, na.rm = TRUE),
    Slope_lwr = s_ci[1], Slope_med = s_ci[2], Slope_upr = s_ci[3],
    P_slope_gt0 = mean(h1_slope > 0, na.rm = TRUE),
    Credible95_positive = isTRUE(s_ci[1] > 0),
    Credible95_negative = isTRUE(s_ci[3] < 0),
    Verdict = if (isTRUE(s_ci[1] > 0)) "Supported (95% CI > 0)"
              else if (isTRUE(s_ci[3] < 0)) "Rejected (95% CI < 0)"
              else "Not supported (CI spans 0)",
    N_draws = sum(is.finite(h1_slope)))

  h3_rows <- lapply(seq_len(nrow(pairs)), function(pp) {
    cc <- h3_coef[, pp]; ci <- q3(cc)
    pgt <- mean(cc > 0, na.rm = TRUE)
    tibble::tibble(
      Centrality = pairs$Centrality[pp], Resilience = pairs$Resilience[pp],
      Coef_mean = mean(cc, na.rm = TRUE),
      Coef_lwr = ci[1], Coef_med = ci[2], Coef_upr = ci[3],
      P_gt0 = pgt, P_direction = max(pgt, 1 - pgt, na.rm = TRUE),
      Credible95 = isTRUE(ci[1] > 0) || isTRUE(ci[3] < 0),
      N_draws = sum(is.finite(cc)))
  })
  h3 <- dplyr::bind_rows(h3_rows)
  h3$Caveat <- "draws-mode: centrality & resilience re-estimated per draw; generated-regressor uncertainty propagated."

  prereg <- as.integer(cfg$PREREGISTERED_EPISODES)
  tol    <- as.integer(cfg$BREAK_ALIGN_TOL %||% 1L)
  epi_lab <- unname(cfg$HISTORICAL_EPISODES[as.character(prereg)])
  finite_by <- byears[!vapply(byears, is.null, logical(1))]
  h4_episode <- tibble::tibble(
    Year = prereg, Episode = epi_lab,
    Post_prob_break = vapply(prereg, function(yr)
      mean(vapply(finite_by, function(b) any(abs(b - yr) <= tol), logical(1))), numeric(1)),
    Tol_years = tol)
  max_b <- max(c(0L, nbreaks), na.rm = TRUE)
  h4_count <- tibble::tibble(
    N_breaks = 0:max_b,
    Posterior_freq = vapply(0:max_b, function(b) mean(nbreaks == b, na.rm = TRUE), numeric(1)))
  all_years <- sort(unique(unlist(finite_by)))
  year_incl <- if (length(all_years)) tibble::tibble(
    Year = all_years,
    Post_prob = vapply(all_years, function(yr)
      mean(vapply(finite_by, function(b) any(b == yr), logical(1))), numeric(1))) else
    tibble::tibble(Year = integer(0), Post_prob = numeric(0))

  list(n_prop = nprop, pairs = pairs,
       tci_bands = tci_bands, tci_draws = tci_mat,
       h1 = h1, h1_slope_draws = h1_slope,
       h3 = h3, h3_coef_draws = h3_coef,
       h4_episode = h4_episode, h4_count = h4_count, h4_year_inclusion = year_incl,
       nbreaks_draws = nbreaks)
}
