build_xy <- function(prep, p) {
  Y <- t(prep$Y)
  k <- nrow(Y); Tlen <- ncol(Y); n <- Tlen - p; nX <- k * p
  Xmat <- matrix(0.0, n, nX)
  for (s in seq_len(n)) {
    t_idx <- s + p; vec <- numeric(0)
    for (lag in seq_len(p)) vec <- c(vec, Y[, t_idx - lag])
    Xmat[s, ] <- vec
  }
  Ydep <- Y[, (p + 1):Tlen, drop = FALSE]
  list(Xmat = Xmat, Ydep = Ydep, n = n, k = k, nX = nX, Tlen = Tlen)
}

ols_var <- function(Xmat, Ydep) {
  k <- nrow(Ydep); nX <- ncol(Xmat); n <- ncol(Ydep)
  XtXi <- solve(crossprod(Xmat) + diag(1e-8, nX))
  B <- matrix(0, k, nX); sigma <- numeric(k); resid <- matrix(0, k, n)
  for (i in 1:k) {
    yi <- Ydep[i, ]
    bi <- as.numeric(XtXi %*% crossprod(Xmat, yi))
    B[i, ] <- bi
    ri <- yi - as.numeric(Xmat %*% bi)
    resid[i, ] <- ri
    sigma[i] <- sqrt(mean(ri^2))
  }
  list(B = B, sigma = sigma, resid = resid)
}

gaussian_crps <- function(y, mu, sdv) {
  sdv <- pmax(sdv, 1e-8)
  z <- (y - mu) / sdv
  sdv * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi))
}

predict_metrics <- function(pred_mean, pred_sd, Ydep) {
  pred_sd <- pmax(pred_sd, 1e-8)
  err <- Ydep - pred_mean
  list(
    rmse = sqrt(mean(err^2)),
    mae  = mean(abs(err)),
    lpds = mean(dnorm(Ydep, pred_mean, pred_sd, log = TRUE)),
    crps = mean(gaussian_crps(Ydep, pred_mean, pred_sd)),
    pit  = as.numeric(pnorm((Ydep - pred_mean) / pred_sd))
  )
}

.spec0 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 8 || stats::sd(x) == 0) return(stats::var(x))
  fit <- tryCatch(stats::ar(x, aic = TRUE, order.max = min(10L, length(x) - 1L)),
                  error = function(e) NULL)
  if (is.null(fit)) return(stats::var(x))
  denom <- (1 - sum(fit$ar))^2
  if (!is.finite(denom) || denom <= 0) return(stats::var(x))
  fit$var.pred / denom
}

geweke_z <- function(x, f1 = 0.1, f2 = 0.5) {
  x <- as.numeric(x); n <- length(x)
  n1 <- floor(f1 * n); n2 <- floor(f2 * n)
  if (n1 < 4 || n2 < 4) return(NA_real_)
  a <- x[1:n1]; b <- x[(n - n2 + 1):n]
  va <- .spec0(a) / n1; vb <- .spec0(b) / n2
  if (!is.finite(va + vb) || (va + vb) <= 0) return(NA_real_)
  (mean(a) - mean(b)) / sqrt(va + vb)
}

mcmc_diag <- function(fit, sectors, cfg = list()) {
  sv <- fit$sv_params
  if (is.null(sv)) return(tibble::tibble(Note = "No SV draws stored in fit"))
  nd <- dim(sv)[1]; k <- dim(sv)[2]
  nchains <- as.integer(fit$diagnostics$n_chains %||% 1L)
  if (nchains < 1L) nchains <- 1L
  kpc <- nd %/% nchains
  if (kpc < 4L) return(tibble::tibble(Note = "Too few draws per chain for diagnostics"))
  rhat_thr <- cfg$RHAT_THRESHOLD %||% 1.01
  ess_thr  <- cfg$ESS_THRESHOLD  %||% 400
  g1 <- cfg$GEWEKE_FRAC1 %||% 0.1; g2 <- cfg$GEWEKE_FRAC2 %||% 0.5
  par_labels <- c("mu", "phi", "sigma2_h")
  rows <- list()
  has_coda <- requireNamespace("coda", quietly = TRUE)
  for (pidx in 1:3) for (i in 1:k) {
    v <- sv[, i, pidx]
    mat <- matrix(v[1:(kpc * nchains)], nrow = kpc, ncol = nchains)
    rh <- tryCatch(posterior::rhat(mat),     error = function(e) NA_real_)
    eb <- tryCatch(posterior::ess_bulk(mat), error = function(e) NA_real_)
    et <- tryCatch(posterior::ess_tail(mat), error = function(e) NA_real_)
    gz <- suppressWarnings(max(abs(apply(mat, 2, function(cc) geweke_z(cc, g1, g2))),
                               na.rm = TRUE))
    if (!is.finite(gz)) gz <- NA_real_
    rl <- NA_real_
    if (has_coda) rl <- tryCatch({
      d <- coda::raftery.diag(as.numeric(mat), q = 0.025, r = 0.005, s = 0.95)
      if (is.matrix(d$resmatrix)) max(d$resmatrix[, "N"], na.rm = TRUE) else NA_real_
    }, error = function(e) NA_real_)
    rows[[length(rows) + 1]] <- tibble::tibble(
      Parameter = sprintf("%s[%s]", par_labels[pidx], sectors[i]),
      Rhat = rh, ESS_bulk = eb, ESS_tail = et, Geweke_absZ = gz, RafteryLewis_N = rl)
  }
  out <- dplyr::bind_rows(rows)
  out$Rhat_pass   <- is.finite(out$Rhat) & out$Rhat < rhat_thr
  out$ESS_pass    <- is.finite(out$ESS_bulk) & out$ESS_bulk >= ess_thr
  out$Geweke_pass <- is.finite(out$Geweke_absZ) & out$Geweke_absZ < 1.96
  out
}

tv_sv_relevance <- function(post, sectors) {
  k <- post$k
  rows <- vector("list", k)
  for (i in 1:k) {
    bi <- post$beta_array[, i, ]
    drift <- mean(apply(bi, 2, stats::sd), na.rm = TRUE)
    level <- mean(abs(colMeans(bi)), na.rm = TRUE)
    si <- post$sigma_array[, i]
    rows[[i]] <- tibble::tibble(
      Sector    = sectors[i],
      TV_drift  = round(drift, 4),
      TV_ratio  = round(drift / (drift + level + 1e-12), 4),
      SV_cv     = round(stats::sd(si) / (mean(si) + 1e-12), 4),
      SV_maxmin = round(max(si) / (min(si) + 1e-12), 3))
  }
  dplyr::bind_rows(rows)
}

model_comparison <- function(prep, post, cfg = list()) {
  xy <- build_xy(prep, post$p); Xmat <- xy$Xmat; Ydep <- xy$Ydep; n <- xy$n; k <- post$k
  ols <- ols_var(Xmat, Ydep)
  tvp_mean <- matrix(0, k, n); const_mean <- matrix(0, k, n)
  for (i in 1:k) {
    tvp_mean[i, ]   <- rowSums(Xmat * post$beta_array[, i, ])
    const_mean[i, ] <- as.numeric(Xmat %*% ols$B[i, ])
  }
  sv_sd <- t(post$sigma_array)
  tvp_resid_sd <- sqrt(rowMeans((Ydep - tvp_mean)^2))
  mk_sd <- function(vec) matrix(rep(vec, n), k, n)
  models <- list(
    "TVP-BVAR-SV (full)"  = list(m = tvp_mean,   s = sv_sd),
    "TVP-BVAR (no SV)"    = list(m = tvp_mean,   s = mk_sd(tvp_resid_sd)),
    "Constant BVAR-SV"    = list(m = const_mean, s = sv_sd),
    "Constant BVAR (std)" = list(m = const_mean, s = mk_sd(ols$sigma)))
  rows <- list(); pit_full <- NULL
  for (nm in names(models)) {
    mm <- predict_metrics(models[[nm]]$m, models[[nm]]$s, Ydep)
    rows[[nm]] <- tibble::tibble(Model = nm,
      RMSE = round(mm$rmse, 4), MAE = round(mm$mae, 4),
      LPDS = round(mm$lpds, 4), CRPS = round(mm$crps, 4))
    if (nm == "TVP-BVAR-SV (full)") pit_full <- mm$pit
  }
  list(table = dplyr::bind_rows(rows), pit = pit_full)
}

waic_dic <- function(prep, post, draws, cfg = list()) {
  if (is.null(draws) || is.null(draws$beta_draws)) return(NULL)
  xy <- build_xy(prep, post$p); Xmat <- xy$Xmat; Ydep <- xy$Ydep; n <- xy$n; k <- post$k
  D <- dim(draws$beta_draws)[1]
  maxd <- min(D, as.integer(cfg$WAIC_MAX_DRAWS %||% 200L))
  sel <- unique(round(seq(1, D, length.out = maxd)))
  maxd <- length(sel)
  n_obs <- k * n
  ll <- matrix(0.0, maxd, n_obs)
  for (di in seq_along(sel)) {
    d <- sel[di]
    for (i in 1:k) {
      mu  <- rowSums(Xmat * draws$beta_draws[d, , i, ])
      sdv <- pmax(draws$sigma_draws[d, , i], 1e-8)
      ll[di, ((i - 1) * n + 1):(i * n)] <- dnorm(Ydep[i, ], mu, sdv, log = TRUE)
    }
  }
  lppd   <- sum(apply(ll, 2, function(c) matrixStats::logSumExp(c) - log(maxd)))
  p_waic <- sum(apply(ll, 2, stats::var))
  waic   <- -2 * (lppd - p_waic)
  Dbar <- -2 * mean(rowSums(ll))
  llhat <- numeric(n_obs)
  for (i in 1:k) {
    mu  <- rowSums(Xmat * post$beta_array[, i, ])
    sdv <- pmax(post$sigma_array[, i], 1e-8)
    llhat[((i - 1) * n + 1):(i * n)] <- dnorm(Ydep[i, ], mu, sdv, log = TRUE)
  }
  Dhat <- -2 * sum(llhat); pD <- Dbar - Dhat; dic <- Dbar + pD
  looic <- NA_real_
  if (requireNamespace("loo", quietly = TRUE))
    looic <- tryCatch(loo::loo(ll)$estimates["looic", "Estimate"],
                      error = function(e) NA_real_)
  list(waic = round(waic, 2), p_waic = round(p_waic, 2), lppd = round(lppd, 2),
       dic = round(dic, 2), pD = round(pD, 2),
       looic = if (is.na(looic)) NA_real_ else round(looic, 2), n_draws = maxd)
}

rolling_vs_tvp_tci <- function(prep, conn, post, cfg = list()) {
  p <- post$p; k <- post$k; nX <- post$nX %||% (k * p)
  xy <- build_xy(prep, p); Xmat <- xy$Xmat; Ydep <- xy$Ydep; n <- xy$n
  W <- as.integer(cfg$ROLLING_WINDOW %||% 40L)
  horizon <- CONFIG$HORIZON_GFEVD
  tvp_tci <- conn$total_conn$TotalConnectedness
  roll <- rep(NA_real_, n)
  if (W < (nX + 2L) || W > n) {
    warning(sprintf("Rolling window W=%d incompatible with n=%d, nX=%d; skipping.", W, n, nX))
  } else {
    for (s in seq.int(W, n)) {
      w <- (s - W + 1):s
      Xw <- Xmat[w, , drop = FALSE]
      XtXi <- tryCatch(solve(crossprod(Xw) + diag(1e-6, nX)), error = function(e) NULL)
      if (is.null(XtXi)) next
      B <- matrix(0, k, nX); sig <- numeric(k)
      for (i in 1:k) {
        yi <- Ydep[i, w]
        bi <- as.numeric(XtXi %*% crossprod(Xw, yi))
        B[i, ] <- bi
        sig[i] <- stats::sd(yi - as.numeric(Xw %*% bi))
      }
      g <- tryCatch(gfevd_at(B, sig, k, p, horizon, CONFIG$EPS), error = function(e) NULL)
      if (is.null(g)) next
      roll[s] <- (sum(g) - sum(diag(g))) / k * 100
    }
  }
  tibble::tibble(TimeIndex = 1:n, Date = as.Date(post$dates),
                 TVP_TCI = tvp_tci, Rolling_TCI = roll)
}

robustness_checklist <- function() {
  tibble::tibble(
    Category = c(
      rep("Evaluation - point forecast", 2),
      rep("Evaluation - density forecast", 3),
      rep("Evaluation - model comparison", 5),
      rep("Evaluation - TVP/SV relevance", 3),
      rep("Robustness - estimation algorithm", 3),
      rep("Robustness - prior sensitivity", 1),
      rep("Robustness - outliers & breaks", 3),
      rep("Robustness - identification", 2),
      rep("Robustness - computational", 1),
      rep("Robustness - empirical", 3)),
    Item = c(
      "RMSE / MAE", "MSE / QLIKE (volatility)",
      "Log predictive density (LPDS)", "CRPS", "PIT calibration",
      "WAIC", "DIC", "LOO-CV", "Marginal likelihood / Bayes factor", "Savage-Dickey ratio",
      "Posterior evidence for time variation", "Forecast gains vs constant BVAR", "Volatility variability (SV)",
      "Convergence: Rhat / ESS / Geweke", "Raftery-Lewis", "Init / burn-in sensitivity",
      "Q dof/scale, SV & Minnesota priors",
      "Stochastic volatility (heteroskedastic shocks)", "Heavy-tailed (Student-t) errors", "Structural-break dating",
      "Variable ordering (Cholesky)", "Sign restrictions / external instruments",
      "FFBS / numerical stability",
      "Rolling-window sub-sample stability", "Comparison with restricted models", "Alternative data transformations"),
    Status = c(
      "Implemented", "Optional",
      "Implemented", "Implemented", "Implemented",
      "Implemented", "Implemented", "Implemented (if loo pkg)", "Not applicable", "Not applicable",
      "Implemented", "Implemented", "Implemented",
      "Implemented", "Implemented (if coda pkg)", "Built-in",
      "Documented",
      "Built-in", "Not applicable", "Implemented",
      "Not applicable", "Not applicable",
      "Built-in",
      "Implemented", "Implemented", "Optional"),
    Note = c(
      "table_33 (one-step panel).", "SV already models conditional variance; QLIKE optional.",
      "table_33.", "table_33 (closed-form Gaussian CRPS).", "figure_eval_01 (uniform => calibrated).",
      "table_34 (from saved per-draw paths).", "table_34.", "table_34 when 'loo' is installed.",
      "No competing point-null of primary interest; descriptive aim.", "Same reason as Bayes factor.",
      "table_32 (time-variation ratio).", "table_33 RMSE/LPDS gain columns.", "table_32 (SV cv & max/min).",
      "table_31 across chains.", "table_31 when 'coda' is installed.", "Multiple over-dispersed chains (N_CHAINS).",
      "Priors reported in fit; full re-estimation sweep is a stated upgrade.",
      "Core model component (KSC mixture SV).", "Gaussian errors; Student-t is a stated extension.", "H4 Bai-Perron stage (07_breaks.R).",
      "Generalized FEVD/GIRF is order-invariant -> N.A.", "Reduced-form descriptive connectedness -> N.A.",
      "Symmetrised covariances + jitter in FFBS (02_estimation.R).",
      "table_35 + figure_eval_02 (rolling vs TVP).", "table_33 (4 nested benchmarks).", "Growth rates used; levels run is optional.")
  )
}

mcmc_diag_derived <- function(draws, cfg, sectors, k, p, nchains) {
  if (is.null(draws$beta_draws) || is.null(draws$sigma_draws))
    return(tibble::tibble(Quantity = NA_character_, Note = "No posterior draws stored"))
  D <- dim(draws$beta_draws)[1]; n_eff <- dim(draws$beta_draws)[2]; nX <- dim(draws$beta_draws)[4]
  nchains <- max(1L, as.integer(nchains))
  kpc <- D %/% nchains
  if (kpc < 2L) return(tibble::tibble(Quantity = NA_character_,
                 Note = sprintf("Too few draws per chain (kpc=%d) for split-Rhat", kpc)))
  rb <- cfg$ROBUSTNESS %||% list()
  rhat_thr <- rb$RHAT_THRESHOLD %||% 1.01
  ess_thr  <- rb$ESS_THRESHOLD %||% 400

  diag_of <- function(mat, label) {
    tibble::tibble(
      Quantity = label,
      Rhat     = tryCatch(posterior::rhat(mat),     error = function(e) NA_real_),
      ESS_bulk = tryCatch(posterior::ess_bulk(mat), error = function(e) NA_real_),
      ESS_tail = tryCatch(posterior::ess_tail(mat), error = function(e) NA_real_))
  }
  scalar_to_chains <- function(v) matrix(v[seq_len(kpc * nchains)], nrow = kpc, ncol = nchains)

  rows <- list()
  ntp <- as.integer(rb$DERIVED_BETA_TIMEPOINTS %||% 3L)
  tps <- unique(round(seq(1, n_eff, length.out = max(2L, ntp))))
  for (i in seq_len(k)) for (t0 in tps) {
    jj <- min(i, nX)
    rows[[length(rows) + 1L]] <- diag_of(scalar_to_chains(draws$beta_draws[, t0, i, jj]),
                                         sprintf("beta[%s, t=%d, own-lag]", sectors[i], t0))
  }
  per_chain <- min(as.integer(rb$DERIVED_TCI_DRAWS_PER_CHAIN %||% 40L), kpc)
  idx_within <- unique(round(seq(1, kpc, length.out = per_chain)))
  per_chain <- length(idx_within)
  tci_scalar <- rep(NA_real_, per_chain * nchains); pos <- 1L
  for (ch in seq_len(nchains)) for (w in idx_within) {
    draw_id <- (ch - 1L) * kpc + w
    g <- tryCatch(suppressWarnings(suppressMessages(
      compute_gfevd_array(draws$beta_draws[draw_id, , , ], draws$sigma_draws[draw_id, , ],
                          k, p, cfg$HORIZON_GFEVD, cfg$EPS))), error = function(e) NULL)
    tci_scalar[pos] <- if (is.null(g)) NA_real_ else mean(compute_tci(g)$TotalConnectedness)
    pos <- pos + 1L
  }
  rows[[length(rows) + 1L]] <- diag_of(matrix(tci_scalar, nrow = per_chain, ncol = nchains),
                                       "mean TCI_t")

  out <- dplyr::bind_rows(rows)
  out$Rhat_pass <- is.finite(out$Rhat) & out$Rhat < rhat_thr
  out$ESS_pass  <- is.finite(out$ESS_bulk) & out$ESS_bulk >= ess_thr
  out$N_chains  <- nchains
  out
}
