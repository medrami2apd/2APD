.nw_vcov <- function(fit, lag = NULL) {
  if (is.null(lag)) sandwich::NeweyWest(fit, prewhite = FALSE, adjust = TRUE)
  else sandwich::NeweyWest(fit, lag = lag, prewhite = FALSE, adjust = TRUE)
}

.cbb_index <- function(n, block_len) {
  if (is.null(block_len)) block_len <- max(1L, floor(n^(1/3)))
  starts <- sample.int(n, size = ceiling(n / block_len), replace = TRUE)
  idx <- unlist(lapply(starts, function(s) ((s - 1 + 0:(block_len - 1)) %% n) + 1L))
  idx[1:n]
}

test_h1 <- function(tci_series, Y = NULL, p = 1L, horizon = 10L,
                    B = 1000L, block_len = NULL, eps = 1e-8) {
  n <- length(tci_series); t_idx <- seq_len(n)
  trend_fit <- stats::lm(tci_series ~ t_idx)
  V <- .nw_vcov(trend_fit)
  ct <- lmtest::coeftest(trend_fit, vcov. = V)
  slope <- ct["t_idx", 1]; slope_p <- ct["t_idx", 4]
  tau <- suppressWarnings(stats::cor(t_idx, tci_series, method = "kendall"))

  mean_fit <- stats::lm(tci_series ~ 1)
  mp <- lmtest::coeftest(mean_fit, vcov. = .nw_vcov(mean_fit))["(Intercept)", 4]

  boot_p <- NA_real_
  if (!is.null(Y) && requireNamespace("vars", quietly = TRUE)) {
    k <- ncol(Y); obs_mean <- mean(tci_series)
    null_means <- rep(NA_real_, B)
    for (b in 1:B) {
      Yb <- vapply(1:k, function(j) Y[.cbb_index(nrow(Y), block_len), j], numeric(nrow(Y)))
      vb <- tryCatch(vars::VAR(as.data.frame(Yb), p = p, type = "const"), error = function(e) NULL)
      if (is.null(vb)) next
      Bcoef <- t(sapply(vb$varresult, function(r) stats::coef(r)[1:(k*p)]))
      Sig <- diag(sqrt(diag(stats::cov(stats::residuals(vb)))), k)
      gg <- gfevd_at(matrix(Bcoef, k, k*p), diag(Sig), k, p, horizon, eps)
      null_means[b] <- (sum(gg) - sum(diag(gg))) / k * 100
    }
    boot_p <- mean(null_means >= obs_mean, na.rm = TRUE)
  }

  tbl <- tibble::tibble(
    Hypothesis = "H1: Densification",
    Metric = c("Kendall tau", "HAC trend slope", "HAC trend p-value",
               "HAC mean p-value (vs 0)", "Bootstrap independence p-value"),
    Value = c(tau, slope, slope_p, mp, boot_p),
    Result = c(if (!is.na(tau) && tau > 0) "Positive" else "Non-positive",
               sprintf("%.4f", slope),
               if (!is.na(slope_p) && slope_p < 0.05) "Significant trend" else "No sig. trend",
               if (!is.na(mp) && mp < 0.05) "> independence" else "n.s.",
               if (!is.na(boot_p) && boot_p < 0.05) "> independence (boot)" else "n.s./skipped"))

  sig_trend <- is.finite(slope_p) && slope_p < 0.05
  indep     <- (is.finite(mp) && mp < 0.05) || (is.finite(boot_p) && boot_p < 0.05)
  verdict <- if (sig_trend && indep) "Supported"
             else if (sig_trend || indep) "Partial"
             else "Not supported"

  list(table = tbl, verdict = verdict,
       slope = slope, slope_p = slope_p, mean_p = mp,
       boot_p = boot_p, kendall_tau = tau,
       sig_trend = sig_trend, exceeds_independence = indep)
}

test_h2 <- function(centrality_tv, neutral_band = 0, B = 1000L) {
  df <- centrality_tv %>%
    dplyr::mutate(State = dplyr::case_when(
      NetDegree >  neutral_band ~ "Transmitter",
      NetDegree < -neutral_band ~ "Absorber",
      TRUE ~ "Neutral")) %>%
    dplyr::arrange(Sector, TimeIndex)

  stay_by_sector <- df %>% dplyr::group_by(Sector) %>%
    dplyr::summarise(p_stay = {
      s <- State; if (length(s) < 2) NA_real_ else mean(s[-1] == s[-length(s)])
    }, n = dplyr::n(), .groups = "drop")
  obs_stay <- mean(stay_by_sector$p_stay, na.rm = TRUE)

  sectors <- unique(df$Sector)
  null_stay <- rep(NA_real_, B)
  for (b in 1:B) {
    ps <- vapply(sectors, function(sc) {
      s <- df$State[df$Sector == sc]
      if (length(s) < 2) return(NA_real_)
      sp <- sample(s); mean(sp[-1] == sp[-length(sp)])
    }, numeric(1))
    null_stay[b] <- mean(ps, na.rm = TRUE)
  }
  perm_p <- mean(null_stay >= obs_stay, na.rm = TRUE)

  tbl <- tibble::tibble(
    Hypothesis = "H2: Role Persistence",
    Metric = c("Mean P(stay)", "Permutation null mean P(stay)", "Permutation p-value"),
    Value = c(obs_stay, mean(null_stay, na.rm = TRUE), perm_p),
    Result = c(if (obs_stay > 0.5) "> 50%" else "<= 50%", "-",
               if (!is.na(perm_p) && perm_p < 0.05) "Persistent (non-random)" else "Random"))

  persistent <- is.finite(perm_p) && perm_p < 0.05
  verdict <- if (persistent) "Supported" else "Not supported"

  list(table = tbl, verdict = verdict,
       mean_p_stay = obs_stay, perm_p = perm_p, persistent = persistent)
}

test_h3 <- function(centrality_tv, resilience_tv,
                    centrality_vars = c("OutDegree","InDegree","Betweenness","Closeness"),
                    resilience_vars = c("Resistance","RecoverySpeed","Reconfiguration"),
                    cluster = "group", padj = "BH") {
  panel <- dplyr::inner_join(centrality_tv, resilience_tv,
                             by = c("Sector","TimeIndex")) %>% as.data.frame()
  zscore <- function(x) { s <- stats::sd(x, na.rm = TRUE); if (is.na(s) || s == 0) x*0 else (x - mean(x, na.rm=TRUE))/s }
  for (v in unique(c(centrality_vars, resilience_vars)))
    if (v %in% names(panel)) panel[[v]] <- zscore(panel[[v]])

  if (nrow(panel) < 30 || length(unique(panel$TimeIndex)) < 3)
    return(tibble::tibble(Hypothesis="H3", Note="Insufficient panel"))
  pdat <- plm::pdata.frame(panel, index = c("Sector","TimeIndex"))

  G <- length(unique(panel$Sector))
  small_g  <- isTRUE(get0("CONFIG", ifnotfound = list())$H3_SMALL_G_CORRECTION %||% TRUE)
  wcb_reps <- as.integer(get0("CONFIG", ifnotfound = list())$H3_WCB_REPS %||% 999L)
  has_wcb  <- small_g && requireNamespace("fixest", quietly = TRUE) &&
              requireNamespace("fwildclusterboot", quietly = TRUE)

  rows <- list()
  for (dv in resilience_vars) for (iv in centrality_vars) {
    if (!all(c(dv, iv) %in% names(panel))) next
    f <- stats::as.formula(sprintf("%s ~ %s", dv, iv))
    mod <- tryCatch(plm::plm(f, data = pdat, model = "within", effect = "twoways"),
                    error = function(e) NULL)
    if (is.null(mod)) next
    V <- tryCatch(plm::vcovHC(mod, type = "HC1", cluster = cluster), error = function(e) NULL)
    ct <- if (is.null(V)) lmtest::coeftest(mod) else lmtest::coeftest(mod, vcov. = V)
    coef_i <- ct[iv, 1]; se_i <- ct[iv, 2]; t_i <- coef_i / se_i
    p_cr1 <- 2 * stats::pt(-abs(t_i), df = max(G - 1L, 1L))
    p_wcb <- NA_real_
    if (has_wcb) p_wcb <- tryCatch({
      fe <- fixest::feols(stats::as.formula(sprintf("%s ~ %s | Sector + TimeIndex", dv, iv)),
                          data = panel)
      bt <- fwildclusterboot::boottest(fe, param = iv, clustid = "Sector",
                                       B = wcb_reps, type = "rademacher")
      fwildclusterboot::pval(bt)
    }, error = function(e) NA_real_)
    p_use <- if (is.finite(p_wcb)) p_wcb else if (small_g) p_cr1 else ct[iv, 4]
    rows[[length(rows)+1L]] <- tibble::tibble(
      Resilience = dv, Centrality = iv, Coef = coef_i, SE = se_i,
      p_HC1 = ct[iv, 4], p_cr1_Gm1 = p_cr1, p_wcb = p_wcb, p_raw = p_use)
  }
  res <- dplyr::bind_rows(rows)
  if (nrow(res) == 0) return(tibble::tibble(Hypothesis="H3", Note="No estimable models"))
  res$p_adj <- stats::p.adjust(res$p_raw, method = padj)
  res$Significant <- res$p_adj < 0.05
  res$Hypothesis <- "H3: Centrality-Resilience"
  res$N_clusters <- G
  res$Inference <- if (has_wcb) sprintf("wild-cluster bootstrap (G=%d)", G)
                   else sprintf("cluster-robust HC1 + t(G-1), G=%d", G)
  res$Caveat <- "single-path point estimate; generated-regressor uncertainty propagated in draws mode (scripts/11_draws_inference.R)"
  res
}

h3_sector_cross_section <- function(centrality_avg, resilience_avg,
                                    centrality_vars = c("OutDegree","InDegree","Betweenness","Closeness"),
                                    resilience_vars = c("Resistance","RecoverySpeed","Reconfiguration")) {
  prof <- dplyr::inner_join(centrality_avg, resilience_avg, by = "Sector") %>%
    as.data.frame()
  cent_keep <- intersect(unique(c(centrality_vars, "NetDegree")), names(prof))
  res_keep  <- intersect(resilience_vars, names(prof))
  prof <- prof[, c("Sector", cent_keep, res_keep), drop = FALSE]

  rank_desc <- function(x) rank(-x, ties.method = "min", na.last = "keep")
  for (v in c(cent_keep, res_keep))
    prof[[paste0("rank_", v)]] <- rank_desc(prof[[v]])

  if (centrality_vars[1] %in% names(prof))
    prof <- prof[order(-prof[[centrality_vars[1]]]), , drop = FALSE]

  rows <- list()
  for (iv in intersect(centrality_vars, names(prof)))
    for (dv in res_keep) {
      x <- prof[[iv]]; y <- prof[[dv]]
      ok <- stats::complete.cases(x, y); n <- sum(ok)
      sp <- if (n >= 4) suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman")) else NA_real_
      pe <- if (n >= 4) suppressWarnings(stats::cor(x[ok], y[ok], method = "pearson")) else NA_real_
      pv <- if (n >= 4) tryCatch(suppressWarnings(
        stats::cor.test(x[ok], y[ok], method = "spearman")$p.value),
        error = function(e) NA_real_) else NA_real_
      rows[[length(rows) + 1L]] <- tibble::tibble(
        Centrality = iv, Resilience = dv, n = n,
        Spearman = sp, Pearson = pe, Spearman_p = pv)
    }
  corr <- dplyr::bind_rows(rows)
  if (nrow(corr) > 0) {
    corr$Spearman_p_adj <- stats::p.adjust(corr$Spearman_p, method = "BH")
    corr$Significant <- !is.na(corr$Spearman_p_adj) & corr$Spearman_p_adj < 0.05
  }
  list(profile = tibble::as_tibble(prof), correlation = corr)
}
