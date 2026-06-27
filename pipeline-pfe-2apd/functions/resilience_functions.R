girf_at <- function(beta_t, sigma_t, k, p, horizon) {
  comp_dim <- k * p
  comp <- matrix(0, comp_dim, comp_dim)
  comp[1:k, 1:comp_dim] <- beta_t
  if (p > 1) comp[(k + 1):comp_dim, 1:(comp_dim - k)] <- diag(comp_dim - k)
  P <- diag(sigma_t, nrow = k)
  girf <- array(0, dim = c(k, k, horizon))
  for (h in 0:(horizon - 1)) {
    psi <- if (h == 0) diag(k) else (comp %^% h)[1:k, 1:k]
    girf[, , h + 1] <- psi %*% P
  }
  girf
}

post_peak_half_life <- function(series, horizon) {
  s <- abs(series)
  pk <- which.max(s); peak <- s[pk]
  if (peak <= 1e-8) return(0)
  post <- s[pk:length(s)]
  hl <- which(post <= 0.5 * peak)[1]
  if (is.na(hl)) horizon - pk + 1 else hl - 1
}

long_run_cumulative_multiplier <- function(series) {
  sum(series)
}
bn_permanent_component <- long_run_cumulative_multiplier

compute_resilience_tv <- function(beta_array, sigma_array, gfevd_array = NULL,
                                  sectors, dates, k, p, horizon,
                                  resist_h = 4L, eps = 1e-2) {
  n_eff <- dim(beta_array)[1]
  if (n_eff < 1) stop("Need >= 1 effective period for resilience.")
  win <- 1:min(horizon, as.integer(resist_h) + 1L)
  rows <- vector("list", 0L)
  for (t in 1:n_eff) {
    girf_t <- girf_at(beta_array[t, , ], sigma_array[t, ], k, p, horizon)
    for (i in 1:k) {
      resistance <- 1 / (max(abs(girf_t[i, , win, drop = FALSE])) + eps)
      hl <- vapply(1:k, function(j) post_peak_half_life(girf_t[i, j, ], horizon), numeric(1))
      recovery <- mean(hl)
      perm <- vapply(1:k, function(j) long_run_cumulative_multiplier(girf_t[i, j, ]), numeric(1))
      reconfig <- sum(abs(perm))
      rows[[length(rows) + 1L]] <- data.frame(
        Sector = sectors[i], TimeIndex = t, Date = dates[t],
        Resistance = resistance, RecoverySpeed = recovery, Reconfiguration = reconfig)
    }
  }
  panel <- dplyr::bind_rows(rows)
  message(sprintf("Resilience panel: %d obs (%d sectors x %d periods)",
                  nrow(panel), k, n_eff))
  panel
}
