vma_theta <- function(beta_t, k, p, horizon) {
  Theta <- array(0, dim = c(horizon, k, k))
  Theta[1, , ] <- diag(k)
  if (horizon >= 2) for (h in 2:horizon) {
    for (lag in 1:p) if (h - lag >= 1) {
      cs <- (lag - 1) * k + 1; ce <- lag * k
      Phi <- matrix(beta_t[, cs:ce], nrow = k, ncol = k)
      Theta[h, , ] <- Theta[h, , ] + Phi %*% Theta[h - lag, , ]
    }
  }
  Theta
}

gfevd_at <- function(beta_t, sigma_t, k, p, horizon, eps = 1e-8) {
  Theta <- vma_theta(beta_t, k, p, horizon)
  P <- diag(sigma_t, nrow = k)
  irf_sq <- matrix(0, k, k)
  for (h in 1:horizon) {
    imp <- Theta[h, , ] %*% P
    irf_sq <- irf_sq + imp^2
  }
  g <- matrix(0, k, k)
  for (i in 1:k) {
    rs <- sum(irf_sq[i, ])
    if (is.finite(rs) && rs > eps) g[i, ] <- irf_sq[i, ] / rs
    else { g[i, ] <- 0; g[i, i] <- 1 }
  }
  g
}

compute_gfevd_array <- function(beta_array, sigma_array, k, p, horizon, eps = 1e-8) {
  n_eff <- dim(beta_array)[1]
  g <- array(0, dim = c(k, k, n_eff))
  for (t in 1:n_eff) g[, , t] <- gfevd_at(beta_array[t, , ], sigma_array[t, ], k, p, horizon, eps)
  message(sprintf("GFEVD: k=%d p=%d horizon=%d periods=%d", k, p, horizon, n_eff))
  g
}

compute_tci <- function(gfevd_array) {
  k <- dim(gfevd_array)[1]; n_eff <- dim(gfevd_array)[3]
  tci <- vapply(1:n_eff, function(t) {
    m <- gfevd_array[, , t]; (sum(m) - sum(diag(m))) / k * 100
  }, numeric(1))
  tibble::tibble(TimeIndex = 1:n_eff, TotalConnectedness = tci)
}

compute_directional <- function(gfevd_array, sectors) {
  k <- dim(gfevd_array)[1]; n_eff <- dim(gfevd_array)[3]
  out <- vector("list", n_eff)
  for (t in 1:n_eff) {
    m <- gfevd_array[, , t]; diag(m) <- 0
    out[[t]] <- tibble::tibble(
      TimeIndex = t, Sector = sectors,
      FromOthers = rowSums(m) * 100,
      ToOthers   = colSums(m) * 100,
      NetSpillover = colSums(m) * 100 - rowSums(m) * 100)
  }
  dplyr::bind_rows(out)
}

compute_spillover_matrix <- function(gfevd_array, sectors, window = NULL) {
  n_eff <- dim(gfevd_array)[3]
  if (is.null(window)) window <- 1:n_eff
  sm <- apply(gfevd_array[, , window, drop = FALSE], c(1, 2), mean) * 100
  diag(sm) <- 0
  dimnames(sm) <- list(sectors, sectors)
  sm
}

companion_spectral_radius <- function(beta_t, k, p) {
  comp_dim <- k * p
  comp <- matrix(0, comp_dim, comp_dim)
  comp[1:k, 1:comp_dim] <- beta_t
  if (p > 1) comp[(k + 1):comp_dim, 1:(comp_dim - k)] <- diag(comp_dim - k)
  max(Mod(eigen(comp, only.values = TRUE)$values))
}
