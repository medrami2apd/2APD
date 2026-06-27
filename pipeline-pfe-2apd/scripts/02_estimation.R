#!/usr/bin/env Rscript
source("scripts/00_setup/setup.R")
stopifnot(exists("start_logging"), exists("end_logging"), exists("%||%"))
start_logging("02_estimation")

if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' is required (mvrnorm).")

cat("\n========================================\n")
cat("TVP-BVAR-SV ESTIMATION (Gibbs + FFBS + KSC-SV)\n")
cat("========================================\n\n")

N_CHAINS   <- CONFIG$N_CHAINS
N_ITER     <- CONFIG$N_ITER
N_BURN     <- CONFIG$N_BURN
THIN       <- CONFIG$THIN
N_TRAIN    <- CONFIG$N_TRAIN
SEED       <- CONFIG$SEED
KAPPA_Q    <- CONFIG$KAPPA_Q
P0_SCALE   <- CONFIG$P0_SCALE
MU_PRIOR_V <- CONFIG$SV_MU_PRIOR_V
PHI_A      <- CONFIG$SV_PHI_A; PHI_B <- CONFIG$SV_PHI_B
S2H_C0     <- CONFIG$SV_S2H_C0; S2H_D0 <- CONFIG$SV_S2H_D0
NU_Q       <- N_TRAIN
set.seed(SEED)

N_CHAINS <- as.integer(N_CHAINS); N_ITER <- as.integer(N_ITER)
N_BURN   <- as.integer(N_BURN);   THIN   <- as.integer(THIN)
if (is.na(N_CHAINS) || N_CHAINS < 1L)
  stop("CONFIG$N_CHAINS must be an integer >= 1 (got ", CONFIG$N_CHAINS, ").")
if (is.na(N_ITER) || N_ITER < 2L)
  stop("CONFIG$N_ITER must be an integer >= 2 (got ", CONFIG$N_ITER, ").")
if (is.na(THIN) || THIN < 1L)
  stop("CONFIG$THIN must be an integer >= 1 (use THIN = 1L for no thinning; got ", CONFIG$THIN, ").")
if (is.na(N_BURN) || N_BURN < 0L)
  stop("CONFIG$N_BURN must be an integer >= 0 (use N_BURN = 0L for no burn-in; got ", CONFIG$N_BURN, ").")
if (N_BURN >= N_ITER)
  stop(sprintf(paste0("CONFIG$N_BURN (%d) must be strictly less than CONFIG$N_ITER (%d), ",
                      "otherwise zero draws are kept. For no burn-in set N_BURN = 0L; ",
                      "a typical choice is N_BURN = N_ITER/2."), N_BURN, N_ITER))

cat("--- Step 1: Loading data ---\n")
if (!file.exists(CONFIG$PATHS$processed)) stop("Processed data not found: ", CONFIG$PATHS$processed)
prep    <- readRDS(CONFIG$PATHS$processed)
Yfull   <- prep$Y
dates   <- prep$dates
sectors <- prep$sectors %||% CONFIG$SECTORS

k    <- CONFIG$K
p    <- as.integer(CONFIG$P_LAGS)
nX   <- CONFIG$NX
Tlen <- nrow(Yfull)
n    <- Tlen - p
if (n < (N_TRAIN + 5)) stop("Too few observations after lagging for N_TRAIN.")

Y <- t(Yfull)
Xmat <- matrix(0.0, nrow = n, ncol = nX)
for (s in seq_len(n)) {
  t_idx <- s + p
  vec <- numeric(0)
  for (lag in seq_len(p)) vec <- c(vec, Y[, t_idx - lag])
  Xmat[s, ] <- vec
}
Ydep <- Y[, (p + 1):Tlen, drop = FALSE]
cat(sprintf("Dimensions: k=%d, T=%d, p=%d, nX=%d, n_eff=%d\n", k, Tlen, p, nX, n))

tr  <- seq_len(N_TRAIN)
Xtr <- Xmat[tr, , drop = FALSE]
XtX_inv <- solve(crossprod(Xtr) + diag(1e-6, nX))
b0_mat   <- matrix(0.0, k, nX)
Vols_mat <- matrix(0.0, k, nX)
h_init   <- numeric(k)
for (i in seq_len(k)) {
  yi <- Ydep[i, tr]
  bi <- XtX_inv %*% crossprod(Xtr, yi)
  b0_mat[i, ] <- as.numeric(bi)
  ri <- yi - as.numeric(Xtr %*% bi)
  s2 <- sum(ri^2) / max(1, (N_TRAIN - nX))
  Vols_mat[i, ] <- pmax(s2 * diag(XtX_inv), 1e-8)
  h_init[i] <- log(max(s2, 1e-6))
}

rinvgamma <- function(shape, rate) 1 / rgamma(1, shape = shape, rate = rate)

KSC_p  <- c(0.00730, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750)
KSC_m  <- c(-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -1.08819) - 1.2704
KSC_v  <- c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261)
KSC_sd <- sqrt(KSC_v)

ffbs_coef <- function(y, Xe, qdiag, rvar, b0, P0) {
  a_filt <- matrix(0.0, n, nX)
  P_filt <- vector("list", n)
  a <- b0; P <- P0; Q <- diag(qdiag, nX)
  for (s in seq_len(n)) {
    P_pred <- P + Q
    xs <- Xe[s, ]
    Px <- P_pred %*% xs
    f  <- as.numeric(crossprod(xs, Px)) + rvar[s]
    K  <- Px / f
    v  <- y[s] - sum(xs * a)
    a  <- a + as.numeric(K) * v
    P  <- P_pred - tcrossprod(as.numeric(K), as.numeric(Px))
    P  <- (P + t(P)) / 2
    a_filt[s, ] <- a; P_filt[[s]] <- P
  }
  B <- matrix(0.0, n, nX)
  B[n, ] <- MASS::mvrnorm(1, a_filt[n, ], P_filt[[n]] + diag(1e-10, nX))
  for (s in (n - 1):1) {
    Pp  <- P_filt[[s]] + Q
    J   <- P_filt[[s]] %*% solve(Pp)
    mks <- a_filt[s, ] + as.numeric(J %*% (B[s + 1, ] - a_filt[s, ]))
    Vks <- P_filt[[s]] - J %*% Pp %*% t(J)
    Vks <- (Vks + t(Vks)) / 2
    B[s, ] <- MASS::mvrnorm(1, mks, Vks + diag(1e-10, nX))
  }
  B
}

sv_update <- function(eps, hprev, mu, phi, s2h) {
  ystar <- log(eps^2 + 1e-7)
  logw <- matrix(0.0, n, 7)
  for (j in 1:7)
    logw[, j] <- log(KSC_p[j]) - log(KSC_sd[j]) -
                 0.5 * ((ystar - hprev - KSC_m[j])^2) / KSC_v[j]
  logw <- logw - apply(logw, 1, max)
  w    <- exp(logw); w <- w / rowSums(w)
  cw   <- t(apply(w, 1, cumsum))
  u    <- runif(n)
  idx  <- max.col(-(cw < u), ties.method = "first")
  mvec <- KSC_m[idx]; Rvec <- KSC_v[idx]
  ytil <- ystar - mvec
  cc <- mu * (1 - phi)
  a_f <- numeric(n); P_f <- numeric(n)
  denom_phi <- max(1 - phi^2, 1e-8)
  a_pred <- mu; P_pred <- s2h / denom_phi
  for (s in seq_len(n)) {
    if (s > 1) { a_pred <- cc + phi * a_f[s - 1]; P_pred <- phi^2 * P_f[s - 1] + s2h }
    f <- P_pred + Rvec[s]; Kf <- P_pred / f
    a_f[s] <- a_pred + Kf * (ytil[s] - a_pred)
    P_f[s] <- (1 - Kf) * P_pred
  }
  h <- numeric(n)
  h[n] <- rnorm(1, a_f[n], sqrt(max(P_f[n], 1e-12)))
  for (s in (n - 1):1) {
    Pp <- phi^2 * P_f[s] + s2h
    J  <- phi * P_f[s] / Pp
    mh <- a_f[s] + J * (h[s + 1] - (cc + phi * a_f[s]))
    Vh <- P_f[s] - J^2 * Pp
    h[s] <- rnorm(1, mh, sqrt(max(Vh, 1e-12)))
  }
  prec <- 1 / MU_PRIOR_V + (1 - phi^2) / s2h + (n - 1) * (1 - phi)^2 / s2h
  num  <- (1 - phi^2) / s2h * h[1] +
          (1 - phi) / s2h * sum(h[2:n] - phi * h[1:(n - 1)])
  mu   <- rnorm(1, num / prec, sqrt(1 / prec))
  lp_phi <- function(ph) {
    if (abs(ph) >= 0.999) return(-Inf)
    pr <- (PHI_A - 1) * log((ph + 1) / 2) + (PHI_B - 1) * log((1 - ph) / 2)
    ll <- 0.5 * log(1 - ph^2) - 0.5 * (1 - ph^2) * (h[1] - mu)^2 / s2h
    res <- (h[2:n] - mu) - ph * (h[1:(n - 1)] - mu)
    ll <- ll - 0.5 * sum(res^2) / s2h
    pr + ll
  }
  phi_prop <- phi + rnorm(1, 0, 0.05)
  if (log(runif(1)) < (lp_phi(phi_prop) - lp_phi(phi))) phi <- phi_prop
  ss  <- (1 - phi^2) * (h[1] - mu)^2 +
         sum(((h[2:n] - mu) - phi * (h[1:(n - 1)] - mu))^2)
  s2h <- rinvgamma(S2H_C0 + n / 2, S2H_D0 + 0.5 * ss)
  list(h = h, mu = mu, phi = phi, s2h = s2h)
}

keep_idx       <- seq(N_BURN + 1, N_ITER, by = THIN)
keep_per_chain <- length(keep_idx)
total_keep     <- keep_per_chain * N_CHAINS
cat(sprintf("\n--- Step 2: Gibbs sampling | %d chains x %d iter (burn %d, thin %d) => %d draws ---\n",
            N_CHAINS, N_ITER, N_BURN, THIN, total_keep))

mon_names <- c(paste0("mu_", seq_len(k)), paste0("phi_", seq_len(k)), paste0("s2h_", seq_len(k)))

run_chain <- function(chain) {
  set.seed(SEED + chain)
  prog_file <- file.path(prog_dir, sprintf("chain_%02d_progress.log", chain))
  cat(sprintf("[chain %d] %d iters (burn %d, thin %d) started %s\n",
              chain, N_ITER, N_BURN, THIN, format(Sys.time())),
      file = prog_file, append = FALSE)
  beta_c  <- array(0.0, dim = c(keep_per_chain, n, k, nX))
  sigma_c <- array(0.0, dim = c(keep_per_chain, n, k))
  sv_c    <- array(0.0, dim = c(keep_per_chain, k, 3))
  mon_c   <- matrix(NA_real_, keep_per_chain, length(mon_names))
  B_cur   <- vector("list", k)
  for (i in seq_len(k)) B_cur[[i]] <- matrix(rep(b0_mat[i, ], each = n), n, nX)
  q_cur   <- matrix(0.0, k, nX); for (i in seq_len(k)) q_cur[i, ] <- (KAPPA_Q^2) * Vols_mat[i, ]
  h_cur   <- matrix(rep(h_init, each = n), n, k)
  mu_cur  <- h_init; phi_cur <- rep(0.95, k); s2h_cur <- rep(0.05, k)
  P0      <- diag(P0_SCALE, nX)
  ck <- 0L; t_chain <- Sys.time()
  for (it in seq_len(N_ITER)) {
    for (i in seq_len(k)) {
      rvar <- exp(h_cur[, i])
      Bi   <- ffbs_coef(Ydep[i, ], Xmat, q_cur[i, ], rvar, b0_mat[i, ], P0)
      B_cur[[i]] <- Bi
      dB <- rbind(Bi[1, ] - b0_mat[i, ], diff(Bi))
      ss_q <- colSums(dB^2)
      for (j in seq_len(nX)) {
        q_cur[i, j] <- rinvgamma(NU_Q / 2 + n / 2,
                                 (NU_Q * (KAPPA_Q^2) * Vols_mat[i, j]) / 2 + ss_q[j] / 2)
      }
      eps  <- Ydep[i, ] - rowSums(Xmat * Bi)
      svr  <- sv_update(eps, h_cur[, i], mu_cur[i], phi_cur[i], s2h_cur[i])
      h_cur[, i] <- svr$h; mu_cur[i] <- svr$mu; phi_cur[i] <- svr$phi; s2h_cur[i] <- svr$s2h
    }
    if (it > N_BURN && ((it - N_BURN) %% THIN == 0)) {
      ck <- ck + 1L
      for (i in seq_len(k)) {
        beta_c[ck, , i, ] <- B_cur[[i]]
        sigma_c[ck, , i]  <- exp(0.5 * h_cur[, i])
      }
      sv_c[ck, , 1] <- mu_cur; sv_c[ck, , 2] <- phi_cur; sv_c[ck, , 3] <- s2h_cur
      mon_c[ck, ]   <- c(mu_cur, phi_cur, s2h_cur)
    }
    if (it == 1L || it %% 100 == 0 || it == N_ITER) {
      phase <- if (it <= N_BURN) "burn-in " else "sampling"
      msg <- sprintf("  chain %d | iter %5d/%d (%s) | kept %d | elapsed %.1fs\n",
                     chain, it, N_ITER, phase, ck,
                     as.numeric(difftime(Sys.time(), t_chain, units = "secs")))
      cat(msg); utils::flush.console()
      cat(msg, file = prog_file, append = TRUE)
    }
  }
  list(beta = beta_c, sigma = sigma_c, sv = sv_c, mon = mon_c)
}

n_cores  <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 1L)
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L
mc_cores <- max(1L, min(N_CHAINS, n_cores))
can_fork <- (.Platform$OS.type != "windows")

prog_dir <- normalizePath(file.path(getwd(), "outputs", "diagnostics"),
                          winslash = "/", mustWork = FALSE)
dir.create(prog_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Live progress -> %s/chain_<id>_progress.log  (open it or `tail -f` to watch)\n",
            prog_dir))

t_start <- Sys.time()
if (mc_cores > 1L && can_fork) {
  cat(sprintf("Running %d chains on %d cores [fork/mclapply].\n", N_CHAINS, mc_cores))
  chain_results <- parallel::mclapply(seq_len(N_CHAINS), run_chain,
                                      mc.cores = mc_cores, mc.preschedule = FALSE,
                                      mc.set.seed = FALSE)
} else if (mc_cores > 1L) {
  cat(sprintf("Running %d chains on %d cores [PSOCK cluster].\n", N_CHAINS, mc_cores))
  cl <- parallel::makeCluster(mc_cores, outfile = "")
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  parallel::clusterExport(cl,
    varlist = c("run_chain", "ffbs_coef", "sv_update", "rinvgamma",
                "b0_mat", "Vols_mat", "h_init", "P0_SCALE", "KAPPA_Q", "NU_Q",
                "N_ITER", "N_BURN", "THIN", "SEED", "keep_per_chain", "prog_dir",
                "k", "n", "nX", "Ydep", "Xmat", "mon_names",
                "KSC_p", "KSC_m", "KSC_v", "KSC_sd",
                "MU_PRIOR_V", "PHI_A", "PHI_B", "S2H_C0", "S2H_D0"),
    envir = environment())
  parallel::clusterEvalQ(cl, requireNamespace("MASS", quietly = TRUE))
  chain_results <- parallel::parLapply(cl, seq_len(N_CHAINS), run_chain)
  parallel::stopCluster(cl)
} else {
  cat(sprintf("Running %d chain(s) sequentially (single core).\n", N_CHAINS))
  chain_results <- lapply(seq_len(N_CHAINS), run_chain)
}
bad <- vapply(chain_results,
              function(x) inherits(x, "try-error") || !is.list(x) || is.null(x$beta),
              logical(1))
if (any(bad)) stop("Gibbs chain(s) failed: ", paste(which(bad), collapse = ", "),
                   " | ", as.character(chain_results[[which(bad)[1]]]))

beta_draws  <- array(0.0, dim = c(total_keep, n, k, nX))
sigma_draws <- array(0.0, dim = c(total_keep, n, k))
sv_par      <- array(0.0, dim = c(total_keep, k, 3),
                     dimnames = list(NULL, NULL, c("mu", "phi", "s2h")))
mon_chain   <- array(NA_real_, dim = c(keep_per_chain, N_CHAINS, length(mon_names)),
                     dimnames = list(NULL, NULL, mon_names))
for (chain in seq_len(N_CHAINS)) {
  cr   <- chain_results[[chain]]
  rows <- ((chain - 1L) * keep_per_chain + 1L):(chain * keep_per_chain)
  beta_draws[rows, , , ] <- cr$beta
  sigma_draws[rows, , ]  <- cr$sigma
  sv_par[rows, , ]       <- cr$sv
  mon_chain[, chain, ]   <- cr$mon
}
cat(sprintf("Sampling complete in %.1f seconds.\n",
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))

cat("\n--- Step 3: Convergence diagnostics ---\n")
rhat_vec <- ess_vec <- setNames(rep(NA_real_, length(mon_names)), mon_names)
if (requireNamespace("posterior", quietly = TRUE) && N_CHAINS >= 2) {
  for (vn in seq_along(mon_names)) {
    draws_mat <- mon_chain[, , vn]
    rv <- tryCatch(posterior::rhat(draws_mat), error = function(e) NA_real_)
    ev <- tryCatch(posterior::ess_bulk(draws_mat), error = function(e) NA_real_)
    rhat_vec[vn] <- rv; ess_vec[vn] <- ev
  }
  max_rhat <- suppressWarnings(max(rhat_vec, na.rm = TRUE))
  min_ess  <- suppressWarnings(min(ess_vec, na.rm = TRUE))
  cat(sprintf("Max R-hat (SV params): %.3f | Min bulk ESS: %.1f\n", max_rhat, min_ess))
  if (is.finite(max_rhat) && max_rhat > 1.01)
    warning("Rhat > 1.01 on SV hyperparameters: inspect before trusting results.")
} else {
  max_rhat <- NA_real_; min_ess <- NA_real_
  cat("posterior package or >=2 chains unavailable; skipped R-hat/ESS.\n")
}
utils::write.csv(
  data.frame(parameter = mon_names, rhat = as.numeric(rhat_vec), ess_bulk = as.numeric(ess_vec)),
  "outputs/diagnostics/convergence_summary.csv", row.names = FALSE)

cat("\n--- Step 4: Assembling posterior-mean arrays ---\n")
beta_array  <- apply(beta_draws, c(2, 3, 4), mean)
sigma_array <- apply(sigma_draws, c(2, 3), mean)
eff_dates   <- dates[(p + 1):Tlen]
n_eff       <- n
if (anyNA(beta_array) || anyNA(sigma_array)) stop("NA in assembled arrays.")
stopifnot(dim(beta_array)[1] == n_eff, dim(beta_array)[2] == k,
          dim(beta_array)[3] == nX, length(eff_dates) == n_eff)
cat(sprintf("beta_array : [%d, %d, %d]  sigma_array: [%d, %d]\n", n_eff, k, nX, n_eff, k))

cat("\n--- Step 5: Saving outputs ---\n")
saveRDS(list(beta_array = beta_array, sigma_array = sigma_array,
             dates = eff_dates, sectors = sectors,
             k = k, p = p, nX = nX, n_eff = n_eff,
             center = prep$center, scale = prep$scale),
        CONFIG$PATHS$posterior)
cat(sprintf("Saved posterior-mean arrays -> %s\n", CONFIG$PATHS$posterior))

if (!is.null(CONFIG$PATHS$posterior_draws)) {
  saveRDS(list(beta_draws = beta_draws,
               sigma_draws = sigma_draws,
               dates = eff_dates, sectors = sectors,
               k = k, p = p, nX = nX, n_eff = n_eff),
          CONFIG$PATHS$posterior_draws)
  cat(sprintf("Saved thinned posterior draws -> %s\n", CONFIG$PATHS$posterior_draws))
}

fit_path <- CONFIG$PATHS$fit %||% CONFIG$PATHS$stan_fit
if (!is.null(fit_path)) {
  saveRDS(list(
    method      = "gibbs_ffbs_tvpbvar_sv",
    beta_array  = beta_array, sigma_array = sigma_array,
    sv_params   = sv_par,
    diagnostics = list(rhat = rhat_vec, ess_bulk = ess_vec,
                       max_rhat = max_rhat, min_ess = min_ess,
                       n_draws = total_keep, n_chains = N_CHAINS,
                       thin = THIN, burn = N_BURN),
    priors = list(KAPPA_Q = KAPPA_Q, NU_Q = NU_Q, PHI_A = PHI_A, PHI_B = PHI_B,
                  S2H_C0 = S2H_C0, S2H_D0 = S2H_D0),
    dates = eff_dates, sectors = sectors, k = k, p = p, nX = nX, n_eff = n_eff),
    fit_path)
  cat(sprintf("Saved fit object -> %s\n", fit_path))
}

cat("\n========================================\n")
cat("Estimation complete (Gibbs/FFBS TVP-BVAR-SV).\n")
cat("Posterior-mean arrays feed connectedness/network/resilience/breaks as-is.\n")
cat("Use posterior_draws.rds to build credible bands (loop GFEVD over draws).\n")
cat("========================================\n")

end_logging()
