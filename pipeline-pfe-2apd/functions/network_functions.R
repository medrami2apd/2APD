build_transmission_adjacency <- function(gfevd_slice, sectors) {
  adj <- t(gfevd_slice)
  diag(adj) <- 0
  dimnames(adj) <- list(sectors, sectors)
  adj
}

.safe_graph <- function(adj) {
  igraph::graph_from_adjacency_matrix(adj, mode = "directed", weighted = TRUE)
}

apply_network_threshold <- function(adj, type = c("none","proportional","absolute"),
                                    value = 0) {
  type <- match.arg(as.character(type), c("none","proportional","absolute"))
  a <- adj; diag(a) <- 0
  if (type == "none" || is.null(value) || value <= 0) return(a)
  pos <- a[a > 0]
  if (length(pos) == 0) return(a)
  if (type == "proportional") {
    thr <- stats::quantile(pos, probs = max(0, 1 - value), na.rm = TRUE)
    a[a < thr] <- 0
  } else {
    a[a < value] <- 0
  }
  a
}

compute_all_centralities <- function(adj, sectors, eps = 1e-8,
                                     threshold_type = "none", threshold_value = 0) {
  out_deg <- rowSums(adj)
  in_deg  <- colSums(adj)
  adj_g <- apply_network_threshold(adj, threshold_type, threshold_value)
  g <- tryCatch(.safe_graph(adj_g), error = function(e) NULL)

  if (!is.null(g) && igraph::ecount(g) > 0) {
    w <- igraph::E(g)$weight
    dist_w <- 1 / (w + eps)
    betw <- tryCatch(igraph::betweenness(g, directed = TRUE, weights = dist_w),
                     error = function(e) rep(NA_real_, length(sectors)))
    clos <- tryCatch(igraph::closeness(g, mode = "out", weights = dist_w, normalized = TRUE),
                     error = function(e) rep(NA_real_, length(sectors)))
    eig  <- tryCatch(igraph::eigen_centrality(g, directed = TRUE, weights = w)$vector,
                     error = function(e) rep(NA_real_, length(sectors)))
    pr   <- tryCatch(igraph::page_rank(g, directed = TRUE, weights = w)$vector,
                     error = function(e) rep(NA_real_, length(sectors)))
  } else {
    betw <- clos <- eig <- pr <- rep(NA_real_, length(sectors))
  }

  tibble::tibble(
    Sector      = sectors,
    OutDegree   = as.numeric(out_deg),
    InDegree    = as.numeric(in_deg),
    NetDegree   = as.numeric(out_deg - in_deg),
    Betweenness = as.numeric(betw),
    Closeness   = as.numeric(clos),
    Eigenvector = as.numeric(eig),
    PageRank    = as.numeric(pr))
}

compute_density <- function(adj) {
  n <- nrow(adj)
  (sum(adj) - sum(diag(adj))) / (n * (n - 1))
}

compute_network_tv <- function(gfevd_array, sectors, dates,
                               threshold_type = "none", threshold_value = 0) {
  n_eff <- dim(gfevd_array)[3]
  cent <- vector("list", n_eff); dens <- numeric(n_eff)
  for (t in 1:n_eff) {
    adj <- build_transmission_adjacency(gfevd_array[, , t], sectors)
    cent[[t]] <- dplyr::mutate(compute_all_centralities(adj, sectors,
                                 threshold_type = threshold_type,
                                 threshold_value = threshold_value),
                               TimeIndex = t, Date = dates[t])
    dens[t] <- compute_density(adj)
  }
  list(centrality_tv = dplyr::bind_rows(cent),
       density_tv = tibble::tibble(TimeIndex = 1:n_eff, Date = dates, Density = dens))
}
