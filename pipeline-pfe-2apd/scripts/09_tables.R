#!/usr/bin/env Rscript
source("scripts/00_setup/setup.R")
start_logging("09_tables")

stopifnot(exists("export_table"))

hyp_path <- file.path(CONFIG$PATHS$hypotheses, "hypothesis_results.rds")
hyp <- readRDS(hyp_path)
brk <- tryCatch(readRDS(CONFIG$PATHS$breaks), error = function(e) NULL)
dri <- tryCatch(readRDS(CONFIG$PATHS$draws_inference), error = function(e) NULL)
draws_mode <- !is.null(dri)

h1_obj <- hyp$h1; h2_obj <- hyp$h2
h1_is_list <- is.list(h1_obj) && !is.data.frame(h1_obj) && !is.null(h1_obj$verdict)
h2_is_list <- is.list(h2_obj) && !is.data.frame(h2_obj) && !is.null(h2_obj$verdict)
h1_tbl <- if (h1_is_list) h1_obj$table else h1_obj
h2_tbl <- if (h2_is_list) h2_obj$table else h2_obj

h1_verdict <- {
  if (draws_mode && !is.null(dri$h1) && "Verdict" %in% names(dri$h1)) {
    dri$h1$Verdict[1]
  } else if (h1_is_list) {
    h1_obj$verdict
  } else {
    sig_trend <- any(grepl("Significant trend", h1_tbl$Result))
    indep     <- any(grepl("> independence", h1_tbl$Result))
    if (sig_trend && indep) "Supported"
    else if (sig_trend || indep) "Partial"
    else "Not supported"
  }
}

h2_verdict <- {
  if (h2_is_list) {
    h2_obj$verdict
  } else {
    pp <- h2_tbl$Value[h2_tbl$Metric == "Permutation p-value"]
    if (length(pp) && !is.na(pp) && pp < 0.05) "Supported" else "Not supported"
  }
}

h3_verdict <- {
  if (draws_mode && !is.null(dri$h3) && "Credible95" %in% names(dri$h3)) {
    nsig <- sum(dri$h3$Credible95, na.rm = TRUE); ntot <- nrow(dri$h3)
    if (nsig > 0)
      sprintf("Partial (%d/%d pairs: 95%% credible interval excludes 0; generated-regressor uncertainty propagated)", nsig, ntot)
    else "Not supported (all 95% credible intervals span 0)"
  } else if (is.data.frame(hyp$h3) && "Significant" %in% names(hyp$h3) && any(hyp$h3$Significant, na.rm = TRUE))
    "Partial (posterior-mean only; generated-regressor uncertainty NOT propagated - run stage 11)"
  else "Not supported"
}

h4_verdict <- {
  if (is.null(brk)) {
    "Not evaluated (breaks output missing - run 07_breaks.R)"
  } else {
    n_break_vec <- vapply(brk, function(x) as.integer(x$n_breaks %||% 0L), integer(1))
    total_breaks <- sum(n_break_vec, na.rm = TRUE)
    aligned <- 0L
    for (nm in names(brk)) {
      a <- brk[[nm]]$alignment
      if (!is.null(a) && nrow(a) > 0 && "Episode" %in% names(a))
        aligned <- aligned + sum(a$Episode != "(no match)", na.rm = TRUE)
    }
    if (total_breaks == 0L) "Not supported"
    else if (aligned > 0L) "Supported (interpretive: breaks align to historical episodes)"
    else "Partial (breaks detected; no episode alignment within tolerance)"
  }
}

h4_detail <- if (!is.null(brk)) {
  paste(vapply(names(brk), function(nm)
    sprintf("%s: %d breaks", brk[[nm]]$label %||% nm,
            as.integer(brk[[nm]]$n_breaks %||% 0L)),
    character(1)), collapse = "; ")
} else NA_character_

verdict <- tibble::tibble(
  Hypothesis = c("H1: Densification", "H2: Role Persistence",
                 "H3: Centrality-Resilience", "H4: Structural Reconfiguration"),
  Evidence   = c(h1_verdict, h2_verdict, h3_verdict, h4_verdict),
  Inference  = c(if (draws_mode) "draws-mode (stage 11)" else "posterior-mean, structured verdict (stage 06)",
                 "within-sector permutation, structured verdict (stage 06)",
                 if (draws_mode) "draws-mode credible intervals (stage 11)" else "posterior-mean (stage 06)",
                 if (draws_mode) "draws-mode break posterior (stage 11) + point-estimate (stage 07)" else "point-estimate (stage 07)"),
  Detail     = c(NA_character_, NA_character_, NA_character_, h4_detail)
)
export_table(verdict, "table_28_hypothesis_verdicts")
print(verdict)
cat(sprintf("\nHypothesis verdicts compiled (H1-H4): %d rows.\n", nrow(verdict)))

end_logging()
