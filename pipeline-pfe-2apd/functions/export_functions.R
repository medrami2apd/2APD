.LOG_ENV <- new.env(parent = emptyenv())

.ensure_dir <- function(path, recursive = TRUE) {
  if (!is.null(path) && nzchar(path) && !dir.exists(path))
    dir.create(path, recursive = recursive, showWarnings = FALSE)
  invisible(path)
}
.tables_dir <- function() {
  d <- tryCatch(CONFIG$PATHS$tables, error = function(e) NULL)
  if (is.null(d) || !nzchar(d)) d <- "outputs/tables"
  .ensure_dir(d); d
}
.logs_dir <- function() {
  d <- tryCatch(CONFIG$PATHS$logs, error = function(e) NULL)
  if (is.null(d) || !nzchar(d)) d <- "outputs/logs"
  .ensure_dir(d); d
}

start_logging <- function(stage) {
  .LOG_ENV$stage <- stage
  .LOG_ENV$t0 <- Sys.time()
  path <- file.path(.logs_dir(), paste0(stage, ".log"))
  .LOG_ENV$con <- tryCatch(file(path, open = "wt"), error = function(e) NULL)
  msg <- sprintf("[%s] START %s", format(.LOG_ENV$t0), stage)
  cat(msg, "\n"); if (!is.null(.LOG_ENV$con)) writeLines(msg, .LOG_ENV$con)
  invisible(TRUE)
}

end_logging <- function() {
  if (is.null(.LOG_ENV$stage)) return(invisible(FALSE))
  dt <- round(as.numeric(difftime(Sys.time(), .LOG_ENV$t0, units = "secs")), 1)
  msg <- sprintf("[%s] END %s (%.1fs)", format(Sys.time()), .LOG_ENV$stage, dt)
  cat(msg, "\n")
  if (!is.null(.LOG_ENV$con)) { writeLines(msg, .LOG_ENV$con); close(.LOG_ENV$con); .LOG_ENV$con <- NULL }
  invisible(TRUE)
}

export_table <- function(df, name) {
  tdir <- .tables_dir()
  csv <- file.path(tdir, paste0(name, ".csv"))
  utils::write.csv(df, csv, row.names = FALSE)
  tryCatch({
    ldir <- file.path(tdir, "latex"); .ensure_dir(ldir)
    tex <- file.path(ldir, paste0(name, ".tex"))
    print(xtable::xtable(df), file = tex, include.rownames = FALSE)
  }, error = function(e) NULL)
  invisible(csv)
}

export_matrix <- function(mat, name) {
  tdir <- .tables_dir()
  df <- as.data.frame(mat); df <- cbind(RowName = rownames(mat), df)
  utils::write.csv(df, file.path(tdir, paste0(name, ".csv")), row.names = FALSE)
}
