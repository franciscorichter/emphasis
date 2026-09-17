# 04b-init-table.R — write results/<tier>/init.csv from the init-arm job files.
#
#   Rscript R/04b-init-table.R [--tier main]
#
# One row per init job (tree x I1/I2/I3): the deficit of the starting point
# handed to MCEM (`start`), the deficit MCEM returns (`after`), whether the
# tree's auto_bounds box contains the exact MLE (`box`), the seconds the
# initialiser itself took (`sec`), and the job outcome (`out`).  A job that
# did not finish has NA in the numeric columns and its outcome in `out`.
# 05-compare.R pairs two of these on tree x cell x cfg.

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, default = NULL) {
  i <- which(args == k); if (length(i)) args[i[1] + 1L] else default
}
TIER <- getarg("--tier", "main")
HERE <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])))
RES  <- file.path(HERE, "..", "results", TIER)
JOBS <- file.path(RES, "jobs")
stopifnot(dir.exists(JOBS))

fs <- list.files(JOBS, pattern = "--I[123]\\.rds$", full.names = TRUE)
g1 <- function(r, k) { v <- r[[k]]; if (is.null(v) || !length(v)) NA else v[1] }
rows <- lapply(fs, function(f) {
  x <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(x)) return(NULL)
  data.frame(tree  = g1(x, "tree_id"),
             cell  = g1(x, "cell"),
             cfg   = g1(x, "config"),
             out   = g1(x, "outcome"),
             n     = g1(x, "n"),
             start = as.numeric(g1(x, "start_delta_ell")),
             after = as.numeric(g1(x, "delta_ell")),
             box   = as.logical(g1(x, "box_contains_mle0")),
             sec   = as.numeric(g1(x, "init_seconds")),
             stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
tab <- tab[order(tab$tree, tab$cfg), ]
out <- file.path(RES, "init.csv")
utils::write.csv(tab, out, row.names = FALSE)
cat(sprintf("[04b-init-table] %d init rows -> %s (ok %d, timeout %d, error %d)\n",
            nrow(tab), out, sum(tab$out == "ok"), sum(tab$out == "timeout"), sum(tab$out == "error")))
