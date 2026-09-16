## H74 remediation harness — containment, tightness and cost of the box
## returned by auto_bounds(), measured on the validation study's own CR trees.
##
##   Rscript H74_fix.R --lib <rlib> --out <csv> [--per-cell 5] [--variant tag]
##                     [--cores 8] [--val <dev/validation>] [--cells a,b,c]
##
## For every tree it records the box, whether the exact MLE (cond = 0, the Nee
## closed form) and the survival-conditioned MLE (DDD::bd_ML, cond = 1, soc = 2)
## lie inside it, the per-axis widths and the area, and the wall clock and
## forward-simulation count auto_bounds spent.  Reference MLEs come from the
## study's reference file when present and are recomputed otherwise.
##
## The draw counter rebinds simulate_div_tree_cpp inside the package namespace,
## so it counts every forward draw auto_bounds makes, retries included.  Run
## with --cores 1 for a clean single-thread cost measurement.

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, d = NULL) { i <- which(args == k); if (length(i)) args[i[1] + 1L] else d }

LIB      <- getarg("--lib", NA_character_)
OUT      <- getarg("--out", "/tmp/h74/out/boxes.csv")
PER_CELL <- as.integer(getarg("--per-cell", "5"))
VARIANT  <- getarg("--variant", "base")
CORES    <- as.integer(getarg("--cores", "1"))
VAL      <- getarg("--val", "/Users/pancho/Code/emphasis/dev/validation")
CELLS    <- getarg("--cells", "")
GAM      <- identical(getarg("--gam", "no"), "yes")

if (!is.na(LIB)) .libPaths(c(LIB, .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape); library(DDD)
                                 library(parallel) })
source(file.path(VAL, "R", "00-common.R"))

td <- readRDS(file.path(VAL, "data", "trees-main.rds"))
rf <- tryCatch(readRDS(file.path(VAL, "data", "reference-main.rds"))$refs,
               error = function(e) NULL)

trees <- Filter(function(x) x$kind == "cr", td$trees)
if (nzchar(CELLS)) {
  keep <- strsplit(CELLS, ",")[[1]]
  trees <- Filter(function(x) x$cell %in% keep, trees)
}
by_cell <- split(trees, vapply(trees, function(x) x$cell, ""))
trees <- unlist(lapply(by_cell, function(g) {
  g <- g[order(vapply(g, function(x) x$index, 1L))]
  g[seq_len(min(PER_CELL, length(g)))]
}), recursive = FALSE)

## --- forward-draw counter --------------------------------------------------
.draws <- new.env(parent = emptyenv()); .draws$n <- 0L
local({
  ns <- asNamespace("emphasis")
  nm <- "simulate_div_tree_cpp"
  orig <- get(nm, envir = ns)
  wrapped <- function(...) { .draws$n <- .draws$n + 1L; orig(...) }
  unlockBinding(nm, ns); assign(nm, wrapped, envir = ns); lockBinding(nm, ns)
})

ref_for <- function(tt) {
  r <- if (!is.null(rf)) rf[[tt$tree_id]] else NULL
  if (!is.null(r) && !is.null(r$mle) && !is.null(r$mle_cond1))
    return(list(mle0 = as.numeric(r$mle), mle1 = as.numeric(r$mle_cond1)))
  m0 <- val_cr_mle(tt$brts, lam_gen = tt$lambda_gen, mu_gen = tt$mu_gen)
  m1 <- tryCatch(.q(DDD::bd_ML(brts = tt$brts,
                               initparsopt = c(max(m0$pars[1], 1e-3), max(m0$pars[2], 1e-3)),
                               idparsopt = 1:2, cond = 1, soc = 2, btorph = 1,
                               verbose = FALSE)),
                 error = function(e) NULL)
  list(mle0 = as.numeric(m0$pars),
       mle1 = if (is.null(m1)) c(NA, NA) else c(m1$lambda0, m1$mu0))
}

one_tree <- function(tt) {
  r <- ref_for(tt)
  .draws$n <- 0L
  t0 <- proc.time()[3]
  ab <- tryCatch(auto_bounds(tt$brts, model = "cr", link = "linear",
                             train_surv_gam = GAM, verbose = FALSE,
                             num_threads = 1L),
                 error = function(e) NULL)
  el <- proc.time()[3] - t0
  lb <- if (is.null(ab)) c(NA, NA) else as.numeric(ab$lower_bound)
  ub <- if (is.null(ab)) c(NA, NA) else as.numeric(ab$upper_bound)
  inbox <- function(th) if (any(is.na(lb)) || any(is.na(th))) NA
                        else all(th >= lb) && all(th <= ub)
  data.frame(
    variant = VARIANT, tree_id = tt$tree_id, cell = tt$cell, n = tt$n,
    eps = tt$eps, T_crown = tt$T_crown,
    lb_lam = lb[1], ub_lam = ub[1], lb_mu = lb[2], ub_mu = ub[2],
    w_lam = ub[1] - lb[1], w_mu = ub[2] - lb[2],
    area = (ub[1] - lb[1]) * (ub[2] - lb[2]),
    mle0_lam = r$mle0[1], mle0_mu = r$mle0[2],
    mle1_lam = r$mle1[1], mle1_mu = r$mle1[2],
    in0 = inbox(r$mle0), in1 = inbox(r$mle1),
    elapsed = el, n_draws = .draws$n, stringsAsFactors = FALSE)
}

res <- if (CORES > 1L) do.call(rbind, mclapply(trees, one_tree, mc.cores = CORES))
       else do.call(rbind, lapply(trees, one_tree))

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.csv(res, OUT, row.names = FALSE)

agg <- aggregate(cbind(in0, in1, area, elapsed, n_draws) ~ eps, data = res, FUN = mean)
cat("\n== by turnover ==\n"); print(agg, digits = 4)
agg2 <- aggregate(cbind(in0, in1) ~ n, data = res, FUN = mean)
cat("\n== by size ==\n"); print(agg2, digits = 4)
cat(sprintf("\nOVERALL variant=%s  trees=%d  in0=%.3f  in1=%.3f  med_area=%.4g  med_s=%.2f  med_draws=%.0f  tot_s=%.0f\n",
            VARIANT, nrow(res), mean(res$in0), mean(res$in1), median(res$area),
            median(res$elapsed), median(res$n_draws), sum(res$elapsed)))
