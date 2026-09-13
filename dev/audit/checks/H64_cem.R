# H64 follow-up: end-to-end CEM with num_threads > maxN (grainsize 0), and
# repeated timing of grainsize 0 vs grainsize maxN.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
hw <- parallel::detectCores()
set.seed(1)
tr   <- ape::rphylo(20, 0.5, 0.2)
brts <- sort(as.numeric(ape::branching.times(tr)), decreasing = TRUE)
pars <- c(0.5, 0, 0, 0, 0.2, 0, 0, 0)

for (rep in 1:3) {
  res <- emphasis:::emphasis_cem(brts, max_iter = 5L, num_points = 20L, max_missing = 1000L,
                          sd_vec = c(0.2, 0, 0, 0, 0.2, 0, 0, 0),
                          lower_bound = c(0.01, 0, 0, 0, 0.001, 0, 0, 0),
                          upper_bound = c(2, 0, 0, 0, 1, 0, 0, 0),
                          maxN = 5L, sample_size = 1L, num_threads = hw,
                          model = c(0L, 0L, 0L), link = 0L)
  if (rep == 1) cat("names:", paste(names(res), collapse = ", "), "\n")
  est <- res$obtained_estim; if (is.null(est)) est <- res$pars
  cat(sprintf("  CEM rep %d: est=(%s) loglik=%s converged=%s iters=%s\n", rep,
              paste(round(est[c(1, 5)], 3), collapse = ", "),
              format(res$loglik, digits = 5), as.character(res$converged),
              as.character(res$iterations %||% NA)))
}

cat("\ntiming: 300 calls, sample_size = maxN = 5\n")
f <- function(nt) emphasis:::augment_trees(brts, pars, 5L, 5L, 10000L, 1e9, as.integer(nt), c(0L,0L,0L), 0L, 1.0)
for (r in 1:3) {
  t0 <- system.time(for (i in 1:300) f(hw))[["elapsed"]]
  t1 <- system.time(for (i in 1:300) f(1L))[["elapsed"]]
  cat(sprintf("  round %d: grain0 (nt=%d) %.2fs | grain5 (nt=1) %.2fs\n", r, hw, t0, t1))
}
