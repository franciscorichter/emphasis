## H24 replication, part 2 (other link/model + MC-scale tol). Run: Rscript dev/audit/checks/H24_replicate2.R
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
set.seed(24)
tr <- NULL
while (is.null(tr)) {
  t0 <- ape::rlineage(0.4, 0.1, Tmax = 8)
  t1 <- tryCatch(ape::drop.fossil(t0), error = function(e) NULL)
  if (!is.null(t1) && Ntip(t1) >= 15 && Ntip(t1) <= 30) tr <- t1
}
cat(sprintf("Tree: %d tips\n", Ntip(tr)))
fit_cem <- function(max_iter, lb, ub, tol = 1e-4, model = "cr", link = "linear") {
  estimate_rates(tr, method = "cem", model = model, link = link,
                 control = list(max_iter = max_iter, num_particles = 20, num_trees = 1,
                                num_threads = 1L, tol = tol, lower_bound = lb, upper_bound = ub,
                                max_time = 60, verbose = FALSE))
}
rep_line <- function(tag, f, tol = 1e-4) {
  d <- f$details
  cat(sprintf("%s: converged=%s k_ran=%d n_improving=%d loglik=%.3f\n", tag, d$converged,
              length(d$best_loglik), sum(c(FALSE, diff(d$best_loglik) >= tol)), f$loglik))
}
for (r in 0) NULL
for (r in 1:2) rep_line(sprintf("P3 dd/linear max_iter=50 rep=%d", r),
                        fit_cem(50L, lb = c(0, -0.03, 0, 0), ub = c(1, 0.01, 0.5, 0.001), model = "dd"))
for (r in 1:3) rep_line(sprintf("P4 cr/linear tol=0.5 max_iter=50 rep=%d", r),
                        fit_cem(50L, lb = c(0, 0), ub = c(1, 0.5), tol = 0.5), tol = 0.5)
