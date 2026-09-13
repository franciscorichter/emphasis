## H24 replication -- independent check. Run: Rscript dev/audit/checks/H24_replicate.R
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })

## Part 1: deterministic test of the arithmetic on the REAL loop.
## tol = 1e6 => no iteration can ever count as improving => a decay after every
## iteration => the floor is first visible at the check of iteration 30, and
## plateau_count >= 5 is satisfied long before. So plateau must fire at k = 30
## exactly when max_iter >= 30 and can never fire when max_iter = 29.
set.seed(101)
tr <- NULL
while (is.null(tr)) {
  t0 <- ape::rlineage(0.5, 0.2, Tmax = 6)
  t1 <- tryCatch(ape::drop.fossil(t0), error = function(e) NULL)
  if (!is.null(t1) && Ntip(t1) >= 20 && Ntip(t1) <= 40) tr <- t1
}
cat(sprintf("Tree: %d tips, crown age %.2f\n", Ntip(tr), max(branching.times(tr))))
lb <- c(0, 0); ub <- c(2, 1)

fit_cem <- function(max_iter, tol = 1e-4, num_trees = 1, num_particles = 20,
                    model = "cr", link = "linear", extra = list()) {
  ctrl <- c(list(max_iter = max_iter, num_particles = num_particles,
                 num_trees = num_trees, num_threads = 1L, tol = tol,
                 lower_bound = lb, upper_bound = ub, max_time = 120, verbose = FALSE),
            extra)
  estimate_rates(tr, method = "cem", model = model, link = link, control = ctrl)
}
for (mi in c(29L, 30L, 35L)) {
  for (r in 1:2) {
    f <- fit_cem(mi, tol = 1e6)
    cat(sprintf("P1 tol=1e6 max_iter=%d rep=%d: converged=%s k_ran=%d\n",
                mi, r, f$details$converged, length(f$details$best_loglik)))
  }
}

## Part 2: pipeline CEM stage with the pipeline's OWN defaults (50 particles x 5 trees, max_iter 20)
for (r in 1:2) {
  tm <- system.time(pp <- emphasis_pipeline(tr, model = "cr", stages = c("bounds", "cem"),
                           control = list(num_threads = 1L, max_time = 300), verbose = FALSE))[3]
  d <- pp$fits$cem$details
  cat(sprintf("P2 pipeline defaults rep=%d: converged=%s k_ran=%d n_improving=%d secs=%.1f\n",
              r, d$converged, length(d$best_loglik), sum(c(FALSE, diff(d$best_loglik) >= 1e-4)), tm))
}

## Part 3: other model / link, estimate_rates default control (max_iter 50), bigger tol
for (spec in list(list(model = "cr", link = "exponential"), list(model = "dd", link = "linear"))) {
  lb <<- if (spec$model == "dd") c(0, 0, -0.05) else c(0, 0)
  ub <<- if (spec$model == "dd") c(2, 1, 0.05) else c(2, 1)
  for (r in 1:2) {
    f <- fit_cem(50L, model = spec$model, link = spec$link)
    d <- f$details
    cat(sprintf("P3 %s/%s max_iter=50 rep=%d: converged=%s k_ran=%d n_improving=%d loglik=%.3f\n",
                spec$model, spec$link, r, d$converged, length(d$best_loglik),
                sum(c(FALSE, diff(d$best_loglik) >= 1e-4)), f$loglik))
  }
}
## Part 4: does a larger tol (MC-scale) let plateau fire within 50? cr/linear
lb <<- c(0, 0); ub <<- c(2, 1)
for (r in 1:3) {
  f <- fit_cem(50L, tol = 0.5)
  d <- f$details
  cat(sprintf("P4 cr/linear tol=0.5 max_iter=50 rep=%d: converged=%s k_ran=%d n_improving=%d\n",
              r, d$converged, length(d$best_loglik), sum(c(FALSE, diff(d$best_loglik) >= 0.5))))
}
