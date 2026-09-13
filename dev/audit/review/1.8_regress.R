# Behaviour outside the rule: max_lambda budget, timing, rejection counters,
# short pars vectors. Run on both builds; counter calls are guarded.
args <- commandArgs(TRUE)
.libPaths(c(args[1], .libPaths()))
suppressMessages(library(emphasis))
ns  <- asNamespace("emphasis")
aug <- get("augment_trees", envir = ns)
has_viol <- exists("thinning_envelope_violations", envir = ns)
cat("build has counter:", has_viol, "\n")

brts7 <- c(6, 4.5, 3.0, 2.0, 1.2, 0.5)

probe <- function(label, pars, model, link, max_lambda, n = 1000L) {
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(aug(brts = brts7, pars = pars, sample_size = n, maxN = 200L * n,
                    max_missing = 500L, max_lambda = max_lambda, num_threads = 1L,
                    model = model, link = link, rho = 1),
                error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(r, "error")) {
    cat(sprintf("%-34s max_lambda=%-8g ERROR (%s)\n", label, max_lambda,
                substr(conditionMessage(r), 1, 60)))
  } else {
    cat(sprintf("%-34s max_lambda=%-8g ok  rej_lam=%-6d %.2fs\n",
                label, max_lambda, r$rejected_lambda, el))
  }
}

cat("\n-- max_lambda budget, D-active model (safety factor doubles the envelope) --\n")
# nh(start) on the first segment is roughly 2*lambda*(1-exp(-mu*T)); scale lambda up
for (ml in c(1e6, 400, 200, 100, 60, 40)) {
  probe("nd lin beta_0=20 beta_D=-0.5", c(20, 0, 0, -0.5, 0.2, 0, 0, 0),
        c(1L, 0L, 1L), 0L, ml, n = 200L)
}

cat("\n-- timing / rejections, CR and dd (constant-rate branch) --\n")
probe("cr lam=.4 mu=.2", c(0.4,0,0,0,0.2,0,0,0), c(0L,0L,0L), 0L, 1e6, n = 4000L)
probe("dd steep",        c(1.2,-0.15,0,0,0.2,0,0,0), c(1L,0L,0L), 0L, 1e6, n = 4000L)

cat("\n-- short pars vector (no 8-slot guard on augment_trees/rcpp_mce) --\n")
for (L in c(4, 7, 8)) {
  p <- c(0.4, 0, 0, 0, 0.2, 0, 0, 0)[seq_len(L)]
  r <- tryCatch(aug(brts = brts7, pars = p, sample_size = 50L, maxN = 10000L,
                    max_missing = 500L, max_lambda = 1e6, num_threads = 1L,
                    model = c(0L,0L,0L), link = 0L, rho = 1),
                error = function(e) e)
  if (inherits(r, "error"))
    cat(sprintf("length(pars)=%d -> ERROR: %s\n", L, substr(conditionMessage(r),1,70)))
  else
    cat(sprintf("length(pars)=%d -> ok, %d trees, mean rows %.1f\n",
                L, length(r$trees), mean(vapply(r$trees, nrow, 0))))
}
