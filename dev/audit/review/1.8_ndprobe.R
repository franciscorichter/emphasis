# (a) locate the max_lambda threshold shift; (b) sampler-vs-charged-density
# self-consistency P(0 missing) vs exp(logg_empty) on the D path.
args <- commandArgs(TRUE)
.libPaths(c(args[1], .libPaths()))
suppressMessages(library(emphasis))
ns <- asNamespace("emphasis"); aug <- get("augment_trees", envir = ns)
has_viol <- exists("thinning_envelope_violations", envir = ns)
viol <- if (has_viol) get("thinning_envelope_violations", envir = ns) else function(...) NA
brts7 <- c(6, 4.5, 3.0, 2.0, 1.2, 0.5)
is_missing_row <- function(df) !(df$t_ext == 0 | df$t_ext == 1e11 | df$t_ext == 5e10)

cat("-- (a) max_lambda threshold, nd linear beta_0=1.5 beta_D=-0.5 --\n")
for (ml in c(200, 100, 50, 30, 20, 15, 12, 10, 8, 6)) {
  r <- tryCatch(aug(brts = brts7, pars = c(1.5, 0, 0, -0.5, 0.3, 0, 0, 0),
                    sample_size = 300L, maxN = 30000L, max_missing = 500L,
                    max_lambda = ml, num_threads = 1L, model = c(1L,0L,1L),
                    link = 0L, rho = 1), error = function(e) e)
  cat(sprintf("  max_lambda=%-5g %s\n", ml,
              if (inherits(r,"error")) "THROWS (all lambda-rejected)"
              else sprintf("ok, rej_lam=%d", r$rejected_lambda)))
}

cat("\n-- (b) P(0 missing) vs exp(logg_empty), D path --\n")
for (bd in c(-1, -0.2, 0, 0.2, 1)) {
  if (has_viol) viol(reset = TRUE)
  n <- 20000L
  r <- aug(brts = brts7, pars = c(0.4, -0.01, 0, bd, 0.2, 0, 0, 0),
           sample_size = n, maxN = 100L * n, max_missing = 500L,
           max_lambda = 1e6, num_threads = 1L, model = c(1L,0L,1L),
           link = 0L, rho = 1)
  nm <- vapply(r$trees, function(df) sum(is_missing_row(df)), numeric(1))
  p0 <- mean(nm == 0); se <- sqrt(p0*(1-p0)/n)
  lge <- r$logg[nm == 0]
  pq <- exp(mean(lge))
  cat(sprintf("  beta_D=%+.1f viol=%-6s p0=%.5f  exp(logg_empty)=%.5f  z=%+.2f  (sd logg_empty=%.2e)\n",
              bd, as.character(viol()), p0, pq, (p0 - pq)/se, sd(lge)))
}
