.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
library(emphasis)
brts20 <- c(5,4.4,3.9,3.6,3.1,2.8,2.5,2.2,2.0,1.7,1.5,1.3,1.1,0.9,0.75,0.6,0.45,0.3,0.15)
say <- function(...) cat("==", ..., "\n")
tryit <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)),
                warning = function(w) paste("WARNING:", conditionMessage(w)))
  cat(sprintf("[%s] ", lab)); if (is.character(r)) cat(r, "\n") else print(r); invisible(r)
}

say("A1: user-supplied init lambda == mu under BDI (guard does not apply)")
set.seed(1)
f <- tryCatch(estimate_rates(brts20, model="cr", method="mcem", init_pars=c(0.5,0.5),
      control=list(lower_bound=c(0,0), upper_bound=c(1,1), sampling="bdi",
                   sample_size=20L, max_iter=5L, num_threads=1L)),
      error=function(e) e)
if (inherits(f,"error")) cat("ERROR:", conditionMessage(f), "\n") else
  cat("pars", f$pars, "loglik", f$loglik, "stop", f$stop_reason, "iters", f$iterations, "\n")

say("A2: default init, near-symmetric box ub=c(1, 1+2e-9) -> lambda-mu = 1e-9")
set.seed(2)
f2 <- tryCatch(estimate_rates(brts20, model="cr", method="mcem",
      control=list(lower_bound=c(0,0), upper_bound=c(1,1+2e-9), sampling="bdi",
                   sample_size=20L, max_iter=5L, num_threads=1L)),
      error=function(e) e)
if (inherits(f2,"error")) cat("ERROR:", conditionMessage(f2), "\n") else
  cat("pars", f2$pars, "loglik", f2$loglik, "stop", f2$stop_reason, "iters", f2$iterations, "\n")

say("A3: maxN < num_trees but sampling = BDI (maxN unused by BDI)")
r3 <- tryCatch(estimate_rates(brts20, model="cr", method="mcem", init_pars=c(1,0.3),
      control=list(lower_bound=c(0,0), upper_bound=c(4,4), sampling="bdi",
                   num_trees=300L, maxN=200L, max_iter=1L, num_threads=1L)),
      error=function(e) conditionMessage(e))
cat(if (is.character(r3)) paste("ERROR:", r3) else "ran fine", "\n")

say("A4: maxN = NA")
r4 <- tryCatch(estimate_rates(brts20, model="cr", method="mcem", init_pars=c(1,0.3),
      control=list(lower_bound=c(0,0), upper_bound=c(4,4), sampling="dynamic_fresh",
                   num_trees=20L, maxN=NA, max_iter=1L, num_threads=1L)),
      error=function(e) conditionMessage(e))
cat(if (is.character(r4)) paste("ERROR:", r4) else "ran fine", "\n")

say("A5: n_pars when everything is fixed (cr, lb == ub)")
r5 <- tryCatch(estimate_rates(brts20, model="cr", method="mcem",
      control=list(lower_bound=c(0.5,0.2), upper_bound=c(0.5,0.2), sampling="bdi",
                   sample_size=20L, max_iter=2L, num_threads=1L)),
      error=function(e) conditionMessage(e))
if (is.character(r5)) cat("ERROR:", r5, "\n") else
  cat("n_pars", r5$n_pars, "class", class(r5$n_pars), "AIC", r5$AIC, "loglik", r5$loglik, "\n")

say("A6: n_pars type + compare_models with a fixed parameter")
set.seed(6)
fa <- suppressWarnings(estimate_rates(brts20, model="cr", method="mcem", init_pars=c(0.6,0.2),
      control=list(lower_bound=c(0,0), upper_bound=c(2,2), sampling="bdi",
                   sample_size=20L, max_iter=2L, num_threads=1L)))
fb <- suppressWarnings(estimate_rates(brts20, model="dd", method="mcem",
      control=list(lower_bound=c(0,-0.1,0,0), upper_bound=c(1,0,0.5,0), sampling="bdi",
                   sample_size=20L, max_iter=2L, num_threads=1L)))
cat("typeof n_pars:", typeof(fa$n_pars), typeof(fb$n_pars), "\n")
print(tryCatch(emphasis:::compare_models(CR=fa, DD=fb), error=function(e) conditionMessage(e)))

say("A7: cem fit -> stop_reason / iterations / print")
set.seed(7)
fc <- tryCatch(suppressWarnings(estimate_rates(brts20, model="cr", method="cem",
      control=list(lower_bound=c(0,0), upper_bound=c(2,2), num_particles=6L,
                   max_iter=2L, sample_size=5L, num_threads=1L, max_time=60))),
      error=function(e) e)
if (inherits(fc,"error")) cat("ERROR:", conditionMessage(fc),"\n") else {
  cat("stop_reason:", fc$stop_reason, " iterations:", fc$iterations,
      " details$converged:", fc$details$converged, "\n")
  print(fc)
}
