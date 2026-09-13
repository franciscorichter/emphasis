.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
suppressMessages(library(emphasis)); ns <- asNamespace("emphasis")
brts1 <- c(5,3.822147,3.059192,1.966150,1.604249,1.358484,0.622046,0.467531,
           0.393533,0.285495,0.276148,0.225482,0.201295,0.154060,0.057677,
           0.044974,0.041380,0.038051,0.013895)
cr <- c(0L,0L,0L); ex <- function(p) ns$.expand_pars(p, cr)
f <- function() { set.seed(20); r <- ns$.mcem_dynamic_fresh(brts1, ex(c(0.5,0.1)),
  sample_size=30L, maxN=3000L, max_missing=1e4, lower_bound=ex(c(0,0)),
  upper_bound=ex(c(2,1)), max_iter=3L, xtol=1e-3, tol=1e-2, patience=3L,
  num_threads=1L, verbose=FALSE, model=cr, link=0L, max_time=120)
  c(r$pars[1], r$pars[5], r$final_IS$fhat) }
for (i in 1:3) cat(sprintf("run %d: %s\n", i, paste(sprintf("%.8f", f()), collapse=" ")))
