.libPaths(c(Sys.getenv("EMPH_LIB"), .libPaths()))
suppressMessages(library(emphasis))
brts1 <- c(5,3.822147,3.059192,1.966150,1.604249,1.358484,0.622046,0.467531,
           0.393533,0.285495,0.276148,0.225482,0.201295,0.154060,0.057677,
           0.044974,0.041380,0.038051,0.013895)
for (nt in c(200L)) {
  t0 <- proc.time()[3]
  f <- estimate_rates(brts1, model="cr", method="mcem", init_pars=c(0.5,0.1),
    control=list(sampling="dynamic_fresh", num_trees=nt, maxN=20000L, max_iter=40L,
      lower_bound=c(0,0), upper_bound=c(2,1), num_threads=1L))
  cat(sprintf("%s num_trees=%d stop=%s iters=%s time=%.1fs pars=%s\n",
    basename(Sys.getenv("EMPH_LIB")), nt,
    f$details$stop_reason,
    if (is.null(f$details$iterations)) nrow(f$details$mcem) else f$details$iterations,
    proc.time()[3]-t0, paste(round(f$pars,4), collapse=",")))
}
