# Phase 2 check (dev/latent-marks/design.md): ML recovery of the latent-rate preset on complete
# trees simulated with marks (cr, lambda0 0.35, mu0 0.07, crown age 5, trend-free ClaDS:
# log_alpha = -sigma^2/2), 200 trees per sigma, joint fit of (lambda0, mu0, sigma).  Run from the
# package root: Rscript dev/latent-marks/recovery-complete.R

suppressMessages(devtools::load_all(quiet = TRUE))
mb <- .resolve_model("cr")
nll <- function(th, trees, la_fixed) {
  sg <- exp(th[3]); mk <- .latent_rate_marks(sg, log_alpha = -sg^2 / 2, kmax = 12)
  -sum(.eval_logf_marks(.expand_pars(th[1:2], mb), trees, mb, mk))
}
for (sg in c(0.5, 0)) {
  set.seed(200 + 10 * sg); mk <- .latent_rate_marks(sg, log_alpha = -sg^2 / 2, kmax = 12)
  tr <- list(); st <- character(0)
  while (length(tr) < 200) { r <- .simulate_marks(c(log(0.35), log(0.07)), "cr", 5, mk, max_N = 20000, max_tries = 1L)
    st <- c(st, r$status); if (r$status == "done" && nrow(r$Ltable) > 2) tr[[length(tr) + 1]] <- r$frame }
  cat(sprintf("sigma %.1f: %s; lineages median %d, max %d\n", sg, paste(names(table(st)), table(st), collapse = " "),
      as.integer(median(sapply(tr, nrow))), max(sapply(tr, nrow))))
  t0 <- Sys.time()
  f <- optim(c(log(0.3), log(0.1), log(0.3)), nll, trees = tr, control = list(maxit = 400))
  cat(sprintf("  fit: lambda0 %.3f (0.35)  mu0 %.3f (0.07)  sigma %.3f (%.2f)  conv %d  %.0f s\n",
      exp(f$par[1]), exp(f$par[2]), exp(f$par[3]), sg, f$convergence, as.numeric(Sys.time() - t0, units = "secs")))
}
