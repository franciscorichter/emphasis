.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
# Each case runs in a child R with a hard 60 s wall clock (perl alarm), since garbage bounds may hang nlopt.
run_case <- function(code, label) {
  cat("=== ", label, "\n")
  pre <- '.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); suppressMessages(library(emphasis)); brts <- c(6,4.5,3.2,2.5,1.8,1.1,0.7,0.3); model <- c(0L,0L,0L); ip <- emphasis:::.expand_pars(c(0.5,0.1), model); '
  f <- tempfile(fileext = ".R"); writeLines(paste0(pre, code), f)
  t0 <- proc.time()[3]
  out <- suppressWarnings(system2("perl", c("-e", shQuote("alarm 60; exec @ARGV"), "Rscript", f), stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status"); cat(out, sep = "\n"); cat(sprintf("[exit status %s, %.1fs]\n", if (is.null(st)) 0 else st, proc.time()[3]-t0))
}
run_case('r <- emphasis:::em_cpp(brts, ip, 20L, 1000L, 1e4, 1e6, rep(0,8), c(2,0,0,0,1,0,0,0), 1e-3, 1L, FALSE, model, 0L, 1.0, NULL); print(c(nlopt=r$nlopt, r$estimates))', "control: 8-element bounds")
for (k in 1:3) run_case('r <- emphasis:::em_cpp(brts, ip, 20L, 1000L, 1e4, 1e6, c(0,0), c(2,1), 1e-3, 1L, FALSE, model, 0L, 1.0, NULL); print(c(nlopt=r$nlopt, r$estimates))', paste("em_cpp: 2-element bounds, 8-element pars, try", k))
run_case('a <- emphasis:::augment_trees(brts, ip, 20L, 1000L, 1e4, 1e6, 1L, model, 0L, 1.0); w <- exp(a$logf - a$logg); w <- w/mean(w); es <- list(trees=a$trees, weights=w, rejected=0L, rejected_overruns=0L, rejected_lambda=0L, rejected_zero_weights=0L, time=0, fhat=0); r <- emphasis:::m_cpp(es, ip, "rpd1", c(0,0), c(2,1), 1e-3, 1L, model, 0L, 1.0, NULL); print(c(nlopt=r$nlopt, r$estimates))', "m_cpp: 2-element bounds")
