lib <- commandArgs(trailingOnly=TRUE)[1]
.libPaths(c(lib, .libPaths()))
library(emphasis)
cat("### lib:", lib, "\n")
brts20 <- c(5,4.4,3.9,3.6,3.1,2.8,2.5,2.2,2.0,1.7,1.5,1.3,1.1,0.9,0.75,0.6,0.45,0.3,0.15)
run <- function(lab, ...) {
  r <- tryCatch(suppressWarnings(estimate_rates(...)), error=function(e) e)
  if (inherits(r,"error")) cat(sprintf("[%s] ERROR: %s\n", lab, conditionMessage(r)))
  else cat(sprintf("[%s] pars=%s loglik=%.4f stop=%s iters=%s\n", lab,
      paste(round(r$pars,5), collapse=","), r$loglik,
      paste0(c(r$stop_reason,"<NULL>")[1]), paste0(c(r$iterations,"<NULL>")[1])))
  invisible(r)
}
say <- function(...) cat("==", ..., "\n")

say("B1 linear cr, bdi, EXPLICIT init lambda==mu=0.5, box [0,1]x[0.5,0.5] (the case the new stop() refuses)")
set.seed(11); run("B1", brts20, model="cr", method="mcem", init_pars=c(0.5,0.5),
  control=list(lower_bound=c(0,0.5), upper_bound=c(1,0.5), sampling="bdi",
               sample_size=20L, max_iter=5L, num_threads=1L))

say("B2 exponential cr, bdi, EXPLICIT init at the symmetric midpoint (-0.5,-0.5)")
set.seed(12); run("B2", brts20, model="cr", method="mcem", link="exponential",
  init_pars=c(-0.5,-0.5),
  control=list(lower_bound=c(-2,-2), upper_bound=c(1,1), sampling="bdi",
               sample_size=20L, max_iter=2L, num_threads=1L))

say("B3 dd, bdi, EXPLICIT init with lambda0==mu0")
set.seed(13); run("B3", brts20, model="dd", method="mcem", init_pars=c(0.5,0,0.5,0),
  control=list(lower_bound=c(0,-0.1,0,-0.1), upper_bound=c(1,0.1,1,0.1),
               sampling="bdi", sample_size=20L, max_iter=3L, num_threads=1L))

say("B4 dd, bdi, DEFAULT init on the same symmetric box (guard fires)")
set.seed(13); run("B4", brts20, model="dd", method="mcem",
  control=list(lower_bound=c(0,-0.1,0,-0.1), upper_bound=c(1,0.1,1,0.1),
               sampling="bdi", sample_size=20L, max_iter=3L, num_threads=1L))

say("B5 thinning, DEFAULT init, symmetric box [0,4]^2 (init silently shifted post-fix)")
set.seed(14); run("B5", brts20, model="cr", method="mcem",
  control=list(lower_bound=c(0,0), upper_bound=c(4,4), sampling="dynamic_fresh",
               sample_size=20L, max_iter=3L, max_time=60, num_threads=1L))

say("B6 nd + gaussian link, DEFAULT init, symmetric box (index/guard sanity)")
set.seed(15); run("B6", brts20, model="nd", method="mcem", link="gaussian",
  control=list(lower_bound=c(0,-1,-1,0,-1,-1), upper_bound=c(2,1,1,2,1,1),
               sample_size=10L, max_iter=1L, max_time=60, num_threads=1L))

say("B7 mu > lambda default init: box lambda [0,1], mu [2,3] (no guard, mu>lambda)")
set.seed(16); run("B7", brts20, model="cr", method="mcem",
  control=list(lower_bound=c(0,2), upper_bound=c(1,3), sampling="bdi",
               sample_size=20L, max_iter=3L, num_threads=1L))
