# (b) contrast: the same weight-convention test on the PRE-FIX build.
# Trees are produced once by the POST-FIX build and reused, so only the
# M-step objective differs between the two runs.
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, .libPaths()))
library(emphasis)
d <- readRDS("/tmp/emphasis_review_estep.rds")
trees <- d$trees; lw <- d$lw; init <- d$init; lb8 <- d$lb8; ub8 <- d$ub8
mb <- c(1L, 0L, 0L); Tc <- d$Tc

mk <- function(w) list(trees = trees, weights = w, rejected = 0L,
                       rejected_overruns = 0L, rejected_lambda = 0L,
                       rejected_zero_weights = 0L, time = 0, fhat = 0)
w_mean1 <- { v <- exp(lw - max(lw)); v / sum(v) * length(v) }
w_max   <- exp(lw - max(lw))
logP <- function(p8) {
  lam <- p8[1]; mu <- p8[5]
  if (lam <= 0) return(-700)
  dd <- lam - mu
  ps <- 1 - (mu * (1 - exp(-dd * Tc))) / (lam - mu * exp(-dd * Tc))
  log(max(min(ps, 1)^2, 1e-300))
}
run <- function(w, cond) {
  r <- emphasis:::m_cpp(e_step = mk(w), init_pars = init, plugin = "rpd1",
                        lower_bound = lb8, upper_bound = ub8, xtol_rel = 1e-6,
                        num_threads = 1L, model = mb, link = 0L, rho = 1,
                        rconditional = cond)
  as.numeric(r$estimates)[c(1, 2, 5)]
}
cat("lib =", lib, "\n")
cat("uncond mean-1    :", run(w_mean1, NULL), "\n")
cat("uncond max-scaled:", run(w_max, NULL), "\n")
c1 <- run(w_mean1, logP); c2 <- run(w_max, logP)
cat("cond   mean-1    :", c1, "\n")
cat("cond   max-scaled:", c2, "\n")
cat("cond  max|diff|  :", max(abs(c1 - c2)), "\n")
