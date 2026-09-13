# (b) end-to-end: a conditioned and an unconditioned MCEM fit on both samplers.
# .run_mcem takes cond_fun directly, so no GAM is needed for the seam test.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
library(emphasis)
set.seed(21)
brts <- sort(ape::branching.times(ape::rphylo(40, birth = 0.6, death = 0.2)), decreasing = TRUE)
Tc <- max(brts)
logP <- function(p8) {                       # log P(both crown clades survive)
  lam <- p8[1]; mu <- p8[5]
  if (lam <= 0) return(-700)
  d <- lam - mu
  ps <- 1 - (mu * (1 - exp(-d * Tc))) / (lam - mu * exp(-d * Tc))
  log(max(min(ps, 1)^2, 1e-300))
}
lb8 <- c(0, 0, 0, 0, 0, 0, 0, 0); ub8 <- c(3, 0, 0, 0, 3, 0, 0, 0)
ip8 <- c(0.6, 0, 0, 0, 0.2, 0, 0, 0)
ctrl <- function(samp) utils::modifyList(emphasis:::estimate_rates_control("mcem"),
  list(sampling = samp, sample_size = 60L, num_trees = 60L, max_iter = 10L,
       num_threads = 1L, verbose = FALSE, lower_bound = lb8, upper_bound = ub8))

out <- list()
for (samp in c("bdi", "dynamic_fresh")) for (cf in c("uncond", "cond")) {
  set.seed(99)
  r <- emphasis:::.run_mcem(brts, ip8, lb8, ub8, ctrl(samp),
                            model = c(0L, 0L, 0L), link = 0L,
                            cond_fun = if (cf == "cond") logP else NULL)
  out[[paste(samp, cf)]] <- c(lambda = r$pars[1], mu = r$pars[5],
                              loglik = r$loglik, iters = r$iterations)
  cat(sprintf("%-16s %-7s lambda=%.5f mu=%.5f loglik=%.4f iters=%d stop=%s\n",
              samp, cf, r$pars[1], r$pars[5], r$loglik, r$iterations, r$stop_reason))
}
d_bdi  <- out[["bdi cond"]] - out[["bdi uncond"]]
d_thin <- out[["dynamic_fresh cond"]] - out[["dynamic_fresh uncond"]]
cat("\nconditioning shift (lambda, mu):\n")
cat("  bdi     :", round(d_bdi[1:2], 5), "\n")
cat("  thinning:", round(d_thin[1:2], 5), "\n")
cat("  same sign on both samplers:",
    all(sign(d_bdi[1:2]) == sign(d_thin[1:2])), "\n")
