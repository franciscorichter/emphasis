args <- commandArgs(trailingOnly = TRUE); lib <- args[1]
.libPaths(c(lib, .libPaths())); suppressPackageStartupMessages(library(emphasis))
cat("lib:", find.package("emphasis"), "\n")
brts12 <- c(6.385824, 2.063997, 1.19255, 0.923743, 0.884126, 0.822976, 0.725585, 0.718801, 0.539214, 0.077287, 0.068123)
# Unit-scale probe: the same tree with branching times x s; MLE rates scale by 1/s.
for (s in c(1, 1000)) {
  set.seed(21)
  t0 <- proc.time()[3]
  fit <- estimate_rates(brts12 * s, model = "cr", method = "mcem", init_pars = c(1.2, 0.9) / s,
                        control = list(lower_bound = c(0, 0), upper_bound = c(3, 3) / s, sampling = "bdi",
                                       sample_size = 200L, num_threads = 1L, max_iter = 25L, tol = 1e-2, patience = 3L))
  m <- fit$details$mcem
  cat(sprintf("s=%g: stop=%s iters=%d pars*s=(%.4f, %.4f) loglik=%.3f  delta_max trace: %s  time=%.0fs\n", s, fit$details$stop_reason, fit$details$iterations, fit$pars[1]*s, fit$pars[2]*s, fit$loglik, paste(signif(m$delta_max[m$m_step %in% TRUE | is.na(m$m_step)], 2), collapse=" "), proc.time()[3]-t0))
}
ml <- tryCatch(DDD::bd_ML(brts12, cond = 0, btorph = 1, soc = 2, verbose = FALSE, initparsopt = c(0.5, 0.2)), error = function(e) NULL)
if (!is.null(ml)) cat(sprintf("DDD::bd_ML (s=1): lambda=%.4f mu=%.4f\n", ml$lambda0, ml$mu0))
