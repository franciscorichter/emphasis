.libPaths(c(commandArgs(TRUE)[1], .libPaths())); library(emphasis)
brts12 <- c(6.385824,2.063997,1.19255,0.923743,0.884126,0.822976,0.725585,0.718801,0.539214,0.077287,0.068123)
fit_box <- function(ub, seed) { set.seed(seed)
  estimate_rates(brts12, model="cr", method="mcem", init_pars=c(1.2,0.9),
    control=list(lower_bound=c(0,0), upper_bound=c(ub,ub), sampling="bdi",
                 sample_size=200L, num_threads=1L, max_iter=25L, tol=1e-2, patience=3L)) }
for (seed in c(21, 7)) {
  P <- t(sapply(c(3,30,300), function(u) { f <- fit_box(u, seed)
        c(unname(f$pars), f$details$iterations, f$loglik) }))
  colnames(P) <- c("lam","mu","iters","loglik")
  cat("seed", seed, "\n"); print(round(P,4))
  cat("  range lam:", round(diff(range(P[,1])),4), " range mu:", round(diff(range(P[,2])),4), "\n")
}
