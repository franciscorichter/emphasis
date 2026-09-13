## H93 (question): simulate_tree() returns two different proposal densities under the one name `log_q`.
## Checks: (a) same tree, same theta: log_q from method="bdi" vs "thinning" — comparable scale?
##         (b) on a FIXED augmented tree z, does logg_thin(z,theta) - logg_bdi(z,theta) depend on theta?
##         (c) does the returned object record which method actually ran (silent fallback at simulate.R:345)?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape)})
set.seed(93)
tr   <- ape::rphylo(20, 0.5, 0.1)
brts <- sort(ape::branching.times(tr), decreasing = TRUE)
tp   <- brts[1]; bt <- sort(tp - brts[-1])
th1 <- c(0.5, 0.3); th2 <- c(0.7, 0.2)
p8  <- function(th) c(th[1], 0,0,0, th[2], 0,0,0)

## (a)
s_b <- simulate_tree(tree = tr, pars = th1, model = "cr", n_trees = 50L, method = "bdi")
s_t <- simulate_tree(tree = tr, pars = th1, model = "cr", n_trees = 50L, method = "thinning")
cat(sprintf("(a) log_q bdi:      mean %.2f sd %.2f  n=%d\n", mean(s_b$log_q), sd(s_b$log_q), length(s_b$log_q)))
cat(sprintf("(a) log_q thinning: mean %.2f sd %.2f  n=%d\n", mean(s_t$log_q), sd(s_t$log_q), length(s_t$log_q)))
cat("(a) names(simulate_tree bdi):", paste(names(s_b), collapse=","), "\n")
cat("(a) names(simulate_tree thinning):", paste(names(s_t), collapse=","), "\n")

## (b) BDI CR density re-implemented so it can be evaluated at any theta on a fixed z
bdi_logg_cr <- function(species, bt, tp, lam0, mu0) {
  births <- sapply(species, `[`, 1); deaths <- sapply(species, `[`, 2)
  ev <- sort(unique(c(0, bt, births, deaths, tp)))
  logg <- 0
  for (j in seq_len(length(ev) - 1)) {
    t0 <- ev[j]; t1 <- ev[j + 1]
    n  <- sum(births <= t0 & deaths > t0)
    k  <- 2L + sum(bt <= t0)
    logg <- logg - emphasis:::.bdi_integral_cr(t0, t1, n, k, lam0, mu0, tp)
  }
  p <- function(t) emphasis:::.bdi_p_cr(t, lam0, mu0, tp)
  logg + sum(log(lam0 * (1 - p(births)))) + sum(log(mu0 / (1 - p(deaths))))
}
res <- replicate(40, {
  a <- emphasis:::.bdi_augment_one(bt, p8(th1), c(0L,0L,0L), 0L, tp)
  if (is.null(a) || length(a$species) == 0) return(rep(NA, 5))
  df <- emphasis:::.bdi_to_tree_df(a$species, bt, tp)
  g1 <- bdi_logg_cr(a$species, bt, tp, th1[1], th1[2])
  g2 <- bdi_logg_cr(a$species, bt, tp, th2[1], th2[2])
  e1 <- emphasis:::eval_logf(p8(th1), list(df), model = c(0L,0L,0L), link = 0L)
  e2 <- emphasis:::eval_logf(p8(th2), list(df), model = c(0L,0L,0L), link = 0L)
  c(check = a$logg - g1, nmiss = length(a$species), d1 = e1$logg - g1, d2 = e2$logg - g2,
    lw_bdi = e1$logf - a$logg)
})
res <- t(res); res <- res[complete.cases(res), , drop = FALSE]
cat(sprintf("(b) reimpl vs package BDI logg: max abs diff %.2e (n=%d draws)\n", max(abs(res[,1])), nrow(res)))
cat(sprintf("(b) logg_thin - logg_bdi at theta1: range [%.2f, %.2f]\n", min(res[,3]), max(res[,3])))
cat(sprintf("(b) change of that difference theta1->theta2: range [%.2f, %.2f]  (0 everywhere iff theta-independent)\n",
            min(res[,4]-res[,3]), max(res[,4]-res[,3])))
cat(sprintf("(b) cor(change, n_missing) = %.3f\n", cor(res[,4]-res[,3], res[,2])))
cat(sprintf("(b) sd(logf - logg_bdi) across draws = %.2e (CR zero-variance)\n", sd(res[,5])))

## (c) silent fallback: gaussian link is not BDI-supported
s_f <- simulate_tree(tree = tr, pars = th1, model = "cr", link = "gaussian", n_trees = 5L, method = "bdi")
cat("(c) method='bdi' + gaussian link: .bdi_supported =", emphasis:::.bdi_supported(c(0L,0L,0L), 2L),
    "; output names:", paste(names(s_f), collapse=","), "; any 'method' field:", "method" %in% names(s_f), "\n")
cat("(c) log_q returned in that case:", round(s_f$log_q, 2), "\n")
