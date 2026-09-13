## H57: BDI logg vs thinning log q on the same augmented tree; is the difference theta-independent?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape)})
set.seed(2)
tr   <- ape::rphylo(20, 0.5, 0.1)
brts <- sort(ape::branching.times(tr), decreasing = TRUE)
tp   <- brts[1]; bt <- sort(tp - brts[-1])

# Re-implementation of the BDI CR logg for a fixed species list at arbitrary theta.
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

th1 <- c(0.5, 0.3); th2 <- c(0.7, 0.2)
p8  <- function(th) c(th[1], 0,0,0, th[2], 0,0,0)

# Draw BDI trees at th1 and sanity-check the re-implementation against the package logg
res <- replicate(30, {
  a <- emphasis:::.bdi_augment_one(bt, p8(th1), c(0L,0L,0L), 0L, tp)
  if (is.null(a) || length(a$species) == 0) return(rep(NA, 6))
  df <- emphasis:::.bdi_to_tree_df(a$species, bt, tp)
  g1 <- bdi_logg_cr(a$species, bt, tp, th1[1], th1[2])
  g2 <- bdi_logg_cr(a$species, bt, tp, th2[1], th2[2])
  e1 <- emphasis:::eval_logf(p8(th1), list(df), model = c(0L,0L,0L), link = 0L)
  e2 <- emphasis:::eval_logf(p8(th2), list(df), model = c(0L,0L,0L), link = 0L)
  c(pkg_minus_reimpl = a$logg - g1, n_missing = length(a$species),
    d1 = e1$logg - g1, d2 = e2$logg - g2,
    lw_bdi = e1$logf - a$logg, lw_thin = e1$logf - e1$logg)
})
res <- t(res); res <- res[complete.cases(res), ]
cat("max |package logg - reimplementation| at theta1:", max(abs(res[, 1])), "\n")
cat("n_missing per draw:", paste(res[, 2], collapse = " "), "\n")
cat("logg_thin - logg_bdi at theta1 (per draw):", round(res[, 3], 3), "\n")
cat("logg_thin - logg_bdi at theta2 (per draw):", round(res[, 4], 3), "\n")
cat("change of that difference between theta1 and theta2:", round(res[, 4] - res[, 3], 3), "\n")
cat(sprintf("sd(lw_bdi) across draws = %.2e  (zero-variance CR property)\n", sd(res[, 5])))
cat(sprintf("sd(lw_thin) across draws = %.2e  (thinning logg paired with BDI draws)\n", sd(res[, 6])))
cat("range lw_bdi:", range(res[, 5]), "\n")
# What simulate_tree returns as log_q under each method
s_b <- simulate_tree(tree = tr, pars = th1, model = "cr", n_trees = 3L, method = "bdi")
s_t <- simulate_tree(tree = tr, pars = th1, model = "cr", n_trees = 3L, method = "thinning")
cat("simulate_tree log_q (bdi):     ", round(s_b$log_q, 3), "\n")
cat("simulate_tree log_q (thinning):", round(s_t$log_q, 3), "\n")
