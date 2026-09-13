## H76 replication: vary model (dd) and link (exponential), smaller n, 3 replicates.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
n <- 400L
run <- function(pars_mat, model, link, k) {
  b <- simulate_tree(pars = pars_mat, max_t = 5, model = model, link = link,
                     max_tries = k, useDDD = FALSE, num_threads = 1L)
  sp <- vapply(b$simulations, `[[`, 0.0, "survival_prob")
  done <- mean(vapply(b$simulations, function(s) s$status == "done", TRUE))
  c(reported = b$survival_prob, done = done, mean_sp = mean(sp),
    n_vals = length(unique(round(sp, 3))))
}
set.seed(1)
# dd, linear link: lambda = 0.6 - 0.01 N, mu = 0.3
dd_mat <- cbind(beta_0 = rep(0.6, n), beta_N = rep(-0.01, n), gamma_0 = rep(0.3, n), gamma_N = rep(0, n))
# cr, exponential link: lambda = exp(-0.7)=0.50, mu = exp(-1.2)=0.30
ex_mat <- cbind(beta_0 = rep(-0.7, n), gamma_0 = rep(-1.2, n))
for (rep in 1:3) for (k in c(0L, 1L, 3L)) {
  r1 <- run(dd_mat, "dd", "linear", k); r2 <- run(ex_mat, "cr", "exponential", k)
  cat(sprintf("rep%d k=%d | dd/linear: reported %.3f done %.3f (mean sp %.3f, %d distinct) | cr/exp: reported %.3f done %.3f (%d distinct)\n",
              rep, k, r1["reported"], r1["done"], r1["mean_sp"], r1["n_vals"], r2["reported"], r2["done"], r2["n_vals"]))
  if (k == 0L) stopifnot(isTRUE(all.equal(r1["reported"], r1["done"], check.attributes = FALSE)),
                         isTRUE(all.equal(r2["reported"], r2["done"], check.attributes = FALSE)))
  else stopifnot(r1["reported"] < r1["done"], r2["reported"] < r2["done"])
}
# single-tree default (max_tries = 1): what values does survival_prob take?
sv <- replicate(200, simulate_tree(pars = c(0.5, 0.3), max_t = 5, model = "cr", useDDD = FALSE)$survival_prob)
cat("single-tree default max_tries=1 values:", paste(sort(unique(sv)), collapse = ","), "\n")
cat("OK\n")
