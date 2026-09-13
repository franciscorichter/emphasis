## H74 replication — variations the verifier did not run:
##  A. moderate turnover (mu/lambda = 0.5, T = 12): is the box fine there?
##  B. high turnover, exponential link: same exclusion on the log scale?
##  C. bd_ML(cond = 0) as well as cond = 1 (emphasis cond = NULL ~ cond = 0)
##  D. hand-supplied wide box + MCEM on a high-turnover tree: can the
##     estimator itself leave the auto box (is the box, not the sampler, the cap)?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape); library(DDD) })

sim_trees <- function(lam, mu, T0, n, nmin = 15, nmax = 60) {
  out <- list(); tries <- 0L
  while (length(out) < n && tries < 600L) {
    tries <- tries + 1L
    s <- simulate_tree(pars = c(lam, mu), max_t = T0, model = "cr", link = "linear",
                       max_tries = 200, max_lin = 5000)
    if (s$status == "done" && !is.null(s$tes)) {
      N <- Ntip(s$tes); if (N >= nmin && N <= nmax) out[[length(out) + 1L]] <- s$tes
    }
  }
  out
}
inside <- function(p, lb, ub) all(p >= lb) && all(p <= ub)
mle <- function(tr, cond) {
  brts <- sort(branching.times(tr), decreasing = TRUE)
  m <- tryCatch(suppressMessages(bd_ML(brts = brts, initparsopt = c(0.5, 0.3), idparsopt = 1:2,
                         cond = cond, soc = 2, btorph = 1, verbose = FALSE)),
                error = function(e) NULL)
  if (is.null(m)) c(NA, NA) else c(m$lambda0, m$mu0)
}

## ---- A. moderate turnover ------------------------------------------------
cat("=== A. moderate turnover: lambda=0.5, mu=0.25, T=12 ===\n")
trA <- sim_trees(0.5, 0.25, 12, 3)
for (tr in trA) {
  ab <- auto_bounds(tr, model = "cr", link = "linear", train_surv_gam = FALSE, verbose = FALSE)
  m1 <- mle(tr, 1); m0 <- mle(tr, 0)
  cat(sprintf("N=%d box lam=[%.3f,%.3f] mu=[%.3f,%.3f] | truth(0.5,0.25) in=%s | bd_ML c1=(%.3f,%.3f) in=%s | c0=(%.3f,%.3f) in=%s\n",
              Ntip(tr), ab$lower_bound[1], ab$upper_bound[1], ab$lower_bound[2], ab$upper_bound[2],
              inside(c(0.5, 0.25), ab$lower_bound, ab$upper_bound),
              m1[1], m1[2], inside(m1, ab$lower_bound, ab$upper_bound),
              m0[1], m0[2], inside(m0, ab$lower_bound, ab$upper_bound)))
}

## ---- B/C. high turnover, both links, both DDD conditionings ---------------
cat("\n=== B/C. high turnover: lambda=1, mu=0.9, T=20 ===\n")
trB <- sim_trees(1, 0.9, 20, 3)
for (tr in trB) {
  abL <- auto_bounds(tr, model = "cr", link = "linear", train_surv_gam = FALSE, verbose = FALSE)
  abE <- auto_bounds(tr, model = "cr", link = "exponential", train_surv_gam = FALSE, verbose = FALSE)
  m1 <- mle(tr, 1); m0 <- mle(tr, 0)
  cat(sprintf("N=%d linear box lam=[%.3f,%.3f] mu=[%.3f,%.3f]; exp box (nat. scale) lam=[%.3f,%.3f] mu=[%.4f,%.3f]\n",
              Ntip(tr), abL$lower_bound[1], abL$upper_bound[1], abL$lower_bound[2], abL$upper_bound[2],
              exp(abE$lower_bound[1]), exp(abE$upper_bound[1]), exp(abE$lower_bound[2]), exp(abE$upper_bound[2])))
  cat(sprintf("   truth in linear=%s exp=%s | bd_ML cond=1 (%.3f,%.3f) in linear=%s exp=%s | cond=0 (%.3f,%.3f) in linear=%s exp=%s\n",
              inside(c(1, .9), abL$lower_bound, abL$upper_bound),
              inside(log(c(1, .9)), abE$lower_bound, abE$upper_bound),
              m1[1], m1[2], inside(m1, abL$lower_bound, abL$upper_bound), inside(log(m1), abE$lower_bound, abE$upper_bound),
              m0[1], m0[2], inside(m0, abL$lower_bound, abL$upper_bound), inside(log(m0), abE$lower_bound, abE$upper_bound)))
}

## ---- D. hand-supplied wide box: can MCEM leave the auto box? --------------
cat("\n=== D. MCEM with hand box [0,3]x[0,3] on high-turnover tree 1 ===\n")
tr <- trB[[1]]
abL <- auto_bounds(tr, model = "cr", link = "linear", train_surv_gam = FALSE, verbose = FALSE)
m0 <- mle(tr, 0)
fits <- lapply(1:2, function(r) tryCatch(
  estimate_rates(tr, model = "cr", link = "linear", method = "mcem",
                 init_pars = c(0.6, 0.4), cond = NULL,
                 control = list(lower_bound = c(0, 0), upper_bound = c(3, 3),
                                max_iter = 10, max_time = 60, sample_size = 100,
                                num_threads = 1, verbose = FALSE)),
  error = function(e) { cat("err:", conditionMessage(e), "\n"); NULL }))
for (f in fits) if (!is.null(f))
  cat(sprintf("hand-box fit: (%.3f, %.3f) | auto ub=(%.3f,%.3f) -> exceeds auto box: %s | bd_ML cond=0 = (%.3f,%.3f)\n",
              f$pars[1], f$pars[2], abL$upper_bound[1], abL$upper_bound[2],
              !inside(f$pars[1:2], abL$lower_bound, abL$upper_bound), m0[1], m0[2]))
