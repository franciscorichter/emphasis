## H12 replication — vary what the verifier did not:
##  1. exponential link (log-rate parametrisation), init with mu > lam
##  2. different trees: 30-tip rphylo, 8-tip tree, SHORT crown age (few
##     immigration events -> can the E-step at mu>lam ever "succeed"?)
##  3. MCEM from mu < lam on trees whose DDD::bd_ML(cond=0) MLE has mu > lam
##     (search random small old clades for such a tree)
##  4. sampler exactness check with fix, on a different tree and rates
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages({library(emphasis); library(ape); library(DDD)})
ns <- asNamespace("emphasis")
int_cr <- ns$.bdi_integral_cr

mk <- function(n, seed, age = 5, b = 0.5, d = 0.2) {
  set.seed(seed); tr <- ape::rphylo(n, b, d)
  tr$edge.length <- tr$edge.length / max(ape::branching.times(tr)) * age; tr
}
aug_try <- function(tr, pars, reps = 5, ss = 5L, link = 0L) {
  out <- replicate(reps, tryCatch({
    a <- ns$.augment_tree_bdi(tr, pars = pars, sample_size = ss, link = link)
    nm <- vapply(a$trees, function(d) sum(d$t_ext == 0), 1L)
    sprintf("ok: %d trees, n_missing=%s, sd(lw)=%.1e, fhat=%.4f", length(a$trees),
            paste(nm, collapse = "/"), sd(a$weights), a$fhat)
  }, error = function(e) paste("ERROR:", conditionMessage(e))))
  print(table(out))
}

cat("=== 1. exponential link, log-rates: pars = log(c(0.3, 0.5)) ===\n")
tr15 <- mk(15, 121)
suppressWarnings(aug_try(tr15, log(c(0.3, 0.5)), link = 1L))
cat("control exponential link log(c(0.5,0.3)):\n"); aug_try(tr15, log(c(0.5, 0.3)), link = 1L)
f_exp <- suppressWarnings(estimate_rates(tr15, model = "cr", method = "mcem", link = "exponential",
                         init_pars = log(c(0.3, 0.5)),
                         control = list(lower_bound = log(c(0.01, 0.01)), upper_bound = log(c(5, 5)),
                                        sample_size = 10L, max_iter = 10L, max_time = 60, num_threads = 1L)))
cat("exp-link MCEM from mu>lam: pars", f_exp$pars, " loglik", f_exp$loglik,
    " stop", f_exp$details$stop_reason, " iters", f_exp$details$iterations, "\n")
f_exp2 <- estimate_rates(tr15, model = "cr", method = "mcem", link = "exponential",
                         init_pars = log(c(0.5, 0.3)),
                         control = list(lower_bound = log(c(0.01, 0.01)), upper_bound = log(c(5, 5)),
                                        sample_size = 10L, max_iter = 10L, max_time = 60, num_threads = 1L))
cat("exp-link MCEM from mu<lam: pars", exp(f_exp2$pars), " stop", f_exp2$details$stop_reason, "\n")

cat("\n=== 2. other trees ===\n")
tr30 <- mk(30, 7, age = 8, b = 0.6, d = 0.4)
cat("-- 30 tips, age 8, (0.3,0.5) --\n"); suppressWarnings(aug_try(tr30, c(0.3, 0.5)))
cat("-- 30 tips, (0.5,0.3) control --\n"); aug_try(tr30, c(0.5, 0.3))
tr8 <- mk(8, 3, age = 5)
cat("-- 8 tips, age 5, (0.3,0.5) --\n"); suppressWarnings(aug_try(tr8, c(0.3, 0.5)))
cat("-- 8 tips, age 0.3 (short: few immigration events), (0.3,0.5), ss=3 --\n")
tr8s <- mk(8, 3, age = 0.3)
suppressWarnings(aug_try(tr8s, c(0.3, 0.5), reps = 10, ss = 3L))
cat("-- same short tree, (0.05, 0.10), ss=3 --\n")
suppressWarnings(aug_try(tr8s, c(0.05, 0.10), reps = 10, ss = 3L))
cat("-- boundary: (0.5, 0.4999) and (0.5, 0.5000001) on tr15 --\n")
aug_try(tr15, c(0.5, 0.4999), reps = 3); suppressWarnings(aug_try(tr15, c(0.5, 0.5000001), reps = 3))

cat("\n=== 3. does any small old clade have bd_ML(cond=0) MLE with mu > lam? ===\n")
found <- NULL
for (s in 1:60) {
  set.seed(1000 + s)
  n <- sample(6:12, 1)
  tr <- tryCatch(ape::rphylo(n, 0.3, 0.25), error = function(e) NULL)
  if (is.null(tr)) next
  tr$edge.length <- tr$edge.length / max(ape::branching.times(tr)) * 10
  brts <- sort(ape::branching.times(tr), decreasing = TRUE)
  ml <- tryCatch(suppressWarnings(DDD::bd_ML(brts, cond = 0, btorph = 0, soc = 2, verbose = FALSE,
                                             initparsopt = c(0.3, 0.2))),
                 error = function(e) NULL)
  if (is.null(ml)) next
  cat(sprintf("seed %d n=%d: bd_ML lam=%.4f mu=%.4f ratio=%.3f\n", 1000 + s, n, ml$lambda0, ml$mu0, ml$mu0 / ml$lambda0))
  if (ml$mu0 > ml$lambda0 && ml$lambda0 > 0.01) { found <- list(tr = tr, ml = ml, seed = 1000 + s); break }
}
if (!is.null(found)) {
  tr <- found$tr; ml <- found$ml
  cat(sprintf("FOUND tree seed %d: bd_ML lam=%.4f mu=%.4f\n", found$seed, ml$lambda0, ml$mu0))
  msgs <- character(0)
  fit <- withCallingHandlers(
    suppressWarnings(estimate_rates(tr, model = "cr", method = "mcem", init_pars = c(ml$lambda0, 0.8 * ml$lambda0),
                   control = list(lower_bound = c(0, 0), upper_bound = c(5, 5),
                                  sample_size = 20L, max_iter = 30L, max_time = 120, num_threads = 1L, verbose = TRUE))),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  cat("UNFIXED BDI MCEM from (lam_hat, 0.8 lam_hat): pars", fit$pars, " stop", fit$details$stop_reason,
      " iters", fit$details$iterations, "\n")
  cat("E-step failure messages:", sum(grepl("E-step failed", msgs)), "\n")
  print(round(fit$details$mcem[, 1:2], 4))
  fit_t <- estimate_rates(tr, model = "cr", method = "mcem", init_pars = c(ml$lambda0, 0.8 * ml$lambda0),
                   control = list(lower_bound = c(0, 0), upper_bound = c(5, 5), sampling = "dynamic_fresh",
                                  sample_size = 50L, max_iter = 30L, max_time = 120, num_threads = 1L))
  cat("THINNING MCEM same start: pars", fit_t$pars, " stop", fit_t$details$stop_reason, "\n")
} else cat("no such tree found in 60 tries\n")

cat("\n=== 4. in-session fix on a different tree/rates: exactness + offset ===\n")
fixed <- function(t1, t2, n, k, lam0, mu0, tp) {
  if (t2 - t1 < 1e-15) return(0)
  d <- lam0 - mu0
  if (abs(d) < 1e-15) {
    a1 <- 1 + lam0 * (tp - t1); a2 <- 1 + lam0 * (tp - t2)
    I_lam <- log(a1 / a2); I_mu <- lam0 * (t2 - t1) + log(a2 / a1)
    return((n + 2L * k) * I_lam + n * I_mu)
  }
  E1 <- exp(-d * (tp - t1)); E2 <- exp(-d * (tp - t2))
  I_lam <- mu0 * (t2 - t1) + log((lam0 - mu0 * E2) / (lam0 - mu0 * E1))
  if (n > 0L) {
    oE1 <- 1 - E1; oE2 <- 1 - E2
    if (abs(oE2) < 1e-300) return(Inf)
    I_mu <- lam0 * (t2 - t1) + log(oE1 / oE2)
  } else I_mu <- 0
  (n + 2L * k) * I_lam + n * I_mu
}
environment(fixed) <- ns
assignInNamespace(".bdi_integral_cr", fixed, ns)
brts30 <- sort(ape::branching.times(tr30), decreasing = TRUE)
for (pr in list(c(0.5, 0.3), c(0.3, 0.5), c(0.2, 0.6), c(1, 1.3))) {
  a <- ns$.augment_tree_bdi(tr30, pars = pr, sample_size = 20L)
  bl <- DDD::bd_loglik(pars1 = c(pr, 0, 0), pars2 = c(0, 0, 0, 0, 2), brts = brts30, missnumspec = 0)
  cat(sprintf("fixed (%.1f,%.1f) 30 tips: n=%d range(lw)=[%.6f,%.6f] fhat-bd_loglik=%.6f  (-log(29!)=%.6f)\n",
              pr[1], pr[2], length(a$trees), min(a$weights), max(a$weights), a$fhat - bl, -lgamma(30)))
}
if (!is.null(found)) {
  tr <- found$tr; ml <- found$ml
  fitf <- estimate_rates(tr, model = "cr", method = "mcem", init_pars = c(ml$lambda0, 0.8 * ml$lambda0),
                   control = list(lower_bound = c(0, 0), upper_bound = c(5, 5),
                                  sample_size = 20L, max_iter = 30L, max_time = 120, num_threads = 1L))
  cat("FIXED BDI MCEM same start: pars", fitf$pars, " stop", fitf$details$stop_reason, " (bd_ML", ml$lambda0, ml$mu0, ")\n")
}
