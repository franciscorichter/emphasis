## H12 — BDI CR sampler in the mu0 > lam0 region.
##
## Claim: for mu0 > lam0, `.bdi_integral_cr` (R/bdi.R:56-60) returns Inf for
## every segment with n > 0 because `1 - E2 < 0 < 1e-300`; uniroot then errors;
## `.mcem_bdi` counts an E-step failure and perturbs toward the box centre, so
## MCEM cannot visit mu >= lam.  Also: last-segment `t_hi = tp - 1e-14`.
##
## Method
##   A. Compare `.bdi_integral_cr` against a numerical integral of the BDI
##      total rate (n+2k)*lam*(1-p) + n*mu/(1-p) for mu<lam (control), mu>lam,
##      mu=lam.  The closed form is mathematically finite for mu>lam.
##   B. Direct call of `.bdi_find_event_time_cr` with mu>lam, n=1.
##   C. `.augment_tree_bdi` / `simulate_tree(method="bdi")` with mu>lam.
##   D. `estimate_rates(method="mcem")` (default sampling = "bdi") with bounds
##      whose centre has mu>lam, and on a pull-of-the-present tree whose bd_ML
##      MLE has mu ~ lam, compared with the thinning sampler and DDD::bd_ML.
##   E. Last-segment guard: H(tp - 1e-14) and P(U > H).
##   F. One-line fix (guard on |1-E2|) applied in-session via
##      assignInNamespace; re-run A-C and check zero-variance weights.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages({library(emphasis); library(ape); library(DDD)})
set.seed(12)
ns <- asNamespace("emphasis")
int_cr  <- ns$.bdi_integral_cr
p_cr    <- ns$.bdi_p_cr
find_ev <- ns$.bdi_find_event_time_cr

cat("=== A. .bdi_integral_cr vs numerical integral ===\n")
num_int <- function(t1, t2, n, k, lam, mu, tp) {
  pfun <- if (lam == mu) function(t) lam * (tp - t) / (1 + lam * (tp - t)) else function(t) p_cr(t, lam, mu, tp)
  rate <- function(t) { p <- pfun(t); (n + 2*k) * lam * (1 - p) + n * mu / (1 - p) }
  integrate(rate, t1, t2, rel.tol = 1e-10)$value
}
tp <- 5
for (case in list(c(0.5, 0.3), c(0.3, 0.5), c(0.4, 0.4), c(0.3, 0.30000001))) {
  lam <- case[1]; mu <- case[2]
  for (n in c(0L, 1L, 3L)) {
    cf <- int_cr(1, 2, n, 2L, lam, mu, tp)
    nm <- num_int(1, 2, n, 2L, lam, mu, tp)
    cat(sprintf("lam=%.8g mu=%.8g n=%d k=2 [1,2]: closed=%s numeric=%.10f  p(1)=%.6f 1-p(1)=%.6f\n",
                lam, mu, n, format(cf, digits = 10), nm,
                p_cr(1, lam, mu, tp), 1 - p_cr(1, lam, mu, tp)))   # note: .bdi_p_cr itself at lam==mu
  }
}

cat("\n=== B. .bdi_find_event_time_cr with mu>lam, n=1 ===\n")
r <- tryCatch(find_ev(1, 0.7, 1L, 2L, 0.3, 0.5, tp, 2), error = function(e) conditionMessage(e))
print(r)
r0 <- tryCatch(find_ev(1, 0.7, 0L, 2L, 0.3, 0.5, tp, 2), error = function(e) conditionMessage(e))
cat("n=0 (no missing lineage yet):", format(r0), "\n")
r1 <- tryCatch(find_ev(1, 0.7, 1L, 2L, 0.5, 0.3, tp, 2), error = function(e) conditionMessage(e))
cat("control mu<lam n=1:", format(r1), "\n")

cat("\n=== C. .augment_tree_bdi / simulate_tree with mu>lam ===\n")
set.seed(121)
tr <- ape::rphylo(15, 0.5, 0.2)       # extant tree, 15 tips
tr$edge.length <- tr$edge.length / max(ape::branching.times(tr)) * 5   # crown age 5
brts <- sort(ape::branching.times(tr), decreasing = TRUE)
cat("crown age", brts[1], " n_tips", Ntip(tr), "\n")
aug_try <- function(pars, reps = 10) {
  out <- replicate(reps, tryCatch({
    a <- ns$.augment_tree_bdi(tr, pars = pars, sample_size = 5L)
    sprintf("ok: %d trees, sd(lw)=%.2e, fhat=%.4f", length(a$trees), sd(a$weights), a$fhat)
  }, error = function(e) paste("ERROR:", conditionMessage(e))))
  print(table(out))
}
cat("-- pars (0.5, 0.3) control --\n"); aug_try(c(0.5, 0.3))
cat("-- pars (0.3, 0.5) mu>lam --\n"); aug_try(c(0.3, 0.5))
cat("-- pars (0.4, 0.4) critical --\n"); aug_try(c(0.4, 0.4))
cat("-- pars (0.4, 0.41) mu slightly > lam --\n"); aug_try(c(0.4, 0.41))

cat("-- simulate_tree(tree, pars=c(0.3,0.5), method='bdi') --\n")
s_bdi <- tryCatch(simulate_tree(tr, pars = c(0.3, 0.5), method = "bdi"),
                  error = function(e) paste("ERROR:", conditionMessage(e)))
str(s_bdi, max.level = 1)
cat("-- simulate_tree(tree, pars=c(0.3,0.5), method='thinning') --\n")
s_th <- tryCatch(simulate_tree(tr, pars = c(0.3, 0.5), method = "thinning"),
                 error = function(e) paste("ERROR:", conditionMessage(e)))
str(s_th, max.level = 1)

cat("\n=== D. estimate_rates CR MCEM (default sampling = bdi) ===\n")
cat("-- D1: box centre has mu>lam: lower=c(0,0), upper=c(1,1.6), init = centre (0.5,0.8) --\n")
t0 <- proc.time()[3]
f1 <- estimate_rates(tr, model = "cr", method = "mcem",
                     control = list(lower_bound = c(0, 0), upper_bound = c(1, 1.6),
                                    sample_size = 20L, max_iter = 20L, max_time = 60,
                                    num_threads = 1L))
cat(sprintf("time %.1fs\n", proc.time()[3] - t0))
cat("pars:", f1$pars, " loglik:", f1$loglik, "\n")
cat("stop_reason:", f1$details$stop_reason, " iterations:", f1$details$iterations, "\n")

cat("-- D2: same box, init_pars = c(0.5, 0.2) (mu<lam), does MCEM cross to mu>lam? --\n")
t0 <- proc.time()[3]
f2 <- estimate_rates(tr, model = "cr", method = "mcem", init_pars = c(0.5, 0.2),
                     control = list(lower_bound = c(0, 0), upper_bound = c(1, 1.6),
                                    sample_size = 20L, max_iter = 25L, max_time = 90,
                                    num_threads = 1L))
cat(sprintf("time %.1fs\n", proc.time()[3] - t0))
cat("pars:", f2$pars, " stop_reason:", f2$details$stop_reason, " iters:", f2$details$iterations, "\n")
print(round(f2$details$mcem[, 1:4], 4))
ml <- DDD::bd_ML(brts, cond = 0, btorph = 0, soc = 2, verbose = FALSE)
cat("DDD::bd_ML (cond=0):", ml$lambda0, ml$mu0, " loglik", ml$loglik, "\n")

cat("-- D3: pull-of-the-present tree (rcoal), bd_ML has mu ~ lam --\n")
set.seed(7)
trc <- ape::rcoal(15); trc$edge.length <- trc$edge.length / max(ape::branching.times(trc)) * 5
brtsc <- sort(ape::branching.times(trc), decreasing = TRUE)
mlc <- DDD::bd_ML(brtsc, cond = 0, btorph = 0, soc = 2, verbose = FALSE)
cat("DDD::bd_ML (cond=0):", mlc$lambda0, mlc$mu0, " loglik", mlc$loglik, "\n")
t0 <- proc.time()[3]
f3 <- estimate_rates(trc, model = "cr", method = "mcem", init_pars = c(1, 0.5),
                     control = list(lower_bound = c(0, 0), upper_bound = c(5, 5),
                                    sample_size = 20L, max_iter = 25L, max_time = 90,
                                    num_threads = 1L))
cat(sprintf("time %.1fs\n", proc.time()[3] - t0))
cat("BDI  pars:", f3$pars, " stop_reason:", f3$details$stop_reason, " iters:", f3$details$iterations, "\n")
print(round(f3$details$mcem[, 1:4], 4))
t0 <- proc.time()[3]
f3t <- estimate_rates(trc, model = "cr", method = "mcem", init_pars = c(1, 0.5),
                      control = list(lower_bound = c(0, 0), upper_bound = c(5, 5),
                                     sampling = "dynamic_fresh",
                                     sample_size = 50L, max_iter = 25L, max_time = 90,
                                     num_threads = 1L))
cat(sprintf("time %.1fs\n", proc.time()[3] - t0))
cat("THIN pars:", f3t$pars, " stop_reason:", f3t$details$stop_reason, " iters:", f3t$details$iterations, "\n")

cat("\n=== E. last-segment guard t_hi = tp - 1e-14 ===\n")
for (n in 1:3) {
  H <- int_cr(4.5, tp - 1e-14, n, 2L, 0.5, 0.3, tp)
  cat(sprintf("n=%d: H(tp-1e-14) = %.3f, P(U>H) = exp(-H) = %.2e; 1-exp(-d*1e-14) = %.6e (d*1e-14 = %.6e)\n",
              n, H, exp(-H), 1 - exp(-0.2 * 1e-14), 0.2e-14))
}

cat("\n=== F. in-session fix: guard on abs(1 - E2) ===\n")
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
    if (abs(oE2) < 1e-300) return(Inf)          # <-- only change
    I_mu <- lam0 * (t2 - t1) + log(oE1 / oE2)
  } else I_mu <- 0
  (n + 2L * k) * I_lam + n * I_mu
}
environment(fixed) <- ns
for (n in c(1L, 3L)) cat(sprintf("fixed lam=0.3 mu=0.5 n=%d: %.10f vs numeric %.10f\n",
                                  n, fixed(1, 2, n, 2L, 0.3, 0.5, tp), num_int(1, 2, n, 2L, 0.3, 0.5, tp)))
assignInNamespace(".bdi_integral_cr", fixed, ns)
cat("-- after fix: pars (0.3, 0.5) --\n"); aug_try(c(0.3, 0.5))
cat("-- after fix: pars (0.4, 0.41) --\n"); aug_try(c(0.4, 0.41))
a <- ns$.augment_tree_bdi(tr, pars = c(0.3, 0.5), sample_size = 30L)
cat(sprintf("fix: n_trees=%d  range(lw)=[%.6f, %.6f]  n_missing per tree: %s\n",
            length(a$trees), min(a$weights), max(a$weights),
            paste(head(vapply(a$trees, function(d) sum(d$t_ext == 0), 1L), 10), collapse = ",")))
bl_lo <- DDD::bd_loglik(pars1 = c(0.3, 0.5, 0, 0), pars2 = c(0, 0, 0, 0, 2), brts = brts, missnumspec = 0)
bl_hi <- DDD::bd_loglik(pars1 = c(0.5, 0.3, 0, 0), pars2 = c(0, 0, 0, 0, 2), brts = brts, missnumspec = 0)
b <- ns$.augment_tree_bdi(tr, pars = c(0.5, 0.3), sample_size = 30L)
cat(sprintf("fhat(0.3,0.5) - bd_loglik = %.6f ;  fhat(0.5,0.3) - bd_loglik = %.6f  (offset should match if fix is consistent)\n",
            a$fhat - bl_lo, b$fhat - bl_hi))
cat("-- after fix: D1 rerun (centre init (0.5,0.8)) --\n")
t0 <- proc.time()[3]
f1b <- estimate_rates(tr, model = "cr", method = "mcem",
                      control = list(lower_bound = c(0, 0), upper_bound = c(1, 1.6),
                                     sample_size = 20L, max_iter = 20L, max_time = 60,
                                     num_threads = 1L))
cat(sprintf("time %.1fs\n", proc.time()[3] - t0))
cat("pars:", f1b$pars, " stop_reason:", f1b$details$stop_reason, " iters:", f1b$details$iterations, "\n")
cat("-- after fix: D3 rerun (rcoal tree) --\n")
t0 <- proc.time()[3]
f3b <- estimate_rates(trc, model = "cr", method = "mcem", init_pars = c(1, 0.5),
                      control = list(lower_bound = c(0, 0), upper_bound = c(5, 5),
                                     sample_size = 20L, max_iter = 25L, max_time = 90,
                                     num_threads = 1L))
cat(sprintf("time %.1fs\n", proc.time()[3] - t0))
cat("BDI-fixed pars:", f3b$pars, " stop_reason:", f3b$details$stop_reason, " iters:", f3b$details$iterations, "\n")
print(round(f3b$details$mcem[, 1:4], 4))
cat("bd_ML:", mlc$lambda0, mlc$mu0, "\n")
