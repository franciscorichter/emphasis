# Item 1.3 review: numerical attack on .bdi_p_cr / .bdi_integral_cr (POST-FIX)
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
suppressWarnings(library(emphasis))
p_cr   <- emphasis:::.bdi_p_cr
int_cr <- emphasis:::.bdi_integral_cr

# Independent reference: p and rates written from scratch (no package code).
p_ref <- function(t, lam, mu, tp) {
  tau <- tp - t
  if (lam == mu) return(1 / (1 + lam * tau))
  d <- lam - mu
  if (d > 0) d / (d - mu * expm1(-d * tau))          # stable near d = 0
  else { e <- exp(d * tau); (d * e) / (lam * e - mu) } # stable for large |d| tau
}
rate_ref <- function(t, n, k, lam, mu, tp) {
  p <- vapply(t, p_ref, numeric(1), lam = lam, mu = mu, tp = tp)
  (n + 2 * k) * lam * (1 - p) + n * mu / (1 - p)
}
num_ref <- function(t1, t2, n, k, lam, mu, tp)
  stats::integrate(rate_ref, t1, t2, n = n, k = k, lam = lam, mu = mu, tp = tp,
                   rel.tol = 1e-12, subdivisions = 2000L)$value

cat("== A. wide grid closed form vs independent integrate ==\n")
set.seed(1)
grid <- expand.grid(lam = c(1e-4, 0.05, 0.5, 3, 20),
                    mu  = c(0, 1e-4, 0.05, 0.5, 3, 20),
                    tp  = c(0.5, 5, 60))
worst <- 0; worst_case <- NULL; nan_cases <- list()
for (i in seq_len(nrow(grid))) {
  lam <- grid$lam[i]; mu <- grid$mu[i]; tp <- grid$tp[i]
  for (seg in list(c(0, tp * 0.3), c(tp * 0.4, tp * 0.98), c(tp * 0.9, tp * 0.999)))
    for (n in c(0L, 1L, 7L)) for (k in c(2L, 40L)) {
      cf <- int_cr(seg[1], seg[2], n, k, lam, mu, tp)
      if (!is.finite(cf)) { nan_cases[[length(nan_cases) + 1]] <- c(lam, mu, tp, seg, n, k, cf); next }
      nm <- tryCatch(num_ref(seg[1], seg[2], n, k, lam, mu, tp), error = function(e) NA)
      if (is.na(nm)) next
      err <- abs(cf - nm) / max(1, abs(nm))
      if (err > worst) { worst <- err; worst_case <- c(lam = lam, mu = mu, tp = tp, t1 = seg[1], t2 = seg[2], n = n, k = k, cf = cf, nm = nm) }
    }
}
cat("worst rel err:", worst, "\n"); print(worst_case)
cat("non-finite closed forms:", length(nan_cases), "\n")
if (length(nan_cases)) print(do.call(rbind, nan_cases))

cat("\n== B. accuracy band around the relative switch (lam = 0.5, tp = 5, seg [1,3], n=2,k=3) ==\n")
# Stable closed-form reference via expm1 (exact in the whole band, no integrate noise):
#   lam - mu*E = d - mu*expm1(-d*tau);  1 - E = -expm1(-d*tau)
ref_stable <- function(t1, t2, n, k, lam, mu, tp) {
  d <- lam - mu; dt <- t2 - t1
  g1 <- expm1(-d * (tp - t1)); g2 <- expm1(-d * (tp - t2))
  I_lam <- mu * dt + log((d - mu * g2) / (d - mu * g1))
  I_mu  <- lam * dt + log(g1 / g2)
  (n + 2 * k) * I_lam + n * I_mu
}
crit <- function(t1, t2, n, k, lam, tp) {
  a1 <- 1 + lam * (tp - t1); a2 <- 1 + lam * (tp - t2)
  (n + 2 * k) * (lam * (t2 - t1) - log(a1 / a2)) + n * (lam * (t2 - t1) + log((tp - t1) / (tp - t2)))
}
lam <- 0.5; tp <- 5
cat(sprintf("critical limit = %.12f ; expm1 ref at d=0 -> NaN expected: %s\n", crit(1,3,2L,3L,lam,tp), ref_stable(1,3,2L,3L,lam,lam,tp)))
for (rd in c(1e-15, 5e-13, 2e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-4)) for (s in c(1, -1)) {
  mu <- lam - s * rd * lam
  v  <- int_cr(1, 3, 2L, 3L, lam, mu, tp)
  nm <- ref_stable(1, 3, 2L, 3L, lam, mu, tp)
  cat(sprintf("d/lam=%+.0e  int_cr=%.12f  expm1_ref=%.12f  |diff|=%.2e  branch=%s\n",
              s * rd, v, nm, abs(v - nm),
              if (abs(lam - mu) <= 1e-12 * max(lam, mu)) "critical" else "general"))
}
cat("-- same band, long tree tp = 200, seg [10, 190], n = 3, k = 20, lam = 0.05 --\n")
lam <- 0.05; tp <- 200
for (rd in c(5e-13, 2e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-6)) for (s in c(1, -1)) {
  mu <- lam - s * rd * lam
  v  <- int_cr(10, 190, 3L, 20L, lam, mu, tp)
  nm <- ref_stable(10, 190, 3L, 20L, lam, mu, tp)
  cat(sprintf("d/lam=%+.0e  int_cr=%.10f  expm1_ref=%.10f  |diff|=%.2e\n", s * rd, v, nm, abs(v - nm)))
}
cat("critical-branch value there:", crit(10, 190, 3L, 20L, lam, tp), "\n")

cat("\n== C. p_cr near the switch ==\n")
for (rd in c(1e-13, 1e-12, 2e-12, 1e-11, 1e-10)) for (s in c(1, -1)) {
  mu <- lam - s * rd * lam
  cat(sprintf("d/lam=%+.0e  p_cr=%.15f  p_ref=%.15f  diff=%.1e\n", s * rd,
              p_cr(1, lam, mu, tp), p_ref(1, lam, mu, tp), p_cr(1, lam, mu, tp) - p_ref(1, lam, mu, tp)))
}

cat("\n== D. overflow probing: (mu - lam) * (tp - t) large (mu > lam) ==\n")
for (cs in list(c(0.1, 8, 100), c(0.5, 15, 50), c(0.5, 150, 5), c(0.5, 200, 5), c(1, 0, 5), c(0, 1, 5), c(0, 0, 5), c(1e-300, 0, 5), c(0, 1e-300, 5))) {
  lam <- cs[1]; mu <- cs[2]; tp <- cs[3]
  cat(sprintf("lam=%g mu=%g tp=%g: p(0)=%.6g p(tp*0.5)=%.6g  int[0,0.5tp] n=1 k=2 = %.6g   n=0: %.6g   int[0.5tp,0.9tp] n=1: %.6g\n",
      lam, mu, tp, p_cr(0, lam, mu, tp), p_cr(tp/2, lam, mu, tp),
      int_cr(0, tp/2, 1L, 2L, lam, mu, tp), int_cr(0, tp/2, 0L, 2L, lam, mu, tp),
      int_cr(tp/2, 0.9*tp, 1L, 2L, lam, mu, tp)))
}
cat("  reference for lam=0.1 mu=8 tp=100 int[0,50] n=1 k=2:", num_ref(0, 50, 1L, 2L, 0.1, 8, 100), "\n")
cat("  reference for lam=0.5 mu=150 tp=5 int[0,2.5] n=1 k=2:", num_ref(0, 2.5, 1L, 2L, 0.5, 150, 5), "\n")
cat("  reference for lam=0.5 mu=200 tp=5 int[0,2.5] n=1 k=2:", num_ref(0, 2.5, 1L, 2L, 0.5, 200, 5), "\n")

cat("\n== E. find_event_time at overflow point ==\n")
fe <- emphasis:::.bdi_find_event_time_cr
set.seed(2)
r <- tryCatch(fe(0, rexp(1), 1L, 2L, 0.5, 200, 5, 2.5), error = function(e) conditionMessage(e))
print(r)
r <- tryCatch(fe(0, rexp(1), 1L, 2L, 0.1, 8, 100, 30), error = function(e) conditionMessage(e))
print(r)

cat("\n== F. vector t and t == tp in p_cr ==\n")
print(p_cr(c(0, 2.5, 5), 0.5, 0.3, 5)); print(p_cr(c(0, 2.5, 5), 0.3, 0.5, 5)); print(p_cr(c(0, 2.5, 5), 0.4, 0.4, 5))

cat("\n== G. negative rates (linear link clamps, but call directly) ==\n")
print(c(p_cr(1, -0.2, 0.3, 5), int_cr(1, 2, 1L, 2L, -0.2, 0.3, 5)))
