.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
suppressMessages(library(emphasis))
p_cr   <- emphasis:::.bdi_p_cr
int_cr <- emphasis:::.bdi_integral_cr
tp <- 5

num_int <- function(t1,t2,n,k,lam,mu,tp,rel.tol=1e-12){
  rate <- function(t){ p <- vapply(t,function(s) p_cr(s,lam,mu,tp),numeric(1))
    (n+2*k)*lam*(1-p) + n*mu/(1-p) }
  tryCatch(stats::integrate(rate,t1,t2,rel.tol=rel.tol)$value, error=function(e) NA_real_)
}
# independent high-accuracy p using expm1 (no cancellation for small d)
p_ref <- function(t,lam,mu,tp){ d<-lam-mu; tau<-tp-t
  # d/(lam - mu*exp(-d*tau)) = d/(lam-mu + mu*(1-exp(-d tau))) = d/(d - mu*expm1(-d*tau))
  d/(d - mu*expm1(-d*tau)) }

cat("=== A: p_cr vs cancellation-free reference across d ===\n")
lam <- 0.5
for (d in c(1e-1,1e-3,1e-6,1e-9,1e-11,1e-12,5e-13,1e-13,1e-15,0,-1e-13,-1e-11,-1e-6,-1e-1)) {
  mu <- lam - d
  a <- p_cr(1, lam, mu, tp); b <- p_ref(1, lam, mu, tp)
  cat(sprintf("d=%+10.1e  p_cr=%.16g  ref=%.16g  relerr=%.2e  branch=%s\n",
      d, a, b, abs(a-b)/abs(b), if (abs(d) <= 1e-12*max(abs(lam),abs(mu))) "crit" else "gen"))
}

cat("\n=== B: integral closed-form vs integrate(), grid incl mu>>lam ===\n")
grid <- list(c(0.5,0.3),c(0.3,0.5),c(0.4,0.4),c(0.01,2),c(2,0.01),c(1e-4,1e-4),
             c(5,5),c(0,0.7),c(0.7,0),c(0,0),c(20,25),c(1e-3,50))
for (g in grid) for (n in c(0L,1L,3L)) {
  cf <- tryCatch(int_cr(0.5,3,n,2L,g[1],g[2],tp), error=function(e) paste("ERR",conditionMessage(e)))
  nm <- num_int(0.5,3,n,2L,g[1],g[2],tp)
  cat(sprintf("lam=%-8g mu=%-8g n=%d  closed=%-20.12g numeric=%-20.12g  diff=%.2e\n",
      g[1],g[2],n, suppressWarnings(as.numeric(cf)), nm,
      suppressWarnings(abs(as.numeric(cf)-nm))))
}

cat("\n=== C: overflow for mu >> lam (|d|*(tp-t1) > 709) ===\n")
for (dd in c(100,140,142,200,500)) {
  lam <- 0.1; mu <- lam + dd
  v0 <- int_cr(0, 3, 0L, 2L, lam, mu, tp)
  v1 <- int_cr(0, 3, 1L, 2L, lam, mu, tp)
  v2 <- int_cr(0, 4.999, 1L, 2L, lam, mu, tp)
  cat(sprintf("mu-lam=%6g  |d|*(tp-t1)=%7.1f  n=0:%-14g n=1:%-14g n=1(t2~tp):%-14g\n",
              dd, dd*(tp-0), v0, v1, v2))
}
cat("\n  same with a deep-time tree (tp=100, moderate rates):\n")
for (dd in c(5,7,7.1,10)) {
  lam<-0.1; mu<-lam+dd
  v <- int_cr(0, 50, 1L, 2L, lam, mu, 100)
  cat(sprintf("  tp=100 mu-lam=%g |d|*tp=%.0f -> %g\n", dd, dd*100, v))
}

cat("\n=== D: accuracy at the switch boundary d = 1e-12*lam (largest d using the limit) ===\n")
for (lam in c(0.01, 0.5, 5, 100)) {
  d <- 1e-12*lam; mu <- lam - d
  a <- p_cr(1, lam, mu, tp); b <- p_ref(1, lam, mu, tp)
  ci <- int_cr(1,3,2L,3L,lam,mu,tp); ni <- num_int(1,3,2L,3L,lam,mu,tp,rel.tol=1e-10)
  cat(sprintf("lam=%-7g d=%8.1e  p relerr=%.2e   int closed=%.12g vs numeric=%.12g (rel %.2e)\n",
      lam,d,abs(a-b)/abs(b), ci, ni, abs(ci-ni)/abs(ni)))
}
cat("\n=== E: t2 == tp exactly and beyond; n=0 vs n>0, both sides ===\n")
for (g in list(c(0.5,0.3),c(0.3,0.5),c(0.4,0.4),c(0,0))) {
  cat(sprintf("lam=%g mu=%g  n=1,t2=tp -> %-8g | n=0,t2=tp -> %-12g | t2=tp+1 -> %g\n",
    g[1],g[2], int_cr(1,tp,1L,2L,g[1],g[2],tp), int_cr(1,tp,0L,2L,g[1],g[2],tp),
    suppressWarnings(int_cr(1,tp+1,1L,2L,g[1],g[2],tp))))
}
