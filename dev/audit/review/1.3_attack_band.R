.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
suppressMessages(library(emphasis))
int_cr <- emphasis:::.bdi_integral_cr; p_cr <- emphasis:::.bdi_p_cr
aug <- emphasis:::.augment_tree_bdi
tp <- 5
brts11 <- c(5,4.745201,4.530461,4.067871,3.622029,2.929002,1.468698,1.386875,1.302139,0.044729)

# cancellation-free reference implementation (expm1 + sign factorisation)
int_ref <- function(t1,t2,n,k,lam,mu,tp){
  if (t2-t1 < 1e-15) return(0)
  d <- lam-mu; tau1 <- tp-t1; tau2 <- tp-t2
  if (d == 0) { a1<-1+lam*tau1; a2<-1+lam*tau2
    I_lam <- lam*(t2-t1) - log(a1/a2)
    I_mu  <- if (n>0) lam*(t2-t1) + log(tau1/tau2) else 0
    return((n+2*k)*I_lam + n*I_mu) }
  if (d > 0) {
    # lam - mu*E = d - mu*expm1(-d tau);  1-E = -expm1(-d tau)
    I_lam <- mu*(t2-t1) + log((d - mu*expm1(-d*tau2))/(d - mu*expm1(-d*tau1)))
    I_mu  <- if (n>0) lam*(t2-t1) + log(expm1(-d*tau1)/expm1(-d*tau2)) else 0
  } else {
    a <- -d  # >0 ; E_i = exp(a tau_i) may overflow -> factor it out
    # lam - mu E_i = -E_i (mu - lam exp(-a tau_i)) ; log ratio = a(tau2-tau1) + log(..)
    num <- mu - lam*exp(-a*tau2); den <- mu - lam*exp(-a*tau1)
    I_lam <- mu*(t2-t1) + a*(tau2-tau1) + log(num/den)
    # 1-E_i = -E_i(1-exp(-a tau_i)) = E_i*expm1(-a tau_i)
    I_mu <- if (n>0) lam*(t2-t1) + a*(tau1-tau2) + log(expm1(-a*tau1)/expm1(-a*tau2)) else 0
  }
  (n+2*k)*I_lam + n*I_mu
}
cat("== which branch at the switch boundary, lam=0.5 ==\n")
for (mu in c(0.5-5e-13, 0.5+5e-13, 0.5-1e-12, 0.5+1e-12)) {
  d <- 0.5-mu; thr <- 1e-12*max(0.5,mu)
  cat(sprintf("mu-0.5=%+.3e  |d|=%.6e thr=%.6e -> %s\n", mu-0.5, abs(d), thr,
      if (abs(d)<=thr) "CRITICAL" else "general"))
}
cat("\n== accuracy of shipped closed form vs cancellation-free ref across the near-critical band ==\n")
lam <- 0.5
for (e in c(-14,-13,-12.5,-12,-11,-10,-9,-8,-7,-6,-4)) for (s in c(1,-1)) {
  d <- s*10^e; mu <- lam-d
  a <- int_cr(1,3,2L,3L,lam,mu,tp); b <- int_ref(1,3,2L,3L,lam,mu,tp)
  cat(sprintf("d=%+9.1e shipped=%.12f ref=%.12f relerr=%.2e %s\n", d,a,b,abs(a-b)/abs(b),
      if (abs(d)<=1e-12*max(lam,mu)) "[crit]" else ""))
}
cat("\n== weight spread (should be 0 for exact CR) across the band, 11-tip tree ==\n")
for (e in c(-13,-12,-11,-10,-9,-8,-7,-6,-4,-2)) for (s in c(1,-1)) {
  d <- s*10^e; mu <- lam - d
  set.seed(11)
  a <- tryCatch(aug(brts11,pars=c(lam,mu),model_bin=c(0L,0L,0L),sample_size=20L,link=0L,rho=1),
                error=function(e) NULL)
  bl <- DDD::bd_loglik(pars1=c(lam,mu,0,0),pars2=c(0,0,1,0,2),brts=brts11,missnumspec=0)
  if (is.null(a)) { cat(sprintf("d=%+9.1e ERROR\n", d)); next }
  cat(sprintf("d=%+9.1e rng(lw)=%.3e  fhat-bd=%+.3e %s\n", d,
      diff(range(a$weights)), a$fhat-bl, if (abs(d)<=1e-12*max(lam,mu)) "[crit]" else ""))
}
cat("\n== overflow zone: shipped vs cancellation-free ref ==\n")
for (dd in c(100,142,200,500,1e4)) {
  lam<-0.1; mu<-lam+dd
  cat(sprintf("mu-lam=%-7g shipped n=1: %-14g ref: %-14g\n", dd,
      int_cr(0,3,1L,2L,lam,mu,tp), int_ref(0,3,1L,2L,lam,mu,tp)))
}
