args <- commandArgs(TRUE); lib <- if (length(args) && args[1] == "pre") "rlib" else "rlib-wave1"
.libPaths(c(file.path("/Users/pancho/.claude/jobs/867af780/tmp", lib), .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(DDD)})
cat("BUILD:", lib, "\n")
cnt <- function(reset=FALSE) if (exists("thinning_envelope_violations", asNamespace("emphasis"))) emphasis:::thinning_envelope_violations(reset) else NA
aug <- function(brts, pars, n, model=c(0L,0L,0L), link=0L, rho=1, threads=1L, maxN=50L*n)
  emphasis:::augment_trees(brts=brts, pars=pars, sample_size=n, maxN=maxN, max_missing=500L,
                           max_lambda=1e6, num_threads=threads, model=model, link=link, rho=rho)
fhat <- function(raw) { lw <- raw$logf - raw$logg; m <- max(lw); w <- exp(lw-m)
  c(fhat = log(mean(w)) + m - log(1 + raw$rejected_zero_weights/length(w)), se = sd(w)/(mean(w)*sqrt(length(w)))) }
ismiss <- function(df) !(df$t_ext == 0 | df$t_ext == 1e11 | df$t_ext == 5e10)
haz <- function(brts, lam, mu, rho=1) { T <- brts[1]; s <- sort(T - brts[-1]); knots <- c(0, s, T); nvec <- 2 + seq_along(knots) - 1; h <- 0
  for (i in seq_len(length(knots)-1)) { a <- knots[i]; b <- knots[i+1]; mue <- max(mu, 1e-10)
    h <- h + nvec[i]*lam*((b-a) - rho*(1/mue)*(exp(-mue*(T-b)) - exp(-mue*(T-a)))) }; h }
brts7 <- c(6, 4.5, 3.0, 2.0, 1.2, 0.5)
cr_pars <- function(lam, mu, link) switch(as.character(link),
  "0" = c(lam,0,0,0,mu,0,0,0), "1" = c(log(lam),0,0,0,log(mu),0,0,0), "2" = c(lam/exp(-0.5),0,0,0,mu/exp(-0.5),0,0,0))
ref <- function(brts, lam, mu) DDD::bd_loglik(pars1=c(lam,mu,0,0), pars2=c(0,0,1,0,2), brts=brts, missnumspec=0)
report <- function(label, brts, lam, mu, link=0L, n=4000L, rho=1, threads=1L, maxN=50L*n) {
  cnt(TRUE); t0 <- Sys.time()
  raw <- aug(brts, cr_pars(lam,mu,link), n, link=link, rho=rho, threads=threads, maxN=maxN)
  f <- fhat(raw); nm <- vapply(raw$trees, function(d) sum(ismiss(d)), 1); p0 <- mean(nm==0); se0 <- sqrt(p0*(1-p0)/n)
  p0c <- exp(-haz(brts, lam, mu, rho)); r <- if (rho==1) ref(brts, lam, mu) else NA
  cat(sprintf("%-28s link=%d rho=%.2f thr=%d | ntrees=%d rej(ov/lam/zw)=%d/%d/%d | P0 %.4f vs %.4f z=%+.1f | fhat %.3f se %.3f ref %.3f gap %+.3f (z %+.1f) | mean#miss %.2f | viol=%s | %.1fs\n",
    label, link, rho, threads, length(raw$trees), raw$rejected_overruns, raw$rejected_lambda, raw$rejected_zero_weights,
    p0, p0c, (p0-p0c)/se0, f["fhat"], f["se"], r, f["fhat"]-r, (f["fhat"]-r)/f["se"], mean(nm), cnt(), as.numeric(Sys.time()-t0, units="secs")))
  invisible(raw)
}
for (lk in 0:2) report("7tip (0.4,0.2)", brts7, 0.4, 0.2, link=lk)
report("7tip mu=lam (0.3,0.3)", brts7, 0.3, 0.3)
report("7tip mu>lam (0.2,0.4)", brts7, 0.2, 0.4)
report("7tip mu=0 (0.4,0)", brts7, 0.4, 0)
report("7tip mu=0 exp-link", brts7, 0.4, 1e-12, link=1L)
report("7tip tiny rates (1e-3,1e-3)", brts7, 1e-3, 1e-3)
report("7tip big rates (2,1.5)", brts7, 2, 1.5, n=2000L)
report("2tip brts=5 (0.4,0.2)", c(5), 0.4, 0.2)
report("2tip mu>lam", c(5), 0.2, 0.5)
report("7tip rho=0.5", brts7, 0.4, 0.2, rho=0.5)
report("7tip rho=0.2", brts7, 0.4, 0.2, rho=0.2)
report("7tip threads=4", brts7, 0.4, 0.2, threads=4L, n=8000L)
set.seed(11); tr <- ape::rphylo(80, 0.5, 0.2); b80 <- sort(ape::branching.times(tr), decreasing=TRUE)
cat("80-tip crown age", b80[1], "\n")
report("80tip (0.5,0.2)", b80, 0.5, 0.2, n=1500L)
report("80tip (0.5,0.45)", b80, 0.5, 0.45, n=1500L)
# sample_size 1, maxN = sample_size
cnt(TRUE); r1 <- aug(brts7, cr_pars(0.4,0.2,0), 1L, maxN=1L); cat("sample_size=1 maxN=1: ntrees", length(r1$trees), "rej", r1$rejected_overruns, r1$rejected_lambda, r1$rejected_zero_weights, "viol", cnt(), "\n")
cnt(TRUE); r1 <- aug(brts7, cr_pars(0.4,0.2,0), 50L, maxN=50L); cat("sample_size=50 maxN=50: ntrees", length(r1$trees), "rej", r1$rejected_overruns, r1$rejected_lambda, r1$rejected_zero_weights, "viol", cnt(), "\n")
# dd model vs DDD::dd_loglik: gap should be constant across the grid
cat("\n-- dd model (K = 40): fhat - dd_loglik(cond=0,btorph=1,soc=2)\n")
for (g in list(c(0.5,0.1), c(0.5,0.3), c(0.8,0.2), c(0.4,0.35))) {
  lam <- g[1]; mu <- g[2]; K <- 40; bN <- -lam/K
  cnt(TRUE); raw <- aug(brts7, c(lam, bN, 0,0, mu,0,0,0), 4000L, model=c(1L,0L,0L))
  f <- fhat(raw); r <- DDD::dd_loglik(pars1=c(lam,mu,K), pars2=c(100,1,0,1,0,2), brts=brts7, missnumspec=0)
  cat(sprintf("  dd (%.2f,%.2f,K=%d): fhat %.3f se %.3f  dd_loglik %.3f  gap %+.3f  rej %d/%d/%d viol=%s\n", lam, mu, K, f["fhat"], f["se"], r, f["fhat"]-r, raw$rejected_overruns, raw$rejected_lambda, raw$rejected_zero_weights, cnt()))
}
