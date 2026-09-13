lib <- commandArgs(TRUE)[1]; .libPaths(c(lib, .libPaths())); library(emphasis)
d <- readRDS("/tmp/1.2_trees.rds"); trees <- d$trees; w <- d$w; T0 <- d$brts[1]
log_psurv <- function(lam, mu, T) { if (abs(lam-mu)<1e-8) mu <- lam-1e-8
  r<-lam-mu; p0<-mu*(1-exp(-r*T))/(lam-mu*exp(-r*T)); 2*log(max(1-p0,1e-300)) }
cond_cr <- function(p) log_psurv(p[1], p[5], T0)
ms <- function(w, init, lb, ub, mb, link, cond=NULL, nt=1L) {
  es <- list(trees=trees, weights=w, rejected=0L, rejected_overruns=0L,
             rejected_lambda=0L, rejected_zero_weights=0L, time=0, fhat=0)
  r <- emphasis:::m_cpp(e_step=es, init_pars=init, plugin="rpd1", lower_bound=lb,
        upper_bound=ub, xtol_rel=1e-6, num_threads=nt, model=mb, link=link, rho=1,
        rconditional=cond); c(as.numeric(r$estimates), opt=r$nlopt) }
res <- list()
# CR layout: only slots 1 and 5 free
lb <- c(1e-3,0,0,0,1e-3,0,0,0); ub <- c(5,0,0,0,5,0,0,0); init <- c(1,0,0,0,0.5,0,0,0)
for (link in 0:2) for (cn in c("null","cond")) {
  cd <- if (cn=="cond") cond_cr else NULL
  base <- ms(w, init, lb, ub, c(0L,0L,0L), link, cd)
  sc <- sapply(c(1/1000,1/137,1/40,1,7,1000), function(k) ms(w*k, init, lb, ub, c(0L,0L,0L), link, cd)[c(1,5)])
  res[[sprintf("cr_link%d_%s",link,cn)]] <- list(base=base[c(1,5)], dev=max(apply(sc,1,function(z) diff(range(z)))), sc=sc)
}
# DD layout: slots 1,2,5 free
lbd <- c(1e-3,-0.5,0,0,1e-3,0,0,0); ubd <- c(5,0,0,0,5,0,0,0); initd <- c(1,-0.05,0,0,0.5,0,0,0)
for (link in 0:1) for (cn in c("null","cond")) {
  cd <- if (cn=="cond") cond_cr else NULL
  sc <- sapply(c(1/137,1,137), function(k) ms(w*k, initd, lbd, ubd, c(1L,0L,0L), link, cd)[c(1,2,5)])
  res[[sprintf("dd_link%d_%s",link,cn)]] <- list(base=sc[,2], dev=max(apply(sc,1,function(z) diff(range(z)))), sc=sc)
}
for (k in names(res)) cat(sprintf("%-18s dev=%-12.4g base=%s\n", k, res[[k]]$dev, paste(signif(res[[k]]$base,7),collapse=", ")))
