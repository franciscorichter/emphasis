lib <- commandArgs(TRUE)[1]; .libPaths(c(lib, .libPaths())); library(emphasis)
d <- readRDS("/tmp/1.2_trees.rds"); trees <- d$trees; w <- d$w; T0 <- d$brts[1]
lb <- c(1e-3,0,0,0,1e-3,0,0,0); ub <- c(5,0,0,0,5,0,0,0); init <- c(1,0,0,0,0.5,0,0,0)
ms <- function(cond) { es <- list(trees=trees, weights=w, rejected=0L, rejected_overruns=0L,
    rejected_lambda=0L, rejected_zero_weights=0L, time=0, fhat=0)
  r <- emphasis:::m_cpp(e_step=es, init_pars=init, plugin="rpd1", lower_bound=lb, upper_bound=ub,
    xtol_rel=1e-6, num_threads=1L, model=c(0L,0L,0L), link=0L, rho=1, rconditional=cond)
  sprintf("est=(%.7g, %.7g) frozen_at_init=%s nlopt=%d",
          r$estimates[1], r$estimates[5], identical(as.numeric(r$estimates), init), r$nlopt) }
lp <- function(lam,mu,T){ if(abs(lam-mu)<1e-8) mu<-lam-1e-8; r<-lam-mu
  p0<-mu*(1-exp(-r*T))/(lam-mu*exp(-r*T)); 2*log(max(1-p0,1e-300)) }
cat("good cond      :", ms(function(p) lp(p[1],p[5],T0)), "\n")
cat("cond -> NA     :", ms(function(p) NA_real_), "\n")
cat("cond -> NaN    :", ms(function(p) NaN), "\n")
cat("cond -> -Inf   :", ms(function(p) -Inf), "\n")
cat("cond NA in box :", ms(function(p) if (p[1] > 1.05) NA_real_ else lp(p[1],p[5],T0)), "\n")
# amplification of the GAM clamp log(max(p,1e-300)) = -690.7755
Q <- function(th) sum(w * emphasis:::eval_logf(th, trees, model=c(0L,0L,0L), link=0L, rho=1)$logf)
th <- c(1.039945,0,0,0,0.5282315,0,0,0)
cat(sprintf("sum_w=%.1f  Q(mle)=%.1f  -Q=%.1f   penalty_at_clamp: PRE=%.1f POST=%.1f\n",
    sum(w), Q(th), -Q(th), -690.7755, sum(w)*-690.7755))
