lib <- commandArgs(TRUE)[1]; .libPaths(c(lib, .libPaths())); library(emphasis)
d <- readRDS("/tmp/1.2_trees.rds"); trees <- d$trees; w <- d$w
lb <- c(1e-3,0,0,0,1e-3,0,0,0); ub <- c(5,0,0,0,5,0,0,0); init <- c(1,0,0,0,0.5,0,0,0)
ms <- function(cond, w) { es <- list(trees=trees, weights=w, rejected=0L, rejected_overruns=0L,
    rejected_lambda=0L, rejected_zero_weights=0L, time=0, fhat=0)
  r <- emphasis:::m_cpp(e_step=es, init_pars=init, plugin="rpd1", lower_bound=lb, upper_bound=ub,
    xtol_rel=1e-6, num_threads=1L, model=c(0L,0L,0L), link=0L, rho=1, rconditional=cond)
  sprintf("(%.6g, %.6g)", r$estimates[1], r$estimates[5]) }
# GAM-shaped conditional: binomial response, survival collapses for small lambda
gamlike <- function(p) log(max(stats::plogis(6*(p[1]-0.8) - 2*p[5]), 1e-300))
cat("gam-like cond, sum_w=40 :", ms(gamlike, w), "\n")
cat("gam-like cond, sum_w=400:", ms(gamlike, w*10), "\n")
cat("gam-like cond, sum_w=4  :", ms(gamlike, w/10), "\n")
