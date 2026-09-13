lib <- commandArgs(TRUE)[1]; .libPaths(c(lib, .libPaths())); library(emphasis)
mb <- 7 + (1:6)/10; me <- 8 + (1:6)/10
tree_bad <- data.frame(brts=c(1:7, mb, me, 8.8, 9), n=c(2:8, 8:13, 14:9, 8, 9),
  t_ext=c(rep(1e11,7), me, rep(0,6), 1e11, 1e11), pd=0,
  tip_start=c(rep(0,7), mb, mb, 0, 0), id=c(0:6, 8:13, 8:13, 7L, -1L),
  parent_id=c(rep(-1L,7), rep(0L,12), -1L, -1L))
dd <- c(1L,0L,0L)
pars8 <- c(1.5,-0.12,0,0,0.4,0,0,0)
cat("logf(tree_bad) =", emphasis:::eval_logf(pars8, list(tree_bad), model=dd, link=0L, rho=1)$logf, "\n")
ms <- function(trees, w, init, lb, ub, cond=NULL) {
  es <- list(trees=trees, weights=w, rejected=0L, rejected_overruns=0L,
             rejected_lambda=0L, rejected_zero_weights=0L, time=0, fhat=0)
  r <- emphasis:::m_cpp(e_step=es, init_pars=init, plugin="rpd1", lower_bound=lb,
        upper_bound=ub, xtol_rel=1e-6, num_threads=1L, model=dd, link=0L, rho=1,
        rconditional=cond); list(est=as.numeric(r$estimates), opt=r$nlopt) }
# box where lambda(13) = b0 + 13*bN <= 1 - 1.3 < 0 everywhere -> logf = -Inf at every point
lb <- c(0.01,-1,0,0,0.001,0,0,0); ub <- c(1,-0.1,0,0,2,0,0,0)
init <- c(0.8,-0.5,0,0,0.4,0,0,0)
lf <- sapply(list(init, lb+1e-6, ub-1e-6, c(1,-0.1,0,0,2,0,0,0)),
  function(p) emphasis:::eval_logf(p, list(tree_bad), model=dd, link=0L, rho=1)$logf)
cat("logf at 4 box points:", lf, "\n")
r <- ms(list(tree_bad), 1, init, lb, ub)
cat("ALL-INFEASIBLE  est =", signif(r$est[c(1,2,5)],8), " nlopt =", r$opt,
    " moved =", !identical(r$est, init), "\n")
cat("logf at returned est =", emphasis:::eval_logf(r$est, list(tree_bad), model=dd, link=0L, rho=1)$logf, "\n")
