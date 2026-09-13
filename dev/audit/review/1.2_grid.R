lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, .libPaths()))
library(emphasis)
d <- readRDS("/tmp/1.2_trees.rds"); trees <- d$trees; w <- d$w; T0 <- d$brts[1]

log_psurv <- function(lam, mu, T) {
  if (abs(lam-mu) < 1e-8) mu <- lam - 1e-8
  r <- lam-mu; p0 <- mu*(1-exp(-r*T))/(lam-mu*exp(-r*T)); 2*log(max(1-p0,1e-300))
}
cond_cr <- function(p) log_psurv(p[1], p[5], T0)

mstep <- function(w, init, lb, ub, mb, link, cond=NULL, nt=1L) {
  es <- list(trees=trees, weights=w, rejected=0L, rejected_overruns=0L,
             rejected_lambda=0L, rejected_zero_weights=0L, time=0, fhat=0)
  r <- emphasis:::m_cpp(e_step=es, init_pars=init, plugin="rpd1", lower_bound=lb,
             upper_bound=ub, xtol_rel=1e-6, num_threads=nt, model=mb,
             link=link, rho=1, rconditional=cond)
  list(est=as.numeric(r$estimates), opt=r$nlopt)
}

models <- list(cr=c(0L,0L,0L), dd=c(1L,0L,0L), d=c(0L,0L,1L), nd=c(1L,0L,1L))
out <- list()
for (mn in names(models)) for (link in 0:2) {
  mb <- models[[mn]]
  lb <- c(1e-3,-1,-1,-1,1e-3,-1,-1,-1); ub <- c(5,1,1,1,5,1,1,1)
  init <- c(1,0.01,0,0.01,0.5,0.01,0,0.01)
  # zero out inactive slots to keep it a legal model
  key <- sprintf("%s_link%d", mn, link)
  a <- tryCatch(mstep(w, init, lb, ub, mb, link), error=function(e) paste("ERR", conditionMessage(e)))
  b <- tryCatch(mstep(w, init, lb, ub, mb, link, cond_cr), error=function(e) paste("ERR", conditionMessage(e)))
  c1 <- tryCatch(mstep(w/137, init, lb, ub, mb, link, cond_cr), error=function(e) paste("ERR", conditionMessage(e)))
  out[[key]] <- list(uncond=a, cond=b, cond_scaled=c1)
}
# boundary / degenerate weight cases, cr only
lb <- c(1e-3,0,0,0,1e-3,0,0,0); ub <- c(5,0,0,0,5,0,0,0); init <- c(1,0,0,0,0.5,0,0,0)
mb <- c(0L,0L,0L)
bd <- list()
bd$allzero_uncond <- tryCatch(mstep(rep(0,length(w)), init, lb, ub, mb, 0L), error=function(e) paste("ERR",conditionMessage(e)))
bd$allzero_cond   <- tryCatch(mstep(rep(0,length(w)), init, lb, ub, mb, 0L, cond_cr), error=function(e) paste("ERR",conditionMessage(e)))
wn <- w; wn[3] <- NaN
bd$nan_w_uncond <- tryCatch(mstep(wn, init, lb, ub, mb, 0L), error=function(e) paste("ERR",conditionMessage(e)))
bd$nan_w_cond   <- tryCatch(mstep(wn, init, lb, ub, mb, 0L, cond_cr), error=function(e) paste("ERR",conditionMessage(e)))
wi <- w; wi[3] <- Inf
bd$inf_w_cond   <- tryCatch(mstep(wi, init, lb, ub, mb, 0L, cond_cr), error=function(e) paste("ERR",conditionMessage(e)))
wneg <- w; wneg[3] <- -w[3]
bd$neg_w_cond   <- tryCatch(mstep(wneg, init, lb, ub, mb, 0L, cond_cr), error=function(e) paste("ERR",conditionMessage(e)))
# single tree, sample_size 1
es1 <- trees[1]
bd$one_tree <- tryCatch({
  es <- list(trees=es1, weights=1, rejected=0L, rejected_overruns=0L,
             rejected_lambda=0L, rejected_zero_weights=0L, time=0, fhat=0)
  r <- emphasis:::m_cpp(e_step=es, init_pars=init, plugin="rpd1", lower_bound=lb, upper_bound=ub,
             xtol_rel=1e-6, num_threads=1L, model=mb, link=0L, rho=1, rconditional=cond_cr)
  list(est=as.numeric(r$estimates), opt=r$nlopt)}, error=function(e) paste("ERR",conditionMessage(e)))
# mu > lambda start, mu == lambda start
bd$mu_gt_lam <- tryCatch(mstep(w, c(0.3,0,0,0,0.9,0,0,0), lb, ub, mb, 0L, cond_cr), error=function(e) paste("ERR",conditionMessage(e)))
bd$mu_eq_lam <- tryCatch(mstep(w, c(0.5,0,0,0,0.5,0,0,0), lb, ub, mb, 0L, cond_cr), error=function(e) paste("ERR",conditionMessage(e)))
bd$near0 <- tryCatch(mstep(w, c(1e-3,0,0,0,1e-3,0,0,0), lb, ub, mb, 0L, cond_cr), error=function(e) paste("ERR",conditionMessage(e)))
# threads
bd$threads4_uncond <- tryCatch(mstep(w, init, lb, ub, mb, 0L, NULL, 4L), error=function(e) paste("ERR",conditionMessage(e)))
bd$threads4_cond   <- tryCatch(mstep(w, init, lb, ub, mb, 0L, cond_cr, 4L), error=function(e) paste("ERR",conditionMessage(e)))
# conditional returning non-finite / NA
bd$cond_inf <- tryCatch(mstep(w, init, lb, ub, mb, 0L, function(p) -Inf), error=function(e) paste("ERR",conditionMessage(e)))
bd$cond_na  <- tryCatch(mstep(w, init, lb, ub, mb, 0L, function(p) NA_real_), error=function(e) paste("ERR",conditionMessage(e)))
# count conditional calls
ncalls <- 0L
cc <- function(p) { ncalls <<- ncalls + 1L; cond_cr(p) }
invisible(tryCatch(mstep(w, init, lb, ub, mb, 0L, cc), error=function(e) NULL))
bd$cond_calls <- ncalls
saveRDS(list(grid=out, bd=bd), commandArgs(TRUE)[2])
cat("done\n")
