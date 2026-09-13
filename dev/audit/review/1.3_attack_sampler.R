lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, .libPaths()))
suppressMessages(library(emphasis))
aug <- emphasis:::.augment_tree_bdi
tp <- 5
brts11 <- c(5,4.745201,4.530461,4.067871,3.622029,2.929002,1.468698,1.386875,1.302139,0.044729)
set.seed(1); brts30 <- sort(runif(29,0,5),decreasing=TRUE); brts30[1] <- 5
brts3  <- c(5, 2.5)

run <- function(brts, pars, link=0L, N=30L, seed=11, model=c(0L,0L,0L), rho=1, maxmiss=1e4L) {
  set.seed(seed)
  a <- tryCatch(aug(brts, pars=pars, model_bin=model, sample_size=N, link=link, rho=rho,
                    max_missing=maxmiss),
                error=function(e) paste("ERROR:", conditionMessage(e)),
                warning=function(w) paste("WARN:", conditionMessage(w)))
  if (is.character(a)) return(list(msg=a))
  list(msg=NULL, ntree=length(a$trees), rng=if(length(a$weights)) diff(range(a$weights)) else NA,
       fhat=a$fhat, w1=if(length(a$weights)) a$weights[1] else NA)
}
bd <- function(brts,pars) tryCatch(DDD::bd_loglik(pars1=c(pars,0,0),pars2=c(0,0,1,0,2),brts=brts,missnumspec=0), error=function(e) NA)

cat("== linear link, both sides, several trees ==\n")
cases <- list(list("11tip",brts11), list("30tip",brts30), list("3tip",brts3))
pl <- list(c(0.5,0.3),c(0.3,0.5),c(0.4,0.4),c(0.5,0.5-5e-13),c(0.5,0.5+5e-13),
           c(0.5,0.5-1e-11),c(2,2.5),c(0.01,0.02),c(0,0),c(0.3,0),c(0,0.3),c(1e-6,1e-6))
for (cs in cases) for (p in pl) {
  r <- run(cs[[2]], p)
  ref <- suppressWarnings(bd(cs[[2]], p))
  cat(sprintf("%-6s lam=%-10g mu=%-14g %s\n", cs[[1]], p[1], p[2],
    if (!is.null(r$msg)) r$msg else sprintf("ntree=%2d rng(lw)=%.2e fhat=%-14.8g bd=%-14.8g diff=%.2e",
      r$ntree, r$rng, r$fhat, ref, abs(r$fhat-ref))))
}
cat("\n== exponential link (link=1), mu>lam, wide range ==\n")
for (p in list(log(c(0.3,0.5)), log(c(0.5,0.3)), log(c(0.4,0.4)), c(-2,-1.5), c(1,1.2), c(2,2.2))) {
  r <- run(brts11, p, link=1L); ref <- suppressWarnings(bd(brts11, exp(p)))
  cat(sprintf("pars=(%g,%g)->lam=%g mu=%g  %s\n", p[1],p[2],exp(p[1]),exp(p[2]),
    if (!is.null(r$msg)) r$msg else sprintf("ntree=%2d rng=%.2e fhat=%.8g bd=%.8g diff=%.2e",
      r$ntree,r$rng,r$fhat,ref,abs(r$fhat-ref))))
}
cat("\n== exponential link, LARGE mu (overflow zone |d|*tp>709) ==\n")
for (p in list(c(0,5), c(0,5.5), c(0,6), c(log(0.1), log(200)))) {
  r <- run(brts11, p, link=1L, N=5L)
  cat(sprintf("pars=(%g,%g) lam=%.4g mu=%.4g |d|*tp=%.0f  %s\n", p[1],p[2],exp(p[1]),exp(p[2]),
      abs(exp(p[1])-exp(p[2]))*tp,
      if (!is.null(r$msg)) r$msg else sprintf("ntree=%d rng=%.2e fhat=%g", r$ntree,r$rng,r$fhat)))
}
cat("\n== sample_size = 1, and dd model (unaffected path) ==\n")
r <- run(brts11, c(0.3,0.5), N=1L); cat("N=1 mu>lam:", if(!is.null(r$msg)) r$msg else sprintf("ntree=%d fhat=%g",r$ntree,r$fhat), "\n")
r <- run(brts11, c(0.3,0.5,0), model=c(1L,0L,0L), N=10L)
cat("dd (model 100) pars=(0.3,0.5,0):", if(!is.null(r$msg)) r$msg else sprintf("ntree=%d rng=%.3g fhat=%g",r$ntree,r$rng,r$fhat), "\n")
r <- run(brts11, c(0.5,-0.001,0.3), model=c(1L,0L,0L), N=10L)
cat("dd (model 100) pars=(0.5,-0.001,0.3):", if(!is.null(r$msg)) r$msg else sprintf("ntree=%d rng=%.3g fhat=%g",r$ntree,r$rng,r$fhat), "\n")
cat("\n== rho < 1 with mu > lam ==\n")
r <- run(brts11, c(0.3,0.5), rho=0.7); cat("rho=0.7:", if(!is.null(r$msg)) r$msg else sprintf("ntree=%d rng=%.2e fhat=%g",r$ntree,r$rng,r$fhat),"\n")
