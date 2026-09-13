lib <- commandArgs(trailingOnly=TRUE)[1]
.libPaths(c(lib, .libPaths())); library(emphasis)
aug <- emphasis:::augment_trees
p8 <- function(...) { v <- c(...); c(v, rep(0, 8-length(v))) }
b4  <- c(4, 2.5, 1.2, 0.6)
b12 <- c(10,8.3,7.1,6.2,5.5,4.4,3.9,3.1,2.2,1.5,0.9,0.3)
say <- function(tag, r) {
  if (inherits(r,"error")) { cat(sprintf("%-38s ERROR: %s\n", tag, sub("\n.*","",conditionMessage(r)))); return(invisible()) }
  lw <- r$logf - r$logg
  cat(sprintf("%-38s M=%3d nt=%s zw=%4d nf=%s ovr=%d lam=%d rej=%d fhat=%.4f\n",
      tag, length(r$logf), if(is.null(r$num_trees)) "NA" else r$num_trees,
      r$rejected_zero_weights, if(is.null(r$rejected_nonfinite)) "NA" else r$rejected_nonfinite,
      r$rejected_overruns, r$rejected_lambda, r$rejected,
      log(sum(exp(lw-max(lw))))+max(lw)-log(length(lw)+r$rejected_zero_weights)))
}
go <- function(tag, brts, pars, N, maxN, model=c(0L,0L,0L), link=0L, nt=1L, rho=1) {
  r <- tryCatch(aug(brts, pars, as.integer(N), as.integer(maxN), 10000L, 1e6,
                    as.integer(nt), as.integer(model), as.integer(link), rho), error=function(e) e)
  say(tag, r)
}
cat("== links, cr, 1 thread ==\n")
for (lk in 0:2) go(sprintf("cr link=%d", lk), b12, p8(0.5,0,0,0,0.1), 40, 2000, c(0L,0L,0L), lk)
cat("== dd / d models ==\n")
go("dd link0", b12, p8(0.6,0,-0.02,0,0.1), 40, 2000, c(1L,0L,0L), 0L)
go("dd link1", b12, p8(log(0.6),0,-0.02,0,log(0.1)), 40, 2000, c(1L,0L,0L), 1L)
go("d  link0", b12, p8(0.6,0,0,-0.01,0.1), 40, 2000, c(0L,1L,0L), 0L)
go("d  link1", b12, p8(log(0.6),0,0,-0.01,log(0.1)), 40, 2000, c(0L,1L,0L), 1L)
cat("== boundaries ==\n")
go("mu = lambda", b12, p8(0.4,0,0,0,0.4), 40, 3000)
go("mu > lambda", b12, p8(0.2,0,0,0,0.6), 40, 3000)
go("rates near 0", b12, p8(1e-6,0,0,0,1e-7), 20, 2000)
go("lambda = 0", b4, p8(0,0,0,0,0.1), 10, 40)
go("mu = 0", b12, p8(0.5,0,0,0,0), 40, 2000)
go("sample_size 1", b4, p8(0.5,0,0,0,0.1), 1, 10)
go("maxN == sample_size", b4, p8(0.5,0,0,0,0.1), 20, 20)
go("maxN == ss, nt=10", b4, p8(0.5,0,0,0,0.1), 1, 1, nt=10L)
go("2-tip tree", c(3,1), p8(0.5,0,0,0,0.1), 20, 500)
go("rho=0.5", b12, p8(0.5,0,0,0,0.1), 40, 2000, rho=0.5)
cat("== dd with lambda hitting 0 (beta_N very negative) ==\n")
go("dd betaN=-0.5", b12, p8(0.6,0,-0.5,0,0.1), 20, 500, c(1L,0L,0L), 0L)
go("dd betaN=-0.05 gauss", b12, p8(0.6,0,-0.05,0,0.1), 20, 500, c(1L,0L,0L), 2L)
cat("== threads on the 12-tip tree ==\n")
for (nt in c(1L,2L,4L,8L)) { s <- replicate(6, { r <- tryCatch(aug(b12,p8(0.5,0,0,0,0.1),40L,3000L,10000L,1e6,nt,c(0L,0L,0L),0L,1), error=function(e) e)
  if (inherits(s0 <- r,"error")) NA_real_ else length(r$logf) }); cat(sprintf("nt=%d  M: %s\n", nt, paste(s, collapse=" "))) }
