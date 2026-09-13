lib <- commandArgs(trailingOnly=TRUE)[1]
.libPaths(c(lib, .libPaths())); library(emphasis)
aug <- emphasis:::augment_trees
p8 <- function(...) { v <- c(...); c(v, rep(0, 8-length(v))) }
b12 <- c(10,8.3,7.1,6.2,5.5,4.4,3.9,3.1,2.2,1.5,0.9,0.3)
f <- function(r) { lw <- r$logf - r$logg; m <- max(lw)
  c(M=length(lw), zw=r$rejected_zero_weights,
    nf=if(is.null(r$rejected_nonfinite)) NA else r$rejected_nonfinite,
    fhat=log(sum(exp(lw-m)))+m-log(length(lw)+r$rejected_zero_weights)) }
run <- function(bN, N=20, maxN=4000, nt=1L)
  tryCatch(f(aug(b12, p8(0.6,0,bN,0,0.1), as.integer(N), as.integer(maxN), 10000L, 1e6,
                 nt, c(1L,0L,0L), 0L, 1)), error=function(e) c(M=NA,zw=NA,nf=NA,fhat=NA))
for (bN in c(-0.10,-0.15,-0.20,-0.30)) {
  m <- t(replicate(8, run(bN)))
  cat(sprintf("betaN=%.2f  M=%s zw=%s nf=%s  fhat mean=%.4f sd=%.4f\n", bN,
    paste(m[,"M"],collapse=","), paste(m[,"zw"],collapse=","), paste(m[,"nf"],collapse=","),
    mean(m[,"fhat"],na.rm=TRUE), sd(m[,"fhat"],na.rm=TRUE)))
}
cat("-- contention with rejections, nt=8 --\n")
m <- t(replicate(15, run(-0.20, N=20, maxN=4000, nt=8L)))
cat(sprintf("nt=8 M=%s\n", paste(m[,"M"],collapse=",")))
cat(sprintf("nt=8 fhat mean=%.4f sd=%.4f\n", mean(m[,"fhat"],na.rm=TRUE), sd(m[,"fhat"],na.rm=TRUE)))
cat("-- benign cr link0, 1 thread, 25 reps --\n")
g <- function() { r <- aug(b12, p8(0.5,0,0,0,0.1), 40L, 3000L, 10000L, 1e6, 1L, c(0L,0L,0L), 0L, 1); f(r) }
m <- t(replicate(25, g()))
cat(sprintf("fhat mean=%.4f sd=%.4f se=%.4f zw=%s\n", mean(m[,"fhat"]), sd(m[,"fhat"]), sd(m[,"fhat"])/5, paste(unique(m[,"zw"]),collapse=",")))
cat("-- cr link1 (exp), 1 thread, 15 reps: H1 cut in natural use --\n")
h <- function() { r <- aug(b12, p8(0.5,0,0,0,0.1), 40L, 3000L, 10000L, 1e6, 1L, c(0L,0L,0L), 1L, 1); f(r) }
m <- t(replicate(15, h()))
cat(sprintf("fhat mean=%.4f sd=%.4f  zw total=%d  nf total=%s\n", mean(m[,"fhat"]), sd(m[,"fhat"]), sum(m[,"zw"]), sum(m[,"nf"])))
