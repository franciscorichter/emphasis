lib <- commandArgs(trailingOnly=TRUE)[1]
.libPaths(c(lib, .libPaths())); library(emphasis)
aug <- emphasis:::augment_trees
p8 <- function(...) { v <- c(...); c(v, rep(0, 8-length(v))) }
b12 <- c(10,8.3,7.1,6.2,5.5,4.4,3.9,3.1,2.2,1.5,0.9,0.3)
one <- function(bN, N, maxN) {
  r <- tryCatch(aug(b12, p8(0.6,0,bN,0,0.1), as.integer(N), as.integer(maxN), 10000L, 1e6,
                    1L, c(1L,0L,0L), 0L, 1), error=function(e) e)
  if (inherits(r,"error")) return(c(ok=0,M=NA,zw=NA,nf=NA,fhat=NA))
  lw <- r$logf - r$logg; m <- max(lw)
  c(ok=1, M=length(lw), zw=r$rejected_zero_weights,
    nf=if(is.null(r$rejected_nonfinite)) NA else r$rejected_nonfinite,
    fhat=log(sum(exp(lw-m)))+m-log(length(lw)+r$rejected_zero_weights))
}
for (bN in c(-0.06,-0.07,-0.08,-0.09)) {
  m <- t(replicate(6, one(bN, 5L, 3000L)))
  cat(sprintf("betaN=%.2f ok=%s zw=%s nf=%s fhat=%.3f\n", bN, paste(m[,"ok"],collapse=","),
      paste(m[,"zw"],collapse=","), paste(m[,"nf"],collapse=","), mean(m[,"fhat"],na.rm=TRUE)))
}
