lib <- commandArgs(trailingOnly=TRUE)[1]
.libPaths(c(lib, .libPaths())); library(emphasis)
aug <- emphasis:::augment_trees
p8 <- function(...) { v <- c(...); c(v, rep(0, 8-length(v))) }
b12 <- c(10,8.3,7.1,6.2,5.5,4.4,3.9,3.1,2.2,1.5,0.9,0.3)
one <- function(nt) {
  r <- tryCatch(aug(b12, p8(0.6,0,-0.07,0,0.1), 10L, 3000L, 10000L, 1e6,
                    as.integer(nt), c(1L,0L,0L), 0L, 1), error=function(e) e)
  if (inherits(r,"error")) return(c(M=NA,zw=NA,fhat=NA))
  lw <- r$logf - r$logg; m <- max(lw)
  c(M=length(lw), zw=r$rejected_zero_weights,
    fhat=log(sum(exp(lw-m)))+m-log(length(lw)+r$rejected_zero_weights))
}
for (nt in c(1L,4L,8L)) { m <- t(replicate(12, one(nt)))
  cat(sprintf("nt=%d M=%s\n   zw=%s\n   fhat mean=%.3f sd=%.3f\n", nt,
      paste(m[,"M"],collapse=","), paste(m[,"zw"],collapse=","),
      mean(m[,"fhat"],na.rm=TRUE), sd(m[,"fhat"],na.rm=TRUE))) }
