lib <- commandArgs(trailingOnly=TRUE)[1]
.libPaths(c(lib, .libPaths())); library(emphasis)
aug <- emphasis:::augment_trees
p8 <- function(...) { v <- c(...); c(v, rep(0, 8-length(v))) }
b12 <- c(10,8.3,7.1,6.2,5.5,4.4,3.9,3.1,2.2,1.5,0.9,0.3)
one <- function(bN, N, maxN) {
  r <- tryCatch(aug(b12, p8(0.6,0,bN,0,0.1), as.integer(N), as.integer(maxN), 10000L, 1e6,
                    1L, c(1L,0L,0L), 0L, 1), error=function(e) e)
  if (inherits(r,"error")) { m <- conditionMessage(r)
    zw <- as.integer(sub(".*; ([0-9]+) zero weights.*","\\1",m))
    nf <- if (grepl("non-finite", m)) as.integer(sub(".*; ([0-9]+) non-finite.*","\\1",m)) else NA
    c(ok=0, M=NA, zw=zw, nf=nf) } else
    c(ok=1, M=length(r$logf), zw=r$rejected_zero_weights,
      nf=if(is.null(r$rejected_nonfinite)) NA else r$rejected_nonfinite)
}
for (bN in c(-0.40,-0.50,-0.60)) for (N in c(2L,5L)) {
  m <- t(replicate(5, one(bN, N, 400L)))
  cat(sprintf("betaN=%.2f N=%d ok=%s M=%s zw=%s nf=%s\n", bN, N,
    paste(m[,"ok"],collapse=","), paste(m[,"M"],collapse=","),
    paste(m[,"zw"],collapse=","), paste(m[,"nf"],collapse=",")))
}
