lib <- commandArgs(trailingOnly=TRUE)[1]
.libPaths(c(lib, .libPaths())); library(emphasis)
aug <- emphasis:::augment_trees; ev <- emphasis:::eval_logf
p8 <- function(...) { v<-c(...); c(v, rep(0,8-length(v))) }
b12 <- c(10,8.3,7.1,6.2,5.5,4.4,3.9,3.1,2.2,1.5,0.9,0.3)
r <- aug(b12, p8(0.6,0,-0.06,0,0.1), 8L, 3000L, 10000L, 1e6, 1L, c(1L,0L,0L), 0L, 1)
for (bN in c(-0.06,-0.5,-2)) {
  z <- ev(p8(0.6,0,bN,0,0.1), r$trees, model=c(1L,0L,0L), link=0L, rho=1)
  cat(sprintf("betaN=%-5g  logf: %s\n            logg: %s\n            lw  : %s\n", bN,
      paste(format(head(z$logf,4), digits=4), collapse=" "),
      paste(format(head(z$logg,4), digits=4), collapse=" "),
      paste(format(head(z$logf-z$logg,4), digits=4), collapse=" ")))
}
