lib <- commandArgs(trailingOnly=TRUE)[1]; case <- commandArgs(trailingOnly=TRUE)[2]
.libPaths(c(lib, .libPaths())); library(emphasis)
aug <- emphasis:::augment_trees
p8 <- c(0.5,0,0,0,0.1,0,0,0); b4 <- c(4,2.5,1.2,0.6)
r <- switch(case,
  "N0"    = aug(b4, p8, 0L, 20L, 10000L, 1e6, 1L, c(0L,0L,0L), 0L, 1),
  "Nneg"  = aug(b4, p8, -1L, 20L, 10000L, 1e6, 1L, c(0L,0L,0L), 0L, 1),
  "maxN0" = aug(b4, p8, 5L, 0L, 10000L, 1e6, 1L, c(0L,0L,0L), 0L, 1),
  "NgtmaxN" = aug(b4, p8, 50L, 5L, 10000L, 1e6, 1L, c(0L,0L,0L), 0L, 1))
cat(case, " M=", length(r$logf), " num_trees=", if (is.null(r$num_trees)) NA else r$num_trees,
    " fhat_ok\n", sep="")
