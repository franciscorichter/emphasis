lib <- commandArgs(trailingOnly=TRUE)[1]
.libPaths(c(lib, .libPaths())); library(emphasis)
aug <- emphasis:::augment_trees
p8 <- c(0.5,0,0,0,0.1,0,0,0)
b12 <- c(10,8.3,7.1,6.2,5.5,4.4,3.9,3.1,2.2,1.5,0.9,0.3)
f <- function(r){lw<-r$logf-r$logg;m<-max(lw);log(sum(exp(lw-m)))+m-log(length(lw)+r$rejected_zero_weights)}
v <- replicate(30, f(aug(b12,p8,40L,3000L,10000L,1e6,1L,c(0L,0L,0L),0L,1)))
cat(sprintf("cr link0 nt=1: mean=%.4f sd=%.4f se=%.4f\n", mean(v), sd(v), sd(v)/sqrt(30)))
w <- replicate(30, f(aug(b12,p8,40L,3000L,10000L,1e6,1L,c(1L,0L,0L),0L,1)))  # dd model, betaN=0
cat(sprintf("dd(bN=0) nt=1: mean=%.4f sd=%.4f se=%.4f\n", mean(w), sd(w), sd(w)/sqrt(30)))
