.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths())); library(emphasis)
aug <- emphasis:::augment_trees
b60 <- c(20.2571,18.9106,18.8223,14.7584,13.9755,11.9092,11.3709,11.2664,11.0612,9.3217,9.2315,8.7633,7.4891,6.9616,6.339,6.2372,6.2025,5.3442,5.3046,5.1328,4.7774,4.7511,3.8254,3.8123,3.3753,3.2648,2.7234,2.7159,2.6661,2.3809,2.2352,2.1715,2.0224,1.9892,1.9652,1.9154,1.7791,1.7301,1.6915,1.529,1.4641,1.3952,1.3293,1.0584,1.0575,1.0053,0.8894,0.8093,0.7871,0.6553,0.6152,0.4845,0.4118,0.3098,0.2323,0.2081,0.1329,0.0911,0.0353)
p8 <- function(l,m) c(l,0,0,0,m,0,0,0)
f <- function(r){lw<-r$logf-r$logg;m<-max(lw);log(sum(exp(lw-m)))+m-log(length(lw)+r$rejected_zero_weights)}
one <- function(s) f(aug(b60*s, p8(0.2/s,0.05/s), 100L,1000L,10000L,1e6,1L,c(0L,0L,0L),0L,1))
a <- replicate(12, one(1)); b <- replicate(12, one(1e6))
n <- 60; pred <- -(n-2)*log(1e6)
cat(sprintf("fhat s=1   : mean=%.3f sd=%.3f\nfhat s=1e6 : mean=%.3f sd=%.3f\n", mean(a), sd(a), mean(b), sd(b)))
cat(sprintf("observed shift (12-rep means) = %.4f  predicted = %.4f  |diff| = %.4f (test tol 1.5, 4-rep means)\n",
            mean(b)-mean(a), pred, abs(mean(b)-mean(a)-pred)))
cat(sprintf("sd of a 4-rep-mean difference ~ %.3f\n", sqrt(sd(a)^2/4 + sd(b)^2/4)))
