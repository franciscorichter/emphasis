lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, .libPaths())); library(emphasis)
A <- emphasis:::.augment_tree_bdi
brts <- c(6, 4.27428208320388, 4.1258975168475, 3.72207978339226,
  3.64993282835291, 3.17502914508657, 2.81163149722479, 2.48700670286076,
  2.38274529774947, 1.94421117175661, 1.86664932095232, 1.81164887493753,
  1.38503498765192, 1.06976980095552, 0.460758774034712, 0.417936203054209,
  0.363178692045526, 0.253829586012209, 0.0517782331486725)
grid <- expand.grid(K = c(20, 50, 1e4), i = 1:3)
lm3 <- rbind(c(0.8,0.3), c(0.6,0.1), c(1.0,0.8))
grid$l0 <- lm3[grid$i,1]; grid$m0 <- lm3[grid$i,2]
for (s in c(101, 202, 303)) {
  set.seed(s)
  g <- sapply(seq_len(nrow(grid)), function(j) {
    l0<-grid$l0[j]; m0<-grid$m0[j]; K<-grid$K[j]
    ref <- DDD::dd_loglik(c(l0,m0,K), c(300,1,0,1,0,2), brts, 0)
    e <- A(brts, c(l0, -(l0-m0)/K, m0, 0), c(1L,0L,0L), sample_size=500L,
           max_missing=1e4L, link=0L, rho=1)
    c(e$fhat - ref, if (is.null(e$acc)) NA else log(e$acc))
  })
  cat(sprintf("seed %d: max|gap|=%.3f range(gap)=%.3f  min_acc_gap=%s\n",
      s, max(abs(g[1,])), diff(range(g[1,])),
      if (all(is.na(g[2,]))) "n/a" else sprintf("%.3f", max(abs(g[1,]-g[2,])))))
}
