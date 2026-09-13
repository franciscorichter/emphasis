source("/Users/pancho/Code/emphasis/dev/audit/review/1.6_common.R")
# does set.seed control em_cpp?
set.seed(2222); a <- drv(brts22, c(1,0.3), c(0,0), c(4,4), max_iter = 2L)$mcem$fhat
set.seed(2222); b <- drv(brts22, c(1,0.3), c(0,0), c(4,4), max_iter = 2L)$mcem$fhat
cat("same seed reproducible:", isTRUE(all.equal(a, b)), "\n")
res <- t(sapply(1:10, function(i) {
  d <- suppressWarnings(drv(brts22, c(1,0.3), c(0,0), c(4,4), sample_size = 300L, maxN = 200L, max_iter = 40L, max_time = 300))
  p <- emphasis:::.contract_pars(d$pars, c(0L,0L,0L))
  c(n_failed = d$n_failed, maxN = d$maxN, iter = d$iterations, dist = max(abs(p - 2)), mu = p[2], stop = d$stop_reason)
}))
print(res)
