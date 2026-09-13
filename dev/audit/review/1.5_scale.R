.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths())); library(emphasis)
brts12 <- c(6.385824,2.063997,1.19255,0.923743,0.884126,0.822976,0.725585,0.718801,0.539214,0.077287,0.068123)
run <- function(s) { set.seed(99)
  f <- estimate_rates(brts12*s, model="cr", method="mcem", init_pars=c(1.2,0.9)/s,
    control=list(lower_bound=c(0,0), upper_bound=c(3,3)/s, sampling="bdi",
                 sample_size=200L, num_threads=1L, max_iter=60L, tol=1e-2, patience=3L))
  m <- f$details$mcem; em <- m[m$m_step,]
  cat(sprintf("scale %-7g stop=%-10s iters=%-3d pars*s = %.4f %.4f   delta[1:5]= %s\n",
      s, f$details$stop_reason, f$details$iterations, f$pars[1]*s, f$pars[2]*s,
      paste(signif(head(em$delta_max,5),2), collapse=" "))) }
for (s in c(1, 100, 1000)) run(s)
