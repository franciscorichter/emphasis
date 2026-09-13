.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
library(emphasis)
brts20 <- c(5,4.4,3.9,3.6,3.1,2.8,2.5,2.2,2.0,1.7,1.5,1.3,1.1,0.9,0.75,0.6,0.45,0.3,0.15)
set.seed(3)
f <- estimate_rates(brts20, model="cr", method="mcem", init_pars=c(1,0.3),
  control=list(lower_bound=c(0,0), upper_bound=c(4,4), sampling="bdi",
               sample_size=20L, max_iter=3L, num_threads=1L))
cat("bdi trace cols:", paste(names(f$details$mcem), collapse=","), "\n")
cat("m_step:", f$details$mcem$m_step, " sum:", sum(f$details$mcem$m_step),
    " iterations:", f$iterations, " nrow:", nrow(f$details$mcem), "\n")
