# H34: does the per-iteration `rejected` trace column include overruns / lambda / zero-weight rejections?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(4)
tr <- ape::rcoal(20); brts <- sort(ape::branching.times(tr), decreasing = TRUE); brts <- brts/max(brts)*5
lb8 <- emphasis:::.expand_pars(c(0, 0), c(0L,0L,0L)); ub8 <- emphasis:::.expand_pars(c(2, 1), c(0L,0L,0L))
ip8 <- emphasis:::.expand_pars(c(0.6, 0.5), c(0L,0L,0L))
msgs <- character()
res <- withCallingHandlers(
  emphasis:::.mcem_dynamic_fresh(brts, ip8, sample_size = 20L, maxN = 5000L,
      max_missing = 2,          # tiny -> many augmentation overruns
      lower_bound = lb8, upper_bound = ub8, max_iter = 2L, xtol = 1e-3, tol = 1e-3, patience = 3L,
      num_threads = 1L, verbose = TRUE, model = c(0L,0L,0L), link = 0L, max_time = 120),
  message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
cat(msgs, sep = "")
print(res$mcem[, c("fhat","delta_max","rejected","num_trees")])
cat("final_IS$n_rejected (rejected+overruns+lambda):", res$final_IS$n_rejected,
    " zero_weights:", res$final_IS$rejected_zero_weights, "\n")
cat("trace `rejected` in last iter == sum of all counters? ",
    tail(res$mcem$rejected,1) == res$final_IS$n_rejected + res$final_IS$rejected_zero_weights, "\n")
