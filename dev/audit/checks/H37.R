.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(1); tr <- ape::rphylo(20, 0.5, 0.1)
mm_log <- numeric(0)
trace(emphasis:::.simulate_particle, quote(assign("mm_log", c(get("mm_log", .GlobalEnv), max_missing), .GlobalEnv)),
      where = asNamespace("emphasis"), print = FALSE)
fit <- estimate_rates(tr, method = "cem", model = "cr",
  control = list(lower_bound = c(0, 0), upper_bound = c(1.5, 1.0),
                 max_missing = 1, maxN = 10, num_particles = 10, max_iter = 10,
                 num_trees = 1, num_threads = 1, max_time = 120))
untrace(emphasis:::.simulate_particle, where = asNamespace("emphasis"))
cat("stop:", fit$details$converged, " iters:", length(fit$details$best_loglik), "\n")
cat("rej_over per iter:", fit$details$history$rej_overruns, "\n")
cat("distinct max_missing seen (in order):", unique(mm_log), "\n")
cat("final-eval max_missing:", tail(mm_log, 1), "\n")
cat("as.integer(1e4*1.1^130) =", suppressWarnings(as.integer(1e4 * 1.1^130)), "\n")
