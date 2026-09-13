.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(1); tr <- ape::rphylo(20, 0.5, 0.1)
brts <- sort(ape::branching.times(tr), decreasing = TRUE)
lb8 <- emphasis:::.expand_pars(c(0, 0), c(0L,0L,0L)); ub8 <- emphasis:::.expand_pars(c(1.5, 1.0), c(0L,0L,0L))
set.seed(2); pop <- emphasis:::.init_population(50L, lb8, ub8)
for (ss in c(1L, 5L, 20L)) {
  input <- list(brts = brts, sample_size = ss, maxN = 10L * ss, max_missing = 1e4, max_lambda = 1e6,
                lower_bound = lb8, upper_bound = ub8, shared_trees = FALSE, bias_correct = FALSE,
                model = c(0L,0L,0L), link = 0L, rho = 1)
  best <- replicate(10, { r <- emphasis:::.eval_particles(pop, input, 1L); max(r$pop$fhat, na.rm = TRUE) })
  cat(sprintf("sample_size=%2d  best fhat over identical population, 10 re-draws: mean=%.3f sd=%.3f range=%.3f  (tol=1e-4)\n",
              ss, mean(best), sd(best), diff(range(best))))
}
