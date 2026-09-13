.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(1); tr <- ape::rphylo(20, 0.5, 0.1)
out <- capture.output(fit <- estimate_rates(tr, method = "cem", model = "cr",
  control = list(lower_bound = c(0, 0), upper_bound = c(1.5, 1.0),
                 num_particles = 10, max_iter = 4, num_trees = 1, num_threads = 1, verbose = TRUE)))
cat(grep("^Iter", out, value = TRUE), sep = "\n")
# with user-supplied sd_vec (compact) as well
out2 <- capture.output(fit2 <- estimate_rates(tr, method = "cem", model = "dd",
  control = list(lower_bound = c(0, -0.05, 0, -0.01), upper_bound = c(1.5, 0.05, 1.0, 0.01),
                 num_particles = 10, max_iter = 3, num_trees = 1, num_threads = 1, verbose = TRUE)))
cat(grep("^Iter", out2, value = TRUE), sep = "\n")
