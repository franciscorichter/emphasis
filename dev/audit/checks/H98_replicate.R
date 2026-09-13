# H98 replication: does the corrected (pars=) form of the three skipped tests run end-to-end?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
run <- function(method, ctrl) {
  sim <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")
  cat(method, ": sim status", sim$status, "tips", if (inherits(sim$tes,"phylo")) ape::Ntip(sim$tes) else NA, "\n")
  fit <- tryCatch(estimate_rates(sim, method = method, model = "cr", control = ctrl),
                  error = function(e) e)
  if (inherits(fit, "error")) cat("  ERROR:", conditionMessage(fit), "\n") else
    cat("  class:", class(fit)[1], "| pars:", round(fit$pars, 3), "| loglik:", fit$loglik, "| method:", fit$method, "\n")
}
run("mcem", list(lower_bound = c(0, 0), upper_bound = c(2, 1), max_iter = 3, sample_size = 30))
run("cem",  list(lower_bound = c(0, 0), upper_bound = c(2, 1), max_iter = 3, num_particles = 20))
run("gam",  list(lower_bound = c(0.1, 0.01), upper_bound = c(1.5, 0.5), grid_points = 6, sample_size = 30))
# Variation: positional form with a real tree in slot 1 does NOT rescue it either
tr <- ape::rphylo(10, 0.5, 0.1)
r <- tryCatch(simulate_tree(tr, max_t = 5, model = "cr"), error = function(e) conditionMessage(e)); cat("tree-only positional:", r, "\n")
