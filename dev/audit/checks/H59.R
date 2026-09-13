## H59: under CR max_tries = sample_size; overflow (max_missing) returns fewer trees, rejected = 0.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape)})
set.seed(4)
tr   <- ape::rphylo(20, 0.5, 0.1)
brts <- sort(ape::branching.times(tr), decreasing = TRUE)
pars8 <- c(0.5, 0,0,0, 0.4, 0,0,0)   # high turnover -> many missing lineages

# Reference: how many missing lineages does a CR draw typically have?
ref <- emphasis:::.augment_tree_bdi(brts, pars8, c(0L,0L,0L), sample_size = 30L, max_missing = 1e4L)
nm <- sapply(ref$trees, function(d) sum(d$t_ext == 0))
cat("missing lineages per draw (max_missing=1e4): ", paste(nm, collapse = " "), "\n")

for (mm in c(2L, 5L, 10L)) {
  a <- emphasis:::.augment_tree_bdi(brts, pars8, c(0L,0L,0L), sample_size = 50L, max_missing = mm)
  cat(sprintf("max_missing=%2d: requested 50, returned %2d trees; fhat=%.3f; names(out)=%s\n",
              mm, length(a$trees), a$fhat, paste(names(a), collapse = ",")))
}

# What .mcem_bdi records
for (mm in c(2L, 10L)) {
cat("== .mcem_bdi with max_missing =", mm, "\n")
m <- emphasis:::.mcem_bdi(brts, pars = c(0.5, 0.4), sample_size = 50L, max_missing = mm, verbose = TRUE,
                          lower_bound = c(0.01, 0.01), upper_bound = c(3, 3),
                          max_iter = 2L, xtol = 1e-3, tol = 1e-3, patience = 2L,
                          num_threads = 1L, model = c(0L,0L,0L), link = 0L)
print(m$mcem[, c("fhat", "rejected", "num_trees")]); cat("stop_reason:", m$stop_reason, "\n")
cat("final_IS$n_rejected:", m$final_IS$n_rejected, "\n")
}
