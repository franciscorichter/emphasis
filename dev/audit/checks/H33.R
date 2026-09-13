# H33: what does max_lambda bound? per-lineage lambda or the total thinning intensity n*lambda*(...)?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(3)
tr <- ape::rcoal(40); brts <- sort(ape::branching.times(tr), decreasing = TRUE); brts <- brts/max(brts)*0.05
lam <- 30; mu <- 1e5    # per-lineage lambda = 30 << 500, but n*lambda > 500 once n >= 17 lineages
p8 <- emphasis:::.expand_pars(c(lam, mu), c(0L,0L,0L))
for (ml in c(500, 1e6)) {
  r <- tryCatch(emphasis:::augment_trees(brts, p8, sample_size = 50L, maxN = 200L, max_missing = 1e4, max_lambda = ml, num_threads = 1L),
                error = function(e) conditionMessage(e))
  if (is.character(r)) cat(sprintf("max_lambda=%g: E-step ERROR: %s\n", ml, r)) else
  cat(sprintf("max_lambda=%g: trees=%d rejected_lambda=%d rejected_overruns=%d zero_w=%d\n", ml,
              length(r$trees), r$rejected_lambda, r$rejected_overruns, r$rejected_zero_weights))
}
cat("MCEM hard-codes max_lambda=1e6 (R/emphasis.R:50); GAM control default 500 (R/inference.R:140); C++ default", 500, "\n")
cat("mcem control names:", paste(names(emphasis:::estimate_rates_control("mcem")), collapse=","), "\n")
