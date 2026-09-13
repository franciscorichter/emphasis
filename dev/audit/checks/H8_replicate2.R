.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
src <- readLines("/Users/pancho/Code/emphasis/dev/audit/checks/H8_replicate.R")
eval(parse(text = src[1:(grep("^brts <- c\\(15", src) - 1)]))   # reuse helper defs only
brts <- c(15, 12.1, 10.3, 9.0, 7.7, 6.2, 5.1, 4.4, 3.0, 2.2, 1.1)
p_mild <- c(0.25, -0.005, 0, 0.02, 0.05, 0, 0, 0.01)
aug <- emphasis:::augment_trees(brts, p_mild, sample_size = 30, maxN = 4000, max_missing = 300, max_lambda = 1e6,
                                num_threads = 1, model = c(0L,0L,1L), link = 1L, rho = 1)
trees <- aug$trees; nmiss <- sapply(trees, function(t) sum(is_ext(t$t_ext)))
cat(sprintf("exp-sampler trees: %d, missing %s\n", length(trees), paste(range(nmiss), collapse="-")))
cpp <- sapply(trees, function(t) cpp_logf(p_mild, t, 0)); ex <- sapply(trees, function(t) logf_R(p_mild, t, 0, "exact")[["loglik"]])
ok <- is.finite(cpp) & is.finite(ex); d <- (cpp - ex)[ok]
cat(sprintf(" linear eval on exp-sampled trees: finite %d/%d; Delta mean %+.3f sd %.3f range [%+.3f, %+.3f]\n", sum(ok), length(ok), mean(d), sd(d), min(d), max(d)))
cat(sprintf(" exp anchor: %.2e\n", max(abs(sapply(trees, function(t) cpp_logf(p_mild, t, 1) - logf_R(p_mild, t, 1, "exact")[["loglik"]])))))
## weighted Q along beta_D, linear sampler, 12 tips
p_gen <- c(0.4, -0.01, 0, 0.05, 0.1, 0, 0, 0.02)
augL <- emphasis:::augment_trees(brts, p_gen, sample_size = 30, maxN = 4000, max_missing = 300, max_lambda = 1e6,
                                 num_threads = 1, model = c(0L,0L,1L), link = 0L, rho = 1)
lw <- augL$logf - augL$logg; w <- exp(lw - max(lw)); grid <- seq(-0.1, 0.15, by = 0.01)
Q <- function(fun) sapply(grid, function(bD) { p <- p_gen; p[4] <- bD; v <- sapply(augL$trees, function(t) fun(p, t)); if (all(is.finite(v))) sum(w*v) else -Inf })
Qc <- Q(function(p,t) cpp_logf(p,t,0)); Qe <- Q(function(p,t) logf_R(p,t,0,"exact")[["loglik"]])
cat(sprintf(" 12-tip weighted Q argmax beta_D: code %.2f, exact %.2f (gen 0.05); finite grid pts code %d exact %d of %d\n",
            grid[which.max(Qc)], grid[which.max(Qe)], sum(is.finite(Qc)), sum(is.finite(Qe)), length(grid)))
