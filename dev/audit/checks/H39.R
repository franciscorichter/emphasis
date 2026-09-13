.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(1); tr <- ape::rphylo(20, 0.5, 0.1)
pm <- as.matrix(expand.grid(beta_0 = seq(0.2, 1, length.out = 4), gamma_0 = seq(0, 0.4, length.out = 4)))
for (nt in c(1L, 2L)) {
  t0 <- proc.time()[3]
  r <- tryCatch(emphasis:::estimate_likelihood_surface(tr, pm, model = "cr", sample_size = 10L,
        num_threads = nt, verbose = FALSE, max_time = 0.001),
        error = function(e) conditionMessage(e))
  cat(sprintf("num_threads=%d  elapsed=%.2fs  result=%s\n", nt, proc.time()[3] - t0,
      if (is.character(r)) paste("ERROR:", r) else sprintf("surface with %d rows, %d finite", nrow(r), sum(is.finite(r$fhat)))))
}
