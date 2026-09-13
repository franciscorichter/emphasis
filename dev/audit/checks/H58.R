## H58: does .bdi_iterate converge within 20 iterations on a strong-DD tree? (no report in the package)
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape)})
ns <- asNamespace("emphasis")
# Patch a copy of .bdi_iterate so it records delta per iteration in .H58_delta
src <- paste(deparse(get(".bdi_iterate", ns), width.cutoff = 500L), collapse = "\n")
src <- sub("if \\(delta < tol\\)\\s*break", "{ .H58_delta[iter] <<- delta; if (delta < tol) break }", src)
stopifnot(any(grepl(".H58_delta", src, fixed = TRUE)))
iter_patched <- eval(parse(text = src), envir = ns)

run <- function(tr, pars8, label, ...) {
  brts <- sort(ape::branching.times(tr), decreasing = TRUE)
  tp <- brts[1]; bt <- sort(tp - brts[-1])
  .H58_delta <<- rep(NA_real_, 100)
  t0 <- proc.time()[3]
  sol <- iter_patched(pars8, c(1L,0L,0L), 0L, bt, tp, ...)
  d <- .H58_delta[!is.na(.H58_delta)]
  cat(sprintf("%-28s tips=%3d iters=%2d final_delta=%.2e converged(<1e-4)=%s  time=%.1fs  Nhat(tp)=%.1f k=%d\n",
              label, length(tr$tip.label), length(d), tail(d, 1), tail(d, 1) < 1e-4,
              proc.time()[3] - t0, sol$Nhat_fun(tp), length(brts) + 1))
  invisible(d)
}
set.seed(3)
tr20 <- ape::rphylo(20, 0.5, 0.1)
tr40 <- ape::rphylo(40, 0.5, 0.1)
# mild DD
d1 <- run(tr20, c(0.6, -0.01, 0,0, 0.1, 0,0,0), "mild DD (K=60)")
# strong DD: K = 20 ~ tips (rates pinned near zero at the tips)
d2 <- run(tr20, c(1.0, -0.05, 0,0, 0.2, 0,0,0), "strong DD (K=20, tips=20)")
# over-saturated: K below the tip count, lambda clipped to 0 by linear link
d3 <- run(tr20, c(1.0, -0.10, 0,0, 0.3, 0,0,0), "over-saturated (K=10)")
# high turnover strong DD on 40 tips
d4 <- run(tr40, c(2.0, -0.05, 0,0, 1.0, 0,0,0), "high turnover DD (K=40)")
d5 <- run(tr40, c(2.0, -0.05, 0,0, 1.0, 0,0,0), "same, no gaussian closure", use_gaussian_closure = FALSE)
cat("delta trajectory, high-turnover:", format(signif(d4, 3)), "\n")
# K above the tip count (no lambda clipping at the tips) but high turnover
d6 <- run(tr40, c(3.0, -0.05, 0,0, 1.5, 0,0,0), "high turnover DD (K=60)")
d7 <- run(tr40, c(1.5, -0.02, 0,0, 0.9, 0,0,0), "high turnover DD (K=75)")
d8 <- run(tr40, c(3.0, -0.05, 0,0, 1.5, 0,0,0), "K=60, max_iter=100", max_iter = 100)
cat("delta trajectory, K=60 (first 20):", format(signif(d6, 3)), "\n")
