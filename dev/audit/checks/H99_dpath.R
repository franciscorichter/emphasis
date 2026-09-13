.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
g <- function(b0, eta) b0 * exp(-(eta - 1)^2 / 2)
brts <- c(4, 2.5, 1.2, 0.6)
aug <- emphasis:::augment_trees(brts, c(0.6, 0, 0, 0, 0.2, 0, 0, 0), 20L, 300L, 1000L, 100, 1L, c(0L, 0L, 0L), 0L)
print(head(aug$trees[[1]], 8)); cat("cols:", names(aug$trees[[1]]), "\n")
p <- c(0.7, 0.15, 0, 0, 0.3, -0.1, 0, 0)     # beta_D = gamma_D = 0
for (lk in c(2L, 1L, 0L)) {
  a <- emphasis:::eval_logf(p, aug$trees, c(1L, 0L, 0L), lk)$logf
  b <- emphasis:::eval_logf(p, aug$trees, c(1L, 0L, 1L), lk)$logf
  cat(sprintf("link %d: max|noD - D(beta_D=0)| = %.3e   (noD[1] = %.4f, D[1] = %.4f)\n", lk, max(abs(a - b)), a[1], b[1]))
}
# analytic: crown-lineage share of the compensator under gaussian, n0 = first-row n
d <- aug$trees[[1]]
lam <- g(p[1], p[2] * d$n); mu <- g(p[5], p[6] * d$n)
dt <- diff(c(0, d$brts))
crown <- 2 * sum(dt * (lam + mu))
a1 <- emphasis:::eval_logf(p, aug$trees[1], c(1L, 0L, 0L), 2L)$logf
b1 <- emphasis:::eval_logf(p, aug$trees[1], c(1L, 0L, 1L), 2L)$logf
cat(sprintf("tree1 gaussian: D - noD = %.4f ; crown-lineage compensator 2*sum(dt*(lam+mu)) = %.4f\n", b1 - a1, crown))
# per-row lineage accounting: how many rows are alive throughout each segment vs n
alive <- sapply(seq_len(nrow(d)), function(i) { prev <- if (i == 1) 0 else d$brts[i - 1]
  sum(d$t_ext[-i] != 0 & d$brts[-i] <= prev & d$t_ext[-i] >= d$brts[i]) })
print(data.frame(brts = d$brts, n = d$n, rows_alive = alive)[1:6, ])
