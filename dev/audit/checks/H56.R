## H56: BDI trees carry pd == 0 and no focal_tip_start; M/D-aware scoring sees M = 0.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape)})
set.seed(1)
tr   <- ape::rphylo(20, 0.5, 0.1)
brts <- sort(ape::branching.times(tr), decreasing = TRUE)

# BDI draw at CR pars
pars8 <- c(0.5, 0, 0, 0, 0.3, 0, 0, 0)
aug <- emphasis:::.augment_tree_bdi(brts, pars8, model_bin = c(0L,0L,0L), sample_size = 5L)
df  <- aug$trees[[which.max(sapply(aug$trees, nrow))]]
cat("BDI df columns:", paste(names(df), collapse = ","), "\n")
cat("BDI rows:", nrow(df), " all pd == 0:", all(df$pd == 0),
    " has focal_tip_start:", "focal_tip_start" %in% names(df), "\n")

# Thinning draw on the same tree: pd is populated
th <- emphasis:::augment_trees(brts, pars8, sample_size = 5L, maxN = 100L, max_missing = 1e4,
                    max_lambda = 500, num_threads = 1L)
tdf <- th$trees[[1]]
cat("Thinning df columns:", paste(names(tdf), collapse = ","), "\n")
cat("Thinning rows:", nrow(tdf), " max pd:", max(tdf$pd), "\n")

# Recompute pd for the BDI df under the package convention:
# P(t) = sum over alive non-extinction nodes with brts <= t of (t - tip_start)
recompute_pd <- function(df) {
  sapply(seq_len(nrow(df)), function(i) {
    t <- df$brts[i]
    sel <- df$t_ext != 0 & df$brts <= t & (df$t_ext > t)  # alive at t
    sum(t - df$tip_start[sel])
  })
}
df2 <- df; df2$pd <- recompute_pd(df)
cat("Recomputed pd range on BDI df:", range(df2$pd), "\n")

# Score with beta_M != 0 (model c(0,1,0), slot 3) : pd = 0 vs recomputed
pM <- c(0.5, 0, -0.01, 0, 0.3, 0, 0, 0)
lf0 <- emphasis:::eval_logf(pM, list(df),  model = c(0L,1L,0L), link = 0L)$logf
lf1 <- emphasis:::eval_logf(pM, list(df2), model = c(0L,1L,0L), link = 0L)$logf
cat(sprintf("beta_M=-0.01: logf(pd=0) = %.4f  logf(pd recomputed) = %.4f  diff = %.4f\n",
            lf0, lf1, lf1 - lf0))
# And under the gate-reachable config (beta_M = 0) pd is inert
lg0 <- emphasis:::eval_logf(pars8, list(df),  model = c(0L,0L,0L), link = 0L)$logf
lg1 <- emphasis:::eval_logf(pars8, list(df2), model = c(0L,0L,0L), link = 0L)$logf
cat(sprintf("beta_M=0 (reachable): logf(pd=0) = %.6f  logf(pd recomputed) = %.6f\n", lg0, lg1))
cat("bdi_supported(c(0,1,0),0):", emphasis:::.bdi_supported(c(0L,1L,0L), 0L), "\n")
