# H44: .extract_Ltable() returns tree$L (FULL L-table incl. true extinct lineages)
# for a simulate_tree() result, so conditional simulation appends augmented
# lineages to a table that already carries the true extinct ones.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape); library(DDD)})

pars <- c(0.6, 0.3); max_t <- 6
# simulate until a CR tree with 15-40 extant tips and >= 5 true extinct rows
set.seed(1)
for (i in 1:200) {
  sim <- simulate_tree(pars = pars, max_t = max_t, model = "cr")
  if (sim$status != "done" || is.null(sim$tes)) next
  n_ext <- sum(sim$L[, 4] != -1)
  if (Ntip(sim$tes) >= 15 && Ntip(sim$tes) <= 40 && n_ext >= 5) break
}
cat(sprintf("Tree: %d extant tips, L rows = %d, true extinct rows = %d\n",
            Ntip(sim$tes), nrow(sim$L), n_ext))

L_full <- emphasis:::.extract_Ltable(sim)          # what .sim_tree_conditional uses
L_tes  <- DDD::phylo2L(sim$tes)
cat(sprintf(".extract_Ltable(sim): %d rows (== tree$L: %s); phylo2L(tes): %d rows\n",
            nrow(L_full), identical(L_full, sim$L), nrow(L_tes)))

brts <- emphasis:::.extract_brts(sim)

analyse <- function(L_extant, n_draws = 50, method = "thinning", label = "") {
  aug <- if (method == "bdi") {
    emphasis:::.augment_tree_bdi(sim, pars = pars, model_bin = c(0L,0L,0L),
                                 sample_size = n_draws, max_missing = 1e4, link = 0L, rho = 1)
  } else {
    emphasis:::.augment_tree_internal(sim, pars = pars, model_bin = c(0L,0L,0L),
                                      sample_size = n_draws, max_missing = 1e4,
                                      num_threads = 1L, link = 0L, rho = 1)
  }
  res <- t(sapply(aug$trees, function(df) {
    n_aug <- sum(df$t_ext != 0 & df$parent_id != -1L & df$t_ext < 1e11)
    L <- emphasis:::.aug_to_Ltable(df, max_t, brts, L_extant)
    phy <- tryCatch(DDD::L2phylo(L, dropextinct = FALSE), error = function(e) NULL)
    min_edge <- if (is.null(phy)) NA else min(phy$edge.length)
    # parent-dies-before-child: new rows only
    new <- L[seq_len(nrow(L)) > nrow(L_extant), , drop = FALSE]
    bad_parent <- 0L
    if (nrow(new) > 0) for (k in seq_len(nrow(new))) {
      p <- which(L[, 3] == new[k, 2])
      if (length(p) == 1 && L[p, 4] != -1 && L[p, 4] > new[k, 1] + 1e-9) bad_parent <- bad_parent + 1L
    }
    # parent born after child (H45 mechanism, distinct from dead-parent)
    young_parent <- 0L
    if (nrow(new) > 0) for (k in seq_len(nrow(new))) {
      p <- which(L[, 3] == new[k, 2])
      if (length(p) == 1 && L[p, 1] < new[k, 1] - 1e-9) young_parent <- young_parent + 1L
    }
    n_extant_tips <- if (is.null(phy)) NA else Ntip(emphasis:::prune_to_extant(phy))
    c(nrow_L = nrow(L), n_aug = n_aug, expected = nrow(L_extant) + n_aug,
      min_edge = min_edge, neg_edge = as.integer(!is.na(min_edge) && min_edge < -1e-9),
      bad_parent = bad_parent, young_parent = young_parent, n_tips_tas = if (is.null(phy)) NA else Ntip(phy),
      n_extant_tas = n_extant_tips)
  }))
  cat(sprintf("\n[%s / %s] draws = %d\n", label, method, nrow(res)))
  cat(sprintf("  mean nrow(L_aug) = %.1f ; mean n_aug = %.1f ; nrow == nrow(L_extant)+n_aug in %d/%d\n",
              mean(res[,"nrow_L"]), mean(res[,"n_aug"]), sum(res[,"nrow_L"] == res[,"expected"]), nrow(res)))
  cat(sprintf("  draws with a negative edge in tas: %d/%d (min edge overall = %.4f)\n",
              sum(res[,"neg_edge"], na.rm=TRUE), nrow(res), min(res[,"min_edge"], na.rm=TRUE)))
  cat(sprintf("  draws with >=1 augmented lineage attached to a parent that died before its birth: %d/%d (total such rows %d)\n",
              sum(res[,"bad_parent"] > 0), nrow(res), sum(res[,"bad_parent"])))
  cat(sprintf("  draws with >=1 augmented lineage attached to a parent born AFTER it (H45): %d/%d (rows %d)\n",
              sum(res[,"young_parent"] > 0), nrow(res), sum(res[,"young_parent"])))
  cat(sprintf("  extinct tips in tas: mean %.1f (true extinct rows in input table = %d, mean n_aug = %.1f)\n",
              mean(res[,"n_tips_tas"] - res[,"n_extant_tas"], na.rm=TRUE),
              sum(L_extant[,4] != -1), mean(res[,"n_aug"])))
  invisible(res)
}

r1 <- analyse(L_full, method = "thinning", label = "L_extant = tree$L (current code)")
r2 <- analyse(L_tes,  method = "thinning", label = "L_extant = phylo2L(tes) (proposed)")
r3 <- analyse(L_full, method = "bdi",      label = "L_extant = tree$L (current code)")
r4 <- analyse(L_tes,  method = "bdi",      label = "L_extant = phylo2L(tes) (proposed)")

# Public API end-to-end: simulate_tree(tree = sim) vs simulate_tree(tree = sim$tes)
cat("\n[public API simulate_tree(tree=..., n_trees = 30)]\n")
for (m in c("thinning", "bdi")) {
  a <- simulate_tree(tree = sim,     pars = pars, model = "cr", n_trees = 30, method = m)
  b <- simulate_tree(tree = sim$tes, pars = pars, model = "cr", n_trees = 30, method = m)
  f <- function(x) { ok <- !vapply(x$trees, is.null, TRUE)
    e <- vapply(x$trees[ok], function(p) min(p$edge.length), 0)
    nt <- vapply(x$trees[ok], Ntip, 0L)
    sprintf("valid %d/%d, neg-edge draws %d, mean Ntip(tas) %.1f", sum(ok), length(ok), sum(e < -1e-9), mean(nt)) }
  cat(sprintf("  %-8s tree=sim:     %s\n", m, f(a)))
  cat(sprintf("  %-8s tree=sim$tes: %s\n", m, f(b)))
}
