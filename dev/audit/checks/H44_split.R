#!/usr/bin/env Rscript
# H44 companion: separate the two causes of an inconsistent parent in the
# augmented L-table -- "parent born after child" (H45's parent-index mechanism,
# present with either base table) and "parent died before child" (H44's
# mechanism: sim$L rows whose birth age equals an observed branching but which
# went extinct after leaving an extant daughter).
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape); library(DDD) })
pars <- c(0.6, 0.25); max_t <- 6
split_bad <- function(L, n_base) {
  ba <- db <- 0L
  if (nrow(L) > n_base) for (k in seq.int(n_base + 1L, nrow(L))) {
    p <- which(L[, 3] == L[k, 2]); if (length(p) != 1L) next
    if (L[p, 1] < L[k, 1]) ba <- ba + 1L
    if (L[p, 4] != -1 && L[p, 4] > L[k, 1]) db <- db + 1L
  }
  c(born_after = ba, died_before = db)
}
for (rep in 1:8) {
  s2 <- NULL
  for (i in 1:300) { s <- simulate_tree(pars = pars, max_t = max_t, model = "cr", max_tries = 50)
    if (s$status == "done" && !is.null(s$tes) && Ntip(s$tes) >= 12 && Ntip(s$tes) <= 40) { s2 <- s; break } }
  if (is.null(s2)) next
  Lf <- s2$L; Le <- DDD::phylo2L(s2$tes); b2 <- emphasis:::.extract_brts(s2); tp2 <- b2[1]
  # how many observed branchings map (nearest birth age) to an extinct row in sim$L?
  n_dead_targets <- sum(sapply(b2[-1], function(b) Lf[which.min(abs(Lf[,1] - b)), 4] != -1))
  a2 <- emphasis:::.augment_tree_bdi(s2, pars = pars, model_bin = c(0L,0L,0L), sample_size = 40L,
                                     max_missing = 1e4, link = 0L, rho = 1.0)
  r2 <- t(sapply(a2$trees, function(df) {
    La <- emphasis:::.aug_to_Ltable(df, tp2, b2, Lf); Lb <- emphasis:::.aug_to_Ltable(df, tp2, b2, Le)
    ta <- tryCatch(DDD::L2phylo(La, dropextinct = FALSE), error = function(e) NULL)
    tb <- tryCatch(DDD::L2phylo(Lb, dropextinct = FALSE), error = function(e) NULL)
    c(f = split_bad(La, nrow(Lf)), e = split_bad(Lb, nrow(Le)),
      ne_f = if (is.null(ta)) NA else min(ta$edge.length) < 0,
      ne_e = if (is.null(tb)) NA else min(tb$edge.length) < 0,
      n_aug = nrow(Lb) - nrow(Le))
  }))
  cat(sprintf("tree %d: %2d tips, %2d extinct, %d/%d observed branchings map to an EXTINCT row in sim$L | augmented lineages: %4d | parent born-after child: full %3d / extant %3d | parent died-before child: full %3d / extant %3d | neg-edge draws (of 40): full %2d / extant %2d\n",
      rep, Ntip(s2$tes), sum(Lf[,4] != -1), n_dead_targets, length(b2) - 1, sum(r2[,"n_aug"]),
      sum(r2[,"f.born_after"]), sum(r2[,"e.born_after"]), sum(r2[,"f.died_before"]), sum(r2[,"e.died_before"]),
      sum(r2[,"ne_f"], na.rm=TRUE), sum(r2[,"ne_e"], na.rm=TRUE)))
}
