#!/usr/bin/env Rscript
# H44 trace: on one simulated tree, list the observed branchings whose
# nearest-birth-age row in sim$L is a TRUE-EXTINCT lineage, then follow the
# augmented rows that .aug_to_Ltable attaches to those labels and compare each
# child's birth age with the parent's death age.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape); library(DDD) })
pars <- c(0.6, 0.25); max_t <- 6
sim <- NULL
for (i in 1:500) { s <- simulate_tree(pars = pars, max_t = max_t, model = "cr", max_tries = 50)
  if (s$status == "done" && !is.null(s$tes) && Ntip(s$tes) >= 12 && Ntip(s$tes) <= 35) {
    L <- s$L; b <- emphasis:::.extract_brts(s)
    tgt <- sapply(b[-1], function(x) which.min(abs(L[,1] - x)))
    if (any(L[tgt, 4] != -1)) { sim <- s; break } } }
stopifnot(!is.null(sim))
L <- sim$L; brts <- emphasis:::.extract_brts(sim); tp <- brts[1]
tgt <- sapply(brts[-1], function(x) which.min(abs(L[,1] - x)))
dead <- which(L[tgt, 4] != -1)
cat(sprintf("tree: %d tips, %d rows, %d true extinct; branchings mapping to extinct rows:\n", Ntip(sim$tes), nrow(L), sum(L[,4] != -1)))
for (j in dead) cat(sprintf("  brts[%d] = %.5f -> sim$L row %d label %d, born %.5f, died %.5f (next observed brts %.5f)\n",
                            j + 1, brts[j + 1], tgt[j], L[tgt[j], 3], L[tgt[j], 1], L[tgt[j], 4],
                            if (j + 2 <= length(brts)) brts[j + 2] else 0))
aug <- emphasis:::.augment_tree_bdi(sim, pars = pars, model_bin = c(0L,0L,0L), sample_size = 40L,
                                    max_missing = 1e4, link = 0L, rho = 1.0)
dead_lbl <- L[tgt[dead], 3]
n_att <- n_dead_before <- n_negedge <- 0L; n_negedge_ext <- 0L
Le <- DDD::phylo2L(sim$tes)
for (df in aug$trees) {
  La <- emphasis:::.aug_to_Ltable(df, tp, brts, L)
  Lb <- emphasis:::.aug_to_Ltable(df, tp, brts, Le)
  if (nrow(La) > nrow(L)) for (k in seq.int(nrow(L) + 1L, nrow(La))) {
    if (La[k, 2] %in% dead_lbl) {
      n_att <- n_att + 1L
      p <- which(La[, 3] == La[k, 2])
      if (La[p, 4] > La[k, 1]) n_dead_before <- n_dead_before + 1L
    }
  }
  ta <- tryCatch(DDD::L2phylo(La, dropextinct = FALSE), error = function(e) NULL)
  tb <- tryCatch(DDD::L2phylo(Lb, dropextinct = FALSE), error = function(e) NULL)
  if (!is.null(ta) && min(ta$edge.length) < 0) n_negedge <- n_negedge + 1L
  if (!is.null(tb) && min(tb$edge.length) < 0) n_negedge_ext <- n_negedge_ext + 1L
}
cat(sprintf("\n40 BDI draws: %d augmented rows attached to a true-extinct label; %d of them born AFTER that parent's death\n", n_att, n_dead_before))
cat(sprintf("negative-edge tas: full-base %d/40, extant-base %d/40\n", n_negedge, n_negedge_ext))
# Show a concrete example
for (df in aug$trees) {
  La <- emphasis:::.aug_to_Ltable(df, tp, brts, L)
  if (nrow(La) > nrow(L)) for (k in seq.int(nrow(L) + 1L, nrow(La))) if (La[k, 2] %in% dead_lbl) {
    p <- which(La[, 3] == La[k, 2])
    cat(sprintf("example: augmented label %d born at age %.4f, parent label %d (born %.4f, died %.4f) -> child born %s parent's death\n",
                La[k, 3], La[k, 1], La[k, 2], La[p, 1], La[p, 4], if (La[p, 4] > La[k, 1]) "AFTER" else "before"))
    break
  }
  if (exists("p")) break
}

# ---- targeted: find a tree where an extinct target row dies BEFORE the next
#      observed branching (the only configuration in which a child can be
#      attached after its parent's death) ---------------------------------------
cat("\n--- targeted search: extinct target row dying before the next observed branching ---\n")
found <- NULL
for (i in 1:3000) {
  s <- simulate_tree(pars = pars, max_t = max_t, model = "cr", max_tries = 50)
  if (s$status != "done" || is.null(s$tes) || Ntip(s$tes) < 8 || Ntip(s$tes) > 40) next
  L <- s$L; b <- emphasis:::.extract_brts(s)
  tgt <- sapply(b[-1], function(x) which.min(abs(L[,1] - x)))
  nxt <- c(b[-(1:2)], 0)
  gap <- which(L[tgt, 4] != -1 & L[tgt, 4] > nxt)
  if (length(gap)) { found <- s; break }
}
if (is.null(found)) { cat("none found in 3000 trees\n") } else {
  L <- found$L; brts <- emphasis:::.extract_brts(found); tp <- brts[1]
  tgt <- sapply(brts[-1], function(x) which.min(abs(L[,1] - x))); nxt <- c(brts[-(1:2)], 0)
  gap <- which(L[tgt, 4] != -1 & L[tgt, 4] > nxt)
  for (j in gap) cat(sprintf("  brts[%d] = %.5f -> row label %d born %.5f died %.5f; next observed brts %.5f -> children born in (%.4f, %.4f) attach to a dead parent\n",
                             j + 1, brts[j + 1], L[tgt[j], 3], L[tgt[j], 1], L[tgt[j], 4], nxt[j], nxt[j], L[tgt[j], 4]))
  Le <- DDD::phylo2L(found$tes)
  aug <- emphasis:::.augment_tree_bdi(found, pars = pars, model_bin = c(0L,0L,0L), sample_size = 60L,
                                      max_missing = 1e4, link = 0L, rho = 1.0)
  db_f <- db_e <- ne_f <- ne_e <- 0L
  for (df in aug$trees) {
    La <- emphasis:::.aug_to_Ltable(df, tp, brts, L); Lb <- emphasis:::.aug_to_Ltable(df, tp, brts, Le)
    cnt <- function(X, nb) { d <- 0L; if (nrow(X) > nb) for (k in seq.int(nb + 1L, nrow(X))) {
      p <- which(X[,3] == X[k,2]); if (length(p) == 1L && X[p,4] != -1 && X[p,4] > X[k,1]) d <- d + 1L }; d }
    db_f <- db_f + cnt(La, nrow(L)); db_e <- db_e + cnt(Lb, nrow(Le))
    ta <- tryCatch(DDD::L2phylo(La, dropextinct = FALSE), error = function(e) NULL)
    tb <- tryCatch(DDD::L2phylo(Lb, dropextinct = FALSE), error = function(e) NULL)
    if (!is.null(ta) && min(ta$edge.length) < 0) ne_f <- ne_f + 1L
    if (!is.null(tb) && min(tb$edge.length) < 0) ne_e <- ne_e + 1L
  }
  cat(sprintf("  %d tips, %d true extinct; 60 BDI draws: children attached AFTER parent's death: full-base %d, extant-base %d; negative-edge tas: full-base %d/60, extant-base %d/60\n",
              Ntip(found$tes), sum(L[,4] != -1), db_f, db_e, ne_f, ne_e))
}
