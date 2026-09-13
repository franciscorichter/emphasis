# H46: forward simulate_tree() with rho < 1 marks dropped tips L[,4] <- max_t.
# In DDD age units (col 1 = birth age, col 4 = death age, present = 0) that is
# "died at the crown", i.e. before its own birth -> negative edges in tas.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape); library(DDD)})
set.seed(46)  # reaches only the R-side rbinom/sample used for tip dropping

R <- 40
res <- data.frame(rep = seq_len(R), n_extant = NA, n_drop = NA, n_ext_true = NA,
                  min_edge_tas = NA, n_neg_edge = NA, drop_death_age = NA,
                  drop_birth_age_max = NA, tes_ok = NA, tas_ok = NA, both_clades = NA)
for (r in seq_len(R)) {
  s <- simulate_tree(pars = c(0.6, 0.1), max_t = 6, model = "cr", rho = 0.5, max_tries = 50)
  if (s$status != "done" || is.null(s$L)) next
  L <- s$L
  drop <- which(L[, 4] == 6)                     # rows stamped max_t
  res$n_extant[r]   <- sum(L[, 4] == -1) + length(drop)
  res$n_drop[r]     <- length(drop)
  res$n_ext_true[r] <- sum(L[, 4] > 0 & L[, 4] < 6)
  res$drop_death_age[r]     <- if (length(drop)) unique(L[drop, 4]) else NA
  res$drop_birth_age_max[r] <- if (length(drop)) max(L[drop, 1]) else NA
  res$tas_ok[r] <- !is.null(s$tas); res$tes_ok[r] <- !is.null(s$tes)
  if (!is.null(s$tas)) {
    res$min_edge_tas[r] <- min(s$tas$edge.length)
    res$n_neg_edge[r]   <- sum(s$tas$edge.length < 0)
  }
  # survivor per crown clade among *kept* extant tips (labels: sign = crown clade)
  kept <- L[L[, 4] == -1, 3]
  res$both_clades[r] <- any(kept > 0) && any(kept < 0)
}
print(res)
cat("\nSummary over", sum(!is.na(res$n_drop)), "successful sims:\n")
cat("  reps with n_drop > 0            :", sum(res$n_drop > 0, na.rm = TRUE), "\n")
cat("  reps with negative tas edges    :", sum(res$n_neg_edge > 0, na.rm = TRUE), "\n")
cat("  min tas edge overall            :", min(res$min_edge_tas, na.rm = TRUE), "\n")
cat("  max birth age of a dropped tip  :", max(res$drop_birth_age_max, na.rm = TRUE),
    " (death age stamped =", unique(na.omit(res$drop_death_age)), ")\n")
cat("  reps where kept tips lost a crown clade:", sum(!res$both_clades, na.rm = TRUE), "\n")

# --- tes: is it really unaffected? tes = L2phylo(L, dropextinct = TRUE) drops
# every row with L[,4] != -1, so the stamped rows vanish; check tip count & ultrametric
s <- simulate_tree(pars = c(0.6, 0.1), max_t = 6, model = "cr", rho = 0.5, max_tries = 50)
cat("\nSingle sim: n_extant_kept =", sum(s$L[,4] == -1), " Ntip(tes) =", Ntip(s$tes),
    " ultrametric(tes) =", is.ultrametric(s$tes), " Ntip(tas) =", Ntip(s$tas),
    " min edge tas =", min(s$tas$edge.length), "\n")

# --- The proposed fix: L[drop,4] <- 0 (died at present = unsampled extant)
L2 <- s$L; L2[L2[,4] == 6, 4] <- 0
tas2 <- DDD::L2phylo(L2, dropextinct = FALSE)
tes2 <- DDD::L2phylo(L2, dropextinct = TRUE)
cat("With L[drop,4] <- 0: min edge tas =", min(tas2$edge.length),
    " Ntip(tas) =", Ntip(tas2), " Ntip(tes) =", Ntip(tes2),
    " tes identical topology/brts to current tes:",
    isTRUE(all.equal(sort(branching.times(tes2)), sort(branching.times(s$tes)))), "\n")

# --- Inheritance by H44 path: .sim_tree_conditional uses tree$L (the stamped table)
# as L_extant, so the stamped rows go into every augmented tas.
aug <- simulate_tree(s, pars = c(0.6, 0.1), model = "cr", n_trees = 10, method = "thinning")
negs <- vapply(aug$trees, function(t) if (is.null(t)) NA_real_ else min(t$edge.length), 0)
cat("Conditional augmentation on the rho<1 sim object: min edge per draw:\n"); print(round(negs, 3))
# Same augmentation but feeding tes only (no stamped L)
aug2 <- simulate_tree(s$tes, pars = c(0.6, 0.1), model = "cr", n_trees = 10, method = "thinning")
negs2 <- vapply(aug2$trees, function(t) if (is.null(t)) NA_real_ else min(t$edge.length), 0)
cat("Conditional augmentation on tes only: min edge per draw:\n"); print(round(negs2, 3))

# --- Does the forward rho path touch a rho = 1 fit?  Count code paths.
cat("\nrho = 1 sim: any L[,4] == max_t rows?",
    any(simulate_tree(pars = c(0.6, 0.1), max_t = 6, model = "cr", rho = 1, max_tries = 50)$L[,4] == 6), "\n")
