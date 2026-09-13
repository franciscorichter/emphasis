# H44 replicate: is the tree$L (full table) vs phylo2L(tes) (extant-only) defect
# robust across trees / parameters? Two separable sub-claims:
#  (A) deterministic: nrow(L_aug) = nrow(tree$L) + n_aug, so tas carries the true
#      extinct lineages in addition to the augmented ones (always true if tree$L used)
#  (B) tree-dependent: augmented lineage attached to an extinct row -> negative edge
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape); library(DDD)})

par_sets <- list(c(0.6, 0.3), c(0.8, 0.5), c(0.5, 0.4), c(1.0, 0.7))
max_ts   <- c(6, 4, 8, 3)
n_draws  <- 30
out <- list()

for (ps in seq_along(par_sets)) {
  pars <- par_sets[[ps]]; max_t <- max_ts[ps]
  got <- 0
  for (i in 1:400) {
    sim <- simulate_tree(pars = pars, max_t = max_t, model = "cr", max_lin = 200)
    if (sim$status != "done" || is.null(sim$tes)) next
    n_ext <- sum(sim$L[, 4] != -1)
    if (Ntip(sim$tes) < 10 || Ntip(sim$tes) > 40 || n_ext < 3) next
    got <- got + 1
    if (got > 4) break

    brts  <- emphasis:::.extract_brts(sim)
    L_full <- emphasis:::.extract_Ltable(sim)
    L_tes  <- DDD::phylo2L(sim$tes)
    stopifnot(identical(L_full, sim$L))

    # deterministic predicate: which rows does id_to_label map observed
    # branchings (pid = 0..length(brts)-2) to in the full table, and are any
    # of them extinct (col4 != -1)?
    obs_ages <- brts[-1]  # brts[pid+2], pid >= 0
    mapped   <- vapply(obs_ages, function(b) which.min(abs(L_full[, 1] - b)), 1L)
    dead_map <- sum(L_full[mapped, 4] != -1)
    # in the extant-only table by construction none are extinct
    mapped_t <- vapply(obs_ages, function(b) which.min(abs(L_tes[, 1] - b)), 1L)
    stopifnot(all(L_tes[mapped_t, 4] == -1))

    run <- function(L_extant) {
      aug <- emphasis:::.augment_tree_internal(sim, pars = pars, model_bin = c(0L,0L,0L),
                sample_size = n_draws, max_missing = 1e4, num_threads = 1L, link = 0L, rho = 1)
      res <- t(sapply(aug$trees, function(df) {
        n_aug <- sum(df$t_ext != 0 & df$parent_id != -1L & df$t_ext < 1e11)
        L <- emphasis:::.aug_to_Ltable(df, max_t, brts, L_extant)
        phy <- tryCatch(DDD::L2phylo(L, dropextinct = FALSE), error = function(e) NULL)
        new <- L[seq_len(nrow(L)) > nrow(L_extant), , drop = FALSE]
        dead_parent <- 0L
        for (k in seq_len(nrow(new))) {
          p <- which(L[, 3] == new[k, 2])
          if (length(p) == 1 && L[p, 4] != -1 && L[p, 4] > new[k, 1] + 1e-9) dead_parent <- dead_parent + 1L
        }
        c(rows_ok = nrow(L) == nrow(L_extant) + n_aug,
          neg = as.integer(!is.null(phy) && min(phy$edge.length) < -1e-9),
          dead_parent = dead_parent,
          n_ext_tips = if (is.null(phy)) NA else Ntip(phy) - Ntip(emphasis:::prune_to_extant(phy)),
          n_aug = n_aug,
          extant_ok = if (is.null(phy)) NA else Ntip(emphasis:::prune_to_extant(phy)) == Ntip(sim$tes))
      }))
      c(rows_ok = mean(res[, "rows_ok"]), neg_draws = sum(res[, "neg"]),
        dead_rows = sum(res[, "dead_parent"]),
        ext_tips = mean(res[, "n_ext_tips"], na.rm = TRUE), n_aug = mean(res[, "n_aug"]),
        extant_ok = mean(res[, "extant_ok"], na.rm = TRUE))
    }
    a <- run(L_full); b <- run(L_tes)
    out[[length(out) + 1]] <- data.frame(
      pars = paste(pars, collapse = "/"), max_t = max_t, tips = Ntip(sim$tes),
      true_ext = n_ext, dead_mapped = dead_map, n_obs_brts = length(obs_ages),
      full_rows_ok = a["rows_ok"], full_neg = a["neg_draws"], full_dead = a["dead_rows"],
      full_exttips = round(a["ext_tips"], 1), full_naug = round(a["n_aug"], 1), full_extant_ok = a["extant_ok"],
      tes_neg = b["neg_draws"], tes_dead = b["dead_rows"], tes_exttips = round(b["ext_tips"], 1),
      tes_extant_ok = b["extant_ok"])
  }
}
res <- do.call(rbind, out); rownames(res) <- NULL
print(res)
cat("\nTrees with >=1 observed branching mapped to an extinct row:",
    sum(res$dead_mapped > 0), "/", nrow(res), "\n")
cat("Among those, trees with >=1 negative-edge draw (full table):",
    sum(res$full_neg > 0 & res$dead_mapped > 0), "; among the rest:",
    sum(res$full_neg > 0 & res$dead_mapped == 0), "\n")
cat("Negative-edge draws with extant-only table (all trees):", sum(res$tes_neg), "\n")
cat("full-table extinct tips ~ true_ext + n_aug:",
    all(abs(res$full_exttips - (res$true_ext + res$full_naug)) < 0.6), "\n")

# public API, n_trees = 1 and >1, with a sim list vs its $tes
cat("\n[public API, last tree]\n")
a1 <- simulate_tree(tree = sim, pars = pars, model = "cr", n_trees = 1, method = "thinning")
cat("n_trees=1 tree=sim: Ntip(tas) =", Ntip(a1$tas), " extinct tips =",
    Ntip(a1$tas) - Ntip(emphasis:::prune_to_extant(a1$tas)), " true extinct rows =", n_ext,
    " log_q =", a1$log_q, "\n")
b1 <- simulate_tree(tree = sim$tes, pars = pars, model = "cr", n_trees = 1, method = "thinning")
cat("n_trees=1 tree=sim$tes: Ntip(tas) =", Ntip(b1$tas), " extinct tips =",
    Ntip(b1$tas) - Ntip(emphasis:::prune_to_extant(b1$tas)), "\n")
