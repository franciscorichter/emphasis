## H45 — parent assignment in the augmentation samplers and what reaches `tas`.
## Self-contained; ~1 min. Uses the scratch build only.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
set.seed(45)

## one reconstructed tree, 20 tips
tr    <- TreeSim::sim.bd.taxa(n = 20, numbsim = 1, lambda = 0.8, mu = 0.4,
                              complete = FALSE)[[1]]
brts  <- emphasis:::.extract_brts(tr)         # ages, crown first
Tc    <- brts[1]
s0    <- Tc - brts[2]                         # forward time of first post-crown split
bt_s  <- sort(Tc - brts[-1])                  # BDI forward node times (crown excluded)
L_ext <- emphasis:::.extract_Ltable(tr)
pars  <- c(0.8, 0.4); mb <- c(0L, 0L, 0L)
S     <- 200L
cat(sprintf("tree: %d tips, crown age %.3f, first post-crown split at forward %.3f (%.1f%% of T)\n",
            Ntip(tr), Tc, s0, 100 * s0 / Tc))

is_aug <- function(df) df$t_ext != 0 & df$t_ext < 1e11   # augmented speciation rows

## ---------------------------------------------------------------- A: thinning
aug <- emphasis:::.augment_tree_internal(tr, pars, mb, sample_size = S,
                                          max_missing = 1e4, max_lambda = 500,
                                          num_threads = 1L)
dfs <- aug$trees
A <- t(sapply(dfs, function(df) {
  a  <- df[is_aug(df), , drop = FALSE]
  m1 <- a$parent_id == -1L
  ## parent node must be older than child; count nodes eligible vs lineages alive
  older_ok <- TRUE; elig <- integer(0); alive <- integer(0)
  for (k in seq_len(nrow(a))) {
    t <- a$brts[k]
    nodes_alive <- df[df$t_ext != 0 & df$brts < t & df$t_ext > t, , drop = FALSE]
    elig  <- c(elig, nrow(nodes_alive))
    alive <- c(alive, a$n[k])          # n = lineages alive on the segment ending at the node
    if (a$parent_id[k] >= 0L) {
      pb <- df$brts[df$id == a$parent_id[k] & df$t_ext != 0]
      older_ok <- older_ok && length(pb) == 1L && pb < t
    }
  }
  c(n_aug = nrow(a), n_m1 = sum(m1), born_before_s0 = sum(a$brts < s0),
    rule_exact = as.integer(all(m1 == (elig == 0L))),
    m1_all_before_s0 = as.integer(all(a$brts[m1] < s0)),
    parent_older = as.integer(older_ok),
    elig_eq_n_minus2 = as.integer(all(elig == alive - 2)))
}))
cat("\n[A] thinning, 200 draws:\n")
cat(sprintf("  augmented lineages total %d; parent_id == -1: %d (%.1f%%)\n",
            sum(A[, "n_aug"]), sum(A[, "n_m1"]), 100 * sum(A[, "n_m1"]) / sum(A[, "n_aug"])))
cat(sprintf("  draws with >=1 parent -1: %d/%d\n", sum(A[, "n_m1"] > 0), S))
cat(sprintf("  born before first post-crown split: %d; rule '-1 iff no node alive at birth' in all draws: %s; all -1 rows born before s0: %s\n",
            sum(A[, "born_before_s0"]), all(A[, "rule_exact"] == 1), all(A[, "m1_all_before_s0"] == 1)))
cat(sprintf("  parent node older than child in all draws: %s\n", all(A[, "parent_older"] == 1)))
cat(sprintf("  eligible parent nodes == alive lineages - 2 in all draws: %s\n",
            all(A[, "elig_eq_n_minus2"] == 1)))

## ---------------------------------------------------------------- B: tas coverage
B <- t(sapply(dfs, function(df) {
  a <- df[is_aug(df), , drop = FALSE]
  L <- emphasis:::.aug_to_Ltable(df, Tc, brts, L_ext)
  ph <- tryCatch(DDD::L2phylo(L, dropextinct = FALSE), error = function(e) NULL)
  c(n_aug = nrow(a), n_m1 = sum(a$parent_id == -1L),
    rows_added = nrow(L) - nrow(L_ext),
    tas_ok = as.integer(!is.null(ph)),
    tas_tips = if (is.null(ph)) NA_integer_ else Ntip(ph))
}))
cat("\n[B] thinning -> tas (same code path as .sim_tree_conditional):\n")
cat(sprintf("  rows_added == n_aug - n_m1 in all draws: %s\n",
            all(B[, "rows_added"] == B[, "n_aug"] - B[, "n_m1"])))
cat(sprintf("  tas tips == 20 + rows_added in all built draws: %s; L2phylo failures: %d\n",
            all(B[B[, "tas_ok"] == 1, "tas_tips"] == 20 + B[B[, "tas_ok"] == 1, "rows_added"]),
            sum(B[, "tas_ok"] == 0)))
cat(sprintf("  lineages in log_q but absent from tas: %d of %d\n",
            sum(B[, "n_m1"]), sum(B[, "n_aug"])))
sc <- emphasis:::.sim_tree_conditional(tr, pars, mb, n_trees = S, method = "thinning")
cat(sprintf("  .sim_tree_conditional(thinning): %d trees, %d NULL tas, log_q finite: %d\n",
            length(sc$trees), sum(sapply(sc$trees, is.null)), sum(is.finite(sc$log_q))))

## ---------------------------------------------------------------- C: BDI
bd <- emphasis:::.augment_tree_bdi(tr, pars, mb, sample_size = S)
C <- t(sapply(bd$trees, function(df) {
  a <- df[is_aug(df), , drop = FALSE]
  pfwd <- bt_s[a$parent_id + 1L]                       # parent node forward time
  det  <- pmax(0L, sapply(a$brts, function(b) sum(bt_s <= b) - 1L))
  L <- emphasis:::.aug_to_Ltable(df, Tc, brts, L_ext)
  new <- L[seq_len(nrow(L)) > nrow(L_ext), , drop = FALSE]
  par_age <- sapply(new[, 2], function(p) L[match(p, L[, 3]), 1])
  ph <- tryCatch(DDD::L2phylo(L, dropextinct = FALSE), error = function(e) NULL)
  c(n_aug = nrow(a), n_m1 = sum(a$parent_id == -1L),
    born_before_s0 = sum(a$brts < s0),
    parent_younger = sum(pfwd > a$brts),               # parent event after child's birth
    deterministic = as.integer(all(a$parent_id == det)),
    L_child_older_than_parent = sum(new[, 1] > par_age),
    tas_ok = as.integer(!is.null(ph)),
    neg_edges = if (is.null(ph)) NA_integer_ else sum(ph$edge.length < 0),
    min_edge = if (is.null(ph)) NA_real_ else min(ph$edge.length),
    tas_tips = if (is.null(ph)) NA_integer_ else Ntip(ph))
}))
cat("\n[C] BDI, 200 draws:\n")
cat(sprintf("  augmented lineages %d; parent_id == -1: %d; parent == most recent observed node in all draws: %s\n",
            sum(C[, "n_aug"]), sum(C[, "n_m1"]), all(C[, "deterministic"] == 1)))
cat(sprintf("  born before first post-crown split: %d; of which parent node YOUNGER than child: %d\n",
            sum(C[, "born_before_s0"]), sum(C[, "parent_younger"])))
cat(sprintf("  L-table rows with child birth age > parent birth age: %d (in %d draws)\n",
            sum(C[, "L_child_older_than_parent"]), sum(C[, "L_child_older_than_parent"] > 0)))
cat(sprintf("  L2phylo failures: %d; built trees with negative edge lengths: %d (%d edges, min edge %.3f); tas tips == 20 + n_aug: %s\n",
            sum(C[, "tas_ok"] == 0), sum(C[, "neg_edges"] > 0, na.rm = TRUE),
            sum(C[, "neg_edges"], na.rm = TRUE), min(C[, "min_edge"], na.rm = TRUE),
            all(C[C[, "tas_ok"] == 1, "tas_tips"] == 20 + C[C[, "tas_ok"] == 1, "n_aug"])))
scb <- emphasis:::.sim_tree_conditional(tr, pars, mb, n_trees = S, method = "bdi")
cat(sprintf("  .sim_tree_conditional(bdi): %d trees, %d NULL tas\n",
            length(scb$trees), sum(sapply(scb$trees, is.null))))
## ape's own check on a BDI tas with a younger parent
bad <- which(C[, "parent_younger"] > 0 & C[, "tas_ok"] == 1)[1]
if (!is.na(bad)) {
  L  <- emphasis:::.aug_to_Ltable(bd$trees[[bad]], Tc, brts, L_ext)
  ph <- DDD::L2phylo(L, dropextinct = FALSE)
  cat(sprintf("  example draw %d: min edge length %.4f, is.ultrametric(prune_to_extant) %s, prune_to_extant tips %d\n",
              bad, min(ph$edge.length), is.ultrametric(emphasis:::prune_to_extant(ph)),
              Ntip(emphasis:::prune_to_extant(ph))))
  ## is the published tas from .sim_tree_conditional(bdi) affected the same way?
  ne <- sapply(scb$trees, function(ph) if (is.null(ph)) NA else sum(ph$edge.length < 0))
  cat(sprintf("  .sim_tree_conditional(bdi) tas with negative edges: %d/%d\n", sum(ne > 0, na.rm = TRUE), S))
}

## ---------------------------------------------------------------- D: does parent_id reach logf/logg?
strip <- function(df) { df$parent_id <- -1L; df$tip_start <- 0; df$focal_tip_start <- 0; df$pd <- 0; df }
cmp <- function(pars8, model) {
  e1 <- emphasis:::eval_logf(pars8, dfs,            model = model, link = 0L, rho = 1)
  e2 <- emphasis:::eval_logf(pars8, lapply(dfs, strip), model = model, link = 0L, rho = 1)
  c(dlogf = max(abs(e1$logf - e2$logf)), dlogg = max(abs(e1$logg - e2$logg)))
}
cat("\n[D] max |change| in logf/logg after erasing parent_id/tip_start/pd (thinning draws):\n")
cat("  cr  (0,0,0): "); print(cmp(c(0.8, 0, 0, 0, 0.4, 0, 0, 0), c(0L, 0L, 0L)))
cat("  dd  (1,0,0): "); print(cmp(c(0.8, -0.01, 0, 0, 0.4, 0, 0, 0), c(1L, 0L, 0L)))
cat("  M   (0,1,0): "); print(cmp(c(0.8, 0, 0.05, 0, 0.4, 0, 0, 0), c(0L, 1L, 0L)))
cat("  D   (0,0,1): "); print(cmp(c(0.8, 0, 0, 0.05, 0.4, 0, 0, 0), c(0L, 0L, 1L)))

## ---------------------------------------------------------------- E: does any fit driver consume tas?
uses <- function(f) any(grepl("aug_to_Ltable|L2phylo|\\$tas", deparse(f)))
cat("\n[E] fit drivers reference .aug_to_Ltable/L2phylo/$tas:",
    "mcem_bdi =", uses(emphasis:::.mcem_bdi),
    "; mcem_dynamic_fresh =", uses(emphasis:::.mcem_dynamic_fresh),
    "; run_mcem =", uses(emphasis:::.run_mcem),
    "; emphasis_cem =", uses(emphasis:::emphasis_cem), "\n")

## ---------------------------------------------------------------- F: supplements
## F1: parent_id ALONE (derived columns untouched) under the D model
only_pid <- function(df) { df$parent_id[is_aug(df)] <- -1L; df }
pD <- c(0.8, 0, 0, 0.05, 0.4, 0, 0, 0)
e1 <- emphasis:::eval_logf(pD, dfs, model = c(0L, 0L, 1L), link = 0L, rho = 1)
e2 <- emphasis:::eval_logf(pD, lapply(dfs, only_pid), model = c(0L, 0L, 1L), link = 0L, rho = 1)
cat(sprintf("[F1] D model, parent_id -> -1 only: max|dlogf| = %.4f, max|dlogg| = %.4f\n",
            max(abs(e1$logf - e2$logf)), max(abs(e1$logg - e2$logg))))
## F2: is the L_extant[1,3] fallback in .aug_to_Ltable ever reached (thinning)?
n_obs <- length(brts) - 1L
fb <- sum(sapply(dfs, function(df) {
  a <- df[is_aug(df) & df$parent_id >= 0L, , drop = FALSE]
  sum(a$parent_id == n_obs)            # closing sentinel id: the only id mapping to NA
}))
cat(sprintf("[F2] thinning rows whose parent id maps to the NA fallback: %d\n", fb))
## F3: crown-lineage exclusion — the thinning proposal never names a crown lineage as parent
cat(sprintf("[F3] crown lineages are never eligible: eligible == n - 2 held in all draws (see [A]); mean n at augmented births = %.2f, so P(crown lineage parent) = 0 vs 2/n = %.2f under exchangeability\n",
            mean(unlist(lapply(dfs, function(df) df$n[is_aug(df)]))),
            mean(unlist(lapply(dfs, function(df) 2 / df$n[is_aug(df)])))))
