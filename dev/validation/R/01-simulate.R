# ---------------------------------------------------------------------------
# 01-simulate.R — tree generation.
#
#   Rscript 01-simulate.R --tier smoke|main|ext [--lib /path/to/rlib]
#
# Writes data/trees-<tier>.rds: a list of tree records, each
#   tree_id, kind ("cr"/"dd"), cell, n, brts (decreasing, crown first),
#   phylo, crown age T, generating parameters, seed, generator call.
#
# Deterministic: seed = digest2int(tree_id).  Nothing here depends on the
# emphasis build, so trees are shared by every tier and every build.
#
# CR: TreeSim::sim.bd.taxa(n, lambda, mu, complete = FALSE) conditions on n,
#     so the crown age is random and recorded.  That conditioning affects the
#     recovery arm (A4) only and is stated there; the primary estimand is the
#     per-tree exact MLE, which is unaffected.
# DD: DDD::dd_sim(pars, age, ddmodel = 1)$tes, accepted only when
#     n in [0.5 K, 1.5 K] so that K stays identifiable and dd_ML is a usable
#     reference; rejected draws are counted, the seed is advanced.
# ---------------------------------------------------------------------------

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f[1])) else
    "/Users/pancho/Code/emphasis/dev/validation/R"
  source(file.path(d, "00-common.R"), local = FALSE)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, default = NULL) {
  i <- which(args == k); if (length(i)) args[i[1] + 1L] else default
}
TIER <- getarg("--tier", "smoke")
if (!is.null(getarg("--lib"))) options(emphasis.lib = getarg("--lib"))

val_load_refs()
suppressPackageStartupMessages(library(digest))

cat(sprintf("[01-simulate] tier = %s\n", TIER))
t_start <- Sys.time()
trees <- list()

# --- CR --------------------------------------------------------------------
cells <- val_cr_cells(TIER)
for (i in seq_len(nrow(cells))) {
  cc <- cells[i, ]
  lam <- cc$lam; mu <- cc$eps * cc$lam
  cell <- sprintf("cr-n%03d-e%02d-l%03d", cc$n, round(10 * cc$eps),
                  round(100 * cc$lam))
  for (k in seq_len(cc$trees)) {
    tid <- sprintf("%s-t%02d", cell, k)
    set.seed(val_seed(tid))
    tr <- TreeSim::sim.bd.taxa(n = cc$n, numbsim = 1, lambda = lam, mu = mu,
                               complete = FALSE)[[1]]
    brts <- sort(as.numeric(ape::branching.times(tr)), decreasing = TRUE)
    trees[[tid]] <- list(
      tree_id = tid, kind = "cr", cell = cell, n = cc$n, eps = cc$eps,
      lambda_gen = lam, mu_gen = mu, brts = brts, T_crown = brts[1L],
      phylo = tr, seed = val_seed(tid), index = k,
      gen_call = sprintf("TreeSim::sim.bd.taxa(n=%d, lambda=%g, mu=%g, complete=FALSE)",
                         cc$n, lam, mu))
  }
}
cat(sprintf("  cr: %d trees in %d cells\n", length(trees), nrow(cells)))

# --- DD --------------------------------------------------------------------
ddcells <- val_dd_cells(TIER)
n_dd_rejected <- 0L
if (nrow(ddcells)) for (i in seq_len(nrow(ddcells))) {
  cc <- ddcells[i, ]
  cell <- sprintf("dd-%s", cc$regime)
  lo <- VAL_DD_BAND[1] * cc$K
  hi <- min(VAL_DD_BAND[2] * cc$K, VAL_DD_NMAX)
  for (k in seq_len(cc$trees)) {
    tid <- sprintf("%s-t%02d", cell, k)
    accepted <- FALSE
    for (attempt in seq_len(60L)) {
      sd_k <- val_seed(paste0(tid, "-a", attempt))
      set.seed(sd_k)
      sim <- .q(DDD::dd_sim(pars = c(cc$lambda0, cc$mu0, cc$K),
                            age = cc$age, ddmodel = 1))
      tes <- sim$tes
      if (is.null(tes) || !inherits(tes, "phylo")) { n_dd_rejected <- n_dd_rejected + 1L; next }
      nt <- ape::Ntip(tes)
      if (nt < lo || nt > hi) { n_dd_rejected <- n_dd_rejected + 1L; next }
      brts <- sort(as.numeric(ape::branching.times(tes)), decreasing = TRUE)
      trees[[tid]] <- list(
        tree_id = tid, kind = "dd", cell = cell, regime = cc$regime,
        n = nt, lambda0_gen = cc$lambda0, mu0_gen = cc$mu0, K_gen = cc$K,
        age = cc$age, stress = cc$stress,
        pars_gen = as.numeric(val_dd_to_emphasis(cc$lambda0, cc$mu0, cc$K)),
        brts = brts, T_crown = brts[1L], phylo = tes, seed = sd_k,
        index = k, attempts = attempt,
        gen_call = sprintf("DDD::dd_sim(pars=c(%g,%g,%g), age=%g, ddmodel=1)",
                           cc$lambda0, cc$mu0, cc$K, cc$age))
      accepted <- TRUE; break
    }
    if (!accepted)
      warning(sprintf("dd tree %s not accepted in 60 attempts (n outside [%.0f, %.0f])",
                      tid, lo, hi), call. = FALSE)
  }
}
cat(sprintf("  dd: %d trees in %d regimes (%d draws rejected on size)\n",
            sum(vapply(trees, function(x) x$kind == "dd", TRUE)),
            nrow(ddcells), n_dd_rejected))

meta <- list(tier = TIER, created = Sys.time(),
             cr_cells = cells, dd_cells = ddcells,
             dd_rejected = n_dd_rejected,
             r_version = as.character(getRversion()),
             treesim = as.character(utils::packageVersion("TreeSim")),
             ddd = as.character(utils::packageVersion("DDD")))
out <- file.path(VAL_DIR$data, sprintf("trees-%s.rds", TIER))
saveRDS(list(trees = trees, meta = meta), out)

sz <- vapply(trees, function(x) x$n, 1L)
cat(sprintf("  n range %d-%d; wrote %d trees to %s (%.1fs)\n",
            min(sz), max(sz), length(trees), out,
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
