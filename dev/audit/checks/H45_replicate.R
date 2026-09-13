## H45 replication — vary tree size, rates, rho, model, and the public API path.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })

is_aug <- function(df) df$t_ext != 0 & df$t_ext < 1e11

summ_thin <- function(tr, pars, mb, S, rho = 1, link = 0L, tag) {
  brts <- emphasis:::.extract_brts(tr); Tc <- brts[1]; s0 <- Tc - brts[2]
  L_ext <- emphasis:::.extract_Ltable(tr); n_obs <- length(brts) - 1L
  aug <- emphasis:::.augment_tree_internal(tr, pars, mb, sample_size = S,
                                            max_missing = 1e4, max_lambda = 500,
                                            num_threads = 1L, link = link, rho = rho)
  dfs <- aug$trees
  M <- t(sapply(dfs, function(df) {
    a <- df[is_aug(df), , drop = FALSE]
    m1 <- a$parent_id == -1L
    ## rule: -1 iff no non-extinction node alive at birth
    elig <- sapply(a$brts, function(t) sum(df$t_ext != 0 & df$brts < t & df$t_ext > t))
    ## children whose parent is a dropped (-1) augmented lineage -> fallback reached
    m1_ids <- a$id[m1]
    child_of_m1 <- sum(a$parent_id %in% m1_ids)
    L <- emphasis:::.aug_to_Ltable(df, Tc, brts, L_ext)
    ph <- tryCatch(DDD::L2phylo(L, dropextinct = FALSE), error = function(e) NULL)
    c(n_aug = nrow(a), n_m1 = sum(m1), rule = as.integer(all(m1 == (elig == 0))),
      elig_n2 = as.integer(all(elig == a$n - 2)),
      child_of_m1 = child_of_m1,
      rows_added = nrow(L) - nrow(L_ext),
      neg = if (is.null(ph)) NA else sum(ph$edge.length < 0),
      tips = if (is.null(ph)) NA else Ntip(ph), fail = as.integer(is.null(ph)))
  }))
  cat(sprintf("[thin %s] %d draws (%d valid): aug %d, -1: %d (%.1f%%), rule -1<=>no node alive: %s, elig==n-2: %s, children of -1 lineages (fallback L_extant[1,3]): %d, rows_added==n_aug-n_m1: %s, tas neg edges: %d trees, L2phylo fail: %d\n",
              tag, S, length(dfs), sum(M[,"n_aug"]), sum(M[,"n_m1"]), 100*sum(M[,"n_m1"])/sum(M[,"n_aug"]),
              all(M[,"rule"]==1), all(M[,"elig_n2"]==1), sum(M[,"child_of_m1"]),
              all(M[,"rows_added"] == M[,"n_aug"] - M[,"n_m1"]),
              sum(M[,"neg"] > 0, na.rm = TRUE), sum(M[,"fail"])))
  invisible(M)
}

summ_bdi <- function(tr, pars, mb, S, link = 0L, tag) {
  brts <- emphasis:::.extract_brts(tr); Tc <- brts[1]; s0 <- Tc - brts[2]
  bt_s <- sort(Tc - brts[-1]); L_ext <- emphasis:::.extract_Ltable(tr)
  stopifnot(emphasis:::.bdi_supported(mb, link))
  bd <- emphasis:::.augment_tree_bdi(tr, pars, mb, sample_size = S, link = link)
  M <- t(sapply(bd$trees, function(df) {
    a <- df[is_aug(df), , drop = FALSE]
    pf <- bt_s[a$parent_id + 1L]
    L <- emphasis:::.aug_to_Ltable(df, Tc, brts, L_ext)
    ph <- tryCatch(DDD::L2phylo(L, dropextinct = FALSE), error = function(e) NULL)
    c(n_aug = nrow(a), n_m1 = sum(a$parent_id == -1L), before_s0 = sum(a$brts < s0),
      younger = sum(pf > a$brts),
      neg = if (is.null(ph)) NA else sum(ph$edge.length < 0),
      minedge = if (is.null(ph)) NA else min(ph$edge.length), fail = as.integer(is.null(ph)))
  }))
  cat(sprintf("[bdi  %s] %d draws (%d valid): aug %d, -1: %d, born before first split: %d, parent node younger than child: %d, tas with neg edges: %d/%d (min edge %.3f), L2phylo fail: %d\n",
              tag, S, length(bd$trees), sum(M[,"n_aug"]), sum(M[,"n_m1"]), sum(M[,"before_s0"]),
              sum(M[,"younger"]), sum(M[,"neg"] > 0, na.rm = TRUE), length(bd$trees),
              min(M[,"minedge"], na.rm = TRUE), sum(M[,"fail"])))
  invisible(M)
}

set.seed(7)
trA <- TreeSim::sim.bd.taxa(n = 40, numbsim = 1, lambda = 1.0, mu = 0.5, complete = FALSE)[[1]]
trB <- TreeSim::sim.bd.taxa(n = 10, numbsim = 1, lambda = 0.5, mu = 0.1, complete = FALSE)[[1]]
for (nm in c("trA", "trB")) {
  tr <- get(nm); b <- emphasis:::.extract_brts(tr)
  cat(sprintf("%s: %d tips, crown %.3f, first post-crown split at %.1f%% of T\n",
              nm, Ntip(tr), b[1], 100 * (b[1] - b[2]) / b[1]))
}

## --- thinning: different tree/rates/rho/link
summ_thin(trA, c(1.0, 0.5), c(0L,0L,0L), 100, tag = "A cr(1,.5)")
summ_thin(trA, c(1.0, 0.9), c(0L,0L,0L), 60,  tag = "A cr(1,.9) high mu")
summ_thin(trB, c(0.5, 0.1), c(0L,0L,0L), 100, tag = "B cr(.5,.1)")
summ_thin(trB, c(0.5, 0.4), c(0L,0L,0L), 100, rho = 0.7, tag = "B cr rho=.7")
summ_thin(trA, c(1.2, -0.01, 0.5, 0), c(1L,0L,0L), 60, link = 1L, tag = "A dd exp link")

## --- BDI: different tree/rates/model/link
summ_bdi(trA, c(1.0, 0.5), c(0L,0L,0L), 100, tag = "A cr(1,.5)")
summ_bdi(trB, c(0.5, 0.1), c(0L,0L,0L), 100, tag = "B cr(.5,.1)")
summ_bdi(trB, c(0.5, 0.4), c(0L,0L,0L), 100, tag = "B cr(.5,.4)")
summ_bdi(trA, c(1.2, -0.01, 0.5, 0), c(1L,0L,0L), 60, link = 0L, tag = "A dd linear")

## --- public API, default method
pub <- simulate_tree(tree = trA, pars = c(1.0, 0.5), model = "cr", n_trees = 100L)
neg <- sapply(pub$trees, function(p) if (is.null(p)) NA else sum(p$edge.length < 0))
cat(sprintf("[public simulate_tree(tree=trA, default method)] %d trees, NULL %d, with negative edges %d, finite log_q %d\n",
            length(pub$trees), sum(is.na(neg)), sum(neg > 0, na.rm = TRUE), sum(is.finite(pub$log_q))))
pubt <- simulate_tree(tree = trA, pars = c(1.0, 0.5), model = "cr", n_trees = 100L, method = "thinning")
negt <- sapply(pubt$trees, function(p) if (is.null(p)) NA else sum(p$edge.length < 0))
cat(sprintf("[public simulate_tree(tree=trA, thinning)] %d trees, NULL %d, with negative edges %d\n",
            length(pubt$trees), sum(is.na(negt)), sum(negt > 0, na.rm = TRUE)))

## --- is the extant part of a negative-edge BDI tas still the input tree?
bad <- which(neg > 0)[1]
if (!is.na(bad)) {
  pe <- emphasis:::prune_to_extant(pub$trees[[bad]])
  cat(sprintf("  example: min edge %.3f; prune_to_extant tips %d, ultrametric %s, RF dist to input = %d\n",
              min(pub$trees[[bad]]$edge.length), Ntip(pe), is.ultrametric(pe),
              phangorn_free_rf <- tryCatch(ape::dist.topo(unroot(pe), unroot(trA)), error = function(e) NA)))
}
