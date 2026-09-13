## H9_replicate.R -- independent replication of H9 with things the verifier
## did not vary: model `d` (not `nd`), mu > 0 (augmented lineages present),
## exponential-link generated tree, invariance checked on logg as well as logf,
## and the fix-sketch premise that the phylo carries E_focal at observed splits.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis)); suppressMessages(library(ape))
set.seed(91)

aug <- function(brts, pars8, mb, link, N = 1L, max_missing = 0L, rho = 1) {
  a <- emphasis:::augment_trees(brts, pars = as.numeric(pars8), sample_size = as.integer(N),
                                maxN = 500L, max_missing = as.integer(max_missing),
                                max_lambda = 1e6, num_threads = 1L,
                                model = as.integer(mb), link = as.integer(link), rho = rho)
  a$trees
}
ev <- function(pars8, dfs, mb, link, rho = 1)
  emphasis:::eval_logf(as.numeric(pars8), dfs, model = as.integer(mb),
                       link = as.integer(link), rho = rho)

cat("==== V1: `d` model, exp link, mu > 0: is linear-link logf flat in beta_D iff n_aug == 0? ====\n")
## d model compact pars: (b0, bD, g0, gD); exponential link; mu = exp(g0)
pars_d <- c(log(0.8), 0.3, log(0.25), 0)
sim <- NULL
for (k in 1:200) {
  s <- simulate_tree(pars = pars_d, max_t = 4, model = "d", link = "exponential",
                     max_lin = 300, useDDD = FALSE)
  if (s$status == "done" && sum(s$L[, 4] == -1) >= 12 && sum(s$L[, 4] == -1) <= 40 && sum(s$L[, 4] != -1) >= 2) { sim <- s; break }
}
stopifnot(!is.null(sim))
tes1 <- DDD::L2phylo(sim$L, dropextinct = TRUE)
brts <- sort(branching.times(tes1), decreasing = TRUE)
cat(sprintf("tree: %d extant tips (L-table rows %d, extinct %d), crown age %.3f\n",
            Ntip(tes1), nrow(sim$L), sum(sim$L[, 4] != -1), brts[1]))
mb_d <- c(0, 0, 1)
## fit-side pars, linear link (the claim is about the linear link's flatness)
p_lin <- function(bD) emphasis:::.expand_pars(c(0.6, bD, 0.15, 0), mb_d)
trees <- aug(brts, p_lin(0), mb_d, link = 0, N = 60L, max_missing = 30L)
trees <- Filter(function(d) nrow(d) > 0, trees)
naug <- sapply(trees, function(d) sum(d$t_ext == 0))          # extinction nodes = augmented lineages (incl. parent_id == -1 ones)
naug_p <- sapply(trees, function(d) sum(d$parent_id >= 0 & d$t_ext != 0))
cat("draws where an augmented lineage has parent_id == -1 (H48 case):", sum(naug != naug_p), "\n")
cat("n_aug distribution over 60 draws:", paste(names(table(naug)), table(naug), sep = ":", collapse = " "), "\n")
bDs <- c(-0.4, 0, 0.3, 0.7)
rng <- t(sapply(trees, function(d) {
  v <- sapply(bDs, function(b) ev(p_lin(b), list(d), mb_d, 0)$logf)
  c(range_logf = diff(range(v)), any_inf = any(!is.finite(v)))
}))
flat0 <- rng[naug == 0, "range_logf"]; flat1 <- rng[naug > 0, "range_logf"]
cat(sprintf("draws with n_aug = 0: %d, max range of logf over beta_D = %.3g\n", length(flat0), if (length(flat0)) max(flat0) else NA))
cat(sprintf("draws with n_aug > 0: %d, min range of logf over beta_D = %.3g, median = %.3g\n",
            length(flat1), if (length(flat1)) min(flat1) else NA, if (length(flat1)) median(flat1) else NA))
## Also the observed-node covariate on this tree: D from e_s - M must be 0 for all parent_id == -1
d0 <- trees[[which(naug == 0)[1]]]          # a draw with no augmented lineage
obs <- d0[d0$parent_id == -1 & d0$t_ext > 1e10, ]
cat("observed nodes: tip_start all 0:", all(obs$tip_start == 0),
    "; pd == (rank)*brts (crown lineages absent):", isTRUE(all.equal(obs$pd, seq_len(nrow(obs)) * obs$brts, tol = 1e-8)), "\n")

cat("\n==== V2: cr / dd with augmented lineages actually present, rho = 1: logf AND logg invariant to pd/tip_start? ====\n")
for (mbn in c("cr", "dd")) for (link in 0:2) {
  mb <- if (mbn == "cr") c(0, 0, 0) else c(1, 0, 0)
  pc <- switch(paste(mbn, link),
               "cr 0" = c(0.7, 0.35), "cr 1" = c(log(0.7), log(0.35)), "cr 2" = c(0.7, 0.35),
               "dd 0" = c(0.9, -0.01, 0.35, 0), "dd 1" = c(log(0.9), -0.01, log(0.35), 0), "dd 2" = c(0.9, 0.03, 0.35, 0.02))
  p8 <- emphasis:::.expand_pars(pc, mb)
  tr <- Filter(function(d) nrow(d) > 0 && sum(d$parent_id >= 0 & d$t_ext != 0) >= 2,
               aug(brts, p8, mb, link, N = 40L, max_missing = 40L))
  if (!length(tr)) { cat(sprintf("%s link=%d: no draw with >=2 augmented lineages\n", mbn, link)); next }
  d0 <- tr[[1]]; na <- sum(d0$parent_id >= 0 & d0$t_ext != 0)
  d1 <- d0; d1$pd <- 1e3 * runif(nrow(d0)); d1$tip_start <- runif(nrow(d0)); d1$focal_tip_start <- runif(nrow(d0))
  e0 <- ev(p8, list(d0), mb, link); e1 <- ev(p8, list(d1), mb, link)
  cat(sprintf("%s link=%d n_aug=%d: logf %.6f vs %.6f (same=%s); logg %.6f vs %.6f (same=%s)\n",
              mbn, link, na, e0$logf, e1$logf, isTRUE(all.equal(e0$logf, e1$logf, tol = 1e-12)),
              e0$logg, e1$logg, isTRUE(all.equal(e0$logg, e1$logg, tol = 1e-12))))
}

cat("\n==== V3: does the reconstructed phylo carry E_focal at each observed split? (mu = 0 tree) ====\n")
pars_nd <- c(0.9, -0.02, 0.25, 0, 0, 0)
sim2 <- NULL
for (k in 1:100) {
  s <- simulate_tree(pars = pars_nd, max_t = 3, model = "nd", link = "linear", max_lin = 200, useDDD = FALSE)
  if (s$status == "done" && nrow(s$L) >= 8 && nrow(s$L) <= 40) { sim2 <- s; break }
}
L <- sim2$L; Tc <- 3; L_birth <- Tc - L[, 1]; o <- order(L_birth); L <- L[o, ]; L_birth <- L_birth[o]
ts_at <- function(lab, t) max(c(L_birth[L[, 3] == lab], L_birth[L[, 2] == lab & L_birth < t]))
ev_t <- L_birth[-(1:2)]; ev_par <- L[-(1:2), 2]
E_true <- mapply(function(t, p) t - ts_at(p, t), ev_t, ev_par)
## phylo side: for each internal node, the edge above it (0 for the root)
tes <- DDD::L2phylo(sim2$L, dropextinct = TRUE); nt <- Ntip(tes); nh <- node.depth.edgelength(tes)
int_nodes <- (nt + 1):(nt + tes$Nnode)
above <- sapply(int_nodes, function(v) { e <- which(tes$edge[, 2] == v); if (length(e)) tes$edge.length[e] else 0 })
fwd_t <- nh[int_nodes]                     # forward time of each internal node (root = 0)
ord <- order(fwd_t); above <- above[ord]; fwd_t <- fwd_t[ord]
## drop the root (crown, fwd 0): remaining internal nodes are the observed events
cat(sprintf("events: %d (L) vs %d (phylo)\n", length(ev_t), length(fwd_t) - 1))
cat("event times match:", isTRUE(all.equal(sort(ev_t), fwd_t[-1], tol = 1e-5)), "\n")
cat("edge-above-node == E_focal(true) at every observed split:",
    isTRUE(all.equal(above[-1][order(fwd_t[-1])], E_true[order(ev_t)], tol = 1e-5)), "\n")
cat("DONE\n")
