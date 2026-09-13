## H16 — mu = 0 at an extinction node vs lambda = 0 at a speciation node in eval_logf.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
set.seed(3)
phy  <- ape::rlineage(0.5, 0.2, Tmax = 6)
phy  <- ape::drop.fossil(phy)
while (ape::Ntip(phy) < 8 || ape::Ntip(phy) > 40) { phy <- ape::drop.fossil(ape::rlineage(0.5, 0.2, Tmax = 6)) }
brts <- sort(ape::branching.times(phy), decreasing = TRUE)
cat(sprintf("tree: %d tips, crown %.2f\n", ape::Ntip(phy), brts[1]))
pars <- c(0.5, 0, 0, 0, 0.3, 0, 0, 0)
raw <- emphasis:::augment_trees(as.numeric(brts), pars, sample_size = 50L, maxN = 5000L,
                     max_missing = 200L, max_lambda = 1e3, num_threads = 1L,
                     model = c(0L,0L,0L), link = 0L, rho = 1)
n_ext <- vapply(raw$trees, function(tr) sum(tr$t_ext == 0), integer(1))
cat("trees with >=1 extinction node:", sum(n_ext > 0), "of", length(raw$trees), "\n")
tr_ext <- raw$trees[n_ext > 0][1:min(5, sum(n_ext > 0))]
## (a) mu = 0, extinction node present
pm0 <- pars; pm0[5] <- 0
ev <- emphasis:::eval_logf(pm0, tr_ext, model = c(0L,0L,0L), link = 0L, rho = 1)
cat("\n(a) mu = 0 with extinction nodes:\n"); print(data.frame(n_ext = n_ext[n_ext>0][1:length(tr_ext)], logf = ev$logf, logg = ev$logg, lw = ev$logf - ev$logg))
cat(sprintf("    per-extinction log(max(mu,1e-300)) = %.2f; lw finite: %s -> kept by E_step filter (isfinite && exp(lw)>0): %s\n",
            log(1e-300), all(is.finite(ev$logf - ev$logg)), all(is.finite(ev$logf - ev$logg) & exp(ev$logf - ev$logg) > 0)))
## (b) lambda = 0 at every speciation node (beta_0 = 0, linear link)
pl0 <- pars; pl0[1] <- 0
ev2 <- emphasis:::eval_logf(pl0, tr_ext, model = c(0L,0L,0L), link = 0L, rho = 1)
cat("\n(b) lambda = 0:\n"); print(data.frame(logf = ev2$logf, logg = ev2$logg, lw = ev2$logf - ev2$logg))
cat("    kept by filter:", all(is.finite(ev2$logf - ev2$logg)), "\n")
## (c) does mu=0 ever produce a *rejected* tree in the sampler itself? run augment_trees at mu=0
raw0 <- tryCatch(emphasis:::augment_trees(as.numeric(brts), pm0, sample_size = 50L, maxN = 5000L,
                     max_missing = 200L, max_lambda = 1e3, num_threads = 1L,
                     model = c(0L,0L,0L), link = 0L, rho = 1), error = function(e) e)
if (inherits(raw0, "error")) { cat("\n(c) augment_trees at mu=0 errored:", conditionMessage(raw0), "\n") } else {
  ne0 <- vapply(raw0$trees, function(tr) sum(tr$t_ext == 0), integer(1))
  cat(sprintf("\n(c) augment_trees at mu=0: %d trees, %d with extinction nodes, rejected_zero_weights=%d, lw range [%.2f, %.2f]\n",
      length(raw0$trees), sum(ne0 > 0), raw0$rejected_zero_weights, min(raw0$logf-raw0$logg), max(raw0$logf-raw0$logg)))
}
