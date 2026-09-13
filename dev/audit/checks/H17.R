## H17 — .ess_from_lw drops non-finite lw; can -Inf reach it?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
cat("ess_from_lw(c(0,-Inf,-Inf)) =", emphasis:::.ess_from_lw(c(0, -Inf, -Inf)), "\n")
cat("ess_from_lw(c(0,-Inf,-Inf)) w/o filter =", { lw <- c(0,-Inf,-Inf); w <- exp(lw - max(lw)); sum(w)^2/sum(w^2) }, "\n")
cat("ess_from_lw(c(-Inf,-Inf)) =", emphasis:::.ess_from_lw(c(-Inf, -Inf)), "\n")
## thinning path: E_step keeps only finite lw -> ESS input always finite
## BDI path: lw = logf(eval_logf) - logg(R); logf can be -Inf when lambda hits 0 (dd linear).
set.seed(5)
phy <- ape::drop.fossil(ape::rlineage(0.6, 0.1, Tmax = 5))
while (ape::Ntip(phy) < 10 || ape::Ntip(phy) > 30) phy <- ape::drop.fossil(ape::rlineage(0.6, 0.1, Tmax = 5))
brts <- sort(ape::branching.times(phy), decreasing = TRUE)
n <- ape::Ntip(phy)
## dd linear with K just above n so augmented lineages push lambda to 0
lambda0 <- 0.6; K <- n + 2
pars <- c(lambda0, -lambda0 / K, 0, 0, 0.1, 0, 0, 0)
inp <- tryCatch({
  ev <- emphasis:::.augment_tree_bdi
  NULL
}, error = function(e) NULL)
## use the BDI E-step directly if exposed; else emulate via thinning trees scored at these pars
raw <- emphasis:::augment_trees(as.numeric(brts), c(lambda0, 0,0,0, 0.1,0,0,0), sample_size = 200L, maxN = 20000L,
                     max_missing = 100L, max_lambda = 1e3, num_threads = 1L,
                     model = c(1L,0L,0L), link = 0L, rho = 1)
ev <- emphasis:::eval_logf(pars, raw$trees, model = c(1L,0L,0L), link = 0L, rho = 1)
lw <- ev$logf - ev$logg
cat(sprintf("n=%d tips, K=%d: scored %d trees at dd-linear pars; logf=-Inf in %d, lw non-finite in %d\n",
            n, K, length(lw), sum(!is.finite(ev$logf)), sum(!is.finite(lw))))
cat(sprintf("ESS (filter) = %.2f ; ESS (all, -Inf as 0 weight) = %.2f ; n_finite = %d\n",
            emphasis:::.ess_from_lw(lw), { w <- exp(lw - max(lw[is.finite(lw)])); sum(w)^2/sum(w^2) }, sum(is.finite(lw))))
## BDI fhat denominator uses length(weights) (bdi.R:699) i.e. includes -Inf; ESS drops them
cat("BDI fhat convention (bdi.R:696-699): log(sum exp(lw - max)/length(lw)) + max -> -Inf trees count in denominator\n")
