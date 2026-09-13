.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
set.seed(7); phy <- ape::rphylo(24, birth = 0.6, death = 0.2)
for (m in c("bdi", "thinning")) {
  st <- simulate_tree(phy, pars = c(0.6, 0.2), model = "cr", n_trees = 40L, method = m, rho = 0.7)
  tr1 <- st$trees[[1]]; cat(m, ": class of tree =", class(tr1), "\n")
  extant <- sapply(st$trees, function(tr) {
    if (inherits(tr, "phylo")) { d <- ape::node.depth.edgelength(tr); tips <- d[1:ape::Ntip(tr)]; sum(abs(tips - max(d)) < 1e-6) }
    else if (is.data.frame(tr)) sum(tr$t_ext == 5e10) + 24 else NA })
  cat("   extant tips per augmented tree: range", range(extant), "; trees with > 24 extant:", sum(extant > 24), "of", length(extant), "\n")
}
