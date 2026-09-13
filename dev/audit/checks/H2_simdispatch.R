.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
set.seed(7); phy <- ape::rphylo(24, birth = 0.6, death = 0.2)
for (m in c("bdi", "thinning")) {
  st <- simulate_tree(phy, pars = c(0.6, 0.2), model = "cr", n_trees = 30L, method = m, rho = 0.7)
  cat(m, ": names =", paste(names(st), collapse = ","), "\n")
  trs <- if (!is.null(st$trees)) st$trees else st
  cnt <- tryCatch(sapply(trs, function(tr) if (is.data.frame(tr)) sum(tr$t_ext == 5e10) else NA), error = function(e) NA)
  cat("   unsampled-extant nodes per tree (max): ", suppressWarnings(max(cnt, na.rm = TRUE)), "\n")
  print(utils::str(st, max.level = 1))
}
