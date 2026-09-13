.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
tryit <- function(pars, lab) {
  cat("---", lab, "pars =", pars, "\n")
  r <- tryCatch(simulate_tree(pars = pars, max_t = 5, model = "cr", link = "gaussian", max_tries = 20), error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r)) { cat(r, "\n"); return(invisible()) }
  str(r, max.level = 1)
  if (!is.null(r$tree)) cat("Ntip =", ape::Ntip(r$tree), "\n")
  if (!is.null(r$Ltable)) cat("Ltable rows =", nrow(r$Ltable), "\n")
}
tryit(c(0.5, 0.1), "positive intercept (control)")
tryit(c(-0.5, 0.1), "negative lambda intercept, lam+mu<0")
tryit(c(-0.5, 1.0), "negative lambda intercept, lam+mu>0 (bernoulli p=-1)")
tryit(c(0.5, -0.1), "negative mu intercept, lam+mu>0 (bernoulli p>1)")
