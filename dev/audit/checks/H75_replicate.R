# H75 replicate: vary tree source, model, method, useDDD, pars matrix, integer/named input.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
run <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) e)
  cat(sprintf("%-55s -> %s\n", label,
      if (inherits(r, "error")) paste("ERROR:", conditionMessage(r))
      else if (!is.null(r$log_q)) sprintf("OK log_q=%s tas=%s", paste(round(r$log_q,2), collapse=","), if (is.null(r$tas)) "NULL" else class(r$tas)[1])
      else if (is.list(r)) sprintf("OK list of %d", length(r)) else "OK"))
  invisible(r)
}
phy <- ape::rphylo(20, 0.6, 0.15)
brts <- ape::branching.times(phy)                     # named numeric, as ape returns it
run("named numeric brts (20 tips), bdi, cr", simulate_tree(tree = brts, pars = c(0.6, 0.15)))
run("unnamed sorted numeric brts, thinning, cr", simulate_tree(tree = sort(unname(brts), decreasing=TRUE), pars = c(0.6, 0.15), method = "thinning"))
run("numeric brts, model dd (4 pars), useDDD=FALSE", simulate_tree(tree = unname(brts), pars = c(0.6, -0.01, 0.15, 0), model = "dd", useDDD = FALSE))
run("integer brts c(5L,3L,1L)", simulate_tree(tree = c(5L, 3L, 1L), pars = c(0.5, 0.1)))
run("numeric brts, pars MATRIX 2 rows", simulate_tree(tree = unname(brts), pars = rbind(c(0.5,0.1), c(0.6,0.15))))
run("numeric brts, link=power, thinning", simulate_tree(tree = unname(brts), pars = c(0.6, 0.15), method = "thinning", link = "power"))
run("numeric brts, rho=0.8", simulate_tree(tree = unname(brts), pars = c(0.6, 0.15), rho = 0.8))
cat("\n-- controls --\n")
run("phylo, bdi, cr", simulate_tree(tree = phy, pars = c(0.6, 0.15)))
run("phylo, useDDD=FALSE", simulate_tree(tree = phy, pars = c(0.6, 0.15), useDDD = FALSE))
sim <- simulate_tree(pars = c(0.5, 0.1), max_t = 5)
run("simulate_tree result as tree", simulate_tree(tree = sim, pars = c(0.5, 0.1)))
cat("\n-- where exactly does it stop? --\n")
r <- tryCatch(withCallingHandlers(simulate_tree(tree = c(5,3,1), pars = c(0.5,0.1)),
      error = function(e) { cat("call stack:", paste(rev(vapply(sys.calls(), function(cl) as.character(cl[[1]])[1], "")), collapse=" <- "), "\n") }), error=function(e) e)
cat("\n-- estimate_rates on numeric brts (control, short) --\n")
fit <- tryCatch(estimate_rates(tree = unname(brts), model = "cr", init_pars = c(0.5, 0.1),
        control = list(max_iter = 2L, max_time = 20, num_threads = 1L, lower_bound = c(0.01,0.001), upper_bound = c(2,1))), error=function(e) e)
cat(if (inherits(fit,"error")) paste("ERROR:", conditionMessage(fit)) else paste("OK pars", paste(round(fit$pars,3), collapse=",")), "\n")
