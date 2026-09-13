.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(5)
sim <- NULL; while (is.null(sim$tes) || ape::Ntip(sim$tes) < 10 || ape::Ntip(sim$tes) > 40) sim <- simulate_tree(pars = c(0.6,0.1), max_t=5, model="cr")
show <- function(lbl, expr) { r <- tryCatch(withCallingHandlers(expr, warning=function(w) {cat("  [warning]", conditionMessage(w), "\n"); invokeRestart("muffleWarning")}), error=function(e) paste("ERROR:", conditionMessage(e)))
  cat(lbl, "->", if (is.list(r)) sprintf("tas=%s log_q=%s", class(r$tas)[1], format(r$log_q)) else r, "\n") }
show("invalid model_bin c(9,9,9) via .sim_tree_conditional (bdi)", emphasis:::.sim_tree_conditional(sim, c(0.6,0.1), c(9L,9L,9L), 1L, TRUE, 0L, 1e4, 500, NULL, 1L, 1.0, "bdi"))
show("invalid model_bin c(9,9,9) (thinning)", emphasis:::.sim_tree_conditional(sim, c(0.6,0.1), c(9L,9L,9L), 1L, TRUE, 0L, 1e4, 500, NULL, 1L, 1.0, "thinning"))
show("max_missing = -1 (thinning)", simulate_tree(tree=sim, pars=c(0.6,0.1), model="cr", max_missing=-1, method="thinning"))
show("max_lambda = 1e-9 (thinning, forces augmentation_lambda throw)", simulate_tree(tree=sim, pars=c(0.6,0.1), model="cr", max_lambda=1e-9, method="thinning"))
show("negative mu, bdi", simulate_tree(tree=sim, pars=c(0.6,-0.5), model="cr", method="bdi"))
# direct call shows the real error text that tryCatch hides
cat("direct .augment_tree_internal with model_bin c(9,9,9):", tryCatch({emphasis:::.augment_tree_internal(sim, pars=c(0.6,0.1), model_bin=c(9L,9L,9L), sample_size=1L, max_missing=1e4, max_lambda=500, maxN=NULL, num_threads=1L, link=0L, rho=1.0); "no error"}, error=function(e) conditionMessage(e)), "\n")
cat("direct .augment_tree_internal with max_lambda 1e-9:", tryCatch({emphasis:::.augment_tree_internal(sim, pars=c(0.6,0.1), model_bin=c(0L,0L,0L), sample_size=1L, max_missing=1e4, max_lambda=1e-9, maxN=NULL, num_threads=1L, link=0L, rho=1.0); "no error"}, error=function(e) conditionMessage(e)), "\n")
