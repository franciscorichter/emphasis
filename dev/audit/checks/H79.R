.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(1)
sim <- NULL; while (is.null(sim$tes) || ape::Ntip(sim$tes) < 10 || ape::Ntip(sim$tes) > 40)
  sim <- simulate_tree(pars = c(0.6, 0.1), max_t = 5, model = "cr")
cat("tips:", ape::Ntip(sim$tes), "\n")
ctl <- list(lower_bound = c(0.05, 0), upper_bound = c(2, 1), max_iter = 2L, num_trees = 20L, max_time = 60, rho = 80)
for (s in c("bdi", "dynamic_fresh")) {
  r <- tryCatch({f <- estimate_rates(sim, method = "mcem", model = "cr", control = c(ctl, list(sampling = s))); paste("NO ERROR; loglik =", round(f$loglik,2))}, error = function(e) paste("ERROR:", conditionMessage(e)), warning = function(w) paste("WARNING:", conditionMessage(w)))
  cat("mcem", s, "rho=80 ->", r, "\n")
}
r <- tryCatch({f <- estimate_rates(sim, method = "cem", model = "cr", control = list(lower_bound = c(0.05, 0), upper_bound = c(2, 1), max_iter = 2L, num_particles = 10L, max_time = 60, rho = 80)); paste("NO ERROR; loglik =", round(f$loglik,2))}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("cem rho=80 ->", r, "\n")
# rho = -1 and rho = 0.5 vs 1 on eval_logf: does C++ silently reset?
aug <- emphasis:::augment_trees(emphasis:::.extract_brts(sim), c(0.6,0,0,0,0.1,0,0,0), 5L, 200L, 1e4, 500, 1L, c(0L,0L,0L), 0L, 1.0)
for (rho in c(1, 0.5, 80, -1, 0)) {
  ev <- tryCatch(emphasis:::eval_logf(c(0.6,0,0,0,0.1,0,0,0), aug$trees, c(0L,0L,0L), 0L, rho), error=function(e) list(logf=conditionMessage(e)))
  cat(sprintf("eval_logf rho=%g -> logf[1]=%s\n", rho, format(ev$logf[1])))
}
# forward sim with rho = 80
r <- tryCatch({simulate_tree(pars = c(0.6,0.1), max_t=5, model="cr", rho=80); "NO ERROR"}, error=function(e) paste("ERROR:", conditionMessage(e)), warning=function(w) paste("WARNING:", conditionMessage(w)))
cat("simulate_tree forward rho=80 ->", r, "\n")
for (rho in c(-1, 0, 0.5)) {
  r <- tryCatch({s <- simulate_tree(pars = c(0.6,0.1), max_t=5, model="cr", rho=rho, max_tries=20); paste("NO ERROR; status", s$status)}, error=function(e) paste("ERROR:", conditionMessage(e)), warning=function(w) paste("WARNING:", conditionMessage(w)))
  cat("simulate_tree forward rho=", rho, "->", r, "\n")
}
