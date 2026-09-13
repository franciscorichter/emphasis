## H73 replication: vary tree, model, and test the fix-sketch scale.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
set.seed(7)
sim_at <- function(b0, g0, link, max_t, n = 12, model = "cr") {
  pm <- if (model == "cr") cbind(b0, g0) else cbind(b0, 0, g0, 0)
  pm <- pm[rep(1, n), , drop = FALSE]
  s <- simulate_tree(pars = pm, max_t = max_t, model = model, link = link,
                     max_tries = 0, max_lin = 2000L, num_threads = 1L)
  st <- sapply(s$simulations, `[[`, "status")
  nt <- sapply(s$simulations, function(x) if (x$status == "done" && !is.null(x$tes)) Ntip(x$tes) else 0L)
  sprintf("done=%d/%d tips(done)=%s", sum(st=="done"), n, paste(nt[st=="done"], collapse=","))
}
report <- function(tree, label, model) {
  max_t <- max(branching.times(tree)); N <- Ntip(tree)
  r_hat <- log(N/2)/max_t; mu_hat <- 0.2*r_hat; lam_hat <- r_hat+mu_hat
  cat(sprintf("\n=== %s: T=%.3f N=%d lam_hat=%.4f log(lam_hat)=%.4f model=%s\n",
              label, max_t, N, lam_hat, log(lam_hat), model))
  cat(" gaussian at centre beta_0=log(lam_hat):   ", sim_at(log(lam_hat), log(mu_hat), "gaussian", max_t, model=model), "\n")
  cat(" gaussian at beta_0=lam_hat (natural):     ", sim_at(lam_hat, mu_hat, "gaussian", max_t, model=model), "\n")
  cat(" gaussian at beta_0=lam_hat*exp(.5) (fix): ", sim_at(lam_hat*exp(.5), mu_hat*exp(.5), "gaussian", max_t, model=model), "\n")
  ab <- auto_bounds(tree, model = model, link = "gaussian", n_test = 3L, bisect_steps = 3L,
                    train_surv_gam = FALSE, verbose = FALSE, num_threads = 1L)
  cat(" auto_bounds centre:", paste(round(ab$center,4), collapse=", "), "\n")
  cat(" auto_bounds lb:    ", paste(round(ab$lower_bound,4), collapse=", "), "\n")
  cat(" auto_bounds ub:    ", paste(round(ab$upper_bound,4), collapse=", "), "\n")
  cat(sprintf(" ASSERT lb[beta_0] > 0: %s ; centre beta_0 > 0: %s\n",
              ab$lower_bound[1] > 0, ab$center[1] > 0))
  invisible(ab)
}
data(bird.orders)
## (a) same tree, other model
report(bird.orders, "bird.orders", "cr")
## (b) a different tree: simulated exponential-link CR tree, ~30 tips
s <- simulate_tree(pars = matrix(c(log(0.3), log(0.05)), 1), max_t = 15, model = "cr",
                   link = "exponential", max_tries = 50, max_lin = 500L, num_threads = 1L)
tr <- s$simulations[[1]]$tes
cat("\nsimulated tree tips:", Ntip(tr), "\n")
report(tr, "simulated exp-link CR tree", "dd")
## (c) tree where log(lam_hat) > 0: rescale bird.orders to crown age 0.9
tr2 <- bird.orders; tr2$edge.length <- tr2$edge.length / max(branching.times(bird.orders)) * 0.9
report(tr2, "bird.orders rescaled T=0.9", "cr")
