.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths())); library(emphasis)
ns <- asNamespace("emphasis"); real_aug <- ns$.augment_tree_bdi
brts12 <- c(6.385824,2.063997,1.19255,0.923743,0.884126,0.822976,0.725585,0.718801,0.539214,0.077287,0.068123)
cr8 <- function(l,m) c(l,0,0,0,m,0,0,0)
bad <- cr8(2.5, 0.2)   # M-step always proposes this; sampler fails there
run <- function(max_iter) {
  seen <- list()
  testthat::with_mocked_bindings(
    .augment_tree_bdi = function(tree, pars, ...) {
      seen[[length(seen)+1L]] <<- pars[1]
      if (isTRUE(all.equal(as.numeric(pars), bad))) stop("sampler fails at theta*")
      real_aug(tree, pars, ...) },
    m_cpp = function(e_step, init_pars, ...) list(estimates = bad, nlopt = 4L, time = 0),
    .package = "emphasis",
    { set.seed(5)
      r <- ns$.mcem_bdi(brts12, pars = cr8(1.2,0.9), sample_size = 10L, max_missing = 1e4,
             lower_bound = cr8(0,0), upper_bound = cr8(4,4), max_iter = max_iter,
             xtol = 1e-3, tol = 1e-2, patience = 3L, num_threads = 1L,
             model = c(0L,0L,0L), link = 0L, rho = 1)
      cat(sprintf("max_iter=%-3d stop=%-14s iters=%d n_failed=%d pars[1]=%.3g loglik=%s rows=%s final_IS=%s\n",
          max_iter, r$stop_reason, r$iterations, r$n_failed, r$pars[1],
          format(r$loglik), if (is.null(r$mcem)) "0" else nrow(r$mcem), !is.null(r$final_IS)))
      cat("   E-step lambda sequence:", paste(signif(unlist(seen),3), collapse=" "), "\n") })
}
run(6); run(7); run(30)
