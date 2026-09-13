.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths())); library(emphasis)
ns <- asNamespace("emphasis")
brts12 <- c(6.385824,2.063997,1.19255,0.923743,0.884126,0.822976,0.725585,0.718801,0.539214,0.077287,0.068123)
real_aug <- ns$.augment_tree_bdi
k <- 0L
# Fail only on the FINAL E-step (calls 1..2 ok, call 3 = final E-step fails)
testthat::with_mocked_bindings(
  .augment_tree_bdi = function(tree, pars, ...) {
    k <<- k + 1L
    if (k >= 3L) stop("simulated failure at theta_K")
    real_aug(tree, pars, ...)
  },
  .package = "emphasis",
  {
    set.seed(20)
    fit <- estimate_rates(brts12, model="cr", method="mcem", init_pars=c(1.2,0.9),
      control=list(lower_bound=c(0,0), upper_bound=c(3,3), sampling="bdi",
                   sample_size=30L, num_threads=1L, max_iter=2L, tol=1e-2))
    cat("driver loglik (details$loglik):", fit$details$loglik, "\n")
    cat("reported fit$loglik           :", fit$loglik, "\n")
    cat("trace rows:", nrow(fit$details$mcem), " last-row m_step:", tail(fit$details$mcem$m_step,1), "\n")
    cat("final_IS null?", is.null(fit$details$final_IS), " loglik_var:", fit$loglik_var, "\n")
    cat("AIC:", fit$AIC, "\n")
    print(fit)
  })
