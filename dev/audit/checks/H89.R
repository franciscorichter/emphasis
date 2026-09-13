.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
set.seed(2)
pars_mat <- cbind(beta_0 = runif(60, 0.3, 1.2), gamma_0 = runif(60, 0.05, 0.5))
sims <- simulate_tree(pars = pars_mat, max_t = 3, model = "cr",
                      max_tries = 0, useDDD = FALSE)
out1 <- capture.output(g <- train_GAM(sims$simulations, pars_mat, model = "cr"))
cat("train_GAM stdout lines:", length(out1), "\n"); print(out1)
cat("train_GAM has 'verbose' formal:", "verbose" %in% names(formals(train_GAM)), "\n")

surface <- data.frame(beta_0 = pars_mat[, 1], gamma_0 = pars_mat[, 2],
                      fhat = -10 + rnorm(60))
out2 <- capture.output(g2 <- emphasis:::train_likelihood_GAM(surface, par_names = c("beta_0", "gamma_0")))
cat("train_likelihood_GAM stdout lines:", length(out2), "\n"); print(out2)
cat("train_likelihood_GAM has 'verbose' formal:",
    "verbose" %in% names(formals(emphasis:::train_likelihood_GAM)), "\n")
# testthat check
r <- tryCatch({ testthat::expect_silent(train_GAM(sims$simulations, pars_mat, model = "cr")); "silent" },
              error = function(e) paste("expect_silent FAILED:", conditionMessage(e)))
cat(r, "\n")
