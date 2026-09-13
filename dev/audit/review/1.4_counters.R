lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, .libPaths())); library(emphasis)
ns <- asNamespace("emphasis")
brts9 <- c(6, 4.8484929804551, 4.40182149214092, 3.10816371045776,
           3.07391367934619, 2.835023069071, 1.8388281974863, 0.504630059261689)
dd <- c(1L,0L,0L)
pars <- c(1.5, -0.12, 0.4, 0)   # H11 B1: lambda(N)=0 for N >= 12.5

set.seed(2)
e <- ns$.augment_tree_bdi(brts9, pars, dd, sample_size = 200L, max_missing = 30L,
                          link = 0L, rho = 1)
cat(sprintf("E-step direct: n_valid=%d n_nonfinite=%d length(logf)=%d any(!finite logf)=%s\n",
            e$n_valid, e$n_nonfinite, length(e$logf), any(!is.finite(e$logf))))

set.seed(2)
fit <- ns$.mcem_bdi(brts = brts9, pars = pars, model = dd, link = 0L,
                    lower_bound = c(0.01,-1,0.001,0), upper_bound = c(5,0,2,0),
                    max_iter = 3L, xtol = 1e-3, tol = 1e-2, patience = 3L,
                    num_threads = 1L, sample_size = 200L, max_missing = 30L,
                    conditional = NULL, rho = 1)
print(fit$mcem[, c("fhat","rejected","n_nonfinite","num_trees","ESS")])
cat("final_IS$rejected_zero_weights =", fit$final_IS$rejected_zero_weights, "\n")
cat("final_IS$n_rejected            =", fit$final_IS$n_rejected, "\n")
cat("loglik =", fit$loglik, " loglik_var =", fit$loglik_var, "\n")
