lib <- commandArgs(trailingOnly=TRUE)[1]
.libPaths(c(lib, .libPaths())); library(emphasis)
b12 <- c(10,8.3,7.1,6.2,5.5,4.4,3.9,3.1,2.2,1.5,0.9,0.3)
cat("start\n"); flush(stdout())
fit <- try(estimate_rates(tree = b12, method = "mcem", model = "cr",
             init_pars = c(0.5, 0.1),
             control = list(num_trees = 0L, max_iter = 2L, sampling = "dynamic_fresh",
                            lower_bound = c(0.01,0.001), upper_bound = c(5,5))), silent=TRUE)
cat("MSG:", if (inherits(fit,"try-error")) conditionMessage(attr(fit,"condition")) else "OK-returned", "\n")
