lib <- commandArgs(trailingOnly=TRUE)[1]
.libPaths(c(lib, .libPaths())); library(emphasis)
b12 <- c(10,8.3,7.1,6.2,5.5,4.4,3.9,3.1,2.2,1.5,0.9,0.3)
cat("calling estimate_rates with num_trees = 0 ...\n"); flush(stdout())
fit <- try(estimate_rates(brts = b12, init_pars = c(0.5, 0.1),
             model = "cr", sampling = "dynamic_fresh",
             lower_bound = c(0.01, 0.001), upper_bound = c(5, 5),
             control = list(num_trees = 0L, max_iter = 2L, verbose = FALSE)), silent = TRUE)
cat("returned:", class(fit), "\n")
