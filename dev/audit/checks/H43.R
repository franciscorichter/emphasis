.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(1); tr <- ape::rphylo(20, 0.5, 0.1)
chk <- function(m, inactive) { m <- as.matrix(m); cat("  inactive cols exactly 0:", all(m[, inactive] == 0), " max|.|=", max(abs(m[, inactive])), "\n") }
# CR: inactive 8-slots 2,3,4,6,7,8
cem <- estimate_rates(tr, method = "cem", model = "cr",
  control = list(lower_bound = c(0, 0), upper_bound = c(1.5, 1.0), num_particles = 20, max_iter = 5, num_trees = 1, num_threads = 1))
cat("CEM cr final_pop$pars:"); chk(cem$details$final_pop$pars, c(2,3,4,6,7,8))
cat("CEM cr best_pars:");      chk(cem$details$best_pars, c(2,3,4,6,7,8))
cat("CEM cr obtained_estim:"); chk(matrix(cem$details$obtained_estim, 1), c(2,3,4,6,7,8))
mc <- estimate_rates(tr, method = "mcem", model = "cr",
  control = list(lower_bound = c(0, 0), upper_bound = c(1.5, 1.0), sample_size = 30, max_iter = 4, num_threads = 1))
tab <- mc$details$mcem[, grep("^par", names(mc$details$mcem))]
cat("MCEM cr M-step estimates (", nrow(tab), "iters):"); chk(tab, c(2,3,4,6,7,8))
# DD (N active): inactive slots 3,4,7,8
mc2 <- estimate_rates(tr, method = "mcem", model = "dd",
  control = list(lower_bound = c(0, -0.05, 0, -0.01), upper_bound = c(1.5, 0.05, 1.0, 0.01), sample_size = 30, max_iter = 4, num_threads = 1))
tab2 <- mc2$details$mcem[, grep("^par", names(mc2$details$mcem))]
cat("MCEM dd M-step estimates (", nrow(tab2), "iters):"); chk(tab2, c(3,4,7,8))
cem2 <- estimate_rates(tr, method = "cem", model = "dd",
  control = list(lower_bound = c(0, -0.05, 0, -0.01), upper_bound = c(1.5, 0.05, 1.0, 0.01), num_particles = 20, max_iter = 5, num_trees = 1, num_threads = 1))
cat("CEM dd final_pop$pars:"); chk(cem2$details$final_pop$pars, c(3,4,7,8))
