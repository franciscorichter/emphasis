# H36: .validate_linear_init_pars guards only the positive covariate extreme; asymmetric <=0 / <0; beta0=0 -> slope 0
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
v <- emphasis:::.validate_linear_init_pars
set.seed(6)
tr <- ape::rcoal(30); brts <- sort(ape::branching.times(tr), decreasing = TRUE); brts <- brts/max(brts)*10
mb <- c(0L,0L,1L); lb <- c(0, -1, 0, -1); ub <- c(2, 1, 1, 1)
cat("A) beta_D = +0.05, beta0 = 0.5: lambda = 0 at D = -10 (negative extreme). Clamped? ")
a <- v(c(0.5, 0.05, 0.1, 0), mb, brts, lb, ub, verbose = TRUE); print(a); cat("   lambda at D=-10:", max(0, a[1] + a[2]*(-10)), "\n")
cat("B) beta0 = 0, beta_D = -0.05 -> "); b <- v(c(0, -0.05, 0.1, 0), mb, brts, lb, ub, verbose = TRUE); print(b)
cat("   lambda identically", b[1], "+", b[2], "*D = 0 everywhere\n")
cat("C) mu: gamma0 = 0.1, gamma_D = -0.01 gives mu = 0 exactly at D = 10 -> "); cc <- v(c(0.5, 0, 0.1, -0.01), mb, brts, lb, ub, verbose = TRUE); print(cc)
cat("   (lambda analogue beta0=0.1, beta_D=-0.01 -> "); d <- v(c(0.1, -0.01, 0.1, 0), mb, brts, lb, ub, verbose = TRUE); print(d)
# range of D actually realised on this tree at the initial (observed) tree: pendant ages in [0,T], mean in [0,T] -> D in (-T, T)
# Does case A collapse in an E-step? (stochastic; count zero weights)
p8 <- emphasis:::.expand_pars(a, mb)
r <- tryCatch(emphasis:::augment_trees(brts, p8, sample_size = 50L, maxN = 500L, max_missing = 1e4, max_lambda = 1e6, num_threads = 1L, model = mb, link = 0L),
              error = function(e) conditionMessage(e))
if (is.character(r)) cat("E-step error:", r, "\n") else cat(sprintf("trees=%d zero_w=%d overruns=%d lambda=%d\n", length(r$trees), r$rejected_zero_weights, r$rejected_overruns, r$rejected_lambda))
