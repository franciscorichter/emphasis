# The survival probability by direct simulation (survival_mc), the local
# surface of log P(S) (survival_surface) and estimate_rates(cond = "surface").
#
# The reference is the closed form under cr on the exponential link: one crown
# lineage leaves a sampled descendant with
#   p = rho d / (rho lambda + (lambda (1 - rho) - mu) e^{-d T}),  d = lambda - mu,
# and P(S) = p^2 (README, "Conditioning on survival"; the experiments
# repository's sim/survival-exact.R).  Seeds are fixed, so each estimate is one
# draw; the tolerances are four standard errors.

logp_cr <- function(b, T, rho = 1) {
  l <- exp(b[1]); m <- exp(b[2]); d <- l - m
  2 * log(rho * d / (rho * l + (l * (1 - rho) - m) * exp(-d * T)))
}

test_that("survival_mc matches the exact cr survival at rho = 1 and rho = 0.7", {
  th <- log(c(0.5, 0.3)); T <- 4
  for (rho in c(1, 0.7)) {
    s <- survival_mc(th, model = "cr", crown_age = T, rho = rho, n_clades = 20000, seed = 3)
    expect_equal(s[["clades"]], 20000)
    expect_lt(abs(s[["logP"]] - logp_cr(th, T, rho)), 4 * s[["se"]])
  }
  # at a higher extinction rate and a stronger sampling loss
  th2 <- log(c(0.8, 0.6))
  s <- survival_mc(th2, model = "cr", crown_age = 3, rho = 0.4, se = 0.02, seed = 11)
  expect_lt(s[["se"]], 0.021)
  expect_lt(abs(s[["logP"]] - logp_cr(th2, 3, 0.4)), 4 * s[["se"]])
})

test_that("survival_mc draws the same clades for the same seed and leaves the R stream alone", {
  set.seed(5); before <- .Random.seed
  a <- survival_mc(log(c(0.5, 0.3)), "cr", crown_age = 4, n_clades = 500, seed = 9)
  expect_identical(.Random.seed, before)
  b <- survival_mc(log(c(0.5, 0.3)), "cr", crown_age = 4, n_clades = 500, seed = 9)
  expect_identical(a[c("logP", "survivors")], b[c("logP", "survivors")])
  # common random numbers: at a nearby point the same clades are drawn, so the
  # difference is smaller than either value's error
  c1 <- survival_mc(log(c(0.52, 0.3)), "cr", crown_age = 4, n_clades = 500, seed = 9)
  expect_lt(abs(c1[["logP"]] - a[["logP"]]), 2 * a[["se"]])
  expect_error(survival_mc(log(c(0.5, 0.3)), "cr", crown_age = 4, rho = 1.2), "rho")
  expect_error(survival_mc(log(0.5), "cr", crown_age = 4))
})

test_that("a clade past the lineage cap counts as surviving", {
  # at a cap of 3 lineages nearly every clade is stopped there before it can die
  s <- survival_mc(log(c(2, 0.01)), "cr", crown_age = 5, n_clades = 300, max_lin = 3, seed = 2)
  expect_gt(s[["P"]], 0.95)
})

test_that("the surface reproduces a direct estimate at its centre within its error", {
  th <- log(c(0.5, 0.3)); T <- 4
  S <- survival_surface(th, model = "cr", crown_age = T, radius = 0.1, se = 0.02, seed = 4)
  expect_s3_class(S, "emphasis_surface")
  expect_equal(S$points, 13L)   # 2 p + 1, p = 6 coefficients in two parameters
  expect_equal(nrow(S$design), 13L)
  d <- survival_mc(th, "cr", crown_age = T, se = 0.02, seed = 777)
  expect_lt(abs(predict(S, th) - d[["logP"]]), 4 * sqrt(S$design$se[1]^2 + d[["se"]]^2))
  expect_lt(abs(predict(S, th) - logp_cr(th, T)), 4 * S$design$se[1])
  # the design value at the centre is the direct estimate there
  expect_lt(abs(predict(S, th) - S$design$logP[1]), 2 * S$design$se[1])
  # predict_survival and estimate_rates' conditioning read it as P
  expect_equal(predict_survival(S, data.frame(beta_0 = th[1], gamma_0 = th[2])),
               exp(predict(S, th)))
  expect_equal(predict(S, th, type = "P"), exp(predict(S, th)))
})

test_that("the surface's gradient is the derivative of its value, inside and on the tangent plane", {
  S <- survival_surface(log(c(0.5, 0.3)), model = "cr", crown_age = 4, radius = 0.1,
                        n_clades = 800, seed = 4)
  fd <- function(th, h = 1e-5) vapply(seq_along(th), function(j) {
    e <- replace(numeric(length(th)), j, h)
    (predict(S, th + e) - predict(S, th - e)) / (2 * h) }, 0)
  inside <- log(c(0.5, 0.3)) + c(0.03, -0.02)
  expect_equal(as.numeric(predict(S, inside, type = "gradient")), fd(inside), tolerance = 1e-4)
  # outside the box along an axis the surface is the tangent plane at the boundary
  out <- log(c(0.5, 0.3)) + c(0.3, 0)
  edge <- log(c(0.5, 0.3)) + c(0.1, 0)
  g <- predict(S, edge, type = "gradient")
  expect_equal(predict(S, out), min(0, predict(S, edge) + sum(g * c(0.2, 0))), tolerance = 1e-10)
  expect_true(emphasis:::.surface_inside(S, inside))
  expect_false(emphasis:::.surface_inside(S, log(c(0.5, 0.3)) + c(0.06, 0)))
})

test_that("cond = \"surface\" is checked", {
  phy <- ape::read.tree(text = "((a:1,b:1):1,(c:1.5,d:1.5):0.5);")
  expect_error(estimate_rates(phy, method = "cem", cond = "surface",
                              control = list(lower_bound = c(-3, -5), upper_bound = c(1, 0))),
               "mcem")
  expect_error(estimate_rates(phy, cond = "exact",
                              control = list(lower_bound = c(-3, -5), upper_bound = c(1, 0))),
               "surface")
  expect_error(estimate_rates(phy, cond = "surface",
                              control = list(lower_bound = c(-3, -5), upper_bound = c(1, 0),
                                             surface_se = -1)),
               "surface_se")
})

# A 30-tip cr tree at crown age 5, shared by the end-to-end fits.
surface_tree <- function() {
  set.seed(2)
  repeat {
    tr <- simulate_tree(pars = log(c(0.6, 0.15)), max_t = 5, model = "cr")
    if (tr$status == "done" && ape::Ntip(tr$tes) >= 15 && ape::Ntip(tr$tes) <= 30) return(tr$tes)
  }
}

test_that("cond = \"surface\" fits cr and records what the conditioning did", {
  skip_on_cran()
  phy <- surface_tree()
  T <- max(ape::branching.times(phy))
  set.seed(1)
  f <- estimate_rates(phy, model = "cr", cond = "surface",
                      control = list(lower_bound = c(-3, -5), upper_bound = c(1, 0),
                                     num_trees = 50, max_iter = 60, max_time = 120))
  expect_s3_class(f, "emphasis_fit")
  expect_true(f$cond)
  cd <- f$conditioning
  expect_identical(cd$method, "surface")
  expect_gte(cd$surfaces, 1L)
  expect_gt(cd$clades, 0)
  expect_true(cd$start %in% c("unconditioned", "nested"))
  expect_true(is.logical(cd$restart))
  expect_true(all(is.finite(f$pars)))
  # the reported log-likelihood is the fit's minus a direct log P(S) at the estimate
  expect_equal(f$loglik, cd$loglik_uncond - cd$logP)
  expect_lt(cd$logP_se, 0.011)
  expect_lt(abs(cd$logP - logp_cr(f$pars, T)), 4 * cd$logP_se + 0.01)
  expect_equal(f$AIC, -2 * f$loglik + 2 * f$n_pars)
  expect_s3_class(cd$surface, "emphasis_surface")
  expect_false(is.null(cd$runs$unconditioned))
})

test_that("cond = \"surface\" fits a lineage-level model (d) and respects max_time", {
  skip_on_cran()
  phy <- surface_tree()
  set.seed(3)
  # the size note (30 tips) is a message, not what is tested
  f <- suppressMessages(estimate_rates(phy, model = "d", cond = "surface",
                      control = list(lower_bound = c(-3, -1, -5, -1), upper_bound = c(1, 1, 0, 1),
                                     num_trees = 50, max_iter = 30, max_time = 120)))
  expect_true(f$cond)
  expect_identical(f$sampling, "dynamic_fresh")
  expect_true(all(is.finite(f$pars)))
  expect_true(is.finite(f$loglik))
  expect_gte(f$conditioning$surfaces, 1L)
  expect_lte(f$iterations, 30L)
  # a short clock ends the fit on time_budget well within a few multiples of it
  t0 <- proc.time()[[3]]
  g <- estimate_rates(phy, model = "cr", cond = "surface",
                      control = list(lower_bound = c(-3, -5), upper_bound = c(1, 0),
                                     num_trees = 50, max_iter = 200, max_time = 2))
  expect_lt(proc.time()[[3]] - t0, 30)
  expect_true(g$stop_reason %in% c("time_budget", "converged"))
})
