# The mean-field sweep is a Picard iteration x <- F(x).  Where its map expands
# the iteration diverges and a bigger budget makes it worse, so the proposal
# gets built on a blown-up mean field (audit H108).  Under-relaxation takes a
# partial step instead.  "auto" halves the relaxation whenever the residual
# grows, which is a no-op wherever the map already contracts.

# A tree and a parameter set where the plain iteration diverges.  Deterministic:
# the seed fixes the tree and the rates are hard-coded.
diverging <- function() {
  set.seed(3)
  phy  <- ape::rphylo(150, 0.5, 0.2)
  brts <- emphasis:::.extract_brts(phy)
  tp   <- brts[1]
  list(bt = sort(tp - brts[-1]), tp = tp,
       mb = emphasis:::.resolve_model("dd"),
       # b0 = 2.4, turnover 0.5, capacity 120
       p8 = emphasis:::.expand_pars(c(2.4, -(2.4 - 1.2) / 120, 1.2, 0),
                                    emphasis:::.resolve_model("dd")))
}

contracting <- function() {
  set.seed(21)
  phy  <- ape::rcoal(60)
  brts <- emphasis:::.extract_brts(phy)
  tp   <- brts[1]
  list(bt = sort(tp - brts[-1]), tp = tp,
       mb = emphasis:::.resolve_model("dd"),
       p8 = emphasis:::.expand_pars(c(0.7, -0.008, 0.25, 0),
                                    emphasis:::.resolve_model("dd")))
}

it <- function(cs, w, max_iter = 200L) {
  suppressWarnings(emphasis:::.bdi_iterate(cs$p8, cs$mb[1:3], 0L, cs$bt, cs$tp,
                                           rho = 1, max_iter = max_iter,
                                           damping = w))
}

test_that("the damping argument is validated", {
  cs <- contracting()
  expect_error(it(cs, 0), "damping")
  expect_error(it(cs, -0.5), "damping")
  expect_error(it(cs, 1.5), "damping")
  expect_error(it(cs, NA_real_), "damping")
  expect_error(it(cs, "half"), "damping")
})

test_that("the plain iteration diverges where the map expands, and worse with a bigger budget", {
  skip_on_cran()
  cs <- diverging()
  short <- it(cs, 1, max_iter = 20L)
  long  <- it(cs, 1, max_iter = 200L)
  expect_false(short$converged)
  expect_false(long$converged)
  expect_gt(short$delta, 1)
  # the signature of divergence rather than slowness: more sweeps, worse residual
  expect_gt(long$delta, short$delta)
})

test_that("under-relaxation contracts the case the plain iteration diverges on", {
  skip_on_cran()
  cs <- diverging()
  for (w in list(0.5, 0.25, "auto")) {
    r <- it(cs, w)
    expect_true(r$converged, info = paste("omega", w))
    expect_lt(r$delta, 1e-3)
  }
})

test_that("auto costs nothing where the map already contracts", {
  skip_on_cran()
  cs <- contracting()
  one  <- it(cs, 1)
  auto <- it(cs, "auto")
  expect_true(one$converged)
  expect_true(auto$converged)
  # auto only halves when the residual grows, which never happens here, so the
  # two trajectories are the same one
  expect_identical(auto$iterations, one$iterations)
  expect_identical(auto$delta, one$delta)
})

test_that("the default leaves an ordinary draw untouched", {
  skip_on_cran()
  set.seed(21); phy <- ape::rcoal(50)
  brts <- emphasis:::.extract_brts(phy)
  mb <- emphasis:::.resolve_model("dd")
  p8 <- emphasis:::.expand_pars(c(0.7, -0.008, 0.25, 0), mb)
  set.seed(4)
  a <- emphasis:::.augment_tree_bdi(brts, p8, model_bin = mb[1:3],
                                    sample_size = 25L, link = 0L, rho = 1,
                                    damping = 1)
  set.seed(4)
  b <- emphasis:::.augment_tree_bdi(brts, p8, model_bin = mb[1:3],
                                    sample_size = 25L, link = 0L, rho = 1,
                                    damping = "auto")
  expect_identical(a$logg, b$logg)
  expect_identical(a$fhat, b$fhat)
})

# Decision D2: a relaxed or non-converged mean field is recorded in the fit's
# notes, with the relaxation, the sweeps, the final residual and whether it
# converged, and print() shows it.

test_that("the iteration reports the relaxation it ended at", {
  cs <- contracting()
  auto <- it(cs, "auto")
  expect_identical(auto$omega, 1)
  expect_identical(auto$damping, "auto")
  expect_identical(it(cs, 0.5)$omega, 0.5)
  skip_on_cran()
  dv <- it(diverging(), "auto")
  expect_lt(dv$omega, 1)
})

mf_fit <- function(omega, converged, damping = "auto") {
  n <- length(omega)
  structure(list(model = c(1L, 0L, 0L, 0L, 0L), pars = c(0.5, -0.01, 0.1, 0),
                 details = list(damping = damping,
                                mcem = data.frame(mf_converged = converged, mf_delta = rep(2e-5, n),
                                                  mf_sweeps = rep(12L, n), mf_omega = omega,
                                                  mf_budget = rep(20L, n)))),
            class = "emphasis_fit")
}

test_that("the mean-field note appears only when the iteration was relaxed, not converged, or set", {
  expect_null(emphasis:::.mean_field_note(mf_fit(c(1, 1), c(TRUE, TRUE))))
  n <- emphasis:::.mean_field_note(mf_fit(c(1, 0.25), c(TRUE, TRUE)))
  expect_match(n, "took 12 of 20 sweeps")
  expect_match(n, "relaxation 0.25")
  expect_match(n, "residual 2.00e-05: converged")
  expect_match(n, "1 relaxed the iteration below the full step and 0 did not converge")
  n <- emphasis:::.mean_field_note(mf_fit(c(1, 1), c(TRUE, FALSE)))
  expect_match(n, "not converged")
  n <- emphasis:::.mean_field_note(mf_fit(c(0.5, 0.5), c(TRUE, TRUE), damping = 0.5))
  expect_match(n, "damping = 0.5")
  # cr: no mean field is solved, nothing to say
  expect_null(emphasis:::.mean_field_note(mf_fit(c(NA, NA), c(NA, NA))))
  expect_identical(names(emphasis:::.fit_notes(mf_fit(c(1, 0.25), c(TRUE, TRUE)), 100L, 0L)),
                   "mean_field")
})

test_that("a fit records a relaxed mean field in its notes and print() says so", {
  skip_on_cran()
  set.seed(21); phy <- ape::rcoal(50)
  ctrl <- list(lower_bound = c(0.1, -0.05, 0.01, -0.05), upper_bound = c(2, 0.05, 1, 0.05),
               num_trees = 20L, max_iter = 2L)
  fit <- suppressWarnings(estimate_rates(phy, method = "mcem", model = "dd", link = "linear",
                                         init_pars = c(0.7, -0.008, 0.25, 0),
                                         control = c(ctrl, list(damping = 0.5))))
  expect_identical(fit$sampling, "bdi")
  expect_true(all(c("mf_converged", "mf_delta", "mf_sweeps", "mf_omega") %in% names(fit$details$mcem)))
  expect_true("mean_field" %in% names(fit$notes))
  expect_match(fit$notes[["mean_field"]], "damping = 0.5")
  expect_output(print(fit), "mean field")
  plain <- suppressWarnings(estimate_rates(phy, method = "mcem", model = "dd", link = "linear",
                                           init_pars = c(0.7, -0.008, 0.25, 0), control = ctrl))
  if (all(plain$details$mcem$mf_omega == 1) && all(plain$details$mcem$mf_converged))
    expect_false("mean_field" %in% names(plain$notes))
})

# The default relaxation depends on the sampling fraction (experiment E10c):
# "auto" at rho = 1, a fixed 1/8 below it.  An explicit setting wins.
test_that("the default damping resolves by rho, and an explicit setting is kept", {
  expect_identical(estimate_rates_control("mcem")$damping, "default")
  expect_identical(emphasis:::.resolve_damping("default", 1), "auto")
  expect_identical(emphasis:::.resolve_damping(NULL, 1), "auto")
  expect_identical(emphasis:::.resolve_damping("default", 0.2), 0.125)
  expect_identical(emphasis:::.resolve_damping("default", 0.999), 0.125)
  expect_identical(emphasis:::.resolve_damping("auto", 0.2), "auto")
  expect_identical(emphasis:::.resolve_damping(0.5, 1), 0.5)
  expect_identical(emphasis:::.resolve_damping(1, 0.2), 1)
})

test_that("the mean-field note names the resolved default", {
  f <- mf_fit(c(0.125, 0.125), c(TRUE, TRUE), damping = 0.125)
  f$details$damping_default <- TRUE
  expect_match(emphasis:::.mean_field_note(f), "damping = 0.125, the default at rho < 1")
  f <- mf_fit(c(1, 0.5), c(TRUE, TRUE), damping = "auto")
  f$details$damping_default <- TRUE
  expect_match(emphasis:::.mean_field_note(f), "damping = \"auto\", the default at rho = 1")
  # set by the user: no "default" in the note
  f <- mf_fit(c(0.125, 0.125), c(TRUE, TRUE), damping = 0.125)
  expect_false(grepl("the default", emphasis:::.mean_field_note(f)))
})

test_that("a fit under rho < 1 runs the mean field at 1/8 and says so", {
  skip_on_cran()
  set.seed(21); phy <- ape::rcoal(40)
  ctrl <- list(lower_bound = c(0.1, -0.05, 0.01, -0.05), upper_bound = c(2, 0.05, 1, 0.05),
               num_trees = 10L, max_iter = 1L)
  fit <- suppressWarnings(estimate_rates(phy, method = "mcem", model = "dd", link = "linear",
                                         init_pars = c(0.7, -0.008, 0.25, 0),
                                         control = c(ctrl, list(rho = 0.5))))
  expect_identical(fit$sampling, "bdi")
  expect_identical(fit$details$damping, 0.125)
  expect_true(all(fit$details$mcem$mf_omega == 0.125))
  expect_match(fit$notes[["mean_field"]], "damping = 0.125, the default at rho < 1")
  fit1 <- suppressWarnings(estimate_rates(phy, method = "mcem", model = "dd", link = "linear",
                                          init_pars = c(0.7, -0.008, 0.25, 0), control = ctrl))
  expect_identical(fit1$details$damping, "auto")
  fit2 <- suppressWarnings(estimate_rates(phy, method = "mcem", model = "dd", link = "linear",
                                          init_pars = c(0.7, -0.008, 0.25, 0),
                                          control = c(ctrl, list(rho = 0.5, damping = "auto"))))
  expect_identical(fit2$details$damping, "auto")
})
