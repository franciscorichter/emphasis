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
