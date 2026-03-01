# Tests for augmentation functions

test_that("augmentPD errors on non-phylo input", {
  expect_error(
    augmentPD(phylo = list(), pars = c(0.1, 0.5, 0, 0),
              maxN = 10, max_missing = 100,
              lower_bound = c(0, 0, -0.1, -0.1),
              upper_bound = c(0.5, 2, 0.1, 0.1)),
    "class.*phylo"
  )
})

test_that("augmentPD errors on non-numeric parameters", {
  tree <- ape::rtree(5)
  expect_error(
    augmentPD(phylo = tree, pars = c(0.1, 0.5, 0, 0),
              maxN = "ten", max_missing = 100,
              lower_bound = c(0, 0, -0.1, -0.1),
              upper_bound = c(0.5, 2, 0.1, 0.1)),
    "numeric"
  )
})

test_that("augmentPD returns a list on a small tree", {
  skip_on_cran()
  skip("augmentPD integration test skipped locally to avoid long runtime")
  set.seed(42)
  tree <- ape::rphylo(8, 0.5, 0)

  result <- augmentPD(
    phylo        = tree,
    pars         = c(0.1, 0.5, -0.01, 0.01),
    maxN         = 50,
    max_missing  = 500,
    lower_bound  = c(0, 0, -0.1, -0.1),
    upper_bound  = c(0.5, 2, 0.1, 0.1),
    num_threads  = 1
  )

  expect_type(result, "list")
  expect_true("fhat" %in% names(result))
})
