# compare_models() chooses the most parsimonious model within `margin` AIC of
# the minimum, and warns when the chosen fit rests on fewer than 10 effective
# draws.

mk <- function(loglik, k, ess, model = c(0L, 0L, 0L, 0L, 0L)) {
  structure(list(model = model, n_pars = as.integer(k), loglik = loglik, AIC = 2 * k - 2 * loglik,
                 pars = rep(0.1, k), details = list(final_IS = list(ESS = ess, n = 300))),
            class = "emphasis_fit")
}

test_that("a richer model within 2 AIC of the minimum is not chosen", {
  tb <- suppressWarnings(compare_models(cr = mk(-100, 2, 300), dd = mk(-97.5, 4, 290)))
  expect_equal(tb$model[1L], "dd")              # dd has the lower AIC, by 1
  expect_identical(attr(tb, "chosen"), "cr")    # within the margin, so the simpler model
  expect_true(all(tb$within_margin))
})

test_that("a richer model beyond the margin is chosen", {
  tb <- suppressWarnings(compare_models(cr = mk(-100, 2, 300), dd = mk(-96, 4, 290)))
  expect_identical(attr(tb, "chosen"), "dd")    # leads by 4
  expect_identical(tb$within_margin, c(TRUE, FALSE))
})

test_that("the margin is an argument", {
  tb <- suppressWarnings(compare_models(cr = mk(-100, 2, 300), dd = mk(-96, 4, 290), margin = 10))
  expect_identical(attr(tb, "chosen"), "cr")
  expect_error(compare_models(cr = mk(-100, 2, 300), dd = mk(-96, 4, 290), margin = -1), "margin")
})

test_that("a chosen model on fewer than 10 effective draws is warned about", {
  expect_warning(tb <- compare_models(cr = mk(-100, 2, 300), nk = mk(-90, 4, 4)), "effective sample of 4.0")
  expect_true(attr(tb, "chosen_starved"))
  expect_identical(attr(tb, "chosen"), "nk")
})
