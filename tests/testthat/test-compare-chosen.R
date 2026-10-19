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

test_that("among models with the fewest parameters within the margin, the lower AIC is chosen", {
  # Meliphagidae in E28 (177 tips): ned and nd both have six parameters and are
  # both within 2 AIC of the minimum, which ned holds; dd is 2.26 behind and
  # outside the margin. The rule names ned, whatever order the fits come in.
  fits <- list(cr  = mk(-604.480113730276, 2, 300),
               dd  = mk(-594.565125715573, 4, 295.5),
               nd  = mk(-592.02332874601,  6, 24.9),
               ned = mk(-591.434418518233, 6, 44.8),
               nk  = mk(-594.7959621419,   6, 20.8))
  tb <- suppressWarnings(do.call(compare_models, fits))
  expect_equal(tb$AIC[tb$model == "ned"], 1194.86883703647, tolerance = 1e-9)
  expect_identical(tb$model[tb$within_margin], c("ned", "nd"))
  expect_identical(attr(tb, "chosen"), "ned")
  for (ord in list(rev(names(fits)), c("nd", "ned", "nk", "dd", "cr")))
    expect_identical(attr(suppressWarnings(do.call(compare_models, fits[ord])), "chosen"), "ned")
})
