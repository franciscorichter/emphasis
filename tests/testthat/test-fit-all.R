# model = "all": the nested ladder fitted from nested starts and ranked; and
# the three evidence-based notes a fit carries (size, starved, turnover).

test_that("the notes state the measured floors, each only where it applies", {
  mk <- function(model, pars, ess = 100, link = 0L) {
    f <- structure(list(pars = pars, model = model, loglik = -10, AIC = 24, n_pars = length(pars),
                        method = "mcem", details = list(final_IS = list(ESS = ess))), class = "emphasis_fit")
    f
  }
  # a clade-level model on a small tree: no size note
  expect_length(emphasis:::.fit_notes(mk(c(1L, 0L, 0L, 0L, 0L), c(0.5, -0.01, 0.1, 0)), n_tips = 30L, link_int = 0L), 0L)
  # a lineage-level model on a small tree: the size note, with the number
  n <- emphasis:::.fit_notes(mk(c(1L, 0L, 1L, 0L, 0L), c(0.5, -0.01, 0.02, 0.1, 0, 0)), n_tips = 30L, link_int = 0L)
  expect_named(n, "size"); expect_match(n[["size"]], "30 tips"); expect_match(n[["size"]], "50 tips")
  # starved
  n <- emphasis:::.fit_notes(mk(c(0L, 0L, 0L, 0L, 0L), c(0.5, 0.1), ess = 4.2), n_tips = 100L, link_int = 0L)
  expect_named(n, "starved"); expect_match(n[["starved"]], "4.2")
  expect_warning(emphasis:::.fit_notes(mk(c(0L, 0L, 0L, 0L, 0L), c(0.5, 0.1), ess = 4.2), 100L, 0L, emit = TRUE), "sampler-limited")
  # turnover, on both links
  n <- emphasis:::.fit_notes(mk(c(0L, 0L, 0L, 0L, 0L), c(0.5, 0.4)), n_tips = 100L, link_int = 0L)
  expect_named(n, "turnover"); expect_match(n[["turnover"]], "0.80")
  n <- emphasis:::.fit_notes(mk(c(0L, 0L, 0L, 0L, 0L), c(log(0.5), log(0.1))), n_tips = 100L, link_int = 1L)
  expect_length(n, 0L)
  expect_message(emphasis:::.fit_notes(mk(c(0L, 0L, 0L, 0L, 0L), c(0.5, 0.4)), 100L, 0L, emit = TRUE), "turnover")
})

test_that("model = \"all\" fits the ladder from nested starts and ranks it", {
  skip_on_cran()
  set.seed(21)
  phy <- ape::rphylo(30L, 0.5, 0.05)
  ctrl <- list(sample_size = 30L, max_iter = 2L, num_threads = 1L, verbose = FALSE,
               models = c("cr", "dd", "nk"), maxN = 20000L, max_missing = 1e4)
  out <- suppressWarnings(suppressMessages(
    estimate_rates(phy, method = "mcem", model = "all", control = ctrl, link = "exponential")))
  expect_s3_class(out, "emphasis_fits")
  expect_equal(out$ladder, c("cr", "dd", "nk"))
  expect_true(all(vapply(out$fits, inherits, TRUE, "emphasis_fit")))
  expect_s3_class(out$comparison, "data.frame")
  expect_setequal(out$comparison$model, c("cr", "dd", "nk"))
  # the nk fit's N coefficient was started from the dd estimate: both finite, both bounded
  expect_true(all(is.finite(out$fits$nk$pars)))
  expect_equal(names(out$fits$nk$pars), emphasis:::.par_names(emphasis:::.resolve_model("nk")))
  expect_output(print(out), "cr -> dd -> nk")
  # a bare branching-time vector loses the lineage-level models
  out2 <- suppressWarnings(suppressMessages(
    estimate_rates(sort(ape::branching.times(phy), decreasing = TRUE), method = "mcem", model = "all",
                   control = utils::modifyList(ctrl, list(models = c("cr", "nk"))), link = "exponential")))
  expect_equal(out2$ladder, "cr")
})
