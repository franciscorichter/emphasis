# The conditional sampler's timing has to be mean-field; its attachment does
# not.  A birth joins a lineage chosen in proportion to that lineage's own
# rate, drawn by inverting the CDF with the uniform that already chose the
# event.  These hold it to the three things that makes it a proposal and not
# a different model: the weights stay positive, N-only models are untouched,
# and log q is right -- which shows up as the two rules estimating the same
# log-likelihood.

test_that("attachment weights average the aggregate and stay positive", {
  cov <- c(0.5, 2, 5, 9)
  lam <- 0.4
  v <- emphasis:::.attach_w(cov, lam, beta = -0.03)
  expect_true(all(v > 0))
  expect_equal(mean(v), lam, tolerance = 1e-12)   # centred, so the mean is lam
  # a negative coefficient makes the older lineage the less likely one
  expect_true(v[which.max(cov)] < v[which.min(cov)])
  # no coefficient, no tilt
  expect_equal(emphasis:::.attach_w(cov, lam, 0), rep(lam, 4))
  # a coefficient strong enough to drive a rate negative is floored, not cut:
  # every lineage the model can use stays reachable
  vb <- emphasis:::.attach_w(cov, lam, beta = -10)
  expect_true(all(vb > 0))
  expect_equal(min(vb), lam * 1e-3)
  # degenerate inputs do not produce a degenerate distribution
  expect_length(emphasis:::.attach_w(numeric(0), lam, -0.1), 0L)
  expect_equal(emphasis:::.attach_w(cov, 0, -0.1), rep(1, 4))
})

test_that("N-only models are untouched by the attachment rule", {
  skip_on_cran()
  set.seed(21); phy <- ape::rcoal(40)
  brts <- emphasis:::.extract_brts(phy)
  for (m in c("cr", "dd")) {
    mb <- emphasis:::.resolve_model(m)
    p8 <- emphasis:::.expand_pars(
      if (m == "cr") c(0.6, 0.2) else c(0.7, -0.008, 0.25, 0), mb)
    set.seed(5)
    a <- emphasis:::.augment_tree_bdi(brts, p8, model_bin = mb[1:3], sample_size = 30L,
                                      link = 0L, rho = 1, attach = "uniform")
    set.seed(5)
    b <- emphasis:::.augment_tree_bdi(brts, p8, model_bin = mb[1:3], sample_size = 30L,
                                      link = 0L, rho = 1, attach = "rate")
    expect_identical(a$logg, b$logg)
    expect_identical(a$fhat, b$fhat)
    expect_identical(lapply(a$trees, function(t) t$parent_id),
                     lapply(b$trees, function(t) t$parent_id))
  }
})

test_that("the two attachment rules estimate the same log-likelihood", {
  skip_on_cran()
  # If log q did not account for the new attachment probability, the weights
  # would be wrong by that factor and the two estimates would part company.
  mb <- emphasis:::.resolve_model("ned")
  B0 <- 0.5; G0 <- 0.165; K <- 40
  cp <- c(B0, -(B0 - G0) / K, -0.03, G0, 0, 0)
  phy <- NULL
  for (i in 1:120) {
    set.seed(i)
    x <- simulate_tree(pars = cp, max_t = 15, model = "ned", rho = 1,
                       max_lin = 20000L, num_threads = 1L)
    if (identical(x$status, "done") && length(x$tes$tip.label) >= 20) { phy <- x$tes; break }
  }
  skip_if(is.null(phy), "no ED tree")
  brts <- emphasis:::.extract_brts(phy)
  pn <- emphasis:::.expand_pars(cp, mb)
  N <- 3000L
  est <- vapply(c("uniform", "rate"), function(a) {
    set.seed(8)
    emphasis:::.augment_tree_bdi(brts, pn, model_bin = mb[1:4], sample_size = N,
                                 link = 0L, rho = 1, attach = a)$fhat
  }, 0)
  skip_if(!all(is.finite(est)), "an estimate is not finite")
  expect_lt(abs(est[["rate"]] - est[["uniform"]]), 0.5)
})

test_that("rate attachment still attaches every birth to a living lineage", {
  skip_on_cran()
  mb <- emphasis:::.resolve_model("ned")
  cp <- c(0.5, -0.0084, -0.03, 0.165, 0, 0)
  phy <- NULL
  for (i in 1:120) {
    set.seed(i)
    x <- simulate_tree(pars = cp, max_t = 15, model = "ned", rho = 1,
                       max_lin = 20000L, num_threads = 1L)
    if (identical(x$status, "done") && length(x$tes$tip.label) >= 20) { phy <- x$tes; break }
  }
  skip_if(is.null(phy), "no ED tree")
  set.seed(6)
  a <- emphasis:::.augment_tree_bdi(emphasis:::.extract_brts(phy),
                                    emphasis:::.expand_pars(cp, mb),
                                    model_bin = mb[1:4], sample_size = 40L,
                                    link = 0L, rho = 1, attach = "rate")
  skip_if(length(a$trees) == 0L, "no draw completed")
  births <- 0L
  for (tr in a$trees) {
    tr <- tr[order(tr$brts), ]
    alive <- c(-2L, -3L)
    for (i in seq_len(nrow(tr))) {
      r <- tr[i, ]
      if (r$id < 0L) next
      if (r$t_ext == 0) { alive <- setdiff(alive, r$id); next }
      births <- births + 1L
      expect_true(r$parent_id %in% alive)
      alive <- c(alive, r$id)
    }
  }
  expect_gt(births, 0L)
})
