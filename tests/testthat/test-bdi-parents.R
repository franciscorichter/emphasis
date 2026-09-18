# The conditional (BDI) sampler records the real parent of every augmented
# birth, so the tree it returns carries the observed topology and can be
# scored by a per-lineage covariate.  Four things have to hold.

test_that("recording parents does not move the random stream", {
  skip_on_cran()
  set.seed(21); phy <- ape::rcoal(30)
  brts <- emphasis:::.extract_brts(phy)
  mb   <- emphasis:::.resolve_model("dd")
  p8   <- emphasis:::.expand_pars(c(0.7, -0.008, 0.25, 0), mb)
  set.seed(5)
  off <- emphasis:::.augment_tree_bdi(brts, p8, model_bin = mb[1:3], sample_size = 40L,
                                      link = 0L, rho = 1, topology = FALSE)
  set.seed(5)
  on  <- emphasis:::.augment_tree_bdi(brts, p8, model_bin = mb[1:3], sample_size = 40L,
                                      link = 0L, rho = 1, topology = TRUE)
  expect_identical(off$logg, on$logg)
  expect_identical(off$n_valid, on$n_valid)
  expect_identical(off$n_rejected, on$n_rejected)
  expect_identical(lapply(off$trees, function(t) t$brts),
                   lapply(on$trees,  function(t) t$brts))
  expect_identical(lapply(off$trees, function(t) t$t_ext),
                   lapply(on$trees,  function(t) t$t_ext))
})

test_that("every augmented birth attaches to a lineage alive at that time", {
  skip_on_cran()
  set.seed(11); phy <- ape::rcoal(40)
  brts <- emphasis:::.extract_brts(phy)
  mb   <- emphasis:::.resolve_model("dd")
  p8   <- emphasis:::.expand_pars(c(0.8, -0.005, 0.45, 0), mb)
  set.seed(9)
  a <- emphasis:::.augment_tree_bdi(brts, p8, model_bin = mb[1:3], sample_size = 40L,
                                    link = 0L, rho = 1)
  skip_if(length(a$trees) == 0L, "no draw completed")
  births <- 0L
  for (tr in a$trees) {
    tr <- tr[order(tr$brts), ]
    expect_true(all(c("focal_tip_start", "clade") %in% names(tr)))
    alive <- c(-2L, -3L)
    for (i in seq_len(nrow(tr))) {
      r <- tr[i, ]
      if (r$id < 0L) next                                  # the closing marker
      if (r$t_ext == 0) { alive <- setdiff(alive, r$id); next }
      births <- births + 1L
      expect_true(r$parent_id %in% alive)
      alive <- c(alive, r$id)
    }
    # ids are unique, so a parent names one lineage and not two
    ids <- tr$id[tr$id >= 0L & tr$t_ext != 0]
    expect_equal(anyDuplicated(ids), 0L)
  }
  expect_gt(births, 0L)
})

test_that("a conditional draw can be scored by the ED likelihood", {
  skip_on_cran()
  set.seed(13); phy <- ape::rcoal(35)
  brts <- emphasis:::.extract_brts(phy)
  mb_dd  <- emphasis:::.resolve_model("dd")
  mb_ned <- emphasis:::.resolve_model("ned")
  p8 <- emphasis:::.expand_pars(c(0.7, -0.008, 0.25, 0), mb_dd)
  pn <- emphasis:::.expand_pars(c(0.7, -0.008, -0.03, 0.25, 0, 0), mb_ned)
  set.seed(4)
  a <- emphasis:::.augment_tree_bdi(brts, p8, model_bin = mb_dd[1:3], sample_size = 30L,
                                    link = 0L, rho = 1)
  skip_if(length(a$trees) == 0L, "no draw completed")
  ev <- emphasis:::eval_logf(pn, a$trees, model = as.integer(mb_ned[1:4]),
                             link = 0L, rho = 1)
  expect_true(all(is.finite(ev$logf)))
  # the ED term moves the score: it is being read, not ignored
  p0 <- emphasis:::.expand_pars(c(0.7, -0.008, 0, 0.25, 0, 0), mb_ned)
  ev0 <- emphasis:::eval_logf(p0, a$trees, model = as.integer(mb_ned[1:4]),
                              link = 0L, rho = 1)
  expect_false(isTRUE(all.equal(ev$logf, ev0$logf)))
})

test_that("the two proposals estimate the same log-likelihood", {
  skip_on_cran()
  set.seed(2); phy <- ape::rcoal(18)
  brts <- emphasis:::.extract_brts(phy)
  mb   <- emphasis:::.resolve_model("dd")
  cp   <- c(0.6, -0.01, 0.2, 0)
  p8   <- emphasis:::.expand_pars(cp, mb)
  N    <- 4000L
  set.seed(31)
  b <- emphasis:::.augment_tree_bdi(brts, p8, model_bin = mb[1:3], sample_size = N,
                                    link = 0L, rho = 1)
  set.seed(31)
  t <- emphasis:::augment_trees(brts, p8, sample_size = N, maxN = 50L * N,
                                max_missing = 5000L, max_lambda = 1e6,
                                num_threads = 1L, model = as.integer(mb[1:4]),
                                link = 0L, rho = 1.0,
                                parent_tip_start = emphasis:::.pts(brts),
                                parent_id = emphasis:::.pid(brts))
  f_bdi  <- b$fhat
  f_thin <- emphasis:::.is_summary(t$logf - t$logg,
                                   n_zero_weight = t$rejected_zero_weights)$fhat
  skip_if(!is.finite(f_bdi) || !is.finite(f_thin), "an estimate is not finite")
  # both are importance-sampling estimates of the same integral; at this draw
  # count they agree to a few hundredths of a log-likelihood unit
  expect_lt(abs(f_bdi - f_thin), 0.15)
})
