# Candidate addition to tests/testthat/test-logsum.R: pins the "<= 0" half of the
# rule (a strictly negative rate), which no existing expectation covers.
mk_tree <- function(k, tp = k + 1) {
  data.frame(brts = c(seq_len(k), tp), n = c(2:(k + 1), k + 2), t_ext = 1e11,
             pd = 0, tip_start = 0, id = c(seq_len(k) - 1L, -1L), parent_id = -1L)
}
test_that("a strictly negative rate (gaussian link, beta_0 < 0) gives -Inf", {
  # gaussian: lambda = beta_0 * exp(-(eta_cov - 1)^2/2); beta_0 = -0.4 gives
  # lambda = -0.2426 < 0 at every node.  Four scored nodes, so the running
  # product is positive-signed; before the fix this returned +Inf.
  p  <- c(-0.4, 0, 0, 0, 0.2, 0, 0, 0)
  lf <- emphasis:::eval_logf(p, list(mk_tree(4)), model = c(1L, 0L, 0L),
                             link = 2L, rho = 1)$logf
  expect_identical(lf, -Inf)
})
