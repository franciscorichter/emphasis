# The measured gate numbers behind tests/testthat/test-d-compensator.R.
#   Rscript dev/gates.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
e <- new.env()
# Source the test file with the testthat verbs stubbed out, so only its helpers
# and fixtures are defined.
local({
  stub <- function(...) invisible(NULL)
  for (nm in c("test_that", "expect_true", "expect_lt", "expect_gt", "expect_equal",
               "expect_identical", "expect_lte", "skip_on_cran", "skip_if"))
    assign(nm, stub, envir = e)
  sys.source("tests/testthat/test-d-compensator.R", envir = e)
})

hands <- list(e$.hand_tree(), e$.hand_tree2())
w <- 0
for (df in hands) for (p in e$.thetas)
  w <- max(w, abs(e$.cpp_logf(p, df) - e$.logf_ref(p, df)))
cat(sprintf("gate 3  hand trees, closed form vs kink-split integrate():  max |diff| = %.3e\n", w))

set.seed(5)
phy  <- ape::rphylo(12L, 0.8, 0.2)
brts <- emphasis:::.extract_brts(phy)
aug <- augment_trees(as.numeric(brts), c(0.9, 0, 0, 0.10, 0.35, 0, 0, 0.03),
                     10L, 20000L, 300L, 1e6, 1L, model = c(0L, 0L, 1L), link = 0L,
                     rho = 1, parent_tip_start = emphasis:::.pts(brts))
w <- 0
for (df in aug$trees) for (p in e$.thetas[c(3L, 5L, 6L)]) {
  a <- e$.cpp_logf(p, df); b <- e$.logf_ref(p, df)
  if (is.finite(a) && is.finite(b)) w <- max(w, abs(a - b))
}
cat(sprintf("gate 3  %d augmented trees, same comparison:                max |diff| = %.3e\n",
            length(aug$trees), w))

trees <- c(aug$trees, hands)
w <- 0
for (link in c(0L, 1L)) {
  grid <- if (link == 0L)
    list(list(c(0L, 0L, 0L), c(0.90,  0.00, 0, 0, 0.35, 0, 0, 0)),
         list(c(1L, 0L, 0L), c(0.90, -0.02, 0, 0, 0.35, 0, 0, 0)),
         list(c(1L, 0L, 0L), c(1.30, -0.06, 0, 0, 0.10, 0, 0, 0)))
  else
    list(list(c(0L, 0L, 0L), c(-0.10,  0.00, 0, 0, -1.20, 0, 0, 0)),
         list(c(1L, 0L, 0L), c(-0.10, -0.01, 0, 0, -1.20, 0, 0, 0)))
  for (g in grid) {
    base <- eval_logf(g[[2L]], trees, model = g[[1L]], link = link, rho = 1)$logf
    dmod <- eval_logf(g[[2L]], trees, model = c(g[[1L]][1L], 0L, 1L), link = link, rho = 1)$logf
    w <- max(w, max(abs(dmod - base)))
  }
}
cat(sprintf("gate 2  beta_D = gamma_D = 0, links 0 and 1, %d trees:       max |diff| = %.3e\n",
            length(trees), w))

set.seed(5)
phy  <- ape::rphylo(14L, 0.8, 0.2)
brts <- emphasis:::.extract_brts(phy)
p <- c(0.9, 0, 0, 0.10, 0.35, 0, 0, 0.03)
aug <- augment_trees(as.numeric(brts), p, 12L, 20000L, 300L, 1e6, 1L,
                     model = c(0L, 0L, 1L), link = 0L, rho = 1,
                     parent_tip_start = emphasis:::.pts(brts))
wr <- 0; wl <- 0; np <- 0L
for (df in aug$trees) {
  ab <- e$.alive_before(df)
  for (i in seq_len(nrow(df))) {
    prev <- if (i == 1L) 0 else df$brts[i - 1L]
    if (df$brts[i] <= prev) next
    t <- seq(prev, df$brts[i], length.out = 9L)[-c(1L, 9L)]
    got <- eval_nh_rate(p, df, t, model = c(0L, 0L, 1L), link = 0L, rho = 1)$pd
    wr <- max(wr, max(abs(got - (length(ab[[i]]) * t - sum(ab[[i]])))))
    wl <- max(wl, max(abs(got - (df$pd[i] + df$n[i] * (t - df$brts[i])))))
    np <- np + length(t)
  }
}
cat(sprintf("gate 4  %d interior times: P vs replayed alive set           max |diff| = %.3e\n", np, wr))
cat(sprintf("gate 4  same, P vs node.pd + node.n * (t - node.brts):       max |diff| = %.3e\n", wl))
