# A tree that is not ultrametric has branching times after the present; the
# estimator refuses it instead of scoring the present as an observed split.

test_that(".extract_brts refuses negative branching times and names the fix", {
  phy <- ape::read.tree(text = "(((a:5,b:1):1,c:1):1,d:1);")   # tips at uneven depths: ape dates one node at -1
  expect_lt(min(ape::branching.times(phy)), 0)
  expect_error(emphasis:::.extract_brts(phy), "not ultrametric")
  expect_error(emphasis:::.extract_brts(c(3, 1, -0.5)), "not ultrametric")
})

test_that("an ultrametric tree and a tree with a tiny rounding offset pass", {
  phy <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  b <- emphasis:::.extract_brts(phy)
  expect_equal(as.numeric(b), c(2, 1, 1))
  expect_false(is.null(attr(b, "parent_tip_start")))
  expect_silent(emphasis:::.extract_brts(c(2, 1, -1e-12)))
})
