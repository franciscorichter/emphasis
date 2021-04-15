context("locate_plugin")

testthat::test_that("test_locate_plugin", {
  testthat::skip_on_os("solaris") # solaris is often very bitchy, better avoid it.
  
  testthat::expect_warning(locate_plugin("noname"))

 # uncomment code once the library below is available (e.g. on a non-private github)
 # library(remphasis_rpd1)
 # testthat::expect_silent(locate_plugin("rpd1"))
 
 # you can also add tests where you check the outcome is right:
 testthat::expect_equal(1 + 1, 2)
 # or perhaps more interesting:
 set.seed(42)
 phy <- TreeSim::sim.bd.age(age = 5, numbsim = 1, lambda = 1, mu = 0)[[1]]
 lambda_estim <- ape::birthdeath(phy)$para[[2]] # this is the lambda estimate
 testthat::expect_equal(1, lambda_estim, tolerance = 0.5)
 
})
