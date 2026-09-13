lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, .libPaths())); library(emphasis)
A <- emphasis:::.augment_tree_bdi
brts <- c(6, 4.27428208320388, 4.1258975168475, 3.72207978339226, 3.64993282835291,
  3.17502914508657, 2.81163149722479, 2.48700670286076, 2.38274529774947,
  1.94421117175661, 1.86664932095232, 1.81164887493753, 1.38503498765192,
  1.06976980095552, 0.460758774034712, 0.417936203054209, 0.363178692045526,
  0.253829586012209, 0.0517782331486725)
set.seed(77); a <- A(brts, c(0.5,0.1), c(0L,0L,0L), sample_size=25L, max_missing=1e4L, link=0L, rho=1)
cat(sprintf("CR  fhat=%.12f  logg[1]=%.12f  ntrees=%d\n", a$fhat, a$logg[1], length(a$trees)))
set.seed(78); b <- A(brts, c(0.8,-0.025,0.3,0), c(1L,0L,0L), sample_size=60L, max_missing=1e4L, link=0L, rho=1)
cat(sprintf("DD  fhat=%.6f  logg[1]=%.12f  ntrees=%d\n", b$fhat, b$logg[1], length(b$trees)))
# rho < 1: the log(acc) factor is applied even though survivors are not f = 0 there
set.seed(79); c1 <- A(brts, c(0.8,-0.025,0.3,0), c(1L,0L,0L), sample_size=60L, max_missing=1e4L, link=0L, rho=0.5)
cat(sprintf("DD rho=0.5 fhat=%.6f acc=%s\n", c1$fhat, format(if (is.null(c1$acc)) NA else c1$acc)))
st <- tryCatch({ set.seed(80)
  tr <- emphasis::simulate_tree(c(0.5,0.1), age = 5, model = "cr")
  "ok" }, error = function(e) conditionMessage(e))
cat("simulate_tree:", st, "\n")
