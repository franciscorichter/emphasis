## H61: 2-tip tree -> seq(0L, -1L) assigned to a zero-length index in .bdi_to_tree_df.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape)})
# The R idiom itself
x <- integer(3); r <- tryCatch({x[integer(0)] <- seq(0L, -1L); "no error"}, warning = function(w) paste("warning:", conditionMessage(w)), error = function(e) paste("error:", conditionMessage(e)))
cat("x[integer(0)] <- seq(0L,-1L):", r, "; x =", x, "\n")

tr2 <- ape::read.tree(text = "(A:5,B:5);")
brts <- sort(ape::branching.times(tr2), decreasing = TRUE)
cat("brts:", brts, "\n")
pars8 <- c(0.5, 0,0,0, 0.3, 0,0,0)
a <- tryCatch(emphasis:::.augment_tree_bdi(brts, pars8, c(0L,0L,0L), sample_size = 20L),
              error = function(e) {cat("ERROR:", conditionMessage(e), "\n"); NULL})
if (!is.null(a)) {
  cat("n trees:", length(a$trees), " logf-logg:", round(range(a$weights), 6), "\n")
  d <- a$trees[[which.max(sapply(a$trees, nrow))]]
  print(d)
  cat("ids of observed rows (should be none):", d$id[d$t_ext == 1e11 & d$brts < brts[1]], "\n")
}
# Direct df conversion with an empty species list
print(emphasis:::.bdi_to_tree_df(list(), numeric(0), 5))
# DD path on the 2-tip tree
b <- tryCatch(emphasis:::.augment_tree_bdi(brts, c(0.6,-0.01,0,0,0.1,0,0,0), c(1L,0L,0L), sample_size = 10L),
              error = function(e) {cat("DD ERROR:", conditionMessage(e), "\n"); NULL})
if (!is.null(b)) cat("DD: n trees:", length(b$trees), " fhat:", b$fhat, "\n")
# Full MCEM on a 2-tip tree
m <- tryCatch(emphasis:::.mcem_bdi(brts, pars = c(0.5, 0.3), sample_size = 20L, max_missing = 1e4L,
                          lower_bound = c(0.01, 0.01), upper_bound = c(3, 3),
                          max_iter = 3L, xtol = 1e-3, tol = 1e-3, patience = 2L,
                          num_threads = 1L, model = c(0L,0L,0L), link = 0L),
              error = function(e) {cat("MCEM ERROR:", conditionMessage(e), "\n"); NULL})
if (!is.null(m)) { print(m$mcem); cat("stop_reason:", m$stop_reason, "\n") }
st <- tryCatch(simulate_tree(tree = tr2, pars = c(0.5, 0.3), model = "cr", n_trees = 5L, method = "bdi"),
               error = function(e) {cat("simulate_tree ERROR:", conditionMessage(e), "\n"); NULL})
if (!is.null(st)) { cat("simulate_tree(bdi) on 2-tip: tas classes:", paste(sapply(st$trees, function(x) class(x)[1]), collapse=","), " log_q:", round(st$log_q,3), "\n")
  ok <- !sapply(st$trees, is.null); if (any(ok)) { t1 <- st$trees[[which(ok)[1]]]; cat("tas tips:", length(t1$tip.label), " min edge:", min(t1$edge.length), "\n") } }
# Compare against the thinning sampler on the same 2-tip tree
th <- tryCatch(emphasis:::augment_trees(brts, pars8, sample_size = 5L, maxN = 100L, max_missing = 1e4,
                             max_lambda = 500, num_threads = 1L),
               error = function(e) {cat("thinning ERROR:", conditionMessage(e), "\n"); NULL})
if (!is.null(th)) cat("thinning on 2-tip: n trees", length(th$trees), " rows[1]:", nrow(th$trees[[1]]), "\n")
# Is the dangling parent_id / negative-edge tas specific to the 2-tip case? Compare a 20-tip tree.
set.seed(7); tr20 <- ape::rphylo(20, 0.5, 0.1)
s20 <- simulate_tree(tree = tr20, pars = c(0.5, 0.3), model = "cr", n_trees = 5L, method = "bdi")
cat("20-tip bdi: min edge per tas:", round(sapply(s20$trees, function(t) if (is.null(t)) NA else min(t$edge.length)), 3), "\n")
a20 <- emphasis:::.augment_tree_bdi(tree = sort(ape::branching.times(tr20), decreasing = TRUE), c(0.5,0,0,0,0.3,0,0,0), c(0L,0L,0L), sample_size = 3L)
cat("20-tip: observed ids present:", paste(sort(unique(a20$trees[[1]]$id[a20$trees[[1]]$t_ext == 1e11])), collapse=","), "\n")
cat("2-tip: parent_id values of augmented rows:", unique(a$trees[[1]]$parent_id[a$trees[[1]]$id >= 0]), " (observed node ids present: none) -> dangling\n")
