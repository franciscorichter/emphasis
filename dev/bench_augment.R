# Augmentation cost against tree size (audit H53).
#
# Model::nh_rate used to recompute the pendant PD by scanning the event list at
# every candidate time, O(N) per call.  It now reads P off the node that governs
# the segment, O(1), and the sweep that keeps node.pd current is carried forward
# with the sampler rather than re-derived.
#
# Measured here on cr, whose rates never read pd, so both builds draw the same
# trees and the timings compare (min of 3, 10 trees per call, one thread):
#
#   tips   nodes     before     after
#    400     ~725   14.30 ms   11.70 ms
#    800    ~1590   60.30 ms   45.00 ms
#   1600    ~2900  191.20 ms  131.20 ms
#
# Still superlinear: insert_species rewrites n over the dirty range and the
# parent draw scans the alive lineages, both O(N) per accepted candidate.  What
# is gone is the pendant-PD scan, which ran per candidate, accepted or not.
#
#   Rscript dev/bench_augment.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

bench <- function(n_tips, model, pars, reps = 3L) {
  set.seed(n_tips)
  phy  <- ape::rphylo(n_tips, 0.8, 0.2)
  brts <- emphasis:::.extract_brts(phy)
  pts  <- emphasis:::.pts(brts)
  draw <- function() augment_trees(as.numeric(brts), pars, 10L, 500000L, 4000L,
                                   1e6, 1L, model = as.integer(model), link = 0L,
                                   rho = 1, parent_tip_start = pts)
  t <- replicate(reps, system.time(draw())[["elapsed"]])
  a <- draw()
  c(nodes = mean(vapply(a$trees, nrow, 0L)), ms_per_tree = 1000 * min(t) / 10)
}

for (spec in list(list("cr", c(0L, 0L, 0L), c(0.7, 0, 0, 0, 0.25, 0, 0, 0)),
                  list("d ", c(0L, 0L, 1L), c(0.7, 0, 0, 0.02, 0.25, 0, 0, 0.005)))) {
  for (n in c(50L, 100L, 200L, 400L, 800L, 1600L)) {
    r <- bench(n, spec[[2L]], spec[[3L]])
    cat(sprintf("%s  tips %4d  mean nodes %7.1f  %8.3f ms/tree\n",
                spec[[1L]], n, r[["nodes"]], r[["ms_per_tree"]]))
  }
}
