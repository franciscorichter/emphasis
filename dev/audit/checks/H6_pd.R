## H6 (side check): why does pd credited on the finished node differ from the
## pd the sampler used?  Prediction from augment_tree.cpp:50 —
## compute_pendant_pd() resets tip_start to 0 for every non-extinction node
## with parent_id == -1, which includes augmented species born before the
## first observed branching (alive_ids empty -> chosen_parent_id = -1,
## augment_tree.cpp:162-166).  During augmentation their tip_start was t_spec
## (make_node, line 25).  So pd_final(t) - pd_sampler(t) should equal the sum of
## brts over such species alive at t (including the node itself).
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
is_ext <- function(te) te == 0; is_tip <- function(te) te == 10e10
is_mis <- function(te) !(is_ext(te) | is_tip(te) | te == 5e10)
brts <- c(5, 3.6, 2.9, 2.2, 1.5, 0.9, 0.4)
E <- emphasis:::augment_trees(brts = brts, pars = c(log(0.6),0,0,0,log(0.25),0,0,0), sample_size = 200L,
                              maxN = 40000L, max_missing = 10000L, max_lambda = 500, num_threads = 1L,
                              model = c(0L,0L,0L), link = 1L, rho = 1)
resid <- c(); n_rows <- 0; n_orphans <- 0
for (d in E$trees) {
  sp <- !is_ext(d$t_ext); mis <- is_mis(d$t_ext)
  ts_sampler <- ifelse(mis, d$brts, 0)
  orphan <- mis & d$parent_id == -1
  n_orphans <- n_orphans + sum(orphan)
  stopifnot(all(d$tip_start[orphan] == 0))                       # the reset happened
  stopifnot(all(d$tip_start[mis & !orphan] == d$brts[mis & !orphan]))
  for (i in which(mis)) {
    t <- d$brts[i]
    alive <- sp & d$brts <= t & d$t_ext > t
    pd_s <- sum(t - ts_sampler[alive])
    pred <- sum(d$brts[alive & orphan])                          # predicted excess
    resid <- c(resid, d$pd[i] - pd_s - pred); n_rows <- n_rows + 1
  }
}
cat(sprintf("augmented nodes: %d, orphan (parent_id == -1) species: %d\n", n_rows, n_orphans))
cat(sprintf("max |pd_final - pd_sampler - sum(brts of orphan species alive)| = %.2e\n", max(abs(resid))))
