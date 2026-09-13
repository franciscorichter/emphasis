.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages({library(emphasis); library(ape); library(DDD)})
msg <- function(...) cat(sprintf(...), file = stderr())
set.seed(7); sdd <- DDD::dd_sim(c(0.8, 0.1, 30), 8)
brts_dd <- sort(sdd$brts, decreasing = TRUE); Nobs <- length(brts_dd) + 1L
msg("dd tree Nobs=%d. Why are K<Nobs particles NA?\n", Nobs)
for (K in c(6, 10, 15, 18, 20, 22, 23)) for (mN in c(10L, 200L)) {
  b0 <- 0.8; bN <- -b0/K
  r <- tryCatch(emphasis:::augment_trees(brts_dd, c(b0,bN,0,0,0.1,0,0,0), 5L, mN, 1e4L, 1e6, 1L, c(1L,0L,0L), 0L, 1.0), error=function(e) e)
  if (inherits(r,"error")) msg("K=%4.1f maxN=%3d -> ERROR: %s\n", K, mN, conditionMessage(r))
  else { lw <- r$logf - r$logg
    msg("K=%4.1f maxN=%3d -> ok: zw=%d ov=%d  logf=%s  fhat=%.2f\n", K, mN, r$rejected_zero_weights, r$rejected_overruns,
        paste(round(r$logf,1), collapse=","), emphasis:::.is_fhat(r$logf, r$logg, n_zero_weight=r$rejected_zero_weights)) }
}
# Is the true likelihood zero for K < Nobs-1?  Score trees drawn at K=30 under K=20 pars.
r30 <- emphasis:::augment_trees(brts_dd, c(0.8,-0.8/30,0,0,0.1,0,0,0), 5L, 200L, 1e4L, 1e6, 1L, c(1L,0L,0L), 0L, 1.0)
for (K in c(20, 22, 23, 30)) {
  ev <- emphasis:::eval_logf(c(0.8,-0.8/K,0,0,0.1,0,0,0), r30$trees, c(1L,0L,0L), 0L, 1.0)
  msg("logf of K=30-drawn trees scored at K=%d: %s\n", K, paste(round(unlist(ev$logf),1), collapse=", "))
}
msg("\nauto_bounds box for the cr tree (pipeline's CEM box):\n")
set.seed(26); repeat { phy <- rlineage(0.5, 0.2, Tmax = 8); ex <- drop.fossil(phy); if (Ntip(ex) >= 15 && Ntip(ex) <= 30) break }
ab <- tryCatch(auto_bounds(ex, model = "cr", verbose = FALSE), error = function(e) e)
if (inherits(ab, "error")) msg("auto_bounds error: %s\n", conditionMessage(ab)) else {
  msg("lower = %s\nupper = %s\n", paste(round(ab$lower_bound,3), collapse=", "), paste(round(ab$upper_bound,3), collapse=", "))
  brts <- emphasis:::.extract_brts(ex)
  corners <- rbind(c(ab$upper_bound[1], ab$upper_bound[2]), c(ab$upper_bound[1], ab$lower_bound[2]),
                   c(ab$lower_bound[1], ab$upper_bound[2]), c(mean(ab$lower_bound[1], ab$upper_bound[1]), ab$upper_bound[2]))
  for (i in seq_len(nrow(corners))) { lam <- corners[i,1]; mu <- corners[i,2]; t0 <- proc.time()[3]
    o <- replicate(3, { r <- tryCatch(emphasis:::augment_trees(brts, c(lam,0,0,0,mu,0,0,0), 5L, 10L, 1e4L, 1e6, 1L, c(0L,0L,0L), 0L, 1.0), error=function(e) e)
                        if (inherits(r,"error")) conditionMessage(r) else sprintf("ok zw=%d ov=%d", r$rejected_zero_weights, r$rejected_overruns) })
    msg("corner lam=%.3f mu=%.3f (%.1fs/3 calls): %s\n", lam, mu, proc.time()[3]-t0, paste(unique(o), collapse=" | ")) }
}
