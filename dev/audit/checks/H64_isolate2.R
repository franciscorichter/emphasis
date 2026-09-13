.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
hc <- parallel::detectCores()
set.seed(6464)
tr <- ape::rphylo(30, birth = 0.6, death = 0.2)
brts <- sort(as.numeric(ape::branching.times(tr)), decreasing = TRUE)
pars <- c(0.6,0,0,0, 0.1,0,0,0)
aug <- function(nt) emphasis:::augment_trees(brts, pars, 6L, 6L, 10000L, 1e6, as.integer(nt), c(0L,0L,0L), 0L, 0.6)
for (rep in 1:3) {
  lw <- lapply(c(6L, hc), function(nt) unlist(lapply(replicate(400, aug(nt), simplify = FALSE), function(r) r$logf - r$logg)))
  cat(sprintf("rep %d cr_rho06 ss=maxN=6: nt6(grain1) mean=%.3f sd=%.3f n=%d | nt%d(grain0) mean=%.3f sd=%.3f n=%d | t p=%.3f KS p=%.3f\n",
      rep, mean(lw[[1]]), sd(lw[[1]]), length(lw[[1]]), hc, mean(lw[[2]]), sd(lw[[2]]), length(lw[[2]]),
      t.test(lw[[1]], lw[[2]])$p.value, suppressWarnings(ks.test(lw[[1]], lw[[2]])$p.value)))
}
