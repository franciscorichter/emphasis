## H64 isolation: remove the H62 early-stop race (sample_size == maxN, so every
## accepted attempt is kept and completion order cannot select) and compare
## grainsize 0 (nt=10, maxN=6) vs grainsize 1 (nt=6, maxN=6) vs serial (nt=1).
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
hc <- parallel::detectCores()
set.seed(6464)
tr <- ape::rphylo(30, birth = 0.6, death = 0.2)
brts <- sort(as.numeric(ape::branching.times(tr)), decreasing = TRUE)

aug <- function(ss, maxN, nt, pars, model = c(0L,0L,0L), link = 0L, rho = 1)
  emphasis:::augment_trees(brts, pars, as.integer(ss), as.integer(maxN),
                           10000L, 1e6, as.integer(nt), model, link, rho)

cfgs <- list(
  dd_lin   = list(pars = c(0.8,-0.01,0,0, 0.1,0,0,0), model = c(1L,0L,0L), link = 0L, rho = 1),
  cr_rho06 = list(pars = c(0.6,0,0,0, 0.1,0,0,0), model = c(0L,0L,0L), link = 0L, rho = 0.6)
)
k <- 6L
for (nm in names(cfgs)) {
  cf <- cfgs[[nm]]
  res <- lapply(c(1L, 6L, hc), function(nt) {
    rs <- replicate(200, aug(k, k, nt, cf$pars, cf$model, cf$link, cf$rho), simplify = FALSE)
    list(ntrees = vapply(rs, function(r) length(r$logf), 1L),
         lw = unlist(lapply(rs, function(r) r$logf - r$logg)))
  })
  names(res) <- c("nt1", "nt6", paste0("nt", hc))
  for (i in seq_along(res))
    cat(sprintf("%-9s %-4s grain=%d : ntrees=%s  lw mean=%.3f sd=%.3f n=%d\n", nm, names(res)[i],
                k %/% c(1L, 6L, hc)[i],
                paste(names(table(res[[i]]$ntrees)), table(res[[i]]$ntrees), sep=":", collapse=" "),
                mean(res[[i]]$lw), sd(res[[i]]$lw), length(res[[i]]$lw)))
  cat(sprintf("%-9s t.test nt1 vs nt%d (grain 0) p=%.3f | nt6 (grain 1) vs nt%d (grain 0) p=%.3f | nt1 vs nt6 p=%.3f | KS nt6 vs nt%d p=%.3f\n\n",
              nm, hc, t.test(res$nt1$lw, res[[3]]$lw)$p.value,
              hc, t.test(res$nt6$lw, res[[3]]$lw)$p.value,
              t.test(res$nt1$lw, res$nt6$lw)$p.value,
              hc, suppressWarnings(ks.test(res$nt6$lw, res[[3]]$lw)$p.value)))
}
## and the confounded design again (ss=2, maxN=6) but grain 1 (nt=6) vs grain 0 (nt=10): both have the race
for (nm in names(cfgs)) {
  cf <- cfgs[[nm]]
  res <- lapply(c(6L, hc), function(nt) {
    rs <- replicate(300, aug(2L, 6L, nt, cf$pars, cf$model, cf$link, cf$rho), simplify = FALSE)
    list(ntrees = vapply(rs, function(r) length(r$logf), 1L),
         lw = unlist(lapply(rs, function(r) r$logf - r$logg)))
  })
  cat(sprintf("race  %-9s nt=6 grain=1: ntrees=%s lw mean=%.3f | nt=%d grain=0: ntrees=%s lw mean=%.3f | t p=%.3f\n", nm,
              paste(names(table(res[[1]]$ntrees)), table(res[[1]]$ntrees), sep=":", collapse=" "), mean(res[[1]]$lw),
              hc, paste(names(table(res[[2]]$ntrees)), table(res[[2]]$ntrees), sep=":", collapse=" "), mean(res[[2]]$lw),
              t.test(res[[1]]$lw, res[[2]]$lw)$p.value))
}
