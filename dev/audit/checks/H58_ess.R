## H58 follow-up: IS efficiency when .bdi_iterate did not converge (40-tip high-turnover DD)
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape)})
set.seed(3); tr20 <- ape::rphylo(20, 0.5, 0.1); tr40 <- ape::rphylo(40, 0.5, 0.1)
ess <- function(lw) { lw <- lw[is.finite(lw)]; if (!length(lw)) return(NA); w <- exp(lw - max(lw)); sum(w)^2 / sum(w^2) }
for (cfg in list(list(tr20, c(1.0,-0.05,0,0,0.2,0,0,0), "strong DD 20 tips (converged)"),
                 list(tr40, c(2.0,-0.05,0,0,1.0,0,0,0), "high-turnover DD 40 tips (NOT converged)"))) {
  brts <- sort(ape::branching.times(cfg[[1]]), decreasing = TRUE)
  t0 <- proc.time()[3]
  a <- emphasis:::.augment_tree_bdi(brts, cfg[[2]], c(1L,0L,0L), sample_size = 50L)
  cat(sprintf("%-42s trees=%d finite(logf)=%d finite(logg)=%d ESS(finite)=%.1f sd(lw)=%.2f fhat=%.2f time=%.0fs\n", cfg[[3]], length(a$trees),
              sum(is.finite(a$logf)), sum(is.finite(a$logg)), ess(a$weights), sd(a$weights[is.finite(a$weights)]), a$fhat, proc.time()[3]-t0))
  cat("   non-finite logf:", head(a$logf[!is.finite(a$logf)]), " non-finite logg:", head(a$logg[!is.finite(a$logg)]), "\n")
}
