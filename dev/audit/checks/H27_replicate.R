## H27 replicate — independent re-run on a DIFFERENT tree (35 tips,
## lambda=0.8, mu=0.5, seed 2027), 2 pipeline replicates, plus a
## stages-subset run (bounds+gam only) that exercises the "GAM surrogate
## becomes the returned loglik/AIC" fallback path directly.
## Common yardstick: 300 BDI trees (zero-variance under CR) + DDD::bd_loglik.
## Budget: < 5 min single-threaded.

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages({ library(emphasis); library(ape) })

set.seed(2027)
tree <- ape::rphylo(35, birth = 0.8, death = 0.5)
brts <- emphasis:::.extract_brts(tree)
model_bin <- c(0L, 0L, 0L); link <- 0L
cat(sprintf("tree: %d tips, crown age %.3f\n", Ntip(tree), brts[1]))

exact_ll <- function(p, cond = 1L)
  DDD::bd_loglik(pars1 = c(p[1], p[2]), pars2 = c(0, cond, 0, 0, 2),
                 brts = brts, missnumspec = 0)

bdi_eval <- function(pars_compact, cond_fun) {
  p8 <- emphasis:::.expand_pars(as.numeric(pars_compact), model_bin)
  bdi <- emphasis:::.augment_tree_bdi(brts, p8, model_bin, sample_size = 300L,
                                      link = link, rho = 1)
  corr <- if (is.null(cond_fun)) 0 else cond_fun(p8)
  c(bdi = bdi$fhat - corr, sd_lw = stats::sd(bdi$logf - bdi$logg),
    exact_c1 = exact_ll(pars_compact, 1L))
}

rows <- list(); init_rows <- list()
t_all <- proc.time()[3]
for (r in 1:2) {
  t0 <- proc.time()[3]
  pl <- emphasis_pipeline(tree, model = "cr", link = "linear",
                          control = list(num_threads = 1L, max_time = 120,
                                         mcem = list(max_iter = 20L)),
                          verbose = FALSE)
  cat(sprintf("\n[rep %d] %.0fs best_stage=%s\n", r, proc.time()[3] - t0, pl$best_stage))
  print(pl)   # what the user actually sees
  cond_fun <- if (!is.null(pl$bounds$survival_gam))
    emphasis:::.build_cond_fun(pl$bounds$survival_gam, model_bin) else NULL
  for (st in c("gam", "cem", "mcem")) {
    f <- pl$fits[[st]]
    if (is.null(f) || !is.finite(f$loglik)) next
    ce <- bdi_eval(f$pars, cond_fun)
    rows[[length(rows) + 1L]] <- data.frame(rep = r, stage = st,
      lambda = f$pars[1], mu = f$pars[2], logged = f$loglik, t(ce),
      gap = f$loglik - ce[["bdi"]])
  }
  ll_cem <- pl$fits$cem$loglik; ll_gam <- pl$fits$gam$loglik
  init_rows[[r]] <- data.frame(rep = r,
    logged_cem = ll_cem, logged_gam = ll_gam,
    exact_cem = exact_ll(pl$fits$cem$pars), exact_gam = exact_ll(pl$fits$gam$pars))
}
res <- do.call(rbind, rows); rownames(res) <- NULL
cat("\n=== logged vs BDI-300 (same cond correction) vs DDD exact ===\n")
print(res, digits = 5)
for (r in 1:2) {
  s <- res[res$rep == r, ]
  cat(sprintf("rep %d ranking logged: %s | bdi: %s | exact: %s\n", r,
              paste(s$stage[order(-s$logged)], collapse = ">"),
              paste(s$stage[order(-s$bdi)], collapse = ">"),
              paste(s$stage[order(-s$exact_c1)], collapse = ">")))
}
cat("\n=== init rule inputs ===\n"); print(do.call(rbind, init_rows), digits = 5)

## --- fallback path: GAM is the only estimating stage -> its surrogate is the
## returned loglik/AIC. Reuse rep-2 bounds to save time.
cat("\n=== bounds+gam only: returned loglik is the GAM surrogate ===\n")
ab <- pl$bounds
t0 <- proc.time()[3]
pg <- emphasis_pipeline(tree, model = "cr", link = "linear",
                        stages = c("gam"),
                        control = list(num_threads = 1L, max_time = 120,
                                       lower_bound = ab$lower_bound,
                                       upper_bound = ab$upper_bound,
                                       survival_gam = ab$survival_gam),
                        verbose = FALSE)
cat(sprintf("(%.0fs) best_stage=%s returned loglik=%.3f AIC=%.3f\n",
            proc.time()[3] - t0, pg$best_stage, pg$loglik, pg$AIC))
ce <- bdi_eval(pg$pars, cond_fun)
cat(sprintf("BDI at returned pars=%.3f  gap(logged - bdi)=%+.3f  exact_c1=%.3f\n",
            ce[["bdi"]], pg$loglik - ce[["bdi"]], ce[["exact_c1"]]))
cat(sprintf("\ntotal %.0fs\n", proc.time()[3] - t_all))
