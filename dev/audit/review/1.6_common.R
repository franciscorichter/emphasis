.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
brts22 <- c(2.984340, 2.456040, 2.041011, 1.919104, 1.703097, 1.313093, 0.939284, 0.715800, 0.553944, 0.455806, 0.427175, 0.308270, 0.174959, 0.141573, 0.112195, 0.069056, 0.058366, 0.044335, 0.035010)
brts1 <- c(5.000000, 3.822147, 3.059192, 1.966150, 1.604249, 1.358484, 0.622046, 0.467531, 0.393533, 0.285495, 0.276148, 0.225482, 0.201295, 0.154060, 0.057677, 0.044974, 0.041380, 0.038051, 0.013895)
ex <- function(p, mb = c(0L,0L,0L)) emphasis:::.expand_pars(p, mb)
drv <- function(brts, init, lb, ub, mb = c(0L,0L,0L), link = 0L, ...) {
  a <- list(brts = brts, pars = ex(init, mb), sample_size = 20L, maxN = 2000L, max_missing = 1e4,
            lower_bound = ex(lb, mb), upper_bound = ex(ub, mb), max_iter = 10L, xtol = 1e-3,
            tol = 1e-2, patience = 3L, num_threads = 1L, verbose = FALSE, model = mb, link = link, max_time = 120)
  a[names(list(...))] <- list(...)
  do.call(emphasis:::.mcem_dynamic_fresh, a)
}
invariants <- function(d, init8, label) {
  m <- d$mcem; ok <- character(0); bad <- character(0)
  chk <- function(cond, what) if (isTRUE(cond)) ok <<- c(ok, what) else bad <<- c(bad, what)
  if (is.null(m)) { cat(sprintf("[%s] mcem NULL; stop=%s iter=%d n_failed=%d\n", label, d$stop_reason, d$iterations, d$n_failed)); return(invisible()) }
  last <- m[nrow(m), ]
  chk(identical(last$final_estep, d$final_estep), "final_estep flag == last row")
  if (d$final_estep) {
    chk(isTRUE(all.equal(as.numeric(last[, paste0("par",1:8)]), as.numeric(d$pars))), "last row pars == pars")
    chk(nrow(m) == d$iterations + 1L, "nrow == iterations+1")
    chk(isTRUE(all.equal(d$final_IS$fhat, last$fhat)), "final_IS$fhat == last fhat")
    chk(is.na(last$delta_max), "last delta NA")
  } else chk(nrow(m) == d$iterations, "nrow == iterations (no final)")
  chk(isTRUE(all.equal(m$n_rejected, rowSums(m[, c("rejected","rejected_overruns","rejected_lambda")]))), "n_rejected == rowSums")
  chk(last$n_rejected == d$final_IS$n_rejected, "last n_rejected == final_IS")
  mm <- m[!m$final_estep, ]
  P <- as.matrix(mm[, paste0("par",1:8)])
  prevP <- rbind(init8, P[-nrow(P), , drop = FALSE])
  ref <- sapply(seq_len(nrow(P)), function(k) max(abs(P[k,] - prevP[k,]) / pmax(abs(prevP[k,]), 1e-2)))
  chk(isTRUE(all.equal(mm$delta_max, ref)) || d$n_failed > 0, "delta recomputed (or failures present)")
  chk(all(diff(m$maxN) >= 0), "maxN monotone")
  chk(is.finite(d$loglik_var) || length(d$final_IS$logf) < 2, "loglik_var finite")
  cat(sprintf("[%s] stop=%s iter=%d n_failed=%d nrow=%d final=%s maxN=%d pars=(%s) fhat=%.3f var=%.4f ESS=%.1f nrej=%d\n   OK: %d  BAD: %s\n",
      label, d$stop_reason, d$iterations, d$n_failed, nrow(m), d$final_estep, d$maxN,
      paste(signif(emphasis:::.contract_pars(d$pars, attr(d,"mb") %||% c(0L,0L,0L)), 3), collapse=","),
      last$fhat, d$loglik_var, d$final_IS$ESS, d$final_IS$n_rejected, length(ok), if (length(bad)) paste(bad, collapse="; ") else "-"))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
