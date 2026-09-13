# Item 1.2 review — run m_cpp on the fixed inputs; invoked with the library path as arg 1.
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(emphasis))
ns <- asNamespace("emphasis"); m_cpp <- ns$m_cpp; eval_logf <- ns$eval_logf
cfg <- readRDS("/Users/pancho/Code/emphasis/dev/audit/review/1.2_inputs.rds")

log_psurv <- function(lam, mu, T) {
  if (abs(lam - mu) < 1e-8) mu <- lam - 1e-8
  r <- lam - mu; p0 <- mu * (1 - exp(-r * T)) / (lam - mu * exp(-r * T))
  2 * log(max(1 - p0, 1e-300))
}
run <- function(c, w = c$w, trees = c$trees, init = c$init, cond = NULL, nt = 1L, xtol = 1e-6, lb = c$lb, ub = c$ub) {
  es <- list(trees = trees, weights = w, rejected = 0L, rejected_overruns = 0L, rejected_lambda = 0L,
             rejected_zero_weights = 0L, time = 0, fhat = 0)
  r <- tryCatch(m_cpp(e_step = es, init_pars = init, plugin = "rpd1", lower_bound = lb, upper_bound = ub,
                      xtol_rel = xtol, num_threads = nt, model = c$bin, link = c$link, rho = 1,
                      rconditional = cond), error = function(e) list(estimates = rep(NA_real_, 8), nlopt = NA, err = conditionMessage(e)))
  list(est = as.numeric(r$estimates), nlopt = r$nlopt, err = r$err, names = names(r))
}
out <- list()
for (nm in names(cfg)) {
  c <- cfg[[nm]]
  out[[paste0(nm, "|nocond")]] <- run(c)
  if (!is.null(c$Tc)) {
    cond <- function(p8) if (c$link == 1L) log_psurv(exp(p8[1]), exp(p8[5]), c$Tc) else log_psurv(p8[1], p8[5], c$Tc)
    out[[paste0(nm, "|cond")]] <- run(c, cond = cond)
    out[[paste0(nm, "|cond|w/7")]] <- run(c, w = c$w / 7, cond = cond)
    out[[paste0(nm, "|nocond|w/7")]] <- run(c, w = c$w / 7)
  }
}
# threads
out[["cr_bdi|nocond|nt4"]] <- run(cfg$cr_bdi, nt = 4L)
out[["cr_bdi|cond|nt4"]] <- run(cfg$cr_bdi, nt = 4L, cond = function(p8) log_psurv(p8[1], p8[5], cfg$cr_bdi$Tc))
# rates near zero
c0 <- cfg$cr_bdi
out[["cr_bdi|init1e-6|lb1e-8"]] <- run(c0, init = c(1e-6,0,0,0,1e-6,0,0,0), lb = c(1e-8,0,0,0,0,0,0,0))
out[["cr_bdi|init_ub"]] <- run(c0, init = c(5,0,0,0,5,0,0,0))
# weight edge cases
out[["cr_bdi|w_all0|nocond"]] <- run(c0, w = c0$w * 0)
out[["cr_bdi|w_all0|cond"]] <- run(c0, w = c0$w * 0, cond = function(p8) log_psurv(p8[1], p8[5], c0$Tc))
wn <- c0$w; wn[1] <- NaN
out[["cr_bdi|w_NaN"]] <- run(c0, w = wn)
wi <- c0$w; wi[1] <- -Inf
out[["cr_bdi|w_-Inf"]] <- run(c0, w = wi)
wneg <- c0$w; wneg[1] <- -1
out[["cr_bdi|w_neg"]] <- run(c0, w = wneg)
# dd: -Inf tree at w = 0 and w = 1, threads 1 and 4; all-infeasible init
d <- cfg$dd_bdi
for (nt in c(1L, 4L)) {
  out[[sprintf("dd_bdi|fin|nt%d", nt)]] <- run(d, nt = nt, xtol = 1e-4)
  out[[sprintf("dd_bdi|bad_w0|nt%d", nt)]] <- run(d, trees = c(d$trees, list(d$tree_bad)), w = c(d$w, 0), nt = nt, xtol = 1e-4)
  out[[sprintf("dd_bdi|bad_w1|nt%d", nt)]] <- run(d, trees = c(d$trees, list(d$tree_bad)), w = c(d$w, 1), nt = nt, xtol = 1e-4)
}
# tree_bad alone (weight 1), and init where every weighted tree has zero density
out[["dd_bdi|bad_only"]] <- run(d, trees = list(d$tree_bad), w = 1, xtol = 1e-4)
out[["dd_bdi|infeasible_init"]] <- run(d, init = c(0.5, -0.1, 0,0, 0.4, 0,0,0), xtol = 1e-4)
# Q at the estimate on the dd finite set, for cross-build comparison
Qfin <- function(th) sum(d$w * eval_logf(th, d$trees, model = d$bin, link = 0L, rho = 1)$logf)
for (k in c("dd_bdi|fin|nt1", "dd_bdi|bad_w0|nt1", "dd_bdi|bad_w1|nt1", "dd_bdi|fin|nt4", "dd_bdi|bad_w0|nt4"))
  out[[k]]$Q <- Qfin(out[[k]]$est)
# logf of tree_bad at bad_w1 estimate
out[["dd_bdi|bad_w1|nt1"]]$logf_bad_at_est <- eval_logf(out[["dd_bdi|bad_w1|nt1"]]$est, list(d$tree_bad), model = d$bin, link = 0L, rho = 1)$logf
out[["dd_bdi|infeasible_init"]]$n_finite_at_est <- sum(is.finite(eval_logf(out[["dd_bdi|infeasible_init"]]$est, d$trees, model = d$bin, link = 0L, rho = 1)$logf))
saveRDS(out, sub("/$", "", paste0("/Users/pancho/Code/emphasis/dev/audit/review/1.2_out_", basename(lib), ".rds")))
for (k in names(out)) {
  o <- out[[k]]
  cat(sprintf("%-32s nlopt=%-3s est=%s%s%s\n", k, as.character(o$nlopt),
              paste(formatC(o$est[o$est != 0 | seq_along(o$est) %in% c(1,5)], digits = 6, format = "g"), collapse = " "),
              if (!is.null(o$Q)) sprintf("  Q=%.6f", o$Q) else "", if (!is.null(o$err)) paste0("  ERR: ", o$err) else ""))
}
cat("names(m_cpp result):", out[[1]]$names, "\n")
