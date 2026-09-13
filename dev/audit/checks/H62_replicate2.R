## H62 replication part 2: em_cpp (does the M-step receive all M trees; is fhat off by log(M/N))
## and estimate_rates at pipeline level (bdi untouched; dynamic_fresh with num_threads>1 inflated).
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
set.seed(11)
tr12 <- ape::rphylo(12, 0.5, 0.1); b12 <- sort(as.numeric(ape::branching.times(tr12)), decreasing = TRUE)
cr8 <- c(0.5, 0, 0, 0, 0.1, 0, 0, 0)
lb <- c(0.01, 0, 0, 0, 0.001, 0, 0, 0); ub <- c(5, 0, 0, 0, 5, 0, 0, 0)

em1 <- function(nt, brts, copy = FALSE) {
  r <- emphasis:::em_cpp(brts = brts, init_pars = cr8, sample_size = 50L, maxN = 5000L,
                         max_missing = 1000L, max_lambda = 1e6, lower_bound = lb, upper_bound = ub,
                         xtol_rel = 1e-3, num_threads = nt, copy_trees = copy, model = c(0L,0L,0L), link = 0L,
                         rho = 1, rconditional = NULL)
  c(trees = if (copy) length(r$trees) else r$trees, nlogf = length(r$logf), nw = length(r$weights),
    fhat = r$fhat, lam = r$estimates[1], mu = r$estimates[5])
}
e1 <- t(replicate(15, em1(1L, b12))); e8 <- t(replicate(30, em1(8L, b12))); e8c <- t(replicate(10, em1(8L, b12, copy = TRUE)))
cat(sprintf("[em_cpp 12tip] 1thr: trees %d..%d, fhat mean %.3f sd %.3f, lambda %.3f, mu %.3f\n",
            min(e1[,"trees"]), max(e1[,"trees"]), mean(e1[,"fhat"]), sd(e1[,"fhat"]), mean(e1[,"lam"]), mean(e1[,"mu"])))
for (nm in c("e8", "e8c")) {
  e <- get(nm); ov <- e[, "trees"] != 50
  cat(sprintf("[em_cpp 12tip] 8thr (%s): overflow %d/%d, trees %d..%d, trees==nlogf==nweights in all runs: %s; fhat(overflow)-mean fhat(1thr): %s vs log(trees/50): %s; lambda %.3f, mu %.3f\n",
              nm, sum(ov), nrow(e), min(e[,"trees"]), max(e[,"trees"]),
              all(e[,"trees"] == e[,"nlogf"] & e[,"trees"] == e[,"nw"]),
              if (any(ov)) sprintf("%.2f..%.2f", min(e[ov,"fhat"]) - mean(e1[,"fhat"]), max(e[ov,"fhat"]) - mean(e1[,"fhat"])) else "NA",
              if (any(ov)) sprintf("%.2f..%.2f", min(log(e[ov,"trees"]/50)), max(log(e[ov,"trees"]/50))) else "NA",
              mean(e[,"lam"]), mean(e[,"mu"])))
}

fit_one <- function(sampling, nt) {
  t0 <- proc.time()[["elapsed"]]
  f <- estimate_rates(tr12, method = "mcem", model = "cr", init_pars = c(0.5, 0.1),
                      control = list(sampling = sampling, num_threads = nt, num_trees = 50L, maxN = 5000L,
                                     max_iter = 4L, patience = 1L, tol = 1e-6,
                                     lower_bound = c(0.01, 0.001), upper_bound = c(5, 5), max_time = 120))
  c(loglik = unname(f$loglik), AIC = unname(f$AIC %||% NA), tab_num_trees = max(f$details$mcem$num_trees),
    lam = unname(f$pars[1]), sec = proc.time()[["elapsed"]] - t0)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
p_bdi1 <- t(replicate(3, fit_one("bdi", 1L)));  p_bdi8 <- t(replicate(3, fit_one("bdi", 8L)))
p_df1  <- t(replicate(4, fit_one("dynamic_fresh", 1L))); p_df8 <- t(replicate(6, fit_one("dynamic_fresh", 8L)))
show <- function(lbl, p) cat(sprintf("[estimate_rates 12tip] %-18s loglik %s | table num_trees %s | sec %s\n", lbl,
                                     paste(round(p[,"loglik"],2), collapse=","), paste(p[,"tab_num_trees"], collapse=","),
                                     paste(round(p[,"sec"],2), collapse=",")))
show("bdi 1thr", p_bdi1); show("bdi 8thr", p_bdi8); show("dynamic_fresh 1thr", p_df1); show("dynamic_fresh 8thr", p_df8)
