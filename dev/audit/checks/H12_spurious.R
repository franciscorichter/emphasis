## H12 follow-up: in the mu>lam region, does .bdi_find_event_time_cr always
## error, or can floating point (t_cur + 1e-15 == t_cur + <1e-15) make
## f(lo) = -U so that uniroot returns a spurious "root" at t_cur?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages({library(emphasis); library(ape); library(DDD)})
ns <- asNamespace("emphasis")
find_ev <- ns$.bdi_find_event_time_cr; int_cr <- ns$.bdi_integral_cr
set.seed(5)
tc <- sort(runif(400, 0, 4.9))
res <- vapply(tc, function(t) {
  r <- tryCatch(suppressWarnings(find_ev(t, 0.7, 1L, 2L, 0.3, 0.5, 5, 4.95)), error = function(e) NA_real_)
  if (is.na(r)) NA_real_ else r - t
}, 1)
cat("n grid points:", length(tc), " errors:", sum(is.na(res)), " spurious roots:", sum(!is.na(res)),
    " (root - t_cur) range:", range(res, na.rm = TRUE), "\n")
for (lo in c(0, 0.5, 1, 2, 4)) {
  sel <- tc >= lo & tc < 2 * max(lo, 0.25)
  cat(sprintf("  t_cur in [%.2g,%.2g): %d pts, %d spurious\n", lo, 2*max(lo,0.25), sum(sel), sum(!is.na(res[sel]))))
}
cat("f(lo) at t_cur=3: int_cr(3, 3+1e-15, 1, 2, 0.3, 0.5, 5) =", int_cr(3, 3 + 1e-15, 1L, 2L, 0.3, 0.5, 5),
    "  (3+1e-15)-3 =", (3 + 1e-15) - 3, "\n")

## consequence: sampler on a low-rate 5-age tree (few immigrations) at mu>lam
mk <- function(n, seed, age = 5, b = 0.5, d = 0.2) {
  set.seed(seed); tr <- ape::rphylo(n, b, d)
  tr$edge.length <- tr$edge.length / max(ape::branching.times(tr)) * age; tr
}
for (cfg in list(list(tr = mk(8, 3, age = 5), pars = c(0.05, 0.10), lab = "8 tips age 5 (0.05,0.10)"),
                 list(tr = mk(15, 121, age = 5), pars = c(0.1, 0.2), lab = "15 tips age 5 (0.1,0.2)"),
                 list(tr = mk(15, 121, age = 5), pars = c(0.3, 0.5), lab = "15 tips age 5 (0.3,0.5)"))) {
  brts <- sort(ape::branching.times(cfg$tr), decreasing = TRUE)
  bl <- DDD::bd_loglik(pars1 = c(cfg$pars, 0, 0), pars2 = c(0, 0, 0, 0, 2), brts = brts, missnumspec = 0)
  out <- replicate(20, tryCatch({
    a <- ns$.augment_tree_bdi(cfg$tr, pars = cfg$pars, sample_size = 5L)
    nm <- vapply(a$trees, function(d) sum(d$t_ext == 0), 1L)
    sprintf("ok n=%d nmiss=%s sd(lw)=%.2e fhat-bd=%.4f", length(a$trees), paste(nm, collapse="/"), sd(a$weights), a$fhat - bl)
  }, error = function(e) "ERROR"))
  cat("--", cfg$lab, " (true offset -log((n-1)!) =", -lgamma(Ntip(cfg$tr)), ") --\n"); print(table(out))
}
