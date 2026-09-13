## H78 replication — independent variation of dev/audit/checks/H78.R
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
options(width = 120, digits = 10)
ev <- function(p8, tr, model, link) emphasis:::eval_logf(p8, list(tr), model = as.integer(model), link = as.integer(link), rho = 1)$logf

## V1. Tree with parent_id >= 0 speciation nodes and one extinction node, so D = E - M != 0.
## Forward time, crown 0 (lineages 0,1 tip_start 0). n = lineages alive during (prev, brts]; pd = P just before event.
tr <- data.frame(
  brts            = c(1.0, 1.5, 2.0, 3.0),
  n               = c(2,   3,   4,   3),
  t_ext           = c(1e11, 1e11, 0.0, 1e11),      # 3rd node is an extinction event
  pd              = c(2.0, 2.5, 3.0, 5.0),
  tip_start       = c(1.0, 1.5, 1.0, 0.0),         # extinction node: lineage born at 1.0
  focal_tip_start = c(0.0, 0.0, 0.0, 0.0),
  id              = c(2L, 3L, 2L, -1L),
  parent_id       = c(0L, 1L, -1L, -1L))
E_of <- function(tr) ifelse(tr$t_ext == 0, tr$brts - tr$tip_start,
                     ifelse(tr$parent_id >= 0, tr$brts - tr$focal_tip_start, tr$pd / tr$n))
hand <- function(tr, p, covfun) {   # mirrors model.hpp piecewise-constant path, link 0, slot3 on
  M <- tr$pd / tr$n; X <- covfun(tr, M)
  lam <- pmax(0, p[1] + p[2]*tr$n + p[3]*M + p[4]*X); mu <- pmax(0, p[5] + p[6]*tr$n + p[7]*M + p[8]*X)
  dt <- diff(c(0, tr$brts)); ext <- tr$t_ext == 0; last <- seq_len(nrow(tr)) == nrow(tr)
  sum(log(lam[!ext & !last])) + sum(log(mu[ext])) - sum(dt * tr$n * (lam + mu))
}
pD <- c(0.8, 0, 0, 0.5, 0.2, 0, 0, 0.1)
v_code <- ev(pD, tr, c(0,0,1), 0)
h_D <- hand(tr, pD, function(tr, M) E_of(tr) - M)
h_E <- hand(tr, pD, function(tr, M) E_of(tr))
h_0 <- hand(tr, pD, function(tr, M) 0 * M)
cat(sprintf("V1 eval_logf(model=c(0,0,1), link=0) = %.8f\n   hand with X = D = E - M: %.8f (|diff| %.1e)\n   hand with X = raw E:     %.8f (|diff| %.1e)\n   hand with X = 0:         %.8f (|diff| %.1e)\n",
            v_code, h_D, abs(h_D - v_code), h_E, abs(h_E - v_code), h_0, abs(h_0 - v_code)))
tr2 <- tr; tr2$focal_tip_start[2] <- 1.0      # change E of node 2 only (pd, n unchanged)
cat(sprintf("   change focal_tip_start of node 2 (E 1.5 -> 0.5): eval_logf = %.8f (moved by %.3e) -> slot 3 reads E via D\n",
            ev(pD, tr2, c(0,0,1), 0), ev(pD, tr2, c(0,0,1), 0) - v_code))
cat(sprintf("   same pars, model=c(0,0,0): %.8f ; model=c(1,1,0): %.8f  -> slot 3 is the only gate (equal to X=0 hand: %s)\n",
            ev(pD, tr, c(0,0,0), 0), ev(pD, tr, c(1,1,0), 0), isTRUE(all.equal(ev(pD, tr, c(0,0,0), 0), h_0))))

## V2. Slot 2 gates nothing under every link
pM <- c(0.8, 0, 0.3, 0, 0.2, 0, 0.1, 0)
for (lk in 0:2) {
  a <- ev(pM, tr, c(0,0,0), lk); b <- ev(pM, tr, c(0,1,0), lk); d <- ev(pM, tr, c(1,1,1), lk); e <- ev(pM, tr, c(0,0,1), lk)
  cat(sprintf("V2 link=%d: c(0,0,0) %.6f | c(0,1,0) %.6f | c(1,1,1) %.6f | c(0,0,1) %.6f  -> slot2 diff %.1e, slot3 diff %.1e\n",
              lk, a, b, d, e, abs(a - b), abs(a - e)))
}

## V3. link = 3 is accepted silently and behaves as linear (documented set {0,1} < accepted set {0,1,2,3,...})
cat(sprintf("V3 eval_logf link=3: %.8f vs link=0: %.8f (diff %.1e); link=7: %.8f\n",
            ev(pD, tr, c(0,0,1), 3), ev(pD, tr, c(0,0,1), 0), abs(ev(pD, tr, c(0,0,1), 3) - ev(pD, tr, c(0,0,1), 0)), ev(pD, tr, c(0,0,1), 7)))

## V4. em_cpp and m_cpp (not exercised by H78.R) with link = 2
brts <- c(4, 3, 2, 1)
p_g <- c(0.8, 0, 0, 0, 0.2, 0, 0, 0)
lb <- c(0.01, 0, 0, 0, 0.001, 0, 0, 0); ub <- c(3, 0, 0, 0, 1, 0, 0, 0)
for (lk in c(2L, 3L)) {
  r <- tryCatch(emphasis:::em_cpp(brts, p_g, 10L, 500L, 50L, 100, lb, ub, 1e-3, 1L, FALSE, c(0L,0L,0L), lk, 1.0, NULL),
                error = function(e) conditionMessage(e))
  if (is.list(r)) cat(sprintf("V4 em_cpp link=%d: estimates = %s, nlopt = %s, fhat finite = %s\n", lk,
                              paste(signif(r$estimates, 4), collapse = ","), r$nlopt, is.finite(r$fhat)))
  else cat(sprintf("V4 em_cpp link=%d: ERROR %s\n", lk, r))
}
aug <- emphasis:::augment_trees(brts = brts, pars = p_g, sample_size = 10, maxN = 500, max_missing = 50,
                                max_lambda = 100, num_threads = 1, model = c(0L,0L,0L), link = 2L, rho = 1)
w <- exp(aug$logf - aug$logg); es <- list(trees = aug$trees, weights = w / sum(w), rejected = 0L, rejected_overruns = 0L,
                                          rejected_lambda = 0L, rejected_zero_weights = 0L, time = 0, fhat = log(mean(w)))
r <- tryCatch(emphasis:::m_cpp(es, p_g, "rpd5c", lb, ub, 1e-3, 1L, c(0L,0L,0L), 2L, 1.0, NULL), error = function(e) conditionMessage(e))
if (is.list(r)) cat(sprintf("V4 m_cpp link=2: estimates = %s, nlopt = %s\n", paste(signif(r$estimates, 4), collapse = ","), r$nlopt)) else cat("V4 m_cpp link=2: ERROR", r, "\n")

## V5. emphasis_cem (R/de.R:548 also says 0/1) with link = 2
r <- tryCatch(emphasis:::emphasis_cem(brts, max_iter = 2L, num_points = 4L, max_missing = 20L, sd_vec = c(0.2, 0, 0, 0, 0.1, 0, 0, 0),
                                     lower_bound = lb, upper_bound = ub, maxN = 50L, sample_size = 3L, link = 2L),
              error = function(e) conditionMessage(e))
if (is.list(r)) cat(sprintf("V5 emphasis_cem link=2: converged = %s, best pars = %s\n", r$converged, paste(signif(r$obtained_estim, 4), collapse = ","))) else cat("V5 emphasis_cem link=2: ERROR", r, "\n")

## V6. R-level range checks: .resolve_link on numerics
cat(sprintf("V6 .resolve_link(2)=%d  .resolve_link(3)=%d  .resolve_link('gaussian')=%d ; .resolve_model(c(0,1,0)) accepted: %s\n",
            emphasis:::.resolve_link(2), emphasis:::.resolve_link(3), emphasis:::.resolve_link("gaussian"),
            !inherits(try(emphasis:::.resolve_model(c(0,1,0)), silent = TRUE), "try-error")))
