## H8 — D-model + linear link: is the compensator in Model::loglik the
## integral of the model's total intensity Σ_s λ_s(t) + μ_s(t), or the
## approximation n_i · (λ_focal + μ_focal) evaluated at the end node?
##
## Method
##   1. R reference for the complete-data log f with the SAME lineage bookkeeping
##      the C++ exp-link branch uses (2 crown lineages born at 0; every
##      non-extinction, non-closing node adds a lineage born at its brts; an
##      extinction node removes the lineage born at its tip_start), with N and
##      M = pd/n frozen at the end node exactly as the C++ does. The only thing
##      the reference changes is the compensator: it integrates
##      Σ_s rate_s(t) over the segment numerically (integrate()).
##   2. Anchor: under link = exponential the C++ branch is exact per-lineage,
##      so reference == C++ must hold to integrate() precision.
##   3. Under link = linear compare C++ loglik against the reference on
##      (a) a hand-built tree with two lineages of different D,
##      (b) real augmented trees, at the generating theta and at other thetas.
##   4. Sanity: beta_D = gamma_D = 0 -> exact agreement; model_bin[3] = 0 (cr/dd)
##      never reaches the D branch.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
set.seed(8)

T_TIP <- 10e10; T_EXT <- 0; T_UNS <- 5e10
is_ext <- function(t) t == T_EXT
is_tip <- function(t) t == T_TIP

## ---- replicate detail::calculate_pendant_pd (model_helpers.hpp:216-229) ----
pendant_pd <- function(tm, tr) {
  s <- tr$brts <= tm & !is_ext(tr$t_ext) & tr$t_ext > tm
  sum(tm - tr$tip_start[s])
}

## ---- per-lineage rates in the code's parameterisation ----------------------
rate_fun <- function(link) {
  if (link == 0) function(b0, eta) pmax(0, b0 + eta)
  else if (link == 1) function(b0, eta) exp(b0 + eta)
  else function(b0, eta) b0 * exp(-0.5 * (eta - 1)^2)
}

## Node-level rates exactly as speciation_rate_ep / extinction_rate_ep
## (model.hpp:143-212): E = brts - tip_start (extinction), brts - focal_tip_start
## (parent_id >= 0), else M.
node_rates <- function(p, node, link) {
  f <- rate_fun(link)
  M <- if (node$n > 0) node$pd / node$n else 0
  E <- if (is_ext(node$t_ext)) node$brts - node$tip_start else
       if (node$parent_id >= 0) node$brts - node$focal_tip_start else M
  D <- E - M
  c(lambda = f(p[1], p[2] * node$n + p[3] * M + p[4] * D),
    mu     = f(p[5], p[6] * node$n + p[7] * M + p[8] * D))
}

## ∫_{t1}^{t2} max(0, a + b t) dt, vectorised over a (b scalar), exact.
int_relu <- function(a, b, t1, t2) {
  sapply(a, function(ai) {
    if (abs(b) < 1e-14) return(max(0, ai) * (t2 - t1))
    tk <- -ai / b                      # kink
    F  <- function(t) ai * t + 0.5 * b * t^2
    seg <- function(u1, u2) if (u2 <= u1) 0 else {
      mid <- 0.5 * (u1 + u2); if (ai + b * mid > 0) F(u2) - F(u1) else 0 }
    if (tk <= t1 || tk >= t2) seg(t1, t2) else seg(t1, tk) + seg(tk, t2)
  })
}

## Complete-data log f.  compensator = "code" reproduces model.hpp:420-426
## (n * (lambda_focal + mu_focal) * dt); "exact" integrates Σ_s rate_s(t).
logf_R <- function(p, tr, link, compensator = c("code", "exact")) {
  compensator <- match.arg(compensator)
  f <- rate_fun(link)
  nn <- nrow(tr)
  ev <- 0; inte <- 0; prev <- 0
  alive_ts <- c(0, 0)                      # two crown lineages, born at 0
  for (i in seq_len(nn)) {
    node <- tr[i, ]
    dt <- node$brts - prev
    r  <- node_rates(p, node, link)
    if (dt > 0) {
      if (compensator == "code") {
        inte <- inte + dt * node$n * (r[["lambda"]] + r[["mu"]])
      } else {
        M <- if (node$n > 0) node$pd / node$n else 0
        ts <- alive_ts
        if (link == 0) {
          ## per lineage the rate is max(0, a + b t): integrate analytically
          a_l <- p[1] + p[2] * node$n + p[3] * M - p[4] * (ts + M); b_l <- p[4]
          a_m <- p[5] + p[6] * node$n + p[7] * M - p[8] * (ts + M); b_m <- p[8]
          inte <- inte + sum(int_relu(a_l, b_l, prev, node$brts)) +
                         sum(int_relu(a_m, b_m, prev, node$brts))
        } else {
          tot <- function(t) sapply(t, function(tt) {
            Dv <- (tt - ts) - M
            sum(f(p[1], p[2] * node$n + p[3] * M + p[4] * Dv)) +
            sum(f(p[5], p[6] * node$n + p[7] * M + p[8] * Dv))
          })
          inte <- inte + integrate(tot, prev, node$brts, rel.tol = 1e-10,
                                   subdivisions = 2000L)$value
        }
      }
    }
    if (is_ext(node$t_ext)) {
      ev <- ev + log(max(r[["mu"]], 1e-300))
      k <- which(abs(alive_ts - node$tip_start) < 1e-12)[1]
      if (is.na(k)) stop("extinction of unknown lineage at node ", i)
      alive_ts <- alive_ts[-k]
    } else {
      if (i != nn) {
        ev <- ev + log(r[["lambda"]])
        alive_ts <- c(alive_ts, node$brts)   # exp branch: sum += exp(-b*node.brts)
      }
    }
    prev <- node$brts
  }
  c(loglik = ev - inte, events = ev, inte = inte)
}

cpp_logf <- function(p, tr, link, model = c(0L, 0L, 1L))
  emphasis:::eval_logf(p, list(tr), model = as.integer(model), link = as.integer(link), rho = 1)$logf

## =============================================================================
## (a) Hand-built tree with two lineages of different D
## =============================================================================
## Observed: crown at 0 (forward), split at 5, closing at 10 (3 tips).
## Augmented: a missing lineage born at 1 from parent id 0 (crown-derived
## lineage born... parent bookkeeping only affects E of its own node), dies at 4.
cat("=== (a) hand-built tree: crown lineages (ts=0) + missing lineage (ts=1) ===\n")
hand <- data.frame(
  brts      = c(1, 4, 5, 10),
  n         = c(2, 3, 2, 3),
  t_ext     = c(4, T_EXT, T_TIP, T_TIP),
  pd        = 0, tip_start = c(1, 1, 0, 0), focal_tip_start = c(0, 1, 0, 0),
  id        = c(4L, 4L, 0L, 1L), parent_id = c(0L, 0L, -1L, -1L))
hand$pd <- sapply(hand$brts, pendant_pd, tr = hand)
print(hand)

## parameters: beta_0 = 0.4, beta_D = 0.1, gamma_0 = 0.1, gamma_D = 0.05
p_a <- c(0.4, 0, 0, 0.1, 0.1, 0, 0, 0.05)
cat("\nSegment [1,4): n = 3 alive lineages with ts = (0, 0, 1); M_seg = pd/n at end node =",
    hand$pd[2] / hand$n[2], "\n")
M2 <- hand$pd[2] / hand$n[2]
D_focal <- (4 - 1) - M2
cat("  end node = extinction of the ts=1 lineage: D_focal = (4-1) - M =", D_focal, "\n")
cat("  code compensator on this segment = 3 * dt * [max(0,.4+.1*D_f) + max(0,.1+.05*D_f)] =",
    3 * 3 * (max(0, .4 + .1 * D_focal) + max(0, .1 + .05 * D_focal)), "\n")
tot_a <- function(t) sapply(t, function(tt) {
  D <- (tt - c(0, 0, 1)) - M2
  sum(pmax(0, .4 + .1 * D)) + sum(pmax(0, .1 + .05 * D)) })
cat("  exact  ∫_1^4 Σ_s (λ_s + μ_s) dt   =", integrate(tot_a, 1, 4)$value, "\n")
cat("  (Σ_s D_s(t) = 3(t - 1/3) - 3M: not zero because M is frozen at the end node\n",
    "   and, even unfrozen, the focal D of one lineage != the mean over lineages)\n")

for (lk in c(0, 1)) {
  cc <- cpp_logf(p_a, hand, lk)
  rc <- logf_R(p_a, hand, lk, "code")
  re <- logf_R(p_a, hand, lk, "exact")
  cat(sprintf("\n link=%d  C++ loglik = %.8f | R(code compensator) = %.8f | R(exact Σ_s) = %.8f | C++ - exact = %+.6f\n",
              lk, cc, rc[["loglik"]], re[["loglik"]], cc - re[["loglik"]]))
  cat(sprintf("          compensators: code = %.8f   exact = %.8f\n", rc[["inte"]], re[["inte"]]))
}

## =============================================================================
## (b) Real augmented trees under D + linear, and the exp-link anchor
## =============================================================================
cat("\n=== (b) augmented trees (thinning) from brts = c(10,7,4,2) under D+linear ===\n")
brts <- c(10, 7, 4, 2)
p_gen <- c(0.3, 0, 0, 0.05, 0.1, 0, 0, 0.02)
aug <- emphasis:::augment_trees(brts, p_gen, sample_size = 40, maxN = 4000,
                                max_missing = 200, max_lambda = 1e6, num_threads = 1,
                                model = as.integer(c(0, 0, 1)), link = 0L, rho = 1)
trees <- aug$trees
nmiss <- sapply(trees, function(t) sum(is_ext(t$t_ext)))
cat("trees:", length(trees), " missing lineages per tree: ",
    paste(range(nmiss), collapse = "-"), "\n")

## anchor: exp link, exact branch in C++ vs R reference
d_exp <- sapply(trees, function(t) cpp_logf(p_gen, t, 1) - logf_R(p_gen, t, 1, "exact")[["loglik"]])
cat(sprintf("exp-link anchor: max |C++ - R(exact)| over %d trees = %.2e\n", length(trees), max(abs(d_exp))))
d_expc <- sapply(trees, function(t) cpp_logf(p_gen, t, 1) - logf_R(p_gen, t, 1, "code")[["loglik"]])
cat(sprintf("exp-link: C++ - R(code compensator): max |.| = %.3f (the C++ exp branch is NOT n*rate_focal)\n", max(abs(d_expc))))

## linear link: C++ reproduces the "code" compensator exactly; differs from exact
d_linc <- sapply(trees, function(t) cpp_logf(p_gen, t, 0) - logf_R(p_gen, t, 0, "code")[["loglik"]])
d_line <- sapply(trees, function(t) cpp_logf(p_gen, t, 0) - logf_R(p_gen, t, 0, "exact")[["loglik"]])
cat(sprintf("linear link: max |C++ - R(code)| = %.2e   (R reproduces model.hpp:424 exactly)\n", max(abs(d_linc))))
cat(sprintf("linear link: C++ - R(exact Σ_s) at generating theta: mean %+.4f, sd %.4f, range [%+.4f, %+.4f]\n",
            mean(d_line), sd(d_line), min(d_line), max(d_line)))
cat(sprintf("             correlation with #missing lineages: %.2f\n", cor(d_line, nmiss)))

## theta-dependence of the discrepancy (what would move the M-step)
cat("\n theta-dependence: Delta(theta) = C++ - exact, same 40 trees, beta_D varied\n")
for (bD in c(-0.1, 0, 0.05, 0.1, 0.2)) {
  p <- p_gen; p[4] <- bD
  cpp <- sapply(trees, function(t) cpp_logf(p, t, 0))
  ref <- sapply(trees, function(t) logf_R(p, t, 0, "exact")[["loglik"]])
  ok <- is.finite(cpp) & is.finite(ref); d <- (cpp - ref)[ok]
  cat(sprintf("   beta_D = %+.2f : mean Delta = %+.4f  sd = %.4f  (finite in both: %d/%d; C++ non-finite: %d, ref non-finite: %d)\n",
              bD, mean(d), sd(d), sum(ok), length(ok), sum(!is.finite(cpp)), sum(!is.finite(ref))))
}

## weighted M-step objective along beta_D: code f vs exact f, same trees, same weights
lw <- aug$logf - aug$logg
w  <- exp(lw - max(lw))
grid <- seq(-0.15, 0.30, by = 0.01)
Q_code  <- sapply(grid, function(bD) { p <- p_gen; p[4] <- bD
  v <- sapply(trees, function(t) cpp_logf(p, t, 0)); if (all(is.finite(v))) sum(w * v) else -Inf })
Q_exact <- sapply(grid, function(bD) { p <- p_gen; p[4] <- bD
  v <- sapply(trees, function(t) logf_R(p, t, 0, "exact")[["loglik"]]); if (all(is.finite(v))) sum(w * v) else -Inf })
cat(sprintf(" grid points with all-finite Q: code %d/%d, exact %d/%d\n",
            sum(is.finite(Q_code)), length(grid), sum(is.finite(Q_exact)), length(grid)))
cat(sprintf("\n argmax_beta_D of weighted Q on this sample: code f -> %.2f, exact f -> %.2f (generating 0.05)\n",
            grid[which.max(Q_code)], grid[which.max(Q_exact)]))

## =============================================================================
## (c) sanity: no D coefficients -> exact agreement; cr/dd never enter the branch
## =============================================================================
cat("\n=== (c) sanity ===\n")
## gamma_D alone (beta_D = 0) still enters: show it
pg <- c(0.3, 0, 0, 0, 0.1, 0, 0, 0.02)
dg <- sapply(trees, function(t) cpp_logf(pg, t, 0) - logf_R(pg, t, 0, "exact")[["loglik"]])
cat(sprintf("beta_D = 0, gamma_D = 0.02, linear: max |C++ - exact| = %.3f (mu side alone)\n", max(abs(dg))))
p0 <- c(0.3, 0, 0, 0, 0.1, 0, 0, 0)
d0 <- sapply(trees, function(t) cpp_logf(p0, t, 0) - logf_R(p0, t, 0, "exact")[["loglik"]])
cat(sprintf("beta_D = gamma_D = 0, D-model, linear: max |C++ - exact| = %.2e\n", max(abs(d0))))
pdd <- c(0.5, -0.02, 0, 0, 0.1, 0, 0, 0)
ddd <- sapply(trees, function(t) cpp_logf(pdd, t, 0, model = c(1, 0, 0)) - logf_R(pdd, t, 0, "exact")[["loglik"]])
cat(sprintf("dd model (model_bin = 1,0,0), linear: max |C++ - exact| = %.2e (n*lambda == Σ_s lambda_s when rates are lineage-free)\n", max(abs(ddd))))

## =============================================================================
## (d) gaussian link (named in the hypothesis as "non-exact"): the C++ branch
##     integrates per lineage via erf but with a different alive set
##     (model.hpp:405-407: nodes with brts <= prev, t_ext >= brts; crown
##     lineages absent). Compare with the exp-branch bookkeeping reference.
## =============================================================================
cat("\n=== (d) gaussian link, D model ===\n")
pgau <- c(0.5, 0, 0, 0.1, 0.1, 0, 0, 0.05)
dga <- sapply(trees, function(t) cpp_logf(pgau, t, 2) - logf_R(pgau, t, 2, "exact")[["loglik"]])
cat(sprintf("gaussian: C++ - R(exact, exp-branch alive set): mean %+.4f sd %.4f range [%+.4f, %+.4f]\n",
            mean(dga), sd(dga), min(dga), max(dga)))
