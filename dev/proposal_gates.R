# dev/proposal_gates.R — the measurements behind tests/testthat/test-proposal-density.R,
# at sample sizes the test suite cannot afford.
#
#   Rscript dev/proposal_gates.R lifetimes    # KS of the drawn lifetimes vs the charged density
#   Rscript dev/proposal_gates.R ess          # effective sample size of a D-model fit
#   Rscript dev/proposal_gates.R unbiased     # E_q[f/q] against a brute-force marginal likelihood
#   Rscript dev/proposal_gates.R all
#
# The cr/dd bit-for-bit gate lives in dev/crdd_invariance.R and is not repeated here.

suppressMessages(devtools::load_all(".", quiet = TRUE))

T_TIP <- 10e10; T_EXT <- 0; T_UNS <- 5e10
is_mis <- function(te) !(te == T_EXT | te == T_TIP | te == T_UNS)
lnk <- function(link, eta) if (link == 1L) exp(eta) else pmax(0, eta)

# The extinction rate sampling_prob charges a lifetime with: the segment's mu,
# which is D-free and reads M where the segment begins.  This is the R side of
# Model::proposal_rates.
mu_charged <- function(df, i, p, link) {
  prev <- if (i == 1L) 0 else df$brts[i - 1L]
  pd <- df$pd[i] + df$n[i] * (prev - df$brts[i])
  M  <- if (df$n[i] > 0) pd / df$n[i] else 0
  max(lnk(link, p[5] + p[6] * df$n[i] + p[7] * M), 1e-10)
}

ess <- function(lw) { w <- exp(lw - max(lw)); sum(w)^2 / sum(w^2) }

.deep_tree <- function() {
  set.seed(11)
  phy <- ape::rphylo(10L, 0.8, 0.25)
  phy$edge.length <- phy$edge.length * 4      # crown age ~10, so mu * (T - t) is
  phy                                          # large enough for the PIT to bite
}

.settings <- list(
  list(name = "d  linear", mb = c(0L, 0L, 1L), link = 0L,
       p = c( 0.12,  0.000, 0, 0.05,  0.45, 0, 0, 0.40)),
  list(name = "nd linear", mb = c(1L, 0L, 1L), link = 0L,
       p = c( 0.14, -0.004, 0, 0.05,  0.45, 0, 0, 0.40)),
  list(name = "d  exp",    mb = c(0L, 0L, 1L), link = 1L,
       p = c(-2.10,  0.000, 0, 0.05, -0.80, 0, 0, 0.60)),
  list(name = "nd exp",    mb = c(1L, 0L, 1L), link = 1L,
       p = c(-2.05, -0.004, 0, 0.05, -0.80, 0, 0, 0.60))
)

# --------------------------------------------------------------------------- #
#  Gate 2 — the lifetimes are drawn from the density they are charged with      #
# --------------------------------------------------------------------------- #

run_lifetimes <- function(N = 2000L) {
  phy  <- .deep_tree()
  brts <- emphasis:::.extract_brts(phy); pts <- emphasis:::.pts(brts); TT <- brts[1]
  cat(sprintf("crown age %.3f, %d tips\n", TT, ape::Ntip(phy)))
  for (s in .settings) {
    a <- augment_trees(as.numeric(brts), s$p, N, 400000L, 800L, 1e6, 1L,
                       model = s$mb, link = s$link, rho = 1, parent_tip_start = pts)
    mu <- l <- r <- numeric(0)
    for (df in a$trees) for (i in which(is_mis(df$t_ext))) {
      mu <- c(mu, mu_charged(df, i, s$p, s$link))
      l  <- c(l, df$t_ext[i] - df$brts[i]); r <- c(r, TT - df$brts[i])
    }
    u <- (1 - exp(-mu * l)) / (1 - exp(-mu * r))     # PIT of the truncated exponential
    cat(sprintf("  %-10s lifetimes %6d   KS p = %-10.4g  median mu*(T-t) = %.2f\n",
                s$name, length(u), suppressWarnings(ks.test(u, "punif"))$p.value,
                median(mu * r)))
  }
}

# --------------------------------------------------------------------------- #
#  Gate 5 — effective sample size                                              #
# --------------------------------------------------------------------------- #

# A second, milder theta on the same tree: the same four models with the D
# coefficients an order of magnitude smaller, where the IS weights are not
# heavy-tailed and the ESS is a stable number rather than a noisy one.
.mild_settings <- list(
  list(name = "d  linear", mb = c(0L, 0L, 1L), link = 0L,
       p = c( 0.30,  0.000, 0, 0.06,  0.16, 0, 0, 0.07)),
  list(name = "nd linear", mb = c(1L, 0L, 1L), link = 0L,
       p = c( 0.34, -0.006, 0, 0.06,  0.16, 0, 0, 0.07)),
  list(name = "d  exp",    mb = c(0L, 0L, 1L), link = 1L,
       p = c(-1.30,  0.000, 0, 0.06, -1.90, 0, 0, 0.10)),
  list(name = "nd exp",    mb = c(1L, 0L, 1L), link = 1L,
       p = c(-1.25, -0.004, 0, 0.06, -1.90, 0, 0, 0.10))
)

run_ess <- function(N = 2000L, reps = 8L) {
  phy  <- .deep_tree()
  brts <- emphasis:::.extract_brts(phy); pts <- emphasis:::.pts(brts)
  for (nm in c("strong D", "mild D")) {
    cat(sprintf("  -- %s --\n", nm))
    for (s in if (nm == "strong D") .settings else .mild_settings) {
      v <- numeric(reps); miss <- numeric(reps); viol <- 0
      for (r in seq_len(reps)) {
        a <- augment_trees(as.numeric(brts), s$p, N, 400000L, 800L, 1e6, 1L,
                           model = s$mb, link = s$link, rho = 1, parent_tip_start = pts)
        e <- eval_logf(s$p, a$trees, model = s$mb, link = s$link, rho = 1)
        v[r] <- ess(e$logf - a$logg)
        miss[r] <- mean(vapply(a$trees, function(d) sum(is_mis(d$t_ext)), 0))
        viol <- viol + a$envelope_violations
      }
      cat(sprintf("  %-10s ESS of %d over %d reps: min %6.1f  median %6.1f  max %6.1f   mean #missing %5.2f   envelope violations %d\n",
                  s$name, N, reps, min(v), stats::median(v), max(v), mean(miss), viol))
    }
  }
}

# --------------------------------------------------------------------------- #
#  Gate 3 — unbiasedness against a brute-force marginal likelihood             #
# --------------------------------------------------------------------------- #
#
# On a 2-tip tree an augmentation with one missing lineage is a pair (t, d):
# born at t from one of the two crown lineages, dead at d.  The observed tree
# has no internal node, so the sampler's candidate-parent list is empty at
# every t and the parent is a crown lineage — which is also every labelled
# attachment the -log(2 tips + Ne) term counts (2 per observed lineage, and on
# this tree both observed lineages are crown lineages carrying tip_start 0, so
# f is the same for all four).  The stratum is therefore exactly
#
#   L1 = 4 * int_0^T int_t^T f(y, z(t, d)) dd dt
#
# and no parent convention is being approximated.  On a tree with an observed
# split that is no longer true: the sampler draws uniformly over the nodes
# alive, which never include the crown lineages, so some labelled attachments
# have proposal probability zero (H45, open).  That is why the brute force is
# done here and not on a larger tree.

gauss_legendre <- function(n) {                    # Golub-Welsch, mapped to [0,1]
  k <- 1:(n - 1); b <- k / sqrt(4 * k^2 - 1)
  J <- diag(0, n); J[cbind(k, k + 1)] <- b; J[cbind(k + 1, k)] <- b
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = (e$values[o] + 1) / 2, w = e$vectors[1, o]^2)
}

aug1 <- function(t, d, TT) data.frame(
  brts = c(t, d, TT), n = c(2, 3, 2), t_ext = c(d, 0, T_TIP),
  pd = c(2 * t, 3 * d - 2 * t, 2 * TT - t), tip_start = c(t, t, TT),
  focal_tip_start = c(0, t, -1), id = c(1L, 1L, 0L), parent_id = c(-1L, -1L, -1L))

aug0 <- function(TT) data.frame(
  brts = TT, n = 2, t_ext = T_TIP, pd = 2 * TT,
  tip_start = TT, focal_tip_start = -1, id = 0L, parent_id = -1L)

L1 <- function(TT, p, mb, lk, n = 90L) {
  gq <- gauss_legendre(n); tt <- TT * gq$x
  trees <- vector("list", n * n); wgt <- numeric(n * n); k <- 0L
  for (ix in seq_len(n)) {
    t <- tt[ix]; dd <- t + (TT - t) * gq$x
    for (iy in seq_len(n)) {
      k <- k + 1L
      trees[[k]] <- aug1(t, dd[iy], TT)
      wgt[k] <- gq$w[ix] * gq$w[iy] * TT * (TT - t)
    }
  }
  4 * sum(wgt * exp(eval_logf(p, trees, model = mb, link = lk, rho = 1)$logf))
}

gate3 <- function(label, TT, p, mb, lk, nbatch, nper, nq = 90L) {
  l0  <- exp(eval_logf(p, list(aug0(TT)), model = mb, link = lk, rho = 1)$logf)
  l1  <- L1(TT, p, mb, lk, nq)
  l1b <- L1(TT, p, mb, lk, as.integer(nq / 2))
  tot <- ss <- tot1 <- ss1 <- totle <- ssle <- w2 <- 0; n2 <- 0L; N <- 0L; nz <- 0L
  for (b in seq_len(nbatch)) {
    a <- augment_trees(TT, p, as.integer(nper), 100L * as.integer(nper), 200L,
                       1e6, 1L, model = mb, link = lk, rho = 1,
                       parent_tip_start = c(-1))
    stopifnot(a$rejected_overruns == 0, a$rejected_lambda == 0,
              a$rejected_nonfinite == 0)
    # a draw whose f is zero carries weight zero but is still a draw from q:
    # it belongs in the denominator, which is why it is counted and not dropped.
    nz <- nz + a$rejected_zero_weights
    w  <- exp(eval_logf(p, a$trees, model = mb, link = lk, rho = 1)$logf - a$logg)
    nm <- vapply(a$trees, function(d) sum(is_mis(d$t_ext)), 0)
    wle <- ifelse(nm <= 1, w, 0)
    N <- N + length(w); tot <- tot + sum(w); ss <- ss + sum(w^2)
    tot1 <- tot1 + sum(w[nm == 1]); ss1 <- ss1 + sum(w[nm == 1]^2)
    totle <- totle + sum(wle); ssle <- ssle + sum(wle^2)
    n2 <- n2 + sum(nm >= 2); w2 <- w2 + sum(w[nm >= 2])
  }
  M <- N + nz                       # completed augmentations, zero-weight included
  E <- tot / M; se <- sqrt(ss / M - E^2) / sqrt(M)
  E1 <- tot1 / M; se1 <- sqrt(ss1 / M - E1^2) / sqrt(M)
  Ele <- totle / M; sele <- sqrt(ssle / M - Ele^2) / sqrt(M)
  cat(sprintf("\n== %s   (draws = %d, of which %d had f = 0)\n", label, M, nz))
  cat(sprintf("   quadrature L1 %.10f (n = %d) vs %.10f (n = %d), rel diff %.1e\n",
              l1, nq, l1b, as.integer(nq / 2), abs(l1 - l1b) / l1))
  cat(sprintf("   L0 %.8f   L1 %.8f   L0 + L1 %.8f\n", l0, l1, l0 + l1))
  cat(sprintf("   one missing lineage : IS %.8f +- %.8f   z = %+.2f\n", E1, se1, (E1 - l1) / se1))
  cat(sprintf("   at most one         : IS %.8f +- %.8f   z = %+.2f\n", Ele, sele, (Ele - l0 - l1) / sele))
  cat(sprintf("   whole sample        : IS %.8f +- %.8f   z = %+.2f\n", E, se, (E - l0 - l1) / se))
  cat(sprintf("   truncation: %d draws with >= 2 missing (%.4f%%), %.4f%% of the weight; MC se %.4f%%\n",
              n2, 100 * n2 / M, 100 * w2 / tot, 100 * se / E))
}

run_unbiased <- function(nbatch = 15L, nper = 20000L) {
  set.seed(7)
  gate3("d linear, E[#missing] 0.53", 3, c(0.30, 0, 0, 0.15, 0.25, 0, 0, 0.10),
        c(0L, 0L, 1L), 0L, nbatch, nper)
  gate3("d linear, E[#missing] 0.036", 3, c(0.02, 0, 0, 0.15, 0.25, 0, 0, 0.10),
        c(0L, 0L, 1L), 0L, nbatch, nper)
  gate3("d exp, E[#missing] 0.036", 3,
        c(log(0.02), 0, 0, 0.15, log(0.25), 0, 0, 0.10), c(0L, 0L, 1L), 1L, nbatch, nper)
}

# --------------------------------------------------------------------------- #

args <- commandArgs(trailingOnly = TRUE)
what <- if (length(args)) args[[1L]] else "all"
if (what %in% c("lifetimes", "all")) { cat("\n-- gate 2: lifetimes --\n");  run_lifetimes() }
if (what %in% c("ess", "all"))       { cat("\n-- gate 5: ESS --\n");        run_ess() }
if (what %in% c("unbiased", "all"))  { cat("\n-- gate 3: unbiasedness --\n"); run_unbiased() }
