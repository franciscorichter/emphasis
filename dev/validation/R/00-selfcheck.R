# ---------------------------------------------------------------------------
# 00-selfcheck.R — tier 0.  Reference identities and build sentinels.
#
# Runs first, serially.  Aborts on a reference failure; writes results/GATE.ok
# on success.  Nothing downstream may run without it: every performance measure
# in the study is a difference to a DDD-computed exact log-likelihood, so an
# encoding error on the reference side would be reported as estimator error.
#
# Assertions that ABORT
#   (a) pars2 layout guard: the layout used by dev/validation/ref_check.R
#       (verbose = 2, soc = 0) differs from the documented one, and only the
#       documented one makes dd_loglik(K = 1e6) agree with bd_loglik.
#   (b) bd_loglik(btorph = 1) equals the Nee closed form to 1e-6.
#   (c) BDI fhat (cr, linear, rho = 1) equals bd_loglik(btorph = 1) with
#       sd(lw) = 0, at mu < lambda — i.e. the constant c_n is 0, measured.
#   (d) bd_ML(cond = 1, soc = 2) equals ape::birthdeath to 1e-4.
#   (e) dd_loglik is insensitive to lx and to the integrator at the optimum.
#   (g) CR scale invariance of log-likelihood differences.
#
# Sentinels that are RECORDED, not fatal (the study is designed to measure them)
#   (c2) BDI at mu > lambda (wave 1.3)
#   (f)  thinning fhat - bd_loglik = 0 +/- 0.05 (wave 1.8)
#   (h)  dd BDI returns acc / n_rejected / n_nonfinite (wave 1.4)
#
# Usage:  Rscript 00-selfcheck.R [--lib /path/to/rlib]
# ---------------------------------------------------------------------------

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f[1])) else
    "/Users/pancho/Code/emphasis/dev/validation/R"
  source(file.path(d, "00-common.R"), local = FALSE)
})

args <- commandArgs(trailingOnly = TRUE)
if (any(args == "--lib")) options(emphasis.lib = args[which(args == "--lib") + 1L])

val_load_refs()
val_load_emphasis()
suppressPackageStartupMessages(library(numDeriv))

FP <- val_build_fingerprint()
cat("emphasis ", FP$version, "  lib ", FP$lib, "\n", sep = "")
cat("DDD ", FP$ddd_version, "  R ", FP$r_version, "\n\n", sep = "")

FAIL <- character(0)
SENT <- list()
chk <- function(id, ok, msg) {
  cat(sprintf("  [%s] %-4s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", id, msg))
  if (!isTRUE(ok)) FAIL <<- c(FAIL, id)
  invisible(ok)
}
sent <- function(id, ok, msg) {
  SENT[[id]] <<- isTRUE(ok)
  cat(sprintf("  [%s] %-4s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", id, msg))
  invisible(ok)
}

# --- three reference trees -------------------------------------------------
mk <- function(seed, n, lam = 0.6, mu = 0.2) {
  set.seed(seed)
  tr <- TreeSim::sim.bd.taxa(n = n, numbsim = 1, lambda = lam, mu = mu,
                             complete = FALSE)[[1]]
  sort(as.numeric(ape::branching.times(tr)), decreasing = TRUE)
}
TREES <- list(t10 = mk(7, 10), t30 = mk(7, 30), t60 = mk(7, 60),
              t30b = mk(11, 30))
THETA <- rbind(c(0.6, 0.2), c(0.8, 0.5), c(0.5, 0.05), c(1.0, 0.9),
               c(0.3, 0.25), c(0.4, 0.6), c(0.5, 0.5))

cat("== (a) pars2 layout guard ==\n")
# dd_loglik(ddmodel = 1) is -Inf whenever mu0 >= lambda0 (beta_N would be
# non-negative), so the dd half of the guard runs only where lambda > mu.
THETA_DD <- THETA[THETA[, 1] > THETA[, 2], , drop = FALSE]
a_diff <- a_match <- numeric(0)
for (nm in names(TREES)) {
  b <- TREES[[nm]]; n <- length(b) + 1L; lx <- val_lx(n)
  for (i in seq_len(nrow(THETA_DD))) {
    l <- THETA_DD[i, 1]; m <- THETA_DD[i, 2]
    right_bd <- ll_cr_ddd(l, m, b)
    right_dd <- ll_dd_ddd(l, m, 1e6, b, lx = lx)
    wrong_bd <- .q(DDD::bd_loglik(pars1 = c(l, m, 0, 0),
                                  pars2 = c(0, 0, 1, 2, 0), brts = b,
                                  missnumspec = 0))
    wrong_dd <- .q(DDD::dd_loglik(pars1 = c(l, m, 1e6),
                                  pars2 = c(lx, 1, 0, 1, 2, 0), brts = b,
                                  missnumspec = 0))
    a_match <- c(a_match, right_dd - right_bd)
    a_diff  <- c(a_diff,  wrong_dd - wrong_bd)
  }
}
# a1 discriminates the documented pars2 layout from the wrong one.  The wrong
# layout (a2) is off by 1.26 nats at minimum, so the threshold only has to sit
# well below that while clearing DDD's own integrator noise.  That noise is
# platform-dependent: max|dd(K=1e6) - bd| is 3.4e-04 under clang/macOS but
# 3.9e-03 under gcc/Linux, and it does NOT shrink with lx -- doubling lx on the
# n = 60 tree moves it -3.90e-03 -> +5.70e-04 -> +5.08e-03, so dd_loglik is
# limited by ODE tolerance, not by truncation.  A2_MIN_GAP / 60 keeps a 60x
# margin to the signal a1 exists to catch.  The measured value is printed so a
# genuine layout regression cannot hide inside the tolerance.
A1_TOL <- 2e-2
chk("a1", max(abs(a_match)) < A1_TOL,
    sprintf("documented layout: max|dd(K=1e6) - bd| = %.2e (< %.0e; a2 gap >= 1.26)",
            max(abs(a_match)), A1_TOL))
chk("a2", min(abs(a_diff)) > 1,
    sprintf("ref_check.R layout (verbose=2, soc=0): gap %.2f to %.2f, theta-dependent",
            min(a_diff), max(a_diff)))

cat("== (b) bd_loglik == Nee closed form ==\n")
b_gap <- numeric(0)
for (nm in names(TREES)) for (i in seq_len(nrow(THETA)))
  b_gap <- c(b_gap, ll_cr_ddd(THETA[i, 1], THETA[i, 2], TREES[[nm]]) -
                    ll_cr_nee(THETA[i, 1], THETA[i, 2], TREES[[nm]]))
chk("b", max(abs(b_gap)) < 1e-6,
    sprintf("max|bd_loglik(btorph=1) - ll_cr_nee| = %.2e", max(abs(b_gap))))

cat("== (c) BDI fhat == bd_loglik; the constant c_n ==\n")
cn <- data.frame()
for (nm in names(TREES)) {
  b <- TREES[[nm]]; n <- length(b) + 1L
  for (i in seq_len(nrow(THETA))) {
    l <- THETA[i, 1]; m <- THETA[i, 2]
    # On a build where the mu >= lambda branch is not fixed, these draws can
    # take minutes (H12's error path is retried), so each is time-limited.
    r <- tryCatch({
      setTimeLimit(elapsed = 120, transient = TRUE)
      on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
      emphasis:::.augment_tree_bdi(b, c(l, m), model_bin = c(0L, 0L, 0L),
                                   sample_size = 50L)
    }, error = function(e) e)
    if (inherits(r, "error")) {
      cn <- rbind(cn, data.frame(tree = nm, n = n, lambda = l, mu = m,
                                 c_n = NA_real_, sd_lw = NA_real_,
                                 err = conditionMessage(r), stringsAsFactors = FALSE))
    } else {
      cn <- rbind(cn, data.frame(tree = nm, n = n, lambda = l, mu = m,
                                 c_n = r$fhat - ll_cr_nee(l, m, b),
                                 sd_lw = stats::sd(r$weights), err = NA_character_,
                                 stringsAsFactors = FALSE))
    }
  }
}
sub <- cn[cn$mu < cn$lambda, ]
chk("c1", all(is.finite(sub$c_n)) && max(abs(sub$c_n)) < 1e-6,
    sprintf("mu < lambda: max|c_n| = %.2e (c_n must be 0; -log((n-1)!) is the btorph=0 constant)",
            suppressWarnings(max(abs(sub$c_n), na.rm = TRUE))))
chk("c1b", all(is.finite(sub$sd_lw)) && max(sub$sd_lw, na.rm = TRUE) < 1e-8,
    sprintf("mu < lambda: max sd(lw) = %.2e (zero-variance weights)",
            suppressWarnings(max(sub$sd_lw, na.rm = TRUE))))
sup <- cn[cn$mu >= cn$lambda, ]
sent("c2", nrow(sup) > 0 && all(is.finite(sup$c_n)) && max(abs(sup$c_n)) < 1e-6,
     sprintf("mu >= lambda (wave 1.3): %d/%d draws finite, max|c_n| = %s",
             sum(is.finite(sup$c_n)), nrow(sup),
             format(suppressWarnings(max(abs(sup$c_n), na.rm = TRUE)), digits = 3)))

cat("== (d) bd_ML(cond = 1) == ape::birthdeath ==\n")
set.seed(7)
phy <- TreeSim::sim.bd.taxa(n = 30, numbsim = 1, lambda = 0.6, mu = 0.2,
                            complete = FALSE)[[1]]
bd <- ape::birthdeath(phy)
rr <- bd$para["b-d"]; aa <- bd$para["d/b"]
ape_mle <- c(unname(rr / (1 - aa)), unname(aa * rr / (1 - aa)))
f1 <- .q(DDD::bd_ML(brts = TREES$t30, cond = 1, btorph = 1, soc = 2,
                    initparsopt = c(0.6, 0.2), idparsopt = 1:2, tdmodel = 0,
                    verbose = FALSE))
chk("d", max(abs(c(f1$lambda0, f1$mu0) - ape_mle)) < 1e-4,
    sprintf("ape (%.4f, %.4f) vs bd_ML cond=1 (%.4f, %.4f)",
            ape_mle[1], ape_mle[2], f1$lambda0, f1$mu0))

cat("== (e) dd_loglik lx / integrator sensitivity ==\n")
e_gap <- numeric(0)
for (nm in names(TREES)) {
  b <- TREES[[nm]]; n <- length(b) + 1L
  for (i in c(1L, 4L)) {
    l <- THETA[i, 1]; m <- THETA[i, 2]; K <- 3 * n
    v1 <- ll_dd_ddd(l, m, K, b, lx = val_lx(n))
    v2 <- ll_dd_ddd(l, m, K, b, lx = 3 * val_lx(n))
    v3 <- ll_dd_ddd(l, m, K, b, lx = val_lx(n),
                    methode = "odeint::runge_kutta_fehlberg78")
    e_gap <- c(e_gap, v2 - v1, v3 - v1)
  }
}
chk("e", max(abs(e_gap)) < 1e-2,
    sprintf("max|dd_loglik difference over lx x3 and integrator| = %.2e", max(abs(e_gap))))

cat("== (f) thinning fhat sentinel (wave 1.8) ==\n")
b <- TREES$t30
thin_fhat <- function(brts, l, m, N = 4000L) {
  p8 <- c(l, 0, 0, 0, m, 0, 0, 0)
  a <- emphasis:::augment_trees(brts, p8, N, 20L * N, 1e4, 1e6, 1L,
                                model = c(0L, 0L, 0L), link = 0L, rho = 1)
  ev <- emphasis:::eval_logf(p8, a$trees, model = c(0L, 0L, 0L), link = 0L, rho = 1)
  lw <- ev$logf - a$logg
  list(fhat = emphasis:::.is_fhat(ev$logf, a$logg,
                                  n_zero_weight = a$rejected_zero_weights),
       ess = sum(exp(lw - max(lw)))^2 / sum(exp(2 * (lw - max(lw)))),
       rzw = a$rejected_zero_weights)
}
# The thinning sampler is clock-seeded, so this sentinel is stochastic: three
# replicates, and the criterion is the mean against its own replicate SE.
rep_gap <- function(brts, l, m, R = 3L) {
  v <- vapply(seq_len(R), function(i) {
    r <- thin_fhat(brts, l, m); c(r$fhat - ll_cr_nee(l, m, brts), r$ess) }, c(0, 0))
  list(mean = mean(v[1, ]), se = stats::sd(v[1, ]) / sqrt(R), ess = mean(v[2, ]))
}
tf  <- rep_gap(TREES$t30,  0.6, 0.2)    # crown 11.5, low ESS
tf2 <- rep_gap(TREES$t30b, 0.6, 0.2)    # crown 7.0, moderate ESS
gap_f  <- tf$mean
gap_f2 <- tf2$mean
sent("f", abs(gap_f2) < max(0.02 + 2 * tf2$se, 0.05) && tf2$ess >= 50,
     sprintf("thinning fhat - bd_loglik: %+.4f +/- %.4f (ESS %.0f) moderate age; %+.4f +/- %.4f (ESS %.0f) long age",
             gap_f2, tf2$se, tf2$ess, gap_f, tf$se, tf$ess))

cat("== (g) CR scale invariance ==\n")
s <- 3
b3 <- b * s
d0 <- ll_cr_nee(0.6, 0.2, b) - ll_cr_nee(0.9, 0.5, b)
d1 <- ll_cr_nee(0.6 / s, 0.2 / s, b3) - ll_cr_nee(0.9 / s, 0.5 / s, b3)
chk("g", abs(d0 - d1) < 1e-8,
    sprintf("loglik difference invariant under brts*%d, rates/%d: %.2e", s, s, abs(d0 - d1)))

cat("== (h) dd BDI wave-1.4 fields ==\n")
set.seed(11)
ddsim <- .q(DDD::dd_sim(pars = c(0.8, 0.1, 40), age = 10, ddmodel = 1))
ddb <- sort(as.numeric(ape::branching.times(ddsim$tes)), decreasing = TRUE)
ddp <- val_dd_to_emphasis(0.8, 0.1, 40)
rdd <- tryCatch(emphasis:::.augment_tree_bdi(ddb, as.numeric(ddp),
                                             model_bin = c(1L, 0L, 0L),
                                             sample_size = 100L),
                error = function(e) e)
has_fields <- !inherits(rdd, "error") &&
  all(c("acc", "n_rejected", "n_nonfinite") %in% names(rdd))
sent("h", has_fields,
     sprintf("dd BDI returns acc/n_rejected/n_nonfinite: %s (fields: %s)",
             has_fields,
             if (inherits(rdd, "error")) conditionMessage(rdd)
             else paste(names(rdd), collapse = ",")))
if (!inherits(rdd, "error")) {
  ll_dd_ref <- ll_dd_ddd(0.8, 0.1, 40, ddb)
  cat(sprintf("       dd BDI fhat - dd_loglik = %+.4f%s\n",
              rdd$fhat - ll_dd_ref,
              if (has_fields && is.finite(rdd$acc))
                sprintf("   (acc = %.3f, -log acc = %+.4f)", rdd$acc, -log(rdd$acc)) else ""))
}

# --- summary ---------------------------------------------------------------
cat("\n---------------------------------------------------------------\n")
cat("c_n table (emphasis BDI fhat - DDD bd_loglik(btorph = 1)):\n")
print(cn, row.names = FALSE, digits = 4)
cat(sprintf("\nlog((n-1)!) at n = 30 is %.4f — that is the btorph = 0 constant,\n",
            lgamma(30)))
cat("not the emphasis constant.  Reported logliks compare with btorph = 1 directly.\n")

out <- list(fingerprint = FP, c_n = cn, fails = FAIL, sentinels = SENT,
            thinning_gap = c(long_age = gap_f, moderate_age = gap_f2),
            thinning_ess = c(long_age = tf$ess, moderate_age = tf2$ess),
            thinning_se = c(long_age = tf$se, moderate_age = tf2$se),
            time = Sys.time())
saveRDS(out, file.path(VAL_DIR$results, "selfcheck.rds"))

if (length(FAIL)) {
  cat("\nGATE FAILED:", paste(FAIL, collapse = ", "), "\n")
  quit(status = 1L)
}
writeLines(c(paste("build:", val_fingerprint_key(FP)),
             paste("time:", format(Sys.time())),
             paste("sentinels:", paste(sprintf("%s=%s", names(SENT),
                                               unlist(SENT)), collapse = " "))),
           file.path(VAL_DIR$results, "GATE.ok"))
cat("\nGATE OK -> ", file.path(VAL_DIR$results, "GATE.ok"), "\n", sep = "")
