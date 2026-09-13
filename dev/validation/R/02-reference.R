# ---------------------------------------------------------------------------
# 02-reference.R — exact MLEs, observed information, and the fixed-parameter
# theta grids.  One record per tree; computed once, cached, never recomputed.
#
#   Rscript 02-reference.R --tier smoke|main|ext [--workers 10] [--lib PATH]
#
# Reads data/trees-<tier>.rds, writes data/reference-<tier>.rds.
#
# CR references (all with btorph = 1, soc = 2, the convention 00-selfcheck.R
# measured emphasis to share with constant 0):
#   theta_MLE0 : argmax of the unconditioned crown likelihood.  Computed twice,
#                by optim(L-BFGS-B) on the closed form from three starts and by
#                DDD::bd_ML(cond = 0); disagreement > 1e-4 per parameter flags
#                the tree and both are re-run from the other's optimum.
#   theta_MLE1 : bd_ML(cond = 1)  == ape::birthdeath — the survival-conditioned
#                MLE the pipeline's GAM stage approximates.
#   theta_MLE2 : bd_ML(cond = 2)  — conditioned on n, the conditional that
#                matches sim.bd.taxa, used only in the recovery arm.
#   SE         : sqrt(diag(solve(-H))) from numDeriv on the closed form; at
#                mu_hat = 0 the profile SE of lambda, and z_mu is undefined.
#
# DD reference: DDD::dd_ML(ddmodel = 1, cond = 0, btorph = 1, soc = 2,
#   res = max(300, 10(n+1))) from two starts (generating pars and 1.5x), higher
#   loglik kept.  conv != 0 or K_hat > 1e4 sets flag "K_unidentified": the tree
#   is counted, reported, and analysed as cr, never scored as estimator error.
#
# The fixed-theta grid (arm C9/D7) is stored with the tree so the fit jobs and
# the analysis use the same points.
# ---------------------------------------------------------------------------

R_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else
    "/Users/pancho/Code/emphasis/dev/validation/R"
})
source(file.path(R_DIR, "00-common.R"))
COMMON <- file.path(R_DIR, "00-common.R")

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, default = NULL) {
  i <- which(args == k); if (length(i)) args[i[1] + 1L] else default
}
TIER    <- getarg("--tier", "smoke")
WORKERS <- as.integer(getarg("--workers", "10"))
if (!is.null(getarg("--lib"))) options(emphasis.lib = getarg("--lib"))

val_load_refs()
suppressPackageStartupMessages({ library(numDeriv); library(future.apply) })

tr_file <- file.path(VAL_DIR$data, sprintf("trees-%s.rds", TIER))
stopifnot(file.exists(tr_file))
TR <- readRDS(tr_file)
cat(sprintf("[02-reference] tier = %s, %d trees, %d workers\n",
            TIER, length(TR$trees), WORKERS))
t_start <- Sys.time()

# --- CR --------------------------------------------------------------------
ref_cr <- function(tt) {
  b <- tt$brts
  t0 <- Sys.time()
  m <- val_cr_mle(b, lam_gen = tt$lambda_gen, mu_gen = tt$mu_gen)
  num1 <- function(x, d = NA_real_) {
    v <- suppressWarnings(as.numeric(x)[1]); if (length(v) && !is.null(v)) v else d
  }
  bd_pars <- function(f) c(num1(f$lambda0), num1(f$mu0))
  bd_ok <- function(f) {
    if (is.null(f)) return(FALSE)
    p <- bd_pars(f)
    all(is.finite(p)) && p[1] > 0 && isTRUE(num1(f$conv) == 0) &&
      is.finite(num1(f$loglik))
  }
  bd0 <- .q(DDD::bd_ML(brts = b, cond = 0, btorph = 1, soc = 2, tdmodel = 0,
                       initparsopt = c(tt$lambda_gen, max(tt$mu_gen, 0.05)),
                       idparsopt = 1:2, verbose = FALSE))
  mle_ok <- function(x) !is.null(x) && is.finite(x$loglik)
  disagree <- if (bd_ok(bd0) && mle_ok(m))
    max(abs(bd_pars(bd0) - m$pars)) else NA_real_
  if (is.na(disagree) || disagree > 1e-4) {
    # re-run each optimiser from the other's optimum, keep the higher loglik
    st <- if (bd_ok(bd0)) bd_pars(bd0) else c(tt$lambda_gen, 0.5 * tt$lambda_gen)
    m2 <- val_cr_mle(b, starts = rbind(st, c(0.8 * tt$lambda_gen, 0.7 * tt$lambda_gen),
                                       c(1.5 * tt$lambda_gen, 0.1 * tt$lambda_gen)))
    bd0b <- .q(DDD::bd_ML(brts = b, cond = 0, btorph = 1, soc = 2, tdmodel = 0,
                          initparsopt = if (mle_ok(m)) pmax(m$pars, c(1e-3, 1e-6))
                                        else c(tt$lambda_gen, 0.5 * tt$lambda_gen),
                          idparsopt = 1:2, verbose = FALSE))
    cand <- list(m, m2)
    lls  <- vapply(cand, function(x) if (!mle_ok(x)) -Inf else x$loglik, 1)
    if (any(is.finite(lls))) m <- cand[[which.max(lls)]]
    if (bd_ok(bd0b) && (!bd_ok(bd0) || num1(bd0b$loglik) > num1(bd0$loglik))) bd0 <- bd0b
    disagree <- if (bd_ok(bd0) && mle_ok(m))
      max(abs(bd_pars(bd0) - m$pars)) else NA_real_
  }
  if (!mle_ok(m))
    return(list(tree_id = tt$tree_id, kind = "cr", n = tt$n, cell = tt$cell,
                flag = "mle_failed", ddd_status = ddd_status_na <- "unavailable",
                elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  # bd_ML returns pars = -1 and conv = -1 when its optimiser fails.  The closed
  # form is the reference surface (asserted equal to bd_loglik in tier 0), so
  # optim's optimum stands; the DDD cross-check is recorded as unavailable.
  ddd_status <- if (bd_ok(bd0)) "ok" else "bd_ML_failed"
  # interior first-order check on the closed form
  gr <- tryCatch(numDeriv::grad(function(p) ll_cr_nee(p[1L], p[2L], b), m$pars),
                 error = function(e) c(NA_real_, NA_real_))
  # scaled so the criterion is "moving one SE changes the log-likelihood by
  # less than 1e-3 nats", which is comparable across n and across time units
  sc <- ifelse(is.finite(m$se) & m$se > 0, m$se, 1)
  grad_ok <- if (isTRUE(m$boundary_mu0))
    is.finite(gr[1]) && abs(gr[1] * sc[1]) < 1e-3 &&
      (!is.finite(gr[2]) || gr[2] * sc[2] < 1e-3)
  else all(is.finite(gr)) && max(abs(gr * sc)) < 1e-3
  bd1 <- .q(DDD::bd_ML(brts = b, cond = 1, btorph = 1, soc = 2, tdmodel = 0,
                       initparsopt = c(tt$lambda_gen, max(tt$mu_gen, 0.05)),
                       idparsopt = 1:2, verbose = FALSE))
  bd2 <- tryCatch(.q(DDD::bd_ML(brts = b, cond = 2, btorph = 1, soc = 2,
                                tdmodel = 0,
                                initparsopt = c(tt$lambda_gen, max(tt$mu_gen, 0.05)),
                                idparsopt = 1:2, verbose = FALSE)),
                  error = function(e) list(lambda0 = NA, mu0 = NA, loglik = NA,
                                           conv = -99))
  box <- val_box_cr(tt$lambda_gen)
  list(
    tree_id = tt$tree_id, kind = "cr", n = tt$n, cell = tt$cell,
    mle = m$pars, mle_loglik = m$loglik, se = m$se,
    boundary_mu0 = m$boundary_mu0, mle_r_negative = m$mle_r_negative,
    ddd_mle0 = bd_pars(bd0), ddd_loglik0 = num1(bd0$loglik),
    ddd_conv0 = num1(bd0$conv), mle_disagreement = disagree,
    ddd_status = ddd_status, grad = gr, grad_ok = grad_ok,
    mle_cond1 = bd_pars(bd1), loglik_cond1 = num1(bd1$loglik),
    cond1_ok = bd_ok(bd1),
    mle_cond2 = bd_pars(bd2), loglik_cond2 = num1(bd2$loglik),
    cond2_ok = bd_ok(bd2),
    ll_at_gen = ll_cr_nee(tt$lambda_gen, tt$mu_gen, b),
    theta_grid = val_theta_grid_cr(m),
    mle_in_box = all(m$pars >= box$lower) && all(m$pars <= box$upper),
    box = box,
    elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# --- DD --------------------------------------------------------------------
ref_dd <- function(tt) {
  b <- tt$brts; n <- tt$n; lx <- val_lx(n)
  t0 <- Sys.time()
  fit_one <- function(start) tryCatch(
    .q(DDD::dd_ML(brts = b, initparsopt = start, idparsopt = 1:3,
                  ddmodel = 1, cond = 0, btorph = 1, soc = 2,
                  res = lx, methode = "analytical",
                  optimmethod = "subplex", verbose = FALSE)),
    error = function(e) NULL)
  f1 <- fit_one(c(tt$lambda0_gen, tt$mu0_gen, tt$K_gen))
  f2 <- fit_one(1.5 * c(tt$lambda0_gen, tt$mu0_gen, tt$K_gen))
  cands <- Filter(function(x) !is.null(x) && is.finite(x$loglik), list(f1, f2))
  if (!length(cands))
    return(list(tree_id = tt$tree_id, kind = "dd", n = n, cell = tt$cell,
                flag = "dd_ML_failed", lx = lx,
                elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  ll <- vapply(cands, function(x) x$loglik, 1)
  best <- cands[[which.max(ll)]]
  two_start_gap <- if (length(cands) == 2L) abs(diff(ll)) else NA_real_
  K_hat <- best$K
  flag <- if (best$conv != 0 || !is.finite(K_hat) || K_hat > 1e4)
    "K_unidentified" else "ok"
  pars_e <- as.numeric(val_dd_to_emphasis(best$lambda, best$mu, K_hat))
  # SE on the emphasis scale from the closed-form-free numeric Hessian of
  # dd_loglik in (beta_0, beta_N, gamma_0); K's flat direction makes this a
  # scale for the grid, not an inferential SE.
  se <- rep(NA_real_, 4)
  if (flag == "ok") {
    fll <- function(p) {
      m <- val_emphasis_to_dd(c(p, 0))
      if (m$status != "ok") return(NA_real_)
      ll_dd_ddd(m$lambda0, m$mu0, m$K, b, lx = lx)
    }
    # central differences (19 dd_loglik evaluations) rather than numDeriv's
    # Richardson extrapolation (70+): one evaluation here is an ODE solve.
    H <- tryCatch(val_hess_fd(fll, pars_e[1:3]), error = function(e) NULL)
    if (!is.null(H) && all(is.finite(H))) {
      V <- tryCatch(solve(-H), error = function(e) NULL)
      if (!is.null(V) && all(diag(V) > 0)) se[1:3] <- sqrt(diag(V))
    }
  }
  if (!is.finite(se[1])) se[1] <- 0.2 * pars_e[1]
  if (!is.finite(se[2])) se[2] <- 0.3 * abs(pars_e[2])
  if (!is.finite(se[3])) se[3] <- max(0.2 * pars_e[3], 0.02)
  box <- val_box_dd(tt$lambda0_gen, tt$mu0_gen, tt$K_gen)
  list(tree_id = tt$tree_id, kind = "dd", n = n, cell = tt$cell,
       regime = tt$regime, lx = lx, flag = flag,
       mle_ddd = c(lambda0 = best$lambda, mu0 = best$mu, K = K_hat),
       mle = pars_e, mle_loglik = best$loglik, se = se,
       conv = best$conv, two_start_gap = two_start_gap,
       ll_at_gen = ll_dd_ddd(tt$lambda0_gen, tt$mu0_gen, tt$K_gen, b, lx = lx),
       theta_grid = val_theta_grid_dd(pars_e, se),
       mle_in_box = all(pars_e >= box$lower - 1e-12) &&
                    all(pars_e <= box$upper + 1e-12),
       box = box,
       elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

future::plan(future::multisession, workers = WORKERS)
on.exit(future::plan(future::sequential), add = TRUE)

refs <- future.apply::future_lapply(
  TR$trees,
  function(tt, .common, .ref_cr, .ref_dd) {
    source(.common, local = FALSE)          # hidden helpers (.q) included
    val_load_refs()
    suppressPackageStartupMessages(library(numDeriv))
    environment(.ref_cr) <- environment(.ref_dd) <- globalenv()
    if (tt$kind == "cr") .ref_cr(tt) else .ref_dd(tt)
  },
  .common = COMMON, .ref_cr = ref_cr, .ref_dd = ref_dd,
  future.seed = TRUE)
names(refs) <- names(TR$trees)

# --- report ----------------------------------------------------------------
cr <- Filter(function(x) x$kind == "cr", refs)
dd <- Filter(function(x) x$kind == "dd", refs)
if (length(cr)) {
  dis <- vapply(cr, function(x) x$mle_disagreement %||% NA_real_, 1)
  bnd <- vapply(cr, function(x) isTRUE(x$boundary_mu0), TRUE)
  inb <- vapply(cr, function(x) isTRUE(x$mle_in_box), TRUE)
  neg <- vapply(cr, function(x) isTRUE(x$mle_r_negative), TRUE)
  gok <- vapply(cr, function(x) isTRUE(x$grad_ok), TRUE)
  bfail <- vapply(cr, function(x) !identical(x$ddd_status, "ok"), TRUE)
  cat(sprintf("  cr: %d trees | max |optim - bd_ML| = %s over %d comparable | bd_ML failed on %d | first-order check ok on %d | mu_hat = 0 on %d | mu>=lambda on %d | MLE in box %d/%d\n",
              length(cr),
              if (all(is.na(dis))) "NA" else format(max(dis, na.rm = TRUE), digits = 3),
              sum(!is.na(dis)), sum(bfail), sum(gok), sum(bnd), sum(neg),
              sum(inb), length(inb)))
  if (any(!is.na(dis) & dis > 1e-4))
    cat("    trees with optim/bd_ML disagreement > 1e-4: ",
        paste(names(cr)[which(!is.na(dis) & dis > 1e-4)], collapse = ", "), "\n")
}
if (length(dd)) {
  fl <- vapply(dd, function(x) x$flag %||% "NA", "")
  el <- vapply(dd, function(x) x$elapsed, 1)
  cat(sprintf("  dd: %d trees | flags: %s | dd_ML %.0f-%.0fs\n",
              length(dd), paste(names(table(fl)), table(fl), sep = "=",
                                collapse = " "), min(el), max(el)))
}

out <- file.path(VAL_DIR$data, sprintf("reference-%s.rds", TIER))
saveRDS(list(refs = refs, tier = TIER, created = Sys.time(),
             ddd = as.character(utils::packageVersion("DDD"))), out)
cat(sprintf("  wrote %s (%.1fs)\n", out,
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
