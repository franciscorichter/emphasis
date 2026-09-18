# mfc_mesh.R -- does re-evaluating the frozen Gillespie rates on a mesh move
# the estimate, and what does it cost?
#
# The approximate Gillespie branch of the conditional sampler holds its rates
# at the value they had when the step started.  Under diversity dependence the
# rate depends on the lineage count, which the step itself changes, so a long
# step is charged at a rate that no longer holds and the proposal density is
# wrong by whatever that drift is.  The importance weights are supposed to
# absorb that, and do -- but only up to Monte Carlo error, so a systematic
# drift shows as a shifted estimate.
#
# The reference is the thinning proposal, which is a different proposal for
# the same target: both estimate log p(obs | theta), so they must agree.
#
#   Rscript dev/mfc_mesh.R [--draws 4000] [--reps 8]
suppressMessages(devtools::load_all(".", quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
ga <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1L] }
DRAWS <- as.integer(ga("--draws", 4000)); REPS <- as.integer(ga("--reps", 8))
MESHES <- c(Inf, 2000, 500, 100)

mb <- .resolve_model("dd")
rows <- list()
for (tip in c(30L, 60L)) {
  set.seed(100 + tip); phy <- ape::rcoal(tip)
  brts <- .extract_brts(phy)
  p8 <- .expand_pars(c(0.7, -(0.7 - 0.25) / (2 * tip), 0.25, 0), mb)

  ref <- vapply(seq_len(REPS), function(r) {
    set.seed(8000 + r)
    t <- tryCatch(augment_trees(brts, p8, sample_size = DRAWS, maxN = 50L * DRAWS,
                                max_missing = 5000L, max_lambda = 1e6, num_threads = 1L,
                                model = as.integer(mb[1:4]), link = 0L, rho = 1.0,
                                parent_tip_start = .pts(brts), parent_id = .pid(brts)),
                  error = function(e) NULL)
    if (is.null(t)) return(NA_real_)
    .is_summary(t$logf - t$logg, n_zero_weight = t$rejected_zero_weights)$fhat
  }, 1)
  r_mu <- mean(ref, na.rm = TRUE); r_se <- stats::sd(ref, na.rm = TRUE) / sqrt(sum(is.finite(ref)))

  for (mesh in MESHES) {
    v <- vapply(seq_len(REPS), function(r) {
      set.seed(8000 + r)
      a <- .augment_tree_bdi(brts, p8, model_bin = mb[1:3], sample_size = DRAWS,
                             link = 0L, rho = 1.0,
                             mesh = if (is.finite(mesh)) mesh else NULL)
      c(a$fhat, .ess_from_lw(a$weights), stats::sd(a$weights[is.finite(a$weights)]))
    }, numeric(3))
    mu <- mean(v[1, ], na.rm = TRUE); se <- stats::sd(v[1, ], na.rm = TRUE) / sqrt(REPS)
    d_se <- (mu - r_mu) / sqrt(se^2 + r_se^2)
    rows[[length(rows) + 1L]] <- data.frame(
      tips = tip, mesh = mesh, fhat = mu, fhat_se = se,
      ref = r_mu, ref_se = r_se, gap_se = d_se,
      ess = mean(v[2, ], na.rm = TRUE), sd_logw = mean(v[3, ], na.rm = TRUE))
    cat(sprintf("[mesh] %2d tips  mesh %6s  fhat %9.4f (se %.4f)  ref %9.4f  gap %+5.1f se  ESS %7.1f  sd(log w) %.3f\n",
                tip, if (is.finite(mesh)) format(mesh) else "off", mu, se, r_mu, d_se,
                mean(v[2, ], na.rm = TRUE), mean(v[3, ], na.rm = TRUE)))
    utils::flush.console()
  }
}
write.csv(do.call(rbind, rows), "dev/mfc_mesh.csv", row.names = FALSE)
