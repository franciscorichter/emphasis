# mfc_ess.R -- how much of a draw survives the importance weights, by sampler.
#
# The conditional (BDI) sampler is exact only under constant rates.  On a tree
# whose speciation rate reads a per-lineage covariate (ED) it is a mean-field
# approximation: it proposes from rates that see only the lineage count, and
# the importance weights carry the rest.  The question this answers is whether
# that trade is worth making -- whether a wrong-but-cheap conditional proposal
# leaves more effective draws than a thinning proposal that reads ED exactly.
#
# Three arms, the same observed trees, the same parameters, the same draw count:
#   thin-ED   thinning, proposal rate reads ED           (model slot 4 = 1)
#   thin-MF   thinning, proposal rate ED-blind           (model slot 4 = 2)
#   bdi       conditional, ED-blind, weights compensate
#
#   Rscript dev/mfc_ess.R [--draws 200] [--trees 3] [--seeds 2] [--sizes 25,50,100]
suppressMessages(devtools::load_all(".", quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
ga <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1L] }
DRAWS <- as.integer(ga("--draws", 200)); NTREE <- as.integer(ga("--trees", 3))
NSEED <- as.integer(ga("--seeds", 2))
SIZES <- as.integer(strsplit(ga("--sizes", "25,50,100"), ",")[[1]])
OUT   <- ga("--out", "dev/mfc_ess.csv")

LINK <- 0L; TURNOVER <- 0.33; B0 <- 0.5; CROWN <- 15
# The ED effect is declared as a fraction of the speciation rate at the average
# lineage, not as a bare coefficient: fair-proportion ED carries the units of
# time, so a fixed coefficient means a different effect at every clade size and
# at a long enough crown age drives lambda negative.
ED_EFFECT <- -0.15
G0   <- TURNOVER * B0
# K such that the N-only part of the rate vanishes at K lineages
mb_ned <- .resolve_model("ned")
mb_dd  <- .resolve_model("dd")

# At beta_ED = 0 the ED process is the diversity-dependent one, and simulating
# it as "dd" skips a per-event sweep over the whole forest that costs more than
# the rest of the calibration together.
sim_raw <- function(cp, T_use, model, seed) {
  set.seed(seed)
  tryCatch(simulate_tree(pars = cp, max_t = T_use, model = model, rho = 1,
                         max_lin = 20000L, num_threads = 1L),
           error = function(e) NULL)
}
draw_one <- function(cp, T_use, model, seed) {
  s <- sim_raw(cp, T_use, model, seed)
  if (is.null(s) || !identical(s$status, "done") || is.null(s$tes)) return(NA_integer_)
  length(s$tes$tip.label)
}

# The clade's size is set at a fixed crown age by its capacity, so the
# capacity is bisected to the target rather than the age: calibrating the age
# instead drives it up under a negative ED effect, and a long age at a small
# size is a tree that is almost all missing lineage.
calibrate_K <- function(n, seed, bED, model, lo = 0.5 * n, hi = 12 * n) {
  best <- 1.4 * n
  for (it in 1:16) {
    mid <- (lo + hi) / 2
    bN  <- -(B0 - G0) / mid
    cp  <- if (model == "dd") c(B0, bN, G0, 0) else c(B0, bN, bED, G0, 0, 0)
    m <- stats::median(vapply(seq_len(12), function(i)
                                draw_one(cp, CROWN, model, seed + 31 * it + i), 1L),
                       na.rm = TRUE)
    if (is.na(m)) { lo <- mid; next }
    best <- mid
    if (m < n) lo <- mid else hi <- mid
    if (abs(m - n) <= 0.05 * n) break
  }
  best
}

# Mean fair-proportion ED of the clade with no ED effect: summed over the tips
# it is Faith's PD, so the mean is PD per tip.
ed_bar <- function(cp0, seed) {
  v <- vapply(seq_len(12), function(i) {
    s <- sim_raw(cp0, CROWN, "dd", seed + 517 + i)
    if (is.null(s) || !identical(s$status, "done") || is.null(s$tas)) return(NA_real_)
    sum(s$tas$edge.length) / length(s$tas$tip.label)
  }, 1)
  stats::median(v, na.rm = TRUE)
}

sim_one <- function(n, seed) {
  # the same clade with no ED effect fixes the capacity and the ED scale, then
  # the coefficient is set from that scale and the capacity re-calibrated
  K0    <- calibrate_K(n, seed, 0, "dd")
  edbar <- ed_bar(c(B0, -(B0 - G0) / K0, G0, 0), seed)
  if (!is.finite(edbar) || edbar <= 0) return(NULL)
  BED   <- ED_EFFECT * B0 / edbar
  K     <- calibrate_K(n, seed, BED, "ned")
  bN    <- -(B0 - G0) / K
  cp    <- c(B0, bN, BED, G0, 0, 0)
  pn    <- .expand_pars(cp, mb_ned)
  for (try in seq_len(400)) {
    set.seed(seed * 1000 + try)
    s <- tryCatch(simulate_tree(pars = cp, max_t = CROWN, model = "ned", rho = 1,
                                max_lin = 20000L, num_threads = 1L),
                  error = function(e) NULL)
    if (is.null(s) || !identical(s$status, "done") || is.null(s$tes)) next
    nt <- length(s$tes$tip.label)
    if (nt >= 0.6 * n && nt <= 1.6 * n)
      return(list(phy = s$tes, n = nt, pars = pn, bN = bN, T_use = CROWN, K = K,
                  bED = BED, edbar = edbar))
  }
  NULL
}

ess_thin <- function(brts, pn, slot4) {
  mb <- mb_ned; mb[4L] <- slot4
  a <- tryCatch(augment_trees(brts, pn, sample_size = DRAWS,
                              maxN = 50L * DRAWS, max_missing = 5000L,
                              max_lambda = 1e6, num_threads = 1L,
                              model = as.integer(mb[1:4]), link = LINK, rho = 1.0,
                              parent_tip_start = .pts(brts), parent_id = .pid(brts)),
                error = function(e) { message("  thin: ", conditionMessage(e)); NULL })
  if (is.null(a) || !length(a$logf)) return(c(ess = NA_real_, n = 0))
  c(ess = .ess_from_lw(a$logf - a$logg), n = length(a$logf))
}

# The mean-field conditional proposal: the rate every lineage is given is the
# one the clade's average lineage has, and the weights carry the rest.  `fold`
# is the cruder version of the same idea -- the mean ED effect folded into the
# intercept once, instead of following P-hat/N-hat over time.
ess_bdi <- function(brts, pn, bN, bED, edbar, fold = FALSE, mesh = NULL) {
  a <- tryCatch({
    if (fold) {
      p8 <- .expand_pars(c(B0 + bED * edbar, bN, G0, 0), mb_dd)
      .augment_tree_bdi(brts, p8, model_bin = mb_dd[1:3], sample_size = DRAWS,
                        link = LINK, rho = 1.0, mesh = mesh)
    } else {
      .augment_tree_bdi(brts, pn, model_bin = mb_ned[1:4], sample_size = DRAWS,
                        link = LINK, rho = 1.0, mesh = mesh)
    }
  }, error = function(e) { message("  bdi: ", conditionMessage(e)); NULL })
  if (is.null(a) || !length(a$trees)) return(c(ess = NA_real_, n = 0))
  if (!fold) return(c(ess = .ess_from_lw(a$weights), n = length(a$trees)))
  ev <- tryCatch(eval_logf(pn, a$trees, model = as.integer(mb_ned[1:4]),
                           link = LINK, rho = 1.0), error = function(e) NULL)
  if (is.null(ev)) return(c(ess = NA_real_, n = 0))
  c(ess = .ess_from_lw(ev$logf - a$logg), n = length(a$trees))
}

rows <- list()
for (n in SIZES) for (tr in seq_len(NTREE)) {
  cat(sprintf("[mfc] n=%3d tree %d: simulating\n", n, tr)); utils::flush.console()
  s <- sim_one(n, 7 * n + tr)
  if (is.null(s)) { cat(sprintf("[mfc] n=%d tree %d: no tree in range\n", n, tr)); next }
  brts <- .extract_brts(s$phy)
  for (sd in seq_len(NSEED)) {
    tm <- numeric(4)
    tm[1] <- system.time({ set.seed(4242 + sd); a <- ess_thin(brts, s$pars, 1L) })[["elapsed"]]
    tm[2] <- system.time({ set.seed(4242 + sd); b <- ess_thin(brts, s$pars, 2L) })[["elapsed"]]
    tm[3] <- system.time({ set.seed(4242 + sd); d <- ess_bdi(brts, s$pars, s$bN, s$bED, s$edbar) })[["elapsed"]]
    tm[4] <- system.time({ set.seed(4242 + sd); e <- ess_bdi(brts, s$pars, s$bN, s$bED, s$edbar, fold = TRUE) })[["elapsed"]]
    rows[[length(rows) + 1L]] <- data.frame(
      n_target = n, n_tips = s$n, tree = tr, seed = sd,
      crown = s$T_use, K = s$K, b_ED = s$bED, ed_bar = s$edbar,
      ess_thin_ed = a[["ess"]], ess_thin_mf = b[["ess"]], ess_bdi = d[["ess"]],
      ess_bdi_fold = e[["ess"]], nd_bdi_fold = e[["n"]], sec_bdi_fold = tm[4],
      nd_thin_ed = a[["n"]], nd_thin_mf = b[["n"]], nd_bdi = d[["n"]],
      sec_thin_ed = tm[1], sec_thin_mf = tm[2], sec_bdi = tm[3])
    cat(sprintf("[mfc] n=%3d (%3d tips) tree %d seed %d  thin-ED %7.1f (%4.1fs)  thin-MF %7.1f (%4.1fs)  bdi %7.1f (%4.1fs)  bdi-fold %7.1f (%4.1fs)\n",
                n, s$n, tr, sd, a[["ess"]], tm[1], b[["ess"]], tm[2], d[["ess"]], tm[3], e[["ess"]], tm[4]))
    utils::flush.console()
  }
}
d <- do.call(rbind, rows)
write.csv(d, OUT, row.names = FALSE)
cat("\n== median ESS of", DRAWS, "draws ==\n")
agg <- aggregate(cbind(ess_thin_ed, ess_thin_mf, ess_bdi, ess_bdi_fold) ~ n_target, d, median, na.rm = TRUE)
print(agg, row.names = FALSE)
win <- sum(agg$ess_bdi > agg$ess_thin_ed & agg$ess_bdi > agg$ess_thin_mf)
cat(sprintf("\nbdi beats both thinning arms in %d of %d cells\n", win, nrow(agg)))
