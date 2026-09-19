# mfc_fit.R -- does the better proposal give a better estimate?
#
# Effective sample size says how much of a draw survives the weights.  What a
# study needs is the estimate.  Same tree, same start, same budget, the two
# proposals side by side: the error in beta_ED against the value the tree was
# simulated from.
#
#   Rscript dev/mfc_fit.R [--trees 6] [--tips 100] [--draws 200] [--iter 40]
suppressMessages(devtools::load_all(".", quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
ga <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1L] }
NTREE <- as.integer(ga("--trees", 6)); TIPS <- as.integer(ga("--tips", 100))
DRAWS <- as.integer(ga("--draws", 200)); ITER <- as.integer(ga("--iter", 40))

B0 <- 0.5; TURNOVER <- 0.33; G0 <- TURNOVER * B0; CROWN <- 15; ED_EFFECT <- -0.15
mb <- .resolve_model("ned")
sim <- function(cp, seed) { set.seed(seed)
  tryCatch(simulate_tree(pars = cp, max_t = CROWN, model = "ned", rho = 1,
                         max_lin = 20000L, num_threads = 1L), error = function(e) NULL) }
ntips <- function(s) if (is.null(s) || !identical(s$status, "done")) NA_integer_ else length(s$tes$tip.label)

# capacity calibrated at the fixed crown age, ED coefficient from the clade's own ED scale
calK <- function(n, seed, bED, model) { lo <- 0.5 * n; hi <- 12 * n; best <- 1.4 * n
  for (it in 1:16) { mid <- (lo + hi) / 2; bN <- -(B0 - G0) / mid
    cp <- if (model == "dd") c(B0, bN, G0, 0) else c(B0, bN, bED, G0, 0, 0)
    m <- stats::median(vapply(seq_len(12), function(i) {
      set.seed(seed + 31 * it + i)
      s <- tryCatch(simulate_tree(pars = cp, max_t = CROWN, model = model, rho = 1,
                                  max_lin = 20000L, num_threads = 1L), error = function(e) NULL)
      if (is.null(s) || !identical(s$status, "done")) NA_integer_ else length(s$tes$tip.label) }, 1L),
      na.rm = TRUE)
    if (is.na(m)) { lo <- mid; next }
    best <- mid; if (m < n) lo <- mid else hi <- mid
    if (abs(m - n) <= 0.05 * n) break }
  best }

rows <- list()
for (tr in seq_len(NTREE)) {
  seed <- 5000 + 97 * tr
  K0 <- calK(TIPS, seed, 0, "dd"); bN0 <- -(B0 - G0) / K0
  ed <- stats::median(vapply(1:12, function(i) {
    s <- sim(c(B0, bN0, 0, G0, 0, 0), seed + 500 + i)
    if (is.na(ntips(s))) NA_real_ else sum(s$tas$edge.length) / length(s$tas$tip.label) }, 1), na.rm = TRUE)
  bED <- ED_EFFECT * B0 / ed
  K <- calK(TIPS, seed, bED, "ned"); bN <- -(B0 - G0) / K
  cp <- c(B0, bN, bED, G0, 0, 0)
  phy <- NULL
  for (i in 1:400) { s <- sim(cp, seed * 7 + i); n <- ntips(s)
    if (!is.na(n) && n >= 0.6 * TIPS && n <= 1.6 * TIPS) { phy <- s$tes; break } }
  if (is.null(phy)) { cat(sprintf("[fit] tree %d: no tree in range\n", tr)); next }

  sED <- 4 * max(abs(bED), 0.05 * B0); sN <- 5 * abs(bN)
  lb <- c(B0 / 5, -sN, -sED, 0, -sN, -sED); ub <- c(B0 * 3, 0, sED, max(3 * G0, B0), sN, sED)
  init <- (lb + ub) / 2; init[1] <- 1.5 * B0
  for (smp in c("bdi", "dynamic_fresh")) {
    ctrl <- list(lower_bound = lb, upper_bound = ub, sampling = smp, sample_size = DRAWS,
                 maxN = 8000L, max_iter = ITER, max_missing = 1e4, num_threads = 1L,
                 rho = 1, xtol = 1e-3, verbose = FALSE, max_time = 900)
    t0 <- proc.time()[3]
    f <- suppressMessages(tryCatch(estimate_rates(phy, method = "mcem", model = "ned",
             init_pars = init, control = ctrl, link = "linear"), error = function(e) e))
    el <- proc.time()[3] - t0
    ok <- !inherits(f, "error")
    est <- if (ok) as.numeric(f$pars) else rep(NA_real_, 6)
    rows[[length(rows) + 1L]] <- data.frame(tree = tr, tips = ape::Ntip(phy), sampler = smp,
      b0 = est[1], bN = est[2], bED = est[3], g0 = est[4],
      b0_true = B0, bN_true = bN, bED_true = bED, g0_true = G0,
      ESS = if (ok) (f$details$final_IS$ESS %||% NA_real_)[1] else NA_real_,
      sec = el, ok = ok)
    cat(sprintf("[fit] tree %d (%3d tips) %-14s bED %8.4f (true %8.4f)  b0 %6.3f (%.3f)  ESS %6.1f  %5.0fs\n",
                tr, ape::Ntip(phy), smp, est[3], bED, est[1], B0,
                rows[[length(rows)]]$ESS, el))
    utils::flush.console()
  }
}
d <- do.call(rbind, rows); write.csv(d, "dev/mfc_fit.csv", row.names = FALSE)
cat("\n== error in beta_ED, as a fraction of the speciation rate ==\n")
for (smp in unique(d$sampler)) { s <- d[d$sampler == smp & d$ok, ]
  e <- (s$bED - s$bED_true) / B0
  cat(sprintf("  %-14s n=%d  bias %+.4f  rmse %.4f  median ESS %5.1f  median %4.0fs\n",
              smp, nrow(s), mean(e), sqrt(mean(e^2)),
              stats::median(s$ESS, na.rm = TRUE), stats::median(s$sec))) }
