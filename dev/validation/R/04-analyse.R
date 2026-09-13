# ---------------------------------------------------------------------------
# 04-analyse.R — summaries, figures, tables, report.  Reads results/<tier>/jobs
# only; it never runs a fit.  The exact log-likelihood is re-evaluated from the
# stored brts and pars (cheap), so every number in the report is a pure
# function of the job RDS files.
#
#   Rscript 04-analyse.R --tier smoke|main|ext [--lib PATH]
#
# Writes  results/<tier>/fits.csv, fhat.csv, pipeline.csv, cells.csv
#         tables/<tier>-T*.md
#         figures/<tier>-F*.pdf
#         results/<tier>/report.md
#
# Estimands (see 00-design.md §5):
#   delta_ell = ll_exact(theta_hat) - ll_exact(theta_MLE)      <= 0
#   z_j       = (theta_hat_j - theta_MLE_j) / SE_j
#   e         = reported loglik - ll_exact(theta_hat)
#   decomposition per cell: bias B, between-tree SD of the systematic error,
#   and the pooled within-tree Monte Carlo SD.
# ---------------------------------------------------------------------------

R_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else
    "/Users/pancho/Code/emphasis/dev/validation/R"
})
source(file.path(R_DIR, "00-common.R"))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, default = NULL) {
  i <- which(args == k); if (length(i)) args[i[1] + 1L] else default
}
TIER <- getarg("--tier", "smoke")
RES  <- file.path(VAL_DIR$results, TIER)
JOBS <- file.path(RES, "jobs")
stopifnot(dir.exists(JOBS))

files <- list.files(JOBS, pattern = "\\.rds$", full.names = TRUE)
cat(sprintf("[04-analyse] tier = %s, %d job files\n", TIER, length(files)))
rows <- lapply(files, function(f) tryCatch(readRDS(f), error = function(e) NULL))
rows <- Filter(Negate(is.null), rows)

# --- build-fingerprint check ----------------------------------------------
keys <- unique(vapply(rows, function(r)
  if (is.null(r$build)) NA_character_ else val_fingerprint_key(r$build), ""))
keys <- keys[!is.na(keys)]
if (length(keys) > 1L)
  stop("rows from more than one emphasis build in ", JOBS, ":\n  ",
       paste(keys, collapse = "\n  "))
BUILD <- if (length(keys)) keys else "unknown"

g1 <- function(r, k, d = NA) { v <- r[[k]]; if (is.null(v)) d else v[1] }

# --- fits table ------------------------------------------------------------
fit_rows <- Filter(function(r) r$kind %in% c("cr_fit", "dd_fit"), rows)
FITS <- if (length(fit_rows)) do.call(rbind, lapply(fit_rows, function(r) {
  p <- r$pars; np <- length(p %||% numeric(0))
  z <- r$z
  data.frame(
    job_id = r$job_id, kind = r$kind, tree_id = r$tree_id, cell = r$cell,
    n = r$n, config = r$config, sampler = r$sampler, N = r$N, init = r$init,
    rep = r$rep, box_scale = r$box_scale, outcome = r$outcome,
    error = substr(r$error %||% "", 1, 200),
    elapsed = g1(r, "elapsed"),
    stop_reason = as.character(g1(r, "stop_reason")),
    iterations = as.numeric(g1(r, "iterations")),
    at_bound = as.logical(g1(r, "at_bound")),
    p1 = if (np >= 1) p[1] else NA_real_, p2 = if (np >= 2) p[2] else NA_real_,
    p3 = if (np >= 3) p[3] else NA_real_, p4 = if (np >= 4) p[4] else NA_real_,
    mle1 = if (!is.null(r$mle)) r$mle[1] else NA_real_,
    mle2 = if (!is.null(r$mle) && length(r$mle) >= 2) r$mle[2] else NA_real_,
    z1 = if (!is.null(z)) z[1] else NA_real_,
    z2 = if (!is.null(z) && length(z) >= 2) z[2] else NA_real_,
    loglik = g1(r, "loglik"), loglik_var = g1(r, "loglik_var"),
    ll_exact_hat = g1(r, "ll_exact_hat"), ll_exact_mle = g1(r, "ll_exact_mle"),
    delta_ell = g1(r, "delta_ell"),
    e_at_hat = g1(r, "e_at_hat"), e_at_prev = g1(r, "e_at_prev"),
    dist_init_se = g1(r, "dist_init_se"),
    K_hat = g1(r, "K_hat"), dd_map_status = as.character(g1(r, "dd_map_status")),
    ESS = if (!is.null(r$final_IS)) r$final_IS$ESS else NA_real_,
    n_rejected = if (!is.null(r$final_IS)) r$final_IS$n_rejected else NA_real_,
    rzw = if (!is.null(r$final_IS)) r$final_IS$rejected_zero_weights else NA_real_,
    acc = if (!is.null(r$final_IS)) r$final_IS$acc else NA_real_,
    sd_lw = if (!is.null(r$final_IS)) r$final_IS$sd_lw else NA_real_,
    stringsAsFactors = FALSE)
})) else data.frame()

# --- fhat table ------------------------------------------------------------
fh_rows <- Filter(function(r) r$kind == "fhat" && !is.null(r$fhat_rows), rows)
FHAT <- if (length(fh_rows)) do.call(rbind, lapply(fh_rows, function(r) {
  d <- r$fhat_rows
  rownames(d) <- NULL
  data.frame(job_id = r$job_id, tree_id = r$tree_id, cell = r$cell, n = r$n,
             outcome = r$outcome, d,
             lambda_pt = r$grid[d$point, 1], mu_pt = r$grid[d$point, 2],
             se1 = r$se_scale[1], se2 = r$se_scale[2],
             mle1 = r$mle[1], mle2 = if (length(r$mle) >= 2) r$mle[2] else NA,
             row.names = NULL, stringsAsFactors = FALSE)
})) else data.frame()

# --- pipeline table --------------------------------------------------------
pp_rows <- Filter(function(r) r$kind == "pipeline", rows)
PIPE <- if (length(pp_rows)) do.call(rbind, lapply(pp_rows, function(r) {
  data.frame(job_id = r$job_id, tree_id = r$tree_id, cell = r$cell, n = r$n,
             outcome = r$outcome, elapsed = g1(r, "elapsed"),
             error = substr(r$error %||% "", 1, 200),
             best_stage = as.character(g1(r, "best_stage")),
             delta_ell = g1(r, "delta_ell"),
             ll_exact_hat = g1(r, "ll_exact_hat"),
             loglik = g1(r, "loglik"),
             box_contains_mle0 = as.logical(g1(r, "box_contains_mle0")),
             box_contains_mle1 = as.logical(g1(r, "box_contains_mle1")),
             d_to_mle0 = if (!is.null(r$pars) && !is.null(r$mle))
               sqrt(sum((r$pars - r$mle)^2)) else NA_real_,
             d_to_mle1 = if (!is.null(r$pars) && !is.null(r$mle_cond1) &&
                             length(r$mle_cond1) == length(r$pars))
               sqrt(sum((r$pars - r$mle_cond1)^2)) else NA_real_,
             stringsAsFactors = FALSE)
})) else data.frame()

for (nm in c("FITS", "FHAT", "PIPE")) {
  d <- get(nm)
  if (nrow(d)) utils::write.csv(d, file.path(RES, sprintf("%s.csv", tolower(nm))),
                                row.names = FALSE)
}

# ---------------------------------------------------------------------------
# CELL SUMMARIES
# ---------------------------------------------------------------------------
qn <- function(x, p) if (all(is.na(x))) NA_real_ else
  unname(stats::quantile(x, p, na.rm = TRUE))

boot_ci <- function(v, B = 2000L) {
  v <- v[is.finite(v)]
  if (length(v) < 3L) return(c(NA_real_, NA_real_))
  m <- replicate(B, mean(sample(v, length(v), TRUE)))
  unname(stats::quantile(m, c(0.025, 0.975)))
}

CELLS <- data.frame()
if (nrow(FITS)) {
  key <- with(FITS, paste(cell, config, sampler, N, init, box_scale, sep = "|"))
  CELLS <- do.call(rbind, lapply(split(seq_len(nrow(FITS)), key), function(ix) {
    d <- FITS[ix, ]
    ok <- d[d$outcome == "ok", ]
    # per-tree replicate structure
    bt <- split(ok, ok$tree_id)
    bias_z1 <- vapply(bt, function(x) mean(x$z1, na.rm = TRUE), 1)
    bias_z2 <- vapply(bt, function(x) mean(x$z2, na.rm = TRUE), 1)
    mcv_z1  <- vapply(bt, function(x) if (nrow(x) > 1) stats::var(x$z1, na.rm = TRUE) else NA_real_, 1)
    mcv_z2  <- vapply(bt, function(x) if (nrow(x) > 1) stats::var(x$z2, na.rm = TRUE) else NA_real_, 1)
    ci1 <- boot_ci(bias_z1)
    data.frame(
      cell = d$cell[1], config = d$config[1], sampler = d$sampler[1],
      N = d$N[1], init = d$init[1], box_scale = d$box_scale[1],
      n = d$n[1], n_fits = nrow(d), n_trees = length(bt),
      frac_ok = mean(d$outcome == "ok"),
      frac_converged = mean(ok$stop_reason == "converged", na.rm = TRUE),
      frac_timeout = mean(d$outcome == "timeout"),
      med_delta_ell = qn(ok$delta_ell, 0.5), p10_delta_ell = qn(ok$delta_ell, 0.10),
      bias_z1 = mean(bias_z1, na.rm = TRUE), bias_z1_lo = ci1[1], bias_z1_hi = ci1[2],
      bias_z2 = mean(bias_z2, na.rm = TRUE),
      sd_between_z1 = stats::sd(bias_z1, na.rm = TRUE),
      mc_sd_z1 = sqrt(mean(mcv_z1, na.rm = TRUE)),
      mc_sd_z2 = sqrt(mean(mcv_z2, na.rm = TRUE)),
      med_abs_e = qn(abs(ok$e_at_hat), 0.5), p90_abs_e = qn(abs(ok$e_at_hat), 0.90),
      med_abs_e_prev = qn(abs(ok$e_at_prev), 0.5),
      med_ESS = qn(ok$ESS, 0.5), med_iter = qn(ok$iterations, 0.5),
      med_sec = qn(ok$elapsed, 0.5), frac_at_bound = mean(ok$at_bound, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }))
  rownames(CELLS) <- NULL
  utils::write.csv(CELLS, file.path(RES, "cells.csv"), row.names = FALSE)
}

# --- decision rules --------------------------------------------------------
verdicts <- data.frame()
if (nrow(CELLS)) {
  verdicts <- do.call(rbind, lapply(seq_len(nrow(CELLS)), function(i) {
    r <- CELLS[i, ]
    rules <- if (grepl("^dd", r$cell)) VAL_RULES$dd else VAL_RULES$cr
    f <- character(0)
    if (!is.na(r$frac_converged) && r$frac_converged < rules$converged_frac)
      f <- c(f, "i:convergence")
    if (r$frac_timeout > 0) f <- c(f, "i:timeout")
    if (!is.na(r$med_delta_ell) && r$med_delta_ell < -rules$delta_ell_median)
      f <- c(f, "ii:median_deficit")
    if (!is.na(r$p10_delta_ell) && r$p10_delta_ell < -rules$delta_ell_p90)
      f <- c(f, "ii:tail_deficit")
    if (!is.na(r$bias_z1) && abs(r$bias_z1) > rules$bias_se) f <- c(f, "iii:bias")
    if (!is.na(r$mc_sd_z1) && r$mc_sd_z1 > rules$mc_sd_se) f <- c(f, "iv:mc_sd")
    thr_med <- if (r$sampler == "bdi" && !grepl("^dd", r$cell)) rules$e_bdi else rules$e_median
    thr_p90 <- if (r$sampler == "bdi" && !grepl("^dd", r$cell)) rules$e_bdi else rules$e_p90
    if (!is.na(r$med_abs_e) && r$med_abs_e > thr_med) f <- c(f, "v:loglik_median")
    if (!is.na(r$p90_abs_e) && r$p90_abs_e > thr_p90) f <- c(f, "v:loglik_tail")
    data.frame(cell = r$cell, config = r$config, sampler = r$sampler, N = r$N,
               init = r$init,
               verdict = if (!length(f)) "adequate" else
                 if (identical(f, "iv:mc_sd")) "adequate with larger N" else "inadequate",
               failing = paste(f, collapse = ","), stringsAsFactors = FALSE)
  }))
}

# ---------------------------------------------------------------------------
# FIGURES
# ---------------------------------------------------------------------------
fig <- function(name, expr, w = 7, h = 5) {
  f <- file.path(VAL_DIR$figures, sprintf("%s-%s.pdf", TIER, name))
  grDevices::pdf(f, width = w, height = h); on.exit(grDevices::dev.off())
  tryCatch(expr, error = function(e) {
    graphics::plot.new(); graphics::title(paste("no data:", conditionMessage(e)))
  })
  invisible(f)
}
pal <- c(bdi = "#2c7fb8", dynamic_fresh = "#d95f0e")

ok_fits <- if (nrow(FITS)) FITS[FITS$outcome == "ok" & is.finite(FITS$delta_ell), ] else FITS

# F1 deficit by n x cell, sampler, at N = 200
fig("F1-deficit", {
  d <- ok_fits[ok_fits$N == 200 & ok_fits$init == "far" & ok_fits$box_scale == 1, ]
  if (!nrow(d)) stop("no N=200 far fits")
  graphics::par(mar = c(9, 4, 3, 1))
  grp <- paste(d$cell, d$sampler, sep = "\n")
  graphics::boxplot(d$delta_ell ~ grp, las = 2, cex.axis = 0.55,
                    col = pal[d$sampler[match(sort(unique(grp)), grp)]],
                    ylab = expression(Delta*"l  (nats)"), xlab = "",
                    main = "F1  log-likelihood deficit at the exact MLE (N = 200)")
  graphics::abline(h = c(0, -0.05, -0.5), lty = c(1, 3, 2), col = "grey40")
})

# F2 variance components of z1
fig("F2-variance", {
  if (!nrow(CELLS)) stop("no cells")
  d <- CELLS[CELLS$N == 200 & CELLS$init == "far" & CELLS$box_scale == 1, ]
  if (!nrow(d)) stop("no cells")
  m <- rbind(abs(d$bias_z1), d$sd_between_z1, d$mc_sd_z1)
  m[!is.finite(m)] <- 0
  graphics::par(mar = c(9, 4, 3, 1))
  graphics::barplot(m, beside = FALSE, las = 2, cex.names = 0.55,
                    names.arg = paste(d$cell, d$sampler, sep = "\n"),
                    col = c("#31a354", "#addd8e", "#f0f0f0"),
                    ylab = "SE units",
                    main = "F2  |bias| / between-tree systematic / Monte Carlo")
  graphics::legend("topright", c("|bias|", "between-tree", "Monte Carlo"),
                   fill = c("#31a354", "#addd8e", "#f0f0f0"), bty = "n", cex = 0.7)
})

# F3 -delta_ell vs N
fig("F3-Nscaling", {
  d <- ok_fits[ok_fits$init == "far" & ok_fits$box_scale == 1, ]
  a <- stats::aggregate(delta_ell ~ N + sampler, data = d, FUN = stats::median)
  if (nrow(a) < 2) stop("need at least two N levels")
  graphics::plot(a$N, pmax(-a$delta_ell, 1e-6), log = "xy", pch = 19,
                 col = pal[a$sampler], xlab = "Monte Carlo sample size N",
                 ylab = expression(-Delta*"l  (nats, median)"),
                 main = "F3  error vs N (slope -1 = MC-dominated, flat = optimiser floor)")
  for (s in unique(a$sampler)) {
    b <- a[a$sampler == s, ]
    if (nrow(b) > 1) graphics::lines(b$N, pmax(-b$delta_ell, 1e-6), col = pal[s])
  }
  graphics::legend("topright", names(pal), col = pal, lty = 1, pch = 19, bty = "n")
})

# F4 init_far vs init_mle, paired
fig("F4-init", {
  d <- ok_fits[ok_fits$N == 200 & ok_fits$box_scale == 1, ]
  a <- stats::aggregate(delta_ell ~ tree_id + sampler + init, data = d, FUN = mean)
  w <- stats::reshape(a, idvar = c("tree_id", "sampler"), timevar = "init",
                      direction = "wide")
  if (!all(c("delta_ell.far", "delta_ell.mle") %in% names(w))) stop("need both inits")
  graphics::plot(w$delta_ell.mle, w$delta_ell.far, pch = 19, col = pal[w$sampler],
                 xlab = expression(Delta*"l  started at the exact MLE"),
                 ylab = expression(Delta*"l  started far"),
                 main = "F4  cost of starting away from the MLE")
  graphics::abline(0, 1, col = "grey50")
  graphics::legend("bottomright", names(pal), col = pal, pch = 19, bty = "n")
})

# F5 reported-loglik error vs ESS
fig("F5-loglik", {
  d <- ok_fits
  if (!nrow(d)) stop("no fits")
  graphics::par(mfrow = c(1, 2))
  t1 <- d[d$sampler == "dynamic_fresh", ]
  graphics::plot(t1$ESS / t1$N, t1$e_at_hat, pch = 19, col = pal["dynamic_fresh"],
                 xlab = "ESS / N", ylab = "reported - exact (nats)",
                 main = "thinning")
  graphics::abline(h = 0, col = "grey50")
  b1 <- d[d$sampler == "bdi", ]
  v <- abs(b1$e_at_hat); v <- v[is.finite(v)]
  if (length(v)) graphics::hist(log10(pmax(v, 1e-16)), breaks = 20, col = pal["bdi"],
                                xlab = "log10 |reported - exact|", main = "BDI")
  else { graphics::plot.new(); graphics::title("BDI: no finite values") }
})

# F6 g(theta) across the fixed-theta grid
fig("F6-surface", {
  if (!nrow(FHAT)) stop("no fixed-theta rows")
  d <- FHAT[is.finite(FHAT$g), ]
  a <- stats::aggregate(g ~ tree_id + point + sampler, data = d, FUN = mean)
  a$x <- (d$lambda_pt[match(paste(a$tree_id, a$point), paste(d$tree_id, d$point))] -
          d$mle1[match(paste(a$tree_id, a$point), paste(d$tree_id, d$point))]) /
         d$se1[match(paste(a$tree_id, a$point), paste(d$tree_id, d$point))]
  graphics::plot(a$x, a$g, pch = 19, col = pal[a$sampler],
                 xlab = "(lambda - lambda_MLE) / SE", ylab = "g(theta) = fhat - exact (nats)",
                 main = "F6  importance-sampling bias across theta")
  graphics::abline(h = 0, col = "grey50")
  graphics::legend("topright", names(pal), col = pal, pch = 19, bty = "n")
})

# F7 recovery
fig("F7-recovery", {
  d <- ok_fits[ok_fits$kind == "cr_fit" & ok_fits$N == 200 & ok_fits$init == "far", ]
  if (!nrow(d)) stop("no cr fits")
  tt <- readRDS(file.path(VAL_DIR$data, sprintf("trees-%s.rds", TIER)))$trees
  lam_gen <- vapply(d$tree_id, function(i) tt[[i]]$lambda_gen, 1)
  graphics::plot(d$mle1 - lam_gen, d$p1 - lam_gen, pch = 19, col = pal[d$sampler],
                 xlab = "exact MLE - generating lambda",
                 ylab = "emphasis - generating lambda",
                 main = "F7  estimator error vs MLE sampling error")
  graphics::abline(0, 1, col = "grey50")
})

# F8 pipeline
fig("F8-pipeline", {
  if (!nrow(PIPE)) stop("no pipeline rows")
  graphics::par(mfrow = c(1, 2))
  graphics::plot(PIPE$d_to_mle0, PIPE$d_to_mle1, pch = 19,
                 col = ifelse(isTRUE(PIPE$box_contains_mle0), "#2c7fb8", "#d95f0e"),
                 xlab = "distance to bd_ML(cond = 0)", ylab = "distance to bd_ML(cond = 1)",
                 main = "F8a  which target")
  graphics::abline(0, 1, col = "grey50")
  graphics::barplot(c(mle0 = mean(PIPE$box_contains_mle0, na.rm = TRUE),
                      mle1 = mean(PIPE$box_contains_mle1, na.rm = TRUE)),
                    ylim = c(0, 1), ylab = "fraction inside auto_bounds box",
                    main = "F8b  H74")
})

# F9 dd deficits
fig("F9-dd", {
  d <- ok_fits[ok_fits$kind == "dd_fit", ]
  if (!nrow(d)) stop("no dd fits")
  graphics::par(mar = c(8, 4, 3, 1))
  graphics::boxplot(d$delta_ell ~ paste(d$cell, d$sampler, sep = "\n"), las = 2,
                    cex.axis = 0.6, ylab = expression(Delta*"l  (nats)"), xlab = "",
                    main = "F9  dd deficit by regime and sampler")
  graphics::abline(h = c(0, -0.2, -1), lty = c(1, 3, 2), col = "grey40")
})

# F10 stop reasons and iterations
fig("F10-stops", {
  if (!nrow(FITS)) stop("no fits")
  graphics::par(mfrow = c(1, 2), mar = c(8, 4, 3, 1))
  tb <- table(FITS$outcome, FITS$sampler)
  graphics::barplot(tb, legend.text = rownames(tb), las = 2, main = "outcome",
                    args.legend = list(bty = "n", cex = 0.7))
  d <- ok_fits
  graphics::boxplot(d$iterations ~ paste(d$n, d$sampler, sep = "\n"), las = 2,
                    cex.axis = 0.6, xlab = "", ylab = "iterations", main = "iterations")
})

# ---------------------------------------------------------------------------
# TABLES + REPORT
# ---------------------------------------------------------------------------
md_table <- function(d, digits = 3) {
  d <- as.data.frame(d)
  for (j in seq_along(d)) if (is.numeric(d[[j]])) d[[j]] <- round(d[[j]], digits)
  c(paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
    apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
wtab <- function(name, lines) {
  f <- file.path(VAL_DIR$tables, sprintf("%s-%s.md", TIER, name))
  writeLines(lines, f); f
}

if (nrow(CELLS)) wtab("T1-cells", md_table(CELLS[, c(
  "cell", "config", "sampler", "N", "init", "n_trees", "n_fits", "frac_ok",
  "frac_converged", "med_delta_ell", "p10_delta_ell", "bias_z1", "mc_sd_z1",
  "med_abs_e", "p90_abs_e", "med_ESS", "med_iter", "med_sec")]))
sc <- tryCatch(readRDS(file.path(VAL_DIR$results, "selfcheck.rds")), error = function(e) NULL)
if (!is.null(sc)) wtab("T2-references", c(
  "# T2 reference reconciliation (tier 0)", "",
  sprintf("- build: %s", val_fingerprint_key(sc$fingerprint)),
  sprintf("- aborting assertions failed: %s",
          if (!length(sc$fails)) "none" else paste(sc$fails, collapse = ", ")),
  sprintf("- sentinels: %s", paste(sprintf("%s=%s", names(sc$sentinels),
                                           unlist(sc$sentinels)), collapse = " ")),
  sprintf("- max |c_n| (emphasis BDI fhat - DDD bd_loglik, btorph = 1, mu < lambda): %.2e",
          max(abs(sc$c_n$c_n[sc$c_n$mu < sc$c_n$lambda]), na.rm = TRUE)),
  sprintf("- thinning fhat - exact: %s",
          paste(sprintf("%s %+.4f (ESS %.0f)", names(sc$thinning_gap),
                        sc$thinning_gap, sc$thinning_ess), collapse = "; "))))
if (nrow(verdicts)) wtab("T3-verdicts", md_table(verdicts))
if (nrow(PIPE)) wtab("T4-pipeline", md_table(PIPE[, c(
  "tree_id", "n", "outcome", "best_stage", "delta_ell", "box_contains_mle0",
  "box_contains_mle1", "d_to_mle0", "d_to_mle1", "elapsed")]))
if (nrow(FHAT)) {
  a <- stats::aggregate(cbind(g, ess) ~ tree_id + sampler, data = FHAT,
                        FUN = function(x) mean(x, na.rm = TRUE))
  wtab("T5-surface", md_table(a))
}

rep_lines <- c(
  sprintf("# emphasis validation — tier %s", TIER), "",
  sprintf("build: `%s`", BUILD),
  sprintf("generated: %s", format(Sys.time())), "",
  "## Outcomes (every job, no exclusions)", "",
  md_table(as.data.frame(table(
    kind = vapply(rows, function(r) r$kind, ""),
    outcome = vapply(rows, function(r) r$outcome %||% "NA", "")))),
  "", "## Cell summary (T1)", "",
  if (nrow(CELLS)) md_table(CELLS[, c("cell", "config", "sampler", "N", "init",
                                      "n_fits", "med_delta_ell", "bias_z1",
                                      "mc_sd_z1", "med_abs_e", "med_ESS",
                                      "med_iter", "med_sec")]) else "-",
  "", "## Verdicts (T3)", "",
  if (nrow(verdicts)) md_table(verdicts) else "-",
  "", "## Figures", "",
  paste0("- ", basename(list.files(VAL_DIR$figures,
                                   pattern = paste0("^", TIER, "-")))),
  "",
  "Timings and iteration counts on a pre-wave-1 build are not final; the",
  "smoke tier's cost table is what sizes the main tier.")
writeLines(rep_lines, file.path(RES, "report.md"))

cat(sprintf("  fits %d | fhat rows %d | pipeline %d | cells %d\n",
            nrow(FITS), nrow(FHAT), nrow(PIPE), nrow(CELLS)))
cat(sprintf("  wrote %s/report.md, %d figures, %d tables\n", RES,
            length(list.files(VAL_DIR$figures, pattern = paste0("^", TIER, "-"))),
            length(list.files(VAL_DIR$tables, pattern = paste0("^", TIER, "-")))))

# --- smoke-tier cost calibration -------------------------------------------
if (TIER == "smoke" && nrow(FITS)) {
  # every kind, not just the mcem fits: 03-fit.R uses this table to size and
  # order the main tier, and the fixed-theta and pipeline jobs are the long ones
  allc <- do.call(rbind, lapply(rows, function(r) if (!identical(r$outcome, "ok"))
    NULL else data.frame(kind = r$kind,
                         sampler = {s <- as.character(r$sampler %||% NA)
                                    if (is.na(s)) "none" else s},
                         N = r$N, n = r$n, elapsed = r$elapsed,
                         stringsAsFactors = FALSE)))
  cal <- stats::aggregate(elapsed ~ kind + sampler + N + n, data = allc,
                          FUN = stats::median)
  utils::write.csv(cal, file.path(RES, "cost-calibration.csv"), row.names = FALSE)
  cat("\n  measured cost per configuration (seconds, 1 thread):\n")
  print(cal, row.names = FALSE)
  cat("\n  main-tier sizing rule: trees_per_cell = min(20, floor(cell_budget /",
      "cost_per_tree)), never below 10; per-job timeout = 3 x this cost,",
      "minimum 120 s.\n")
  # A configuration that timed out in the smoke tier has NO measured cost; its
  # main-tier timeout must come from val_timeout(), not from this table.
  to <- FITS[FITS$outcome %in% c("timeout", "crash"), ]
  if (nrow(to)) {
    u <- unique(to[, c("kind", "sampler", "N", "n", "config")])
    cat("\n  no measured cost (timed out or crashed in the smoke tier):\n")
    print(u, row.names = FALSE)
    cat("  -> these cells keep the val_timeout() defaults in the main tier",
        "and are re-measured there.\n")
  }
}
