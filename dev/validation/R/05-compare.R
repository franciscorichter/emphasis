# 05-compare.R — pair two completed result directories of the same tier and
# report what moved between them.
#
#   Rscript R/05-compare.R --a results-run2-20260916 --b results-run3-20260917 [--tier main]
#
# Both directories must hold results/<tier>/{pipe,init,cells,fits}.csv as
# 04-analyse.R writes them.  Rows are paired on job_id (pipe, fits), on
# tree x cell x cfg (init) and on cell x config (cells), so a job present in
# one run only is counted, not compared.
#
# What it answers, in order:
#   1. pipeline arm  — does the box contain the exact MLE (cond = 0 / 1), what
#      the fit lands at, and what it costs.  This is the arm auto_bounds feeds,
#      so it is where a bounds change shows.
#   2. init arm      — the three starting points handed to the same MCEM run,
#      before and after.
#   3. fit cells     — C1..C8 and the DD configs use a box scaled around the
#      MLE, not auto_bounds; a change confined to auto_bounds must leave these
#      within Monte Carlo noise.  Cells whose median deficit moved by more than
#      0.05 nats or whose bias moved by more than 0.1 SE are listed.
#   4. outcomes      — ok / timeout / error per job kind.
#
# Output: a markdown report on stdout.  Nothing is written under either
# results directory.

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, default = NULL) {
  i <- which(args == k); if (length(i)) args[i[1] + 1L] else default
}
A    <- getarg("--a"); B <- getarg("--b")
TIER <- getarg("--tier", "main")
if (is.null(A) || is.null(B)) stop("usage: 05-compare.R --a <dir> --b <dir> [--tier main]")

rd <- function(dir, f) {
  p <- file.path(dir, "results", TIER, f)
  if (!file.exists(p)) stop("missing: ", p)
  utils::read.csv(p, stringsAsFactors = FALSE)
}
la <- basename(A); lb <- basename(B)

num <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
pct <- function(x) ifelse(is.na(x), "NA", sprintf("%.0f%%", 100 * x))
row <- function(...) cat("|", paste(..., sep = " | "), "|\n")
hdr <- function(...) { row(...); cat("|", paste(rep("---", length(c(...))), collapse = " | "), "|\n") }
med <- function(x) if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)

cat(sprintf("# %s vs %s — tier %s\n\n", la, lb, TIER))
cat(sprintf("A = `%s`\nB = `%s`\n\n", A, B))
for (d in c(A, B)) {
  bt <- file.path(d, "results", TIER, "build.txt")
  if (file.exists(bt)) cat(sprintf("- `%s` build: `%s`\n", basename(d), trimws(readLines(bt, n = 1))))
}
cat("\n")

# --- 1. pipeline arm --------------------------------------------------------
pa <- rd(A, "pipe.csv"); pb <- rd(B, "pipe.csv")
pj <- merge(pa, pb, by = "job_id", suffixes = c(".a", ".b"))
cat("## 1. Pipeline arm (auto_bounds → GAM → CEM → MCEM)\n\n")
cat(sprintf("%d jobs in A, %d in B, %d paired.\n\n", nrow(pa), nrow(pb), nrow(pj)))

psum <- function(p) c(
  n_ok      = sum(p$outcome == "ok"),
  n_timeout = sum(p$outcome == "timeout"),
  contain0  = mean(p$box_contains_mle0, na.rm = TRUE),
  contain1  = mean(p$box_contains_mle1, na.rm = TRUE),
  med_dl    = med(p$delta_ell),
  p10_dl    = if (all(is.na(p$delta_ell))) NA else unname(stats::quantile(p$delta_ell, 0.1, na.rm = TRUE)),
  within01  = mean(p$delta_ell >= -0.1, na.rm = TRUE),
  med_sec   = med(p$elapsed))
sa <- psum(pa); sb <- psum(pb)
hdr("quantity", la, lb)
row("fits ok / timeout", sprintf("%d / %d", sa["n_ok"], sa["n_timeout"]), sprintf("%d / %d", sb["n_ok"], sb["n_timeout"]))
row("box contains bd_ML(cond = 0)", pct(sa["contain0"]), pct(sb["contain0"]))
row("box contains bd_ML(cond = 1)", pct(sa["contain1"]), pct(sb["contain1"]))
row("median Δℓ", num(sa["med_dl"]), num(sb["med_dl"]))
row("10th pct Δℓ", num(sa["p10_dl"]), num(sb["p10_dl"]))
row("within 0.1 nats", pct(sa["within01"]), pct(sb["within01"]))
row("median wall-clock (s)", num(sa["med_sec"], 1), num(sb["med_sec"], 1))
cat("\n")

cat("### By cell\n\n")
hdr("cell", "n", paste0("contain1 ", la), paste0("contain1 ", lb), paste0("med Δℓ ", la), paste0("med Δℓ ", lb), paste0("med s ", la), paste0("med s ", lb))
for (cl in sort(unique(pj$cell.a))) {
  s <- pj[pj$cell.a == cl, ]
  row(cl, nrow(s),
      pct(mean(s$box_contains_mle1.a, na.rm = TRUE)), pct(mean(s$box_contains_mle1.b, na.rm = TRUE)),
      num(med(s$delta_ell.a)), num(med(s$delta_ell.b)),
      num(med(s$elapsed.a), 0), num(med(s$elapsed.b), 0))
}
cat("\n")

cat("### Paired, split on whether the box contains bd_ML(cond = 1)\n\n")
hdr("containment A → B", "jobs", paste0("med Δℓ ", la), paste0("med Δℓ ", lb), "med Δ(Δℓ) B−A")
ok <- !is.na(pj$box_contains_mle1.a) & !is.na(pj$box_contains_mle1.b)
for (tr in list(c(FALSE, TRUE), c(TRUE, TRUE), c(TRUE, FALSE), c(FALSE, FALSE))) {
  s <- pj[ok & pj$box_contains_mle1.a == tr[1] & pj$box_contains_mle1.b == tr[2], ]
  if (nrow(s) == 0) next
  row(sprintf("%s → %s", if (tr[1]) "yes" else "no", if (tr[2]) "yes" else "no"), nrow(s),
      num(med(s$delta_ell.a)), num(med(s$delta_ell.b)), num(med(s$delta_ell.b - s$delta_ell.a)))
}
cat("\n")
bs <- table(a = pa$best_stage, useNA = "ifany"); bs2 <- table(b = pb$best_stage, useNA = "ifany")
cat(sprintf("Best stage, %s: %s.  %s: %s.\n\n", la,
            paste(names(bs), bs, sep = "=", collapse = ", "), lb,
            paste(names(bs2), bs2, sep = "=", collapse = ", ")))

# Where the time goes: each finished pipeline job carries run_log, a data
# frame with one row per stage (bounds, gam, cem, mcem) and its elapsed
# seconds.  A timed-out job has no run_log, so this is the cost of the jobs
# that finished; the timeout count next to it says how many did not.
stage_times <- function(dir) {
  fs <- list.files(file.path(dir, "results", TIER, "jobs"), pattern = "--C10[.]rds$", full.names = TRUE)
  rows <- lapply(fs, function(f) {
    x <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(x) || !identical(x$outcome, "ok") || !is.data.frame(x$run_log)) return(NULL)
    data.frame(job_id = x$job_id, cell = x$cell, stage = x$run_log$stage,
               elapsed = x$run_log$elapsed, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
sta <- stage_times(A); stb <- stage_times(B)
if (!is.null(sta) && !is.null(stb)) {
  cat("### Wall-clock by pipeline stage (finished jobs only)\n\n")
  stages <- c("bounds", "gam", "cem", "mcem")
  hdr("stage", paste0("median s ", la), paste0("median s ", lb), "ratio B/A", paste0("share of job ", la), paste0("share of job ", lb))
  tot_a <- sum(sta$elapsed); tot_b <- sum(stb$elapsed)
  for (s in stages) {
    a <- sta$elapsed[sta$stage == s]; b <- stb$elapsed[stb$stage == s]
    if (!length(a) && !length(b)) next
    row(s, num(med(a), 1), num(med(b), 1), num(med(b) / med(a), 2),
        pct(sum(a) / tot_a), pct(sum(b) / tot_b))
  }
  cat(sprintf("\nFinished pipeline jobs with a stage log: %d in %s, %d in %s; timed out: %d and %d.\n\n",
              length(unique(sta$job_id)), la, length(unique(stb$job_id)), lb,
              sum(pa$outcome == "timeout"), sum(pb$outcome == "timeout")))
  cat("### GAM stage by cell (median s)\n\n")
  hdr("cell", paste0("gam ", la), paste0("gam ", lb), paste0("whole job ", la), paste0("whole job ", lb))
  for (cl in sort(unique(c(sta$cell, stb$cell)))) {
    ga <- sta$elapsed[sta$cell == cl & sta$stage == "gam"]; gb <- stb$elapsed[stb$cell == cl & stb$stage == "gam"]
    ja <- tapply(sta$elapsed[sta$cell == cl], sta$job_id[sta$cell == cl], sum)
    jb <- tapply(stb$elapsed[stb$cell == cl], stb$job_id[stb$cell == cl], sum)
    row(cl, num(med(ga), 0), num(med(gb), 0), num(med(ja), 0), num(med(jb), 0))
  }
  cat("\n")
}

# --- 2. init arm -------------------------------------------------------------
ia <- rd(A, "init.csv"); ib <- rd(B, "init.csv")
cat("## 2. Initialiser arm (one auto_bounds box per tree, three starting points)\n\n")
cat("`box` is whether that tree's auto_bounds box contains the exact MLE; `start`/`after` are the deficit handed to MCEM and the deficit it returns.\n\n")
ij <- merge(ia, ib, by = c("tree", "cell", "cfg"), suffixes = c(".a", ".b"))
cat(sprintf("%d rows in A, %d in B, %d paired.\n\n", nrow(ia), nrow(ib), nrow(ij)))
hdr("cfg", "n A/B", paste0("box contains MLE ", la), paste0("box contains MLE ", lb),
    paste0("med start Δℓ ", la), paste0("med start Δℓ ", lb),
    paste0("med after Δℓ ", la), paste0("med after Δℓ ", lb),
    paste0("within 0.1 ", la), paste0("within 0.1 ", lb), paste0("med s ", la), paste0("med s ", lb))
for (cf in sort(unique(c(ia$cfg, ib$cfg)))) {
  a <- ia[ia$cfg == cf, ]; b <- ib[ib$cfg == cf, ]
  row(cf, sprintf("%d/%d", nrow(a), nrow(b)),
      pct(mean(a$box, na.rm = TRUE)), pct(mean(b$box, na.rm = TRUE)),
      num(med(a$start)), num(med(b$start)),
      num(med(a$after)), num(med(b$after)),
      pct(mean(a$after >= -0.1, na.rm = TRUE)), pct(mean(b$after >= -0.1, na.rm = TRUE)),
      num(med(a$sec), 0), num(med(b$sec), 0))
}
cat("\n")
cat("### Containment per tree (I1 rows), by cell\n\n")
hdr("cell", "trees", paste0("contained ", la), paste0("contained ", lb), "gained", "lost")
i1 <- ij[ij$cfg == "I1", ]
for (cl in sort(unique(i1$cell))) {
  s <- i1[i1$cell == cl, ]
  row(cl, nrow(s), sprintf("%d", sum(s$box.a, na.rm = TRUE)), sprintf("%d", sum(s$box.b, na.rm = TRUE)),
      sum(!s$box.a & s$box.b, na.rm = TRUE), sum(s$box.a & !s$box.b, na.rm = TRUE))
}
cat("\n### Paired on containment, deficit after MCEM (all cfgs)\n\n")
hdr("containment A → B", "rows", paste0("med after ", la), paste0("med after ", lb), "med Δ(after) B−A", paste0("med s ", la), paste0("med s ", lb))
okb <- !is.na(ij$box.a) & !is.na(ij$box.b)
for (tr in list(c(FALSE, TRUE), c(TRUE, TRUE), c(TRUE, FALSE), c(FALSE, FALSE))) {
  s <- ij[okb & ij$box.a == tr[1] & ij$box.b == tr[2], ]
  if (nrow(s) == 0) next
  row(sprintf("%s → %s", if (tr[1]) "yes" else "no", if (tr[2]) "yes" else "no"), nrow(s),
      num(med(s$after.a)), num(med(s$after.b)), num(med(s$after.b - s$after.a)),
      num(med(s$sec.a), 0), num(med(s$sec.b), 0))
}
cat("\n")

# --- 3. fit cells (regression guard) ----------------------------------------
ca <- rd(A, "cells.csv"); cb <- rd(B, "cells.csv")
key <- function(d) paste(d$cell, d$config, sep = "|")
cj <- merge(transform(ca, k = key(ca)), transform(cb, k = key(cb)), by = "k", suffixes = c(".a", ".b"))
cat("## 3. Fit cells C1–C8 / D1–D4 (boxes scaled around the MLE — a bounds change must not move these)\n\n")
cat(sprintf("%d cells in A, %d in B, %d paired.\n\n", nrow(ca), nrow(cb), nrow(cj)))
d_dl <- cj$med_delta_ell.b - cj$med_delta_ell.a
d_bz <- cj$bias_z1.b - cj$bias_z1.a
d_cv <- cj$frac_converged.b - cj$frac_converged.a
d_sc <- cj$med_sec.b / cj$med_sec.a
cat(sprintf("Median deficit moved by median %s nats (IQR %s to %s); bias_z1 by median %s SE; frac_converged by median %s; wall-clock ratio B/A median %s.\n\n",
            num(med(d_dl)), num(stats::quantile(d_dl, .25, na.rm = TRUE)), num(stats::quantile(d_dl, .75, na.rm = TRUE)),
            num(med(d_bz)), num(med(d_cv)), num(med(d_sc), 2)))
mv <- cj[(!is.na(d_dl) & abs(d_dl) > 0.05) | (!is.na(d_bz) & abs(d_bz) > 0.1), ]
if (nrow(mv) == 0) {
  cat("No cell moved by more than 0.05 nats in median deficit or 0.1 SE in bias.\n\n")
} else {
  cat(sprintf("%d cells moved beyond 0.05 nats or 0.1 SE:\n\n", nrow(mv)))
  hdr("cell", "config", paste0("med Δℓ ", la), paste0("med Δℓ ", lb), paste0("bias_z1 ", la), paste0("bias_z1 ", lb), "n_fits A/B")
  for (i in seq_len(nrow(mv))) {
    m <- mv[i, ]
    row(m$cell.a, m$config.a, num(m$med_delta_ell.a), num(m$med_delta_ell.b),
        num(m$bias_z1.a), num(m$bias_z1.b), sprintf("%d/%d", m$n_fits.a, m$n_fits.b))
  }
  cat("\n")
}
# --- 4. outcomes --------------------------------------------------------------
fa <- rd(A, "fits.csv"); fb <- rd(B, "fits.csv")
cat("## 4. Outcomes by kind\n\n")
oa <- as.data.frame(table(kind = fa$kind, outcome = fa$outcome)); ob <- as.data.frame(table(kind = fb$kind, outcome = fb$outcome))
oj <- merge(oa, ob, by = c("kind", "outcome"), all = TRUE, suffixes = c(".a", ".b"))
oj[is.na(oj)] <- 0
hdr("kind", "outcome", la, lb)
for (i in seq_len(nrow(oj))) row(oj$kind[i], oj$outcome[i], oj$Freq.a[i], oj$Freq.b[i])
cat("\n")
fj <- merge(fa[, c("job_id", "outcome", "elapsed", "at_bound")], fb[, c("job_id", "outcome", "elapsed", "at_bound")], by = "job_id", suffixes = c(".a", ".b"))
cat(sprintf("Paired fit jobs: %d.  at_bound A: %s, B: %s.  Wall-clock ratio B/A median %s.\n",
            nrow(fj), pct(mean(fj$at_bound.a, na.rm = TRUE)), pct(mean(fj$at_bound.b, na.rm = TRUE)),
            num(med(fj$elapsed.b / fj$elapsed.a), 2)))
