# ---------------------------------------------------------------------------
# 03-fit.R — the emphasis job runner.
#
#   Rscript 03-fit.R --tier smoke|main|ext [--workers 10] [--lib PATH]
#                    [--only cr_fit,dd_fit,fhat,pipeline] [--dry-run]
#                    [--force] [--no-gate]
#
# One OS process per fit (callr::r_bg), so a hang or a C++ abort costs one job.
# Every job writes results/<tier>/jobs/<job_id>.rds; a rerun skips job_ids that
# already have a file, so the run is resumable and a crash loses nothing.
#
# Scheduling: jobs are ordered by est_cost_s descending (longest first) and
# dispatched dynamically to `workers` slots.  Each job has a hard timeout;
# emphasis's own max_time is set 30 s below it, so a clean "time_budget" stop
# is recorded before the driver kills the process.
#
# The driver refuses to start unless results/GATE.ok exists (written by
# 00-selfcheck.R) and its build fingerprint matches the library in use.
# ---------------------------------------------------------------------------

R_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else
    "/Users/pancho/Code/emphasis/dev/validation/R"
})
source(file.path(R_DIR, "00-common.R"))
COMMON <- file.path(R_DIR, "00-common.R")
WORKER <- file.path(R_DIR, "03-worker.R")

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, default = NULL) {
  i <- which(args == k); if (length(i)) args[i[1] + 1L] else default
}
TIER    <- getarg("--tier", "smoke")
WORKERS <- as.integer(getarg("--workers", "10"))
LIB     <- getarg("--lib", Sys.getenv("EMPHASIS_LIB", NA_character_))
ONLY    <- if (!is.null(getarg("--only"))) strsplit(getarg("--only"), ",")[[1]] else NULL
JOBS_RE <- getarg("--jobs")   # regex on job_id, applied after --only; for stragglers
DRY     <- any(args == "--dry-run")
FORCE   <- any(args == "--force")
NOGATE  <- any(args == "--no-gate")
options(emphasis.lib = LIB)

suppressPackageStartupMessages({ library(callr); library(digest) })

RES <- file.path(VAL_DIR$results, TIER)
JOBS <- file.path(RES, "jobs")
dir.create(JOBS, recursive = TRUE, showWarnings = FALSE)

trees_file <- file.path(VAL_DIR$data, sprintf("trees-%s.rds", TIER))
refs_file  <- file.path(VAL_DIR$data, sprintf("reference-%s.rds", TIER))
stopifnot(file.exists(trees_file), file.exists(refs_file))
TR  <- readRDS(trees_file)$trees
REF <- readRDS(refs_file)$refs

# --- gate ------------------------------------------------------------------
gate <- file.path(VAL_DIR$results, "GATE.ok")
if (!NOGATE && !file.exists(gate))
  stop("results/GATE.ok is missing: run 00-selfcheck.R first.")

# --- build fingerprint -----------------------------------------------------
val_load_emphasis(LIB)
FP <- val_build_fingerprint()
fp_file <- file.path(RES, "build.txt")
key <- val_fingerprint_key(FP)
if (file.exists(fp_file)) {
  old <- readLines(fp_file)[1]
  if (!identical(old, key) && !FORCE)
    stop("results/", TIER, " holds rows from a different emphasis build:\n  on disk: ",
         old, "\n  in use : ", key,
         "\nUse a fresh results directory, or --force to override.")
} else if (!DRY) writeLines(key, fp_file)
cat(sprintf("[03-fit] tier=%s  workers=%d  build=%s\n", TIER, WORKERS, key))

# ---------------------------------------------------------------------------
# JOB TABLE
# ---------------------------------------------------------------------------
jobs <- list()
addjob <- function(...) jobs[[length(jobs) + 1L]] <<- list(...)

# Costs: the smoke tier's measured table when it exists, otherwise the pre-fix
# formula.  Used for longest-first packing and the budget check only.
CALIB <- if (TIER == "smoke") NULL else val_load_calibration("smoke")
if (!is.null(CALIB)) cat("  cost model: results/smoke/cost-calibration.csv\n") else
  cat("  cost model: pre-fix formula (no smoke calibration on disk)\n")

# Per-job hard timeout: the larger of the kind/size default and 3x the cost
# model, floor 120 s.  With a smoke calibration on disk this is the design's
# "3 x the measured cost" rule; without one it is the pre-fix default.
#
# --timeout-mult <x> scales every job's timeout by x.  It exists for re-running
# jobs that a configuration change made slower than the calibration assumed:
# delete their results/<tier>/jobs/<job_id>.rds, then re-run with --only <kind>
# and a multiplier, and the resume logic picks up only those.  The value used
# is recorded in each job's row as timeout_s.
TMULT <- as.numeric(getarg("--timeout-mult", "1"))
tmo <- function(kind, n, N, sampler = NA) {
  est <- val_est_cost(kind, sampler, n, N, calib = CALIB)
  TMULT * max(120, val_timeout(kind, n, N), 3 * est)
}

cr_ids <- names(TR)[vapply(TR, function(x) x$kind == "cr", TRUE)]
dd_ids <- names(TR)[vapply(TR, function(x) x$kind == "dd", TRUE)]

# which trees are in subset S (the first `val_subset_size` per cell)
subset_of <- function(ids) {
  cells <- vapply(TR[ids], function(x) x$cell, "")
  unlist(lapply(split(ids, cells), function(v) {
    n <- TR[[v[1]]]$n
    head(v[order(vapply(TR[v], function(x) x$index, 1L))],
         val_subset_size(n, TIER))
  }), use.names = FALSE)
}
cr_sub <- if (length(cr_ids)) subset_of(cr_ids) else character(0)
dd_sub <- if (length(dd_ids)) subset_of(dd_ids) else character(0)

# --- mcem fits -------------------------------------------------------------
cfg_cr <- val_cr_configs(TIER)
for (tid in cr_ids) {
  tt <- TR[[tid]]; n <- tt$n
  for (i in seq_len(nrow(cfg_cr))) {
    cf <- cfg_cr[i, ]
    if (cf$subset && !(tid %in% cr_sub)) next
    # the box x 10 sub-arm runs only on the n = 50, eps = 0.6 cell
    if (cf$box_scale != 1 && !(n == 50L && isTRUE(abs(tt$eps - 0.6) < 1e-9))) next
    for (r in seq_len(cf$reps)) {
      jid <- sprintf("%s--%s-r%d", tid, cf$config, r)
      addjob(job_id = jid, kind = "cr_fit", tree_id = tid, cell = tt$cell,
             n = n, config = cf$config, sampler = cf$sampler, N = cf$N,
             init = cf$init, rep = r, box_scale = cf$box_scale,
             seed = val_seed(jid),
             timeout_s = tmo("cr_fit", n, cf$N, cf$sampler),
             est_cost_s = val_est_cost("cr_fit", cf$sampler, n, cf$N, calib = CALIB))
    }
  }
}

cfg_dd <- val_dd_configs(TIER)
for (tid in dd_ids) {
  tt <- TR[[tid]]; n <- tt$n
  if (identical(REF[[tid]]$flag, "dd_ML_failed")) next   # no reference
  for (i in seq_len(nrow(cfg_dd))) {
    cf <- cfg_dd[i, ]
    if (cf$subset && !(tid %in% dd_sub)) next
    for (r in seq_len(cf$reps)) {
      jid <- sprintf("%s--%s-r%d", tid, cf$config, r)
      addjob(job_id = jid, kind = "dd_fit", tree_id = tid, cell = tt$cell,
             n = n, config = cf$config, sampler = cf$sampler, N = cf$N,
             init = cf$init, rep = r, box_scale = 1,
             seed = val_seed(jid),
             timeout_s = tmo("dd_fit", n, cf$N, cf$sampler),
             est_cost_s = val_est_cost("dd_fit", cf$sampler, n, cf$N, calib = CALIB))
    }
  }
}

# --- fixed-theta grids (C9 / D7) -------------------------------------------
for (tid in c(cr_sub, dd_sub)) {
  tt <- TR[[tid]]; n <- tt$n
  if (tt$kind == "dd" && identical(REF[[tid]]$flag, "dd_ML_failed")) next
  # the fixed-theta arm at N = 2000 costs ~9 x 3 E-steps per tree; above
  # n = 50 that is the most expensive job in the study, so N drops there
  Nfh <- if (TIER == "smoke") 500L else if (n <= 50L) 2000L else 1000L
  addjob(job_id = sprintf("%s--C9bdi", tid), kind = "fhat", tree_id = tid,
         cell = tt$cell, n = n, config = "C9", sampler = "bdi", N = Nfh,
         init = NA_character_, rep = 1L, reps = 2L, box_scale = 1,
         seed = val_seed(paste0(tid, "C9b")),
         timeout_s = tmo("fhat", n, Nfh, "bdi"),
         est_cost_s = val_est_cost("fhat", "bdi", n, Nfh, calib = CALIB))
  if (n <= 100L && TIER != "ext")  # thinning at n = 200 exceeds the 120 s C++ E-step cap
    addjob(job_id = sprintf("%s--C9thin", tid), kind = "fhat", tree_id = tid,
           cell = tt$cell, n = n, config = "C9", sampler = "dynamic_fresh",
           N = Nfh, init = NA_character_, rep = 1L, reps = 3L, box_scale = 1,
           seed = val_seed(paste0(tid, "C9t")),
           timeout_s = tmo("fhat", n, Nfh, "dynamic_fresh"),
           est_cost_s = val_est_cost("fhat", "dynamic_fresh", n, Nfh, calib = CALIB))
}

# --- pipeline (C10) --------------------------------------------------------
for (tid in if (TIER == "ext") character(0) else cr_sub) {
  tt <- TR[[tid]]
  addjob(job_id = sprintf("%s--C10", tid), kind = "pipeline", tree_id = tid,
         cell = tt$cell, n = tt$n, config = "C10", sampler = NA_character_,
         N = 200L, init = NA_character_, rep = 1L, box_scale = 1,
         seed = val_seed(paste0(tid, "C10")),
         timeout_s = tmo("pipeline", tt$n, 200L, NA),
         est_cost_s = val_est_cost("pipeline", NA, tt$n, calib = CALIB))
}

# --- initialiser comparison (I1 cem / I2 naive / I3 gam) -------------------
# Same tree, same box, three starting points.  Answers whether the
# cross-entropy search earns the draws it spends.
for (tid in if (TIER == "ext") character(0) else cr_sub) {
  tt <- TR[[tid]]
  for (cfg in c("I1", "I2", "I3")) {
    addjob(job_id = sprintf("%s--%s", tid, cfg), kind = "init", tree_id = tid,
           cell = tt$cell, n = tt$n, config = cfg, sampler = "bdi",
           N = 200L, init = cfg, rep = 1L, box_scale = 1,
           seed = val_seed(paste0(tid, cfg)),
           timeout_s = tmo("pipeline", tt$n, 200L, NA),
           est_cost_s = val_est_cost("pipeline", NA, tt$n, calib = CALIB))
  }
}

JT <- do.call(rbind, lapply(jobs, function(j)
  data.frame(j[c("job_id", "kind", "tree_id", "cell", "n", "config", "sampler",
                 "N", "init", "rep", "box_scale", "seed", "timeout_s",
                 "est_cost_s")], stringsAsFactors = FALSE)))
if (!is.null(ONLY)) { keep <- JT$kind %in% ONLY; JT <- JT[keep, ]; jobs <- jobs[keep] }
if (!is.null(JOBS_RE)) { keep <- grepl(JOBS_RE, JT$job_id); JT <- JT[keep, ]; jobs <- jobs[keep] }
ord <- order(-JT$est_cost_s)
JT <- JT[ord, ]; jobs <- jobs[ord]
rownames(JT) <- NULL
saveRDS(JT, file.path(RES, "job-table.rds"))

cat(sprintf("  %d jobs; estimated %.0f core-seconds (%.2f h at %d workers)\n",
            nrow(JT), sum(JT$est_cost_s),
            sum(JT$est_cost_s) / 3600 / WORKERS, WORKERS))
tb <- do.call(rbind, lapply(split(JT, JT$kind), function(d)
  data.frame(kind = d$kind[1], jobs = nrow(d), core_s = round(sum(d$est_cost_s)))))
print(tb, row.names = FALSE)

if (DRY) { cat("  --dry-run: job table written, nothing executed.\n"); quit(save = "no") }

# ---------------------------------------------------------------------------
# PROCESS POOL
# ---------------------------------------------------------------------------
done_file <- function(jid) file.path(JOBS, paste0(jid, ".rds"))
pending <- which(!file.exists(vapply(JT$job_id, done_file, "")))
cat(sprintf("  %d jobs already done, %d to run\n", nrow(JT) - length(pending),
            length(pending)))

logf <- file.path(RES, "driver.log")
dlog <- function(...) {
  msg <- sprintf("%s %s", format(Sys.time(), "%H:%M:%S"), sprintf(...))
  cat(msg, "\n", sep = ""); cat(msg, "\n", sep = "", file = logf, append = TRUE)
}
dlog("start tier=%s pending=%d workers=%d build=%s", TIER, length(pending),
     WORKERS, key)

wfun <- function(job, common, worker, lib, trees_file, refs_file, out_file) {
  source(common); source(worker)
  val_run_job(job, lib, trees_file, refs_file, out_file)
}

running <- list()
next_i <- 1L
t_run0 <- Sys.time()
n_ok <- n_to <- n_cr <- n_err <- 0L
last_reported <- -1L

launch <- function(i) {
  j <- jobs[[i]]
  of <- done_file(j$job_id)
  p <- callr::r_bg(wfun,
                   args = list(job = j, common = COMMON, worker = WORKER,
                               lib = LIB, trees_file = trees_file,
                               refs_file = refs_file, out_file = of),
                   supervise = TRUE, stderr = "|", stdout = "|")
  list(proc = p, job = j, out = of, t0 = Sys.time())
}

while (next_i <= length(pending) || length(running)) {
  while (length(running) < WORKERS && next_i <= length(pending)) {
    running[[length(running) + 1L]] <- launch(pending[next_i])
    next_i <- next_i + 1L
  }
  Sys.sleep(0.5)
  keep <- logical(length(running))
  for (k in seq_along(running)) {
    rr <- running[[k]]
    el <- as.numeric(difftime(Sys.time(), rr$t0, units = "secs"))
    alive <- rr$proc$is_alive()
    if (alive && el > rr$job$timeout_s) {
      try(rr$proc$kill(), silent = TRUE)
      saveRDS(c(rr$job, list(outcome = "timeout", error = "driver kill",
                             elapsed = el, build = FP)), rr$out)
      dlog("TIMEOUT %s (%.0fs)", rr$job$job_id, el); n_to <- n_to + 1L
      keep[k] <- FALSE; next
    }
    if (!alive) {
      errtxt <- tryCatch(paste(utils::tail(rr$proc$read_all_error_lines(), 5),
                               collapse = " | "), error = function(e) "")
      if (!file.exists(rr$out)) {
        saveRDS(c(rr$job, list(outcome = "crash", error = errtxt, elapsed = el,
                               build = FP)), rr$out)
        dlog("CRASH   %s (%.0fs) %s", rr$job$job_id, el, substr(errtxt, 1, 160))
        n_cr <- n_cr + 1L
      } else {
        r <- tryCatch(readRDS(rr$out), error = function(e) NULL)
        oc <- r$outcome %||% "error"
        if (identical(oc, "ok")) n_ok <- n_ok + 1L else {
          n_err <- n_err + 1L
          dlog("ERROR   %s: %s", rr$job$job_id, substr(r$error %||% "", 1, 160))
        }
      }
      keep[k] <- FALSE; next
    }
    keep[k] <- TRUE
  }
  running <- running[keep]
  ndone <- n_ok + n_to + n_cr + n_err
  if (ndone > 0 && ndone %% 25 == 0 && ndone != last_reported) {
    last_reported <- ndone
    el <- as.numeric(difftime(Sys.time(), t_run0, units = "secs"))
    dlog("progress %d/%d done in %.0fs; projected total %.0fs",
         ndone, length(pending), el, el / ndone * length(pending))
  }
}

el <- as.numeric(difftime(Sys.time(), t_run0, units = "secs"))
dlog("done: ok=%d error=%d timeout=%d crash=%d in %.0fs wall", n_ok, n_err,
     n_to, n_cr, el)
cat(sprintf("[03-fit] %d job files in %s\n",
            length(list.files(JOBS, pattern = "\\.rds$")), JOBS))
