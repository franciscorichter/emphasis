#!/usr/bin/env Rscript
# H63 — is `num_threads` honoured by the TBB code in E_step / augment_trees /
# M_step?  The hypothesis: `tbb::task_arena arena(num_threads)` is built but the
# parallel_for runs outside `arena.execute`, so it runs in the default arena
# (hardware concurrency); `num_threads` only sets the grainsize.
#
# Measure: process CPU time (user+sys) / wall == core-equivalents.  The R-side
# `unpack` of the result runs single-threaded on the main thread after the
# parallel region, so for the E-step the parallel region's core usage is
#     cores_cpp = (cpu_total - unpack_wall) / cpp_wall,
# with cpp_wall = the E-step's own `time` (ms) and unpack_wall = wall - cpp_wall.
# Also reported: OS threads alive in the process (ps -M) and the number of
# trees returned (an overshoot above sample_size is the H62 race).
#
# The tree is large (150 tips, ~400 extinct lineages per augmentation, ~22 ms
# per tree at 1 thread) so that the C++ phase dominates.  maxN is kept small
# (5 x sample_size) so that an H62 overshoot cannot explode the unpack.
#
# Modes (the script re-invokes itself with Rscript so that "fresh process"
# cases are really fresh):
#   estep        E-step at num_threads = 1, 2, 3, 5, 10 in ONE process
#   gen <file>   generate an E-step list (em_cpp copy_trees = TRUE) -> RDS
#   mstep_fresh <file>   m_cpp(num_threads = 1) as the FIRST TBB call in the process
#   mstep_after <file>   E-step (nt = 1) first, then m_cpp(num_threads = 1 / 10)
#   fork         parent warms the pool at nt = 10, then mclapply children at nt = 1 and 4
#   (no arg)     orchestrate all of the above
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "all"
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))

hw <- parallel::detectCores(logical = TRUE)

# ---- fixed tree / model ----------------------------------------------------
set.seed(63)
tr    <- ape::rphylo(150, 0.5, 0.3)
brts  <- sort(as.numeric(ape::branching.times(tr)), decreasing = TRUE)
pars8 <- c(0.5, 0, 0, 0, 0.3, 0, 0, 0)          # cr: beta0, ..., gamma0, ...
model <- c(0L, 0L, 0L)
N     <- 150L
maxN  <- 5L * N

nthreads_alive <- function() {
  out <- tryCatch(system2("ps", c("-M", Sys.getpid()), stdout = TRUE, stderr = FALSE),
                  error = function(e) character())
  max(0L, length(out) - 1L)
}

measure <- function(expr) {
  t0 <- proc.time()
  val <- force(expr)
  dt <- proc.time() - t0
  cpu <- unname(dt["user.self"] + dt["sys.self"])
  wall <- unname(dt["elapsed"])
  list(val = val, cpu = cpu, wall = wall, cores = cpu / max(wall, 1e-9))
}

estep <- function(nt) {
  emphasis:::augment_trees(brts = brts, pars = pars8, sample_size = N, maxN = maxN,
                           max_missing = 2000L, max_lambda = 1e4,
                           num_threads = as.integer(nt),
                           model = model, link = 0L, rho = 1.0)
}

fmt_e <- function(tag, nt, m) {
  r <- m$val
  cpp <- r$time / 1000
  unpack <- max(m$wall - cpp, 0)
  cores_cpp <- (m$cpu - unpack) / cpp
  cat(sprintf("%-22s nt=%2d  cpu=%5.2fs wall=%5.2fs  cpp=%5.2fs unpack=%4.2fs  cores_cpp=%5.2f  n_trees=%4d (N=%d)  threads_alive=%2d\n",
              tag, nt, m$cpu, m$wall, cpp, unpack, cores_cpp, length(r$logf), N, nthreads_alive()))
}
fmt <- function(tag, nt, m) {
  cat(sprintf("%-22s nt=%2d  cpu=%5.2fs wall=%5.2fs  cores=%5.2f  threads_alive=%2d\n",
              tag, nt, m$cpu, m$wall, m$cores, nthreads_alive()))
}

# ---- modes ------------------------------------------------------------------
if (mode == "estep") {
  cat(sprintf("hardware threads = %d ; threads_alive before any TBB call = %d\n", hw, nthreads_alive()))
  for (nt in c(1L, 2L, 3L, 5L, 10L, 1L)) fmt_e("E-step augment_trees", nt, measure(estep(nt)))
}

if (mode == "gen") {
  # E-step via em_cpp so the returned list carries `weights` and `fhat`
  # (what m_cpp expects).  The cr model is N-only, so the 3-column tree
  # layout is enough for the M-step (see H69).
  e <- emphasis:::em_cpp(brts = brts, init_pars = pars8, sample_size = N, maxN = maxN,
                         max_missing = 2000L, max_lambda = 1e4,
                         lower_bound = rep(-1e6, 8), upper_bound = rep(1e6, 8), xtol_rel = 1e-3,
                         num_threads = 1L, copy_trees = TRUE, model = model, link = 0L, rho = 1.0)
  saveRDS(e, args[2])
  cat(sprintf("generated E-step: %d trees, %.0f nodes/tree\n", length(e$trees),
              mean(sapply(e$trees, nrow))))
}

mstep_run <- function(e, nt, reps = 5L) {
  lb <- c(1e-4, 0, 0, 0, 1e-4, 0, 0, 0); ub <- c(1e6, 0, 0, 0, 1e6, 0, 0, 0)   # inactive slots pinned
  measure(for (r in seq_len(reps))
    emphasis:::m_cpp(e_step = e, init_pars = c(0.4, 0, 0, 0, 0.2, 0, 0, 0), plugin = "",
                     lower_bound = lb, upper_bound = ub, xtol_rel = 1e-4,
                     num_threads = as.integer(nt), model = model, link = 0L, rho = 1.0))
}

if (mode == "mstep_fresh") {
  e <- readRDS(args[2])
  cat(sprintf("threads_alive before any TBB call = %d\n", nthreads_alive()))
  fmt("M-step first-in-process", 1L, mstep_run(e, 1L))
  fmt("M-step again", 1L, mstep_run(e, 1L))
  fmt("M-step", 10L, mstep_run(e, 10L))
  fmt("M-step", 1L, mstep_run(e, 1L))
}

if (mode == "mstep_after") {
  e <- readRDS(args[2])
  fmt_e("E-step (warm-up)", 1L, measure(estep(1L)))
  fmt("M-step after E-step", 1L, mstep_run(e, 1L))
  fmt("M-step after E-step", 10L, mstep_run(e, 10L))
}

if (mode == "fork") {
  fmt_e("parent E-step (warm-up)", 10L, measure(estep(10L)))
  child <- function(i, nt) {
    t_alive <- nthreads_alive()
    m <- measure(estep(nt))
    cpp <- m$val$time / 1000; unpack <- max(m$wall - cpp, 0)
    sprintf("child %d nt=%d: threads_alive=%d cpu=%.2fs wall=%.2fs cpp=%.2fs cores_cpp=%.2f n_trees=%d",
            i, nt, t_alive, m$cpu, m$wall, cpp, (m$cpu - unpack) / cpp, length(m$val$logf))
  }
  for (nt in c(1L, 4L)) {
    tw <- proc.time()
    res <- parallel::mclapply(1:4, child, nt = nt, mc.cores = 4L)
    tw <- (proc.time() - tw)["elapsed"]
    ok <- !vapply(res, inherits, logical(1), "try-error")
    cat(sprintf("mclapply(4 children, cpp nt=%d): wall=%.2fs, all returned=%s\n", nt, tw, all(ok)))
    for (r in res) cat("   ", if (inherits(r, "try-error")) paste("ERROR:", r) else r, "\n")
  }
}

if (mode == "all") {
  run <- function(...) {
    cat("\n$ Rscript H63.R", paste(c(...)), "\n")
    out <- system2("Rscript", c(shQuote(self), ...), stdout = TRUE, stderr = TRUE)
    cat(paste(grep("built under R version|^Warning message:$", out, value = TRUE, invert = TRUE),
              collapse = "\n"), "\n")
  }
  f <- tempfile(fileext = ".rds")
  run("estep")
  run("gen", f)
  run("mstep_fresh", f)
  run("mstep_after", f)
  run("fork")
  for (nt in c(1L, 3L, 6L)) run("threads", nt)
}

# ---- threads <nt>: load-independent count of threads that executed chunks ----
# Fresh process, ONE E-step at num_threads = nt, then per-thread cumulative
# user CPU from `ps -M`.  Worker threads exist only because of this one
# parallel_for, so #threads with utime > 0.2 s == #threads that ran chunks.
# If the arena were honoured this is <= nt; with the range split into
# 2^ceil(log2(nt)) leaves in the default arena it can exceed nt.
if (mode == "threads") {
  nt <- as.integer(args[2])
  m <- measure(estep(nt))
  out <- system2("ps", c("-M", Sys.getpid()), stdout = TRUE, stderr = FALSE)
  tm <- regmatches(out, gregexpr("[0-9]+:[0-9]{2}\\.[0-9]{2}", out))
  utime <- vapply(tm, function(v) { if (length(v) < 2) return(NA_real_)
    p <- as.numeric(strsplit(v[2], "[:]")[[1]]); p[1] * 60 + p[2] }, numeric(1))
  utime <- utime[!is.na(utime)]
  busy <- sum(utime > 0.2)
  cat(sprintf("threads mode nt=%2d: cpp=%.2fs n_trees=%d  threads_alive=%d  threads_with_utime>0.2s=%d  [utimes: %s]\n",
              nt, m$val$time / 1000, length(m$val$logf), length(utime), busy,
              paste(sprintf("%.2f", sort(utime, decreasing = TRUE)), collapse = " ")))
}
