#!/usr/bin/env Rscript
# H63 replication (independent). Varies what the verifier did not: a different
# tree (100 tips, rphylo(100, 0.8, 0.5)), nt = 4 (a power of two: prediction
# exactly 4 busy threads) vs nt = 5 (prediction 8), em_cpp(num_threads = 1) as
# a whole (does its M-step spawn a pool?), and forks from a COLD parent (no
# TBB use) and a parent warmed at nt = 1 (market exists, workers never launched).
# Modes: threads <nt> | em1 | fork_cold | fork_warm1 | fork_warm10 | all
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "all"
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))

set.seed(7)
tr    <- ape::rphylo(100, 0.8, 0.5)
brts  <- sort(as.numeric(ape::branching.times(tr)), decreasing = TRUE)
pars8 <- c(0.8, 0, 0, 0, 0.5, 0, 0, 0)
model <- c(0L, 0L, 0L)
N     <- 200L
maxN  <- 5L * N

ps_utimes <- function() {
  out <- system2("ps", c("-M", Sys.getpid()), stdout = TRUE, stderr = FALSE)
  tm <- regmatches(out, gregexpr("[0-9]+:[0-9]{2}\\.[0-9]{2}", out))
  u <- vapply(tm, function(v) { if (length(v) < 2) return(NA_real_)
    p <- as.numeric(strsplit(v[2], "[:]")[[1]]); p[1] * 60 + p[2] }, numeric(1))
  u[!is.na(u)]
}
measure <- function(expr) {
  t0 <- proc.time(); val <- force(expr); dt <- proc.time() - t0
  list(val = val, cpu = unname(dt["user.self"] + dt["sys.self"]), wall = unname(dt["elapsed"]))
}
estep <- function(nt, n = N, mx = maxN) {
  emphasis:::augment_trees(brts = brts, pars = pars8, sample_size = n, maxN = mx,
                           max_missing = 2000L, max_lambda = 1e4, num_threads = as.integer(nt),
                           model = model, link = 0L, rho = 1.0)
}
report <- function(tag, nt, m) { force(m)
  u <- ps_utimes(); cpp <- m$val$time / 1000; unpack <- max(m$wall - cpp, 0)
  cat(sprintf("%-28s nt=%2d cpp=%5.2fs cpu=%5.2fs wall=%5.2fs cores_cpp=%4.2f n_trees=%4d threads_alive=%2d busy(>0.2s)=%d [%s]\n",
              tag, nt, cpp, m$cpu, m$wall, (m$cpu - unpack) / cpp, length(m$val$logf),
              length(u), sum(u > 0.2), paste(sprintf("%.2f", sort(u, decreasing = TRUE)), collapse = " ")))
}

if (mode == "threads") {
  nt <- as.integer(args[2])
  cat(sprintf("fresh process: threads_alive before = %d\n", length(ps_utimes())))
  report("E-step fresh", nt, measure(estep(nt)))
}

if (mode == "em1") {
  cat(sprintf("fresh process: threads_alive before = %d\n", length(ps_utimes())))
  m <- measure(emphasis:::em_cpp(brts = brts, init_pars = pars8, sample_size = 600L, maxN = 3000L,
                                 max_missing = 2000L, max_lambda = 1e4,
                                 lower_bound = c(1e-4,0,0,0,1e-4,0,0,0), upper_bound = c(1e6,0,0,0,1e6,0,0,0),
                                 xtol_rel = 1e-4, num_threads = 1L, copy_trees = FALSE,
                                 model = model, link = 0L, rho = 1.0))
  u <- ps_utimes()
  cat(sprintf("em_cpp(num_threads=1): cpu=%.2fs wall=%.2fs cores=%.2f time(E+M)=%.2fs n_trees=%d threads_alive=%d busy(>0.2s)=%d [%s]\n",
              m$cpu, m$wall, m$cpu / m$wall, m$val$time / 1000, m$val$trees, length(u), sum(u > 0.2),
              paste(sprintf("%.2f", sort(u, decreasing = TRUE)), collapse = " ")))
}

child <- function(i, nt) {
  before <- length(ps_utimes()); m <- measure(estep(nt)); u <- ps_utimes()
  cpp <- m$val$time / 1000; unpack <- max(m$wall - cpp, 0)
  sprintf("child %d nt=%d: threads_before=%d threads_after=%d busy(>0.2s)=%d cpu=%.2fs cpp=%.2fs cores_cpp=%.2f n_trees=%d",
          i, nt, before, length(u), sum(u > 0.2), m$cpu, cpp, (m$cpu - unpack) / cpp, length(m$val$logf))
}
forktest <- function(label) {
  cat(sprintf("parent threads_alive before fork = %d\n", length(ps_utimes())))
  for (nt in c(1L, 4L)) {
    tw <- proc.time()
    res <- parallel::mclapply(1:3, child, nt = nt, mc.cores = 3L)
    tw <- (proc.time() - tw)["elapsed"]
    cat(sprintf("%s: mclapply(3 children, cpp nt=%d) wall=%.2fs\n", label, nt, tw))
    for (r in res) cat("   ", if (inherits(r, "try-error")) paste("ERROR:", r) else r, "\n")
  }
}
if (mode == "fork_cold")   forktest("cold parent (no TBB use)")
if (mode == "fork_warm1")  { report("parent E-step", 1L, measure(estep(1L)));  forktest("parent warmed nt=1") }
if (mode == "fork_warm10") { report("parent E-step", 10L, measure(estep(10L))); forktest("parent warmed nt=10") }

if (mode == "all") {
  run <- function(...) {
    cat("\n$ Rscript H63_replicate.R", paste(c(...)), "\n")
    out <- system2("Rscript", c(shQuote(self), ...), stdout = TRUE, stderr = TRUE)
    cat(paste(grep("built under R version|^Warning message:$", out, value = TRUE, invert = TRUE), collapse = "\n"), "\n")
  }
  for (nt in c(1L, 2L, 4L, 5L, 8L)) run("threads", nt)
  run("em1")
  run("fork_cold"); run("fork_warm1"); run("fork_warm10")
}
