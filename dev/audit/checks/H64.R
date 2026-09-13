## H64: grainsize = maxN / num_threads is 0 when maxN < num_threads.
## Self-contained. Harness mode runs the risky (possibly hanging) calls in a
## child Rscript under `perl alarm` so a TBB hang is reported, not suffered.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
args <- commandArgs(trailingOnly = TRUE)
hc <- parallel::detectCores()   # == std::thread::hardware_concurrency() here

make_brts <- function(n = 12) {
  set.seed(64)
  tr <- ape::rphylo(n, birth = 0.5, death = 0.1)
  sort(as.numeric(ape::branching.times(tr)), decreasing = TRUE)
}
brts <- make_brts()
pars_mu0 <- c(0.5, 0, 0, 0, 0.0, 0, 0, 0)  # lambda=0.5, mu=0: every attempt accepted
pars_mu  <- c(0.5, 0, 0, 0, 0.1, 0, 0, 0)  # lambda=0.5, mu=0.1

aug <- function(ss, maxN, nt, pars = pars_mu0)
  tryCatch(emphasis:::augment_trees(brts, pars, as.integer(ss), as.integer(maxN),
                                    10000L, 1e6, as.integer(nt)),
           error = function(e) structure(list(msg = conditionMessage(e)), class = "augerr"))

if (length(args) && args[1] == "child") {
  cat("hardware_concurrency (detectCores):", hc, "\n")
  ## A. grainsize-0 calls return at all?
  for (maxN in c(1L, 2L, 5L, 9L)) {
    t0 <- Sys.time()
    r <- aug(1L, maxN, hc)
    cat(sprintf("A maxN=%d nt=%d grain=%d -> %s, %.3fs\n", maxN, hc, maxN %/% hc,
                if (inherits(r, "augerr")) paste("ERROR:", r$msg) else paste("ntrees", length(r$logf)),
                as.numeric(Sys.time() - t0, units = "secs")))
  }
  ## B. index completeness: sample_size == maxN == k, mu=0 (no rejection possible),
  ##    grainsize 0 for every k < hc. A skipped index would give < k trees -> error.
  for (nt in c(1L, hc)) for (k in 1:(hc - 1)) {
    res <- replicate(30, { r <- aug(k, k, nt); if (inherits(r, "augerr")) -1L else length(r$logf) })
    cat(sprintf("B nt=%2d k=%d grain=%d : ok=%d/30 (ntrees==k), errors=%d, rej_zero(last)=%s\n",
                nt, k, k %/% nt, sum(res == k), sum(res < 0),
                { r <- aug(k, k, nt); if (inherits(r, "augerr")) "NA" else r$rejected_zero_weights }))
  }
  ## C. distribution of IS log-weights unchanged by grainsize 0 (mu=0.1, ss=3, maxN=9)
  lw <- lapply(c(1L, hc), function(nt) {
    unlist(replicate(150, { r <- aug(3L, 9L, nt, pars_mu); if (inherits(r, "augerr")) NULL else r$logf - r$logg }, simplify = FALSE))
  })
  cat(sprintf("C nt=1 : n=%d mean lw=%.4f sd=%.4f | nt=%d (grain 0): n=%d mean lw=%.4f sd=%.4f | t.test p=%.3f\n",
              length(lw[[1]]), mean(lw[[1]]), sd(lw[[1]]), hc, length(lw[[2]]), mean(lw[[2]]), sd(lw[[2]]),
              t.test(lw[[1]], lw[[2]])$p.value))
  ## E. overshoot of sample_size (H62 race) by configuration: does grainsize 0 change it?
  for (cfg in list(c(9L, 1L), c(9L, hc), c(90L, hc), c(900L, hc))) {
    n <- replicate(300, { r <- aug(3L, cfg[1], cfg[2], pars_mu); if (inherits(r, "augerr")) NA_integer_ else length(r$logf) })
    cat(sprintf("E ss=3 maxN=%d nt=%d grain=%d : table(ntrees)= %s ; errors=%d\n", cfg[1], cfg[2], cfg[1] %/% cfg[2],
                paste(names(table(n)), table(n), sep = ":", collapse = " "), sum(is.na(n))))
  }
  ## D. arithmetic of the defect for the CEM default and whitebox
  for (nt in c(1, 4, 10, 12, 32)) cat(sprintf("D maxN=10 num_threads=%d -> grainsize=%d\n", nt, 10 %/% nt))
  quit(status = 0)
}

## harness: run child under a 240 s alarm
cmd <- sprintf("perl -e 'alarm 240; exec @ARGV' Rscript %s child 2>&1; echo EXIT=$?", shQuote(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))))
out <- system(cmd, intern = TRUE)
cat(out, sep = "\n")
