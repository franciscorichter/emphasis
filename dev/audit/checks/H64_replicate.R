## H64 replication: vary what the verifier did not — 30-tip tree, dd model,
## exponential link, rho = 0.6 (incomplete-sampling sampler branch), and the
## public estimate_rates() path with num_threads > maxN. Runs in a child under
## `perl alarm` so a hang is reported rather than suffered.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
args <- commandArgs(trailingOnly = TRUE)
hc <- parallel::detectCores()

if (length(args) && args[1] == "child") {
  set.seed(6464)
  tr <- ape::rphylo(30, birth = 0.6, death = 0.2)
  brts <- sort(as.numeric(ape::branching.times(tr)), decreasing = TRUE)
  cat("hc =", hc, " tips =", length(brts) + 1, "\n")

  aug <- function(ss, maxN, nt, pars, model = c(0L,0L,0L), link = 0L, rho = 1)
    tryCatch(emphasis:::augment_trees(brts, pars, as.integer(ss), as.integer(maxN),
                                      10000L, 1e6, as.integer(nt), model, link, rho),
             error = function(e) structure(list(msg = conditionMessage(e)), class = "augerr"))

  ## configurations: (name, pars, model, link, rho)
  cfgs <- list(
    cr_mu0     = list(pars = c(0.6,0,0,0, 0.0,0,0,0), model = c(0L,0L,0L), link = 0L, rho = 1),
    dd_lin     = list(pars = c(0.8,-0.01,0,0, 0.1,0,0,0), model = c(1L,0L,0L), link = 0L, rho = 1),
    dd_explink = list(pars = c(log(0.8),-0.01,0,0, log(0.1),0,0,0), model = c(1L,0L,0L), link = 1L, rho = 1),
    cr_rho06   = list(pars = c(0.6,0,0,0, 0.1,0,0,0), model = c(0L,0L,0L), link = 0L, rho = 0.6)
  )

  ## 1. grainsize-0 calls return, and index completeness for mu=0 (ss == maxN == k)
  for (k in c(1L, 2L, 4L, 7L, hc - 1L)) {
    res <- replicate(20, { r <- aug(k, k, hc, cfgs$cr_mu0$pars); if (inherits(r, "augerr")) -1L else length(r$logf) })
    cat(sprintf("1 cr mu=0 k=%d nt=%d grain=%d: ntrees==k in %d/20, errors=%d\n", k, hc, k %/% hc, sum(res == k), sum(res < 0)))
  }

  ## 2. per-configuration: acceptance rate and IS log-weight distribution, nt=1 (grain maxN) vs nt=hc (grain 0)
  for (nm in names(cfgs)) {
    cf <- cfgs[[nm]]
    out <- lapply(c(1L, hc), function(nt) {
      t0 <- Sys.time()
      rs <- replicate(120, aug(2L, 6L, nt, cf$pars, cf$model, cf$link, cf$rho), simplify = FALSE)
      el <- as.numeric(Sys.time() - t0, units = "secs")
      err <- sum(vapply(rs, inherits, TRUE, "augerr"))
      rs <- rs[!vapply(rs, inherits, TRUE, "augerr")]
      nt_ <- vapply(rs, function(r) length(r$logf), 1L)
      rej <- vapply(rs, function(r) r$rejected_zero_weights + r$rejected_overruns + r$rejected_lambda + r$rejected, 1)
      lw  <- unlist(lapply(rs, function(r) r$logf - r$logg))
      list(err = err, ntrees = nt_, rej = rej, lw = lw, el = el)
    })
    p <- tryCatch(t.test(out[[1]]$lw, out[[2]]$lw)$p.value, error = function(e) NA_real_)
    cat(sprintf("2 %-10s nt=1 : err=%d ntrees=%s rej/call=%.2f lw mean=%.3f sd=%.3f n=%d (%.1fs)\n", nm, out[[1]]$err,
                paste(names(table(out[[1]]$ntrees)), table(out[[1]]$ntrees), sep=":", collapse=" "),
                mean(out[[1]]$rej), mean(out[[1]]$lw), sd(out[[1]]$lw), length(out[[1]]$lw), out[[1]]$el))
    cat(sprintf("2 %-10s nt=%d: err=%d ntrees=%s rej/call=%.2f lw mean=%.3f sd=%.3f n=%d (%.1fs)  t.test p=%.3f\n", nm, hc, out[[2]]$err,
                paste(names(table(out[[2]]$ntrees)), table(out[[2]]$ntrees), sep=":", collapse=" "),
                mean(out[[2]]$rej), mean(out[[2]]$lw), sd(out[[2]]$lw), length(out[[2]]$lw), out[[2]]$el, p))
  }

  ## 3. public path: estimate_rates with num_threads = hc > maxN
  run_pub <- function(method, ctrl) {
    t0 <- Sys.time()
    r <- tryCatch(estimate_rates(tr, method = method, model = "cr", control = ctrl),
                  error = function(e) structure(list(msg = conditionMessage(e)), class = "puberr"))
    el <- as.numeric(Sys.time() - t0, units = "secs")
    if (inherits(r, "puberr")) cat(sprintf("3 %s: ERROR %s (%.1fs)\n", method, r$msg, el))
    else cat(sprintf("3 %s: pars=%s loglik=%.3f (%.1fs)\n", method, paste(round(r$pars, 3), collapse=","), r$loglik, el))
  }
  run_pub("cem",  list(lower_bound = c(0.1, 0), upper_bound = c(1.5, 0.5), num_threads = hc, maxN = 4L,
                       max_iter = 3L, num_particles = 20L, max_time = 60))
  run_pub("mcem", list(lower_bound = c(0.1, 0), upper_bound = c(1.5, 0.5), num_threads = hc, maxN = 5L,
                       sample_size = 3L, sampling = "dynamic_fresh", max_iter = 3L, max_time = 60))
  quit(status = 0)
}

cmd <- sprintf("perl -e 'alarm 280; exec @ARGV' Rscript %s child 2>&1; echo EXIT=$?",
               shQuote(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))))
cat(system(cmd, intern = TRUE), sep = "\n")
