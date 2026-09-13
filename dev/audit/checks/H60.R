## H60: DD branch `break` on omp<1e-12 / total<1e-15 — how many rejections come from there vs genuine survivors?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape)})
ns <- asNamespace("emphasis")
src <- paste(deparse(get(".bdi_augment_one", ns), width.cutoff = 500L), collapse = "\n")
src <- sub("if \\(omp < 1e-12\\)\\s*break", "if (omp < 1e-12) { .H60_reason <<- c(.H60_reason, if (n_alive > 0L) 'omp_break_alive' else 'omp_break_empty'); break }", src)
src <- sub("if \\(total < 1e-15\\)\\s*break", "if (total < 1e-15) { .H60_reason <<- c(.H60_reason, if (n_alive > 0L) 'total_break_alive' else 'total_break_empty'); break }", src)
src <- sub("if \\(n_alive > 0L\\)\\s*return\\(NULL\\)", "if (n_alive > 0L) { .H60_reason <<- c(.H60_reason, 'reject_survivor'); return(NULL) }", src)
src <- sub("if \\(n_total > max_missing\\)\\s*return\\(NULL\\)", "if (n_total > max_missing) { .H60_reason <<- c(.H60_reason, 'reject_overflow'); return(NULL) }", src)
stopifnot(lengths(regmatches(src, gregexpr(".H60_reason", src, fixed = TRUE))) == 8)
aug_patched <- eval(parse(text = src), envir = ns)

run <- function(tr, pars8, label, n = 300) {
  brts <- sort(ape::branching.times(tr), decreasing = TRUE)
  tp <- brts[1]; bt <- sort(tp - brts[-1])
  sol <- emphasis:::.bdi_iterate(pars8, c(1L,0L,0L), 0L, bt, tp)
  tab <- table(unlist(lapply(seq_len(n), function(i) {
    .H60_reason <<- character(0)
    a <- aug_patched(bt, pars8, c(1L,0L,0L), 0L, tp, sol$p_fun, sol$Nhat_fun, sol$Phat_fun, sol$Ehat_fun, 1e4L)
    r <- .H60_reason
    if (is.null(a)) r[grepl("reject", r)][1] else if (length(r)) paste0("accepted_after_", paste(unique(r), collapse = "+")) else "accepted_clean"
  })))
  cat("==", label, "\n"); print(tab)
  p <- sol$p_fun(seq(0, tp, length.out = 200))
  cat(sprintf("   p(t) range: [%.3g, %.3g]; min(1-p)=%.3g\n", min(p), max(p), min(1 - p)))
}
set.seed(5)
tr20 <- ape::rphylo(20, 0.5, 0.1)
run(tr20, c(0.6, -0.01, 0,0, 0.1, 0,0,0), "mild DD (K=60)")
run(tr20, c(1.0, -0.05, 0,0, 0.2, 0,0,0), "strong DD (K=20)")
run(tr20, c(1.0, -0.10, 0,0, 0.3, 0,0,0), "over-saturated (K=10, lambda clipped to 0)")
run(tr20, c(1.0, -0.10, 0,0, 0.9, 0,0,0), "over-saturated + high mu")
