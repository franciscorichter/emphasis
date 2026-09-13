# H46 replication: vary model (dd), rho (0.8), max_t (10), link (exponential)
# and check (a) every stamped row gives exactly one negative pendant edge in tas,
# (b) tes untouched, (c) fix L[drop,4] <- 0 yields non-negative edges.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape); library(DDD)})
set.seed(4646)

run <- function(pars, model, rho, max_t, link, R = 20) {
  out <- data.frame(n_drop = NA_integer_, n_neg = NA_integer_, min_edge = NA_real_,
                    neg_are_dropped_tips = NA, tes_ultra = NA, tes_ntip_ok = NA,
                    fix_min_edge = NA_real_, fix_tes_same = NA)[rep(1, R), ]
  for (r in seq_len(R)) {
    s <- simulate_tree(pars = pars, max_t = max_t, model = model, rho = rho,
                       link = link, max_tries = 100, max_lin = 500)
    if (s$status != "done" || is.null(s$L) || is.null(s$tas)) next
    L <- s$L
    drop <- which(L[, 4] == max_t)
    out$n_drop[r]   <- length(drop)
    out$n_neg[r]    <- sum(s$tas$edge.length < 0)
    out$min_edge[r] <- min(s$tas$edge.length)
    # are the negative edges exactly the pendant edges of the dropped-tip labels?
    neg <- s$tas$edge[s$tas$edge.length < 0, , drop = FALSE]
    tipn <- neg[neg[, 2] <= Ntip(s$tas), 2]
    lbls <- as.integer(sub("^t", "", s$tas$tip.label[tipn]))
    out$neg_are_dropped_tips[r] <- setequal(lbls, abs(L[drop, 3])) && nrow(neg) == length(tipn)
    out$tes_ultra[r]   <- is.ultrametric(s$tes)
    out$tes_ntip_ok[r] <- Ntip(s$tes) == sum(L[, 4] == -1)
    L2 <- L; L2[drop, 4] <- 0
    tas2 <- DDD::L2phylo(L2, dropextinct = FALSE)
    tes2 <- DDD::L2phylo(L2, dropextinct = TRUE)
    out$fix_min_edge[r] <- min(tas2$edge.length)
    out$fix_tes_same[r] <- isTRUE(all.equal(sort(branching.times(tes2)),
                                            sort(branching.times(s$tes))))
  }
  out
}

cat("=== dd model, exponential link, rho = 0.8, max_t = 10 ===\n")
a <- run(pars = c(log(0.8), -0.01, log(0.15), 0), model = "dd", rho = 0.8, max_t = 10,
         link = "exponential")
print(a)
cat("\n=== cr model, linear link, rho = 0.9, max_t = 3 ===\n")
b <- run(pars = c(1.0, 0.2), model = "cr", rho = 0.9, max_t = 3, link = "linear")
print(b)

summ <- function(x, nm) {
  ok <- !is.na(x$n_drop)
  cat(nm, ": sims ok =", sum(ok),
      "| n_drop>0 =", sum(x$n_drop[ok] > 0),
      "| n_neg == n_drop in all =", all(x$n_neg[ok] == x$n_drop[ok]),
      "| neg edges are exactly dropped tips' pendant edges =", all(x$neg_are_dropped_tips[ok]),
      "| tes ultrametric =", all(x$tes_ultra[ok]),
      "| tes ntip ok =", all(x$tes_ntip_ok[ok]),
      "| fix min edge >= 0 =", all(x$fix_min_edge[ok] >= 0),
      "| fix leaves tes brts unchanged =", all(x$fix_tes_same[ok]), "\n")
}
summ(a, "dd/exp")
summ(b, "cr/lin")

# Does a rho<1 sim object, fed straight into emphasis(), see the stamped table?
# inference extracts brts from tes; run a 2-iteration fit and inspect the brts used.
s <- simulate_tree(pars = c(0.6, 0.1), max_t = 6, model = "cr", rho = 0.5, max_tries = 100)
brts_tes <- sort(branching.times(s$tes), decreasing = TRUE)
brts_used <- emphasis:::.extract_brts(s)
cat("\nbrts used by inference equal tes branching times:",
    isTRUE(all.equal(unname(brts_used), unname(brts_tes))), "\n")
cat("crown age of tes:", max(brts_tes), " max_t:", 6, "\n")
