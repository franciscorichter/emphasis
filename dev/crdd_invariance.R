# dev/crdd_invariance.R — the cr/dd bit-for-bit gate.
#
# cr (model_bin = c(0,0,0)) and dd (c(1,0,0)) read neither M nor D, so no change
# to the D machinery may move their log f or log q by a single bit.  The check is
# a fixed corpus of augmented trees scored with eval_logf, before and against
# after, under identical().
#
# Augmentation is seeded from the clock inside C++ and is not reproducible, so
# the corpus itself is the thing that has to be carried across the two builds:
# capture it once from the build being compared against, then score it again
# from the build under test.
#
#   Rscript dev/crdd_invariance.R capture   [file]   # writes the corpus + scores
#   Rscript dev/crdd_invariance.R compare   [file]   # rescores and diffs
#
# The default file is dev/.crdd_corpus.rds, which is NOT tracked (see dev/README);
# regenerate it against whichever commit you are comparing to.

suppressMessages(library(methods))

args  <- commandArgs(trailingOnly = TRUE)
mode  <- if (length(args) >= 1L) args[[1L]] else "compare"
store <- if (length(args) >= 2L) args[[2L]] else "dev/.crdd_corpus.rds"

suppressMessages(devtools::load_all(".", quiet = TRUE))

# --------------------------------------------------------------------------- #
#  The corpus: augmented trees over several observed trees and several thetas,  #
#  drawn under cr and under dd, with and without the observed topology.         #
# --------------------------------------------------------------------------- #

.corpus <- function() {
  out <- list()
  for (seed in c(3L, 9L, 21L)) {
    set.seed(seed)
    phy  <- ape::rphylo(14L, 0.7, 0.15)
    brts <- emphasis:::.extract_brts(phy)
    for (with_topology in c(TRUE, FALSE)) {
      pts <- if (with_topology) emphasis:::.pts(brts) else numeric(0)
      for (spec in list(list(model = c(0L, 0L, 0L),
                             pars  = c(0.8, 0, 0, 0, 0.25, 0, 0, 0)),
                        list(model = c(1L, 0L, 0L),
                             pars  = c(0.9, -0.02, 0, 0, 0.30, 0, 0, 0)))) {
        a <- augment_trees(as.numeric(brts), spec$pars,
                           12L, 20000L, 400L, 1e6, 1L,
                           model = spec$model, link = 0L, rho = 1,
                           parent_tip_start = pts)
        if (!length(a$trees)) next
        out[[length(out) + 1L]] <- list(seed = seed, topology = with_topology,
                                        model = spec$model, trees = a$trees)
      }
    }
  }
  out
}

# Every cr/dd (model, link, theta) the scorer can be asked for.
.score_grid <- list(
  list(model = c(0L, 0L, 0L), link = 0L, pars = c(0.80,  0.00, 0, 0, 0.25,  0, 0, 0)),
  list(model = c(0L, 0L, 0L), link = 0L, pars = c(1.40,  0.00, 0, 0, 0.05,  0, 0, 0)),
  list(model = c(1L, 0L, 0L), link = 0L, pars = c(0.90, -0.02, 0, 0, 0.30,  0, 0, 0)),
  list(model = c(1L, 0L, 0L), link = 0L, pars = c(1.20, -0.05, 0, 0, 0.10,  0, 0, 0)),
  list(model = c(0L, 0L, 0L), link = 1L, pars = c(-0.20, 0.00, 0, 0, -1.20, 0, 0, 0)),
  list(model = c(1L, 0L, 0L), link = 1L, pars = c(-0.15,-0.01, 0, 0, -1.00, 0, 0, 0)),
  list(model = c(0L, 0L, 0L), link = 2L, pars = c(0.80,  0.00, 0, 0, 0.25,  0, 0, 0)),
  list(model = c(1L, 0L, 0L), link = 2L, pars = c(0.90,  0.05, 0, 0, 0.30,  0, 0, 0))
)

.score <- function(corpus) {
  lapply(corpus, function(block) {
    lapply(.score_grid, function(g) {
      list(rho1 = eval_logf(g$pars, block$trees, model = g$model,
                            link = g$link, rho = 1.0),
           rho8 = eval_logf(g$pars, block$trees, model = g$model,
                            link = g$link, rho = 0.8))
    })
  })
}

if (mode == "capture") {
  corpus <- .corpus()
  n <- sum(vapply(corpus, function(b) length(b$trees), 0L))
  cat(sprintf("corpus: %d blocks, %d augmented trees, %d nodes total\n",
              length(corpus), n,
              sum(unlist(lapply(corpus, function(b)
                vapply(b$trees, nrow, 0L))))))
  saveRDS(list(corpus = corpus, scores = .score(corpus),
               commit = system("git rev-parse HEAD", intern = TRUE)),
          store)
  cat("written: ", store, "\n", sep = "")
} else if (mode == "compare") {
  if (!file.exists(store))
    stop("no baseline at ", store, " — run `Rscript dev/crdd_invariance.R capture` ",
         "against the commit you are comparing to")
  base <- readRDS(store)
  now  <- .score(base$corpus)
  ok   <- identical(base$scores, now)
  n    <- sum(vapply(base$corpus, function(b) length(b$trees), 0L))
  cat(sprintf("baseline commit : %s\n", base$commit))
  cat(sprintf("trees scored    : %d, over %d (model, link, rho) settings\n",
              n, 2L * length(.score_grid)))
  if (ok) {
    cat("identical()     : TRUE — cr and dd are unchanged, bit for bit\n")
  } else {
    d <- 0
    for (i in seq_along(now)) for (j in seq_along(now[[i]]))
      for (k in c("rho1", "rho8"))
        d <- max(d, max(abs(unlist(now[[i]][[j]][[k]]) -
                            unlist(base$scores[[i]][[j]][[k]]))))
    cat(sprintf("identical()     : FALSE — max |difference| = %.3e\n", d))
  }
  quit(status = if (ok) 0L else 1L)
} else {
  stop("mode must be capture or compare")
}
