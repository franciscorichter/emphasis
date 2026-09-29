#' The per-lineage design of an augmented tree
#'
#' One row per lineage alive on each segment between consecutive events of an
#' augmented tree: the segment's bounds, the covariates the rates read there
#' (\code{N}, the clade mean pendant age \code{M}, the lineage's pendant start
#' \code{ts} so that \code{D = (t - ts) - M}, and its \code{ED} at the segment
#' start, growing at slope one within it), and whether the lineage's own event
#' ends the segment (\code{event}: 0 none, 1 speciation, 2 extinction).  These
#' are the numbers the likelihood is built on; \code{\link{covariate_path}}
#' runs a selection path over them.
#'
#' @param tree An augmented-tree data frame with the topology (the frames
#'   \code{augment_trees} returns when given parent ids), or a \code{phylo}
#'   object, which is taken as fully observed.
#' @return A data frame with columns \code{seg}, \code{lineage}, \code{id},
#'   \code{event}, \code{t0}, \code{t1}, \code{N}, \code{M}, \code{ts},
#'   \code{ed0}.
#' @export
lineage_table <- function(tree) {
  if (inherits(tree, "phylo")) tree <- .observed_frame(tree)
  lineage_table_cpp(tree)
}

# A phylo as the estimator's own frame with no augmented lineage: a constant
# rate with no extinction inserts nothing, so the node list is the observed one.
.observed_frame <- function(phy) {
  brts <- .extract_brts(phy)
  a <- augment_trees(as.numeric(brts), c(0.3, 0, 0, 0, 0, 0, 0, 0),
                     1L, 500L, 100L, 1e6, 1L,
                     model = c(0L, 0L, 0L), link = 0L, rho = 1,
                     parent_tip_start = .pts(brts), parent_id = .pid(brts))
  a$trees[[1L]]
}

#' Selection along a path over the covariates
#'
#' Instead of comparing nested models one pair at a time, this traces the
#' differential-geometric LARS path (Augugliaro, Mineo & Wit 2013, package
#' \pkg{dglars}) over all the covariates a rate can read --- \code{N},
#' \code{D} and \code{ED} --- on the same augmented trees an M-step maximises
#' over.  Under the exponential link the complete-data log-likelihood is a
#' log-linear point process, which on a fine time grid is a Poisson regression
#' with one row per lineage and grid cell; the speciation and extinction
#' halves separate, so each rate gets its own path.  The path's order of
#' entry says which covariates the data ask for first, and a criterion along
#' it (BIC by default) picks a point.
#'
#' Importance weights are folded in by resampling the trees in proportion to
#' them, so the frame is unweighted; the grid's constant exposure is absorbed
#' by the intercept.  An extinction in the first half of a cell has no row for
#' the dying lineage at that cell's midpoint, so one is added at the event
#' time with the cell's exposure; with the default grid that touches under a
#' percent of the compensator.
#'
#' @param trees A list of augmented-tree data frames with the topology, as
#'   \code{augment_trees} returns them with parent ids, or a single
#'   \code{phylo} taken as fully observed.
#' @param log_weights Optional log importance weights, one per tree; trees are
#'   resampled in proportion to them.
#' @param n_trees How many trees to use (resampled with replacement when
#'   \code{log_weights} is given, the first \code{n_trees} otherwise).
#' @param grid Number of equal cells the crown-to-present interval is cut into.
#' @param covariates Which covariates to offer the path; any of \code{"N"},
#'   \code{"D"}, \code{"ED"}, \code{"EDc"} (ED centred on the alive lineages).
#' @param rate \code{"speciation"}, \code{"extinction"} or both.
#' @param criterion \code{"BIC"} or \code{"AIC"} along the path.
#' @param control Passed to \code{dglars::dglars.fit}.
#' @return A list with one element per rate: the \code{dglars} fit, the entry
#'   order of the covariates (the first step at which each is active), the
#'   coefficients at the criterion's optimum with the intercept corrected for
#'   the cell width, the chosen active set, and the frame's size.
#' @export
covariate_path <- function(trees, log_weights = NULL, n_trees = 40L, grid = 500L,
                           covariates = c("N", "D", "ED"),
                           rate = c("speciation", "extinction"),
                           criterion = c("BIC", "AIC"), control = list()) {
  if (!requireNamespace("dglars", quietly = TRUE))
    stop("covariate_path() needs the dglars package: install.packages(\"dglars\")")
  criterion <- match.arg(criterion)
  rate <- match.arg(rate, several.ok = TRUE)
  covariates <- match.arg(covariates, c("N", "D", "ED", "EDc"), several.ok = TRUE)
  if (inherits(trees, "phylo")) trees <- list(.observed_frame(trees))
  if (!is.list(trees) || !length(trees)) stop("'trees' must be a non-empty list of augmented trees")
  # which trees: resampled by weight, or the first n
  if (!is.null(log_weights)) {
    if (length(log_weights) != length(trees)) stop("one log weight per tree")
    w <- exp(log_weights - max(log_weights)); w <- w / sum(w)
    pick <- sample.int(length(trees), n_trees, replace = TRUE, prob = w)
  } else {
    pick <- seq_len(min(n_trees, length(trees)))
  }
  crown <- max(vapply(trees[pick], function(d) max(d$brts), 0))
  delta <- crown / grid
  tm <- (seq_len(grid) - 0.5) * delta
  frame <- do.call(rbind, lapply(pick, function(i) .grid_rows(lineage_table_cpp(trees[[i]]), tm, delta)))
  X <- cbind(N = frame$N, D = frame$D, ED = frame$ED, EDc = frame$EDc)[, covariates, drop = FALSE]
  out <- list()
  for (r in rate) {
    y <- if (r == "speciation") frame$y_spec else frame$y_ext
    fit <- dglars::dglars.fit(X, y, family = stats::poisson("log"), control = control)
    cf <- as.matrix(fit$beta)                    # (p + 1) x steps, intercept first
    rownames(cf) <- c("(Intercept)", covariates)
    active <- cf[-1L, , drop = FALSE] != 0
    entry <- apply(active, 1L, function(a) if (any(a)) which(a)[1L] else NA_integer_)
    crit <- (if (criterion == "BIC") stats::BIC(fit) else stats::AIC(fit))$val
    best <- which.min(crit)
    coef_best <- cf[, best]
    coef_best[1L] <- coef_best[1L] - log(delta)   # the cell width sits in the intercept
    out[[r]] <- list(fit = fit, entry_order = sort(entry[!is.na(entry)]),
                     entry_step = entry, step = best, criterion = criterion,
                     coef = coef_best, active = names(entry)[active[, best]],
                     n_rows = nrow(frame), n_events = sum(y), n_trees = length(pick), grid = grid)
  }
  out
}

# The grid rows of one tree from its lineage table: every (lineage, cell)
# whose midpoint the lineage is alive at, covariates at the midpoint, and the
# events assigned to the cell that contains them.
.grid_rows <- function(tab, tm, delta) {
  if (!nrow(tab)) return(NULL)
  segs <- unique(tab$seg)
  rows <- lapply(segs, function(s) {
    r <- tab[tab$seg == s, ]
    j <- which(tm >= r$t0[1L] & tm < r$t1[1L])
    if (!length(j)) return(NULL)
    n <- nrow(r)
    e <- r[rep(seq_len(n), times = length(j)), ]
    e$cell <- rep(j, each = n)
    e$tm <- tm[e$cell]
    e
  })
  e <- do.call(rbind, rows)
  # events: the row of the event lineage in the cell containing t1
  ev <- tab[tab$event != 0L, ]
  e$y_spec <- 0L; e$y_ext <- 0L
  if (nrow(ev)) {
    ev$cell <- pmin(ceiling(ev$t1 / delta), length(tm))
    key_e <- paste(e$lineage, e$cell)
    key_v <- paste(ev$lineage, ev$cell)
    hit <- match(key_v, key_e)
    for (q in seq_len(nrow(ev))) {
      if (is.na(hit[q])) {
        # no row at that cell's midpoint (the lineage died in the first half
        # of the cell, or the segment holds no midpoint): add one at the event
        add <- ev[q, ]; add$tm <- add$t1; add$y_spec <- 0L; add$y_ext <- 0L
        e <- rbind(e, add); hit[q] <- nrow(e)
      }
      if (ev$event[q] == 1L) e$y_spec[hit[q]] <- 1L else e$y_ext[hit[q]] <- 1L
    }
  }
  e$D  <- (e$tm - e$ts) - e$M
  e$ED <- e$ed0 + (e$tm - e$t0)
  # ED centred on the lineages alive at the same cell
  e$EDc <- e$ED - stats::ave(e$ED, e$cell, FUN = mean)
  e[, c("cell", "tm", "N", "D", "ED", "EDc", "y_spec", "y_ext")]
}
