#' Simulate a phylogenetic tree under the general diversification model
#'
#' Unified wrapper for forward and conditional tree simulation.  When
#' \code{tree = NULL} (default) a complete tree is simulated forward from the
#' crown; when \code{tree} is supplied the observed extant tree is augmented
#' with stochastically drawn extinct lineages.
#'
#' @section Model specification:
#' The diversification model is specified by a 3-element binary integer vector
#' \code{model = c(use_N, use_P, use_E)} that switches covariates on or off:
#' \describe{
#'   \item{\code{use_N}}{N(t): current species richness.}
#'   \item{\code{use_P}}{P(t): total pendant branch-length of alive lineages.}
#'   \item{\code{use_E}}{E(s,t): pendant edge of lineage \eqn{s} (= time since
#'     its last divergence). Makes rates lineage-specific.}
#' }
#' Named string shortcuts are also accepted:
#' \tabular{ll}{
#'   \code{"cr"} \tab \code{c(0,0,0)} constant rate\cr
#'   \code{"dd"} \tab \code{c(1,0,0)} diversity dependence\cr
#'   \code{"pd"} \tab \code{c(0,1,0)} phylogenetic diversity dependence\cr
#'   \code{"ep"} \tab \code{c(0,0,1)} evolutionary pendant\cr
#' }
#' Mixed models (e.g., \code{c(1,1,0)}) are fully supported.
#'
#' @section Parameter vector:
#' \code{pars} always follows the order:
#' \preformatted{c(beta_0, [beta_N], [beta_P], [beta_E],
#'   gamma_0, [gamma_N], [gamma_P], [gamma_E])}
#' Only the parameters for active covariates (1s in \code{model}) are included,
#' giving a vector of length \code{2 + 2 * sum(model)}.  Lambda parameters
#' come first, then mu parameters, both in N-P-E order.
#'
#' The rates are (identity link, clamped to 0):
#' \deqn{\lambda(s,t) = \max(0,\; \beta_0 + \beta_N N + \beta_P P + \beta_E E_s)}
#' \deqn{\mu(s,t)     = \max(0,\; \gamma_0 + \gamma_N N + \gamma_P P + \gamma_E E_s)}
#'
#' @param tree Optional extant tree for conditional simulation.  Accepts
#'   \code{NULL} (forward sim), a \code{\link{simulate_tree}} result list, a
#'   \code{phylo} object, or a numeric vector of branching times (crown age
#'   first, sorted decreasing).
#' @param pars Numeric parameter vector; see \emph{Parameter vector} above.
#' @param max_t Crown age.  Required for forward simulation (\code{tree = NULL});
#'   ignored otherwise.
#' @param model Model specification: a string (\code{"cr"}, \code{"dd"},
#'   \code{"pd"}, \code{"ep"}) or a length-3 binary integer vector
#'   \code{c(use_N, use_P, use_E)}.  Default \code{"cr"}.
#' @param max_lin Maximum number of lineages before the simulation is
#'   considered too large.  Default \code{1e6}.
#' @param max_tries Maximum retries after extinction or overflow.  Default
#'   \code{100}.
#' @param useDDD Logical.  Convert the L-table to \code{phylo} objects via
#'   \pkg{DDD}.  Default \code{TRUE}.
#' @return A named list:
#' \describe{
#'   \item{\code{tes}}{Extant-only phylogeny (\code{phylo}).}
#'   \item{\code{tas}}{Full phylogeny with extinct lineages (\code{phylo} or
#'     \code{NULL}).}
#'   \item{\code{L}}{L-table matrix (DDD format).}
#'   \item{\code{brts}}{Branching times of the extant tree.}
#'   \item{\code{status}}{\code{"done"} or \code{"failed"}.}
#'   \item{\code{model}}{The resolved model vector \code{c(use_N, use_P, use_E)}.}
#'   \item{\code{pars}}{The parameter vector supplied.}
#' }
#' @examples
#' \dontrun{
#' # Constant-rate birth-death (cr = c(0,0,0))
#' tr <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")
#'
#' # Diversity-dependent (dd = c(1,0,0))
#' tr <- simulate_tree(pars = c(0.8, -0.02, 0.2, -0.01), max_t = 5, model = "dd")
#'
#' # Mixed N+P model  c(1,1,0): 6 parameters
#' tr <- simulate_tree(pars = c(0.5, -0.01, 0.005, 0.15, 0.002, -0.001),
#'                     max_t = 5, model = c(1, 1, 0))
#'
#' # Conditional simulation — augment extant tree with extinct lineages
#' aug <- simulate_tree(tree = tr, pars = c(0.5, -0.01, 0.005, 0.15, 0.002, -0.001),
#'                      model = c(1, 1, 0))
#' aug$tas
#' }
#' @export
simulate_tree <- function(tree      = NULL,
                          pars,
                          max_t     = NULL,
                          model     = "cr",
                          max_lin   = 1e6,
                          max_tries = 100,
                          useDDD    = TRUE) {

  model_bin <- .resolve_model(model)

  if (!is.null(tree)) return(.sim_tree_conditional(tree, pars, model_bin))

  if (is.null(max_t)) stop("'max_t' is required for forward simulation (tree = NULL).")

  # Validate parameter length: 2 + 2 * sum(model_bin) active parameters
  expected_n <- 2L + 2L * sum(model_bin)
  if (length(pars) != expected_n) {
    stop(sprintf(
      paste0("model = c(%s) requires %d parameters:\n",
             "  c(beta_0%s, gamma_0%s)\n",
             "  (lambda params first, then mu params, both in N-P-E order)"),
      paste(model_bin, collapse = ", "), expected_n,
      if (model_bin[1]) ", beta_N"  else "",
      if (model_bin[2]) ", beta_P"  else "",
      if (model_bin[3]) ", beta_E"  else "",
      if (model_bin[1]) ", gamma_N" else "",
      if (model_bin[2]) ", gamma_P" else "",
      if (model_bin[3]) ", gamma_E" else ""
    ))
  }

  full_pars <- .expand_pars(pars, model_bin)
  raw <- simulate_div_tree_cpp(full_pars, model_bin, max_t,
                                as.integer(max_lin), as.integer(max_tries))

  tes <- tas <- brts <- NULL
  if (raw$status == "done" && useDDD) {
    tes  <- DDD::L2phylo(raw$Ltable, dropextinct = TRUE)
    tas  <- DDD::L2phylo(raw$Ltable, dropextinct = FALSE)
    brts <- DDD::L2brts(raw$Ltable,  dropextinct = TRUE)
  }

  list(tes    = tes,
       tas    = tas,
       L      = raw$Ltable,
       brts   = brts,
       status = raw$status,
       model  = model_bin,
       pars   = pars)
}


# --------------------------------------------------------------------------- #
#  Internal helpers                                                            #
# --------------------------------------------------------------------------- #

# Parse model argument → integer vector c(use_N, use_P, use_E)
#' @keywords internal
.resolve_model <- function(model) {
  shortcuts <- list(cr = c(0L, 0L, 0L),
                    dd = c(1L, 0L, 0L),
                    pd = c(0L, 1L, 0L),
                    ep = c(0L, 0L, 1L))
  if (is.character(model)) {
    model <- match.arg(model, names(shortcuts))
    return(shortcuts[[model]])
  }
  model <- as.integer(model)
  if (length(model) != 3L || !all(model %in% 0:1))
    stop("'model' must be a string (\"cr\",\"dd\",\"pd\",\"ep\") or a length-3 binary integer vector.")
  model
}

# Expand compact user pars → full 8-element vector for C++
# Full layout: c(beta_0, beta_N, beta_P, beta_E, gamma_0, gamma_N, gamma_P, gamma_E)
#' @keywords internal
.expand_pars <- function(pars, model_bin) {
  active <- which(model_bin == 1L)   # indices 1..3 → N, P, E slots
  n_lam  <- 1L + length(active)
  beta   <- c(pars[1L], rep(0, 3))
  gamma  <- c(pars[n_lam + 1L], rep(0, 3))
  if (length(active) > 0L) {
    beta [active + 1L] <- pars[2L:n_lam]
    gamma[active + 1L] <- pars[(n_lam + 2L):length(pars)]
  }
  c(beta, gamma)
}


# --------------------------------------------------------------------------- #
#  Conditional simulation helpers                                              #
# --------------------------------------------------------------------------- #

#' @keywords internal
.sim_tree_conditional <- function(tree, pars, model_bin) {
  brts  <- .extract_brts(tree)
  max_t <- brts[1]

  # Resolve tes
  tes <- if (inherits(tree, "phylo")) {
    tree
  } else if (is.list(tree) && !inherits(tree, "phylo") && !is.null(tree$tes)) {
    tree$tes
  } else if (is.list(tree) && !inherits(tree, "phylo") && !is.null(tree$tas)) {
    prune_to_extant(tree$tas)
  } else {
    NULL
  }

  L_extant <- if (is.list(tree) && !inherits(tree, "phylo") && !is.null(tree$L)) {
    tree$L
  } else if (inherits(tree, "phylo")) {
    tryCatch(DDD::phylo2L(tree), error = function(e) NULL)
  } else {
    NULL
  }

  aug <- tryCatch(
    .augment_tree_internal(tree,
                           pars        = pars,
                           sample_size = 1L,
                           max_missing = 1e4,
                           max_lambda  = 500,
                           soc         = 2L,
                           num_threads = 1L),
    error = function(e) NULL
  )

  if (is.null(aug) || length(aug$trees) == 0L) {
    return(list(tes = tes, tas = NULL, L = NULL, brts = brts,
                status = "failed", model = model_bin, pars = pars))
  }

  df  <- aug$trees[[1L]]
  L   <- .aug_to_Ltable(df, max_t, brts, L_extant)
  tas <- if (!is.null(L)) {
    tryCatch(DDD::L2phylo(L, dropextinct = FALSE), error = function(e) NULL)
  } else {
    NULL
  }

  list(tes    = tes,
       tas    = tas,
       L      = L,
       brts   = brts,
       status = if (!is.null(tas)) "done" else "failed",
       model  = model_bin,
       pars   = pars)
}


# Convert an augmented-tree data frame (from mc_loglik) into a DDD-compatible
# L-table that includes both the extant tree and the simulated extinct lineages.
#' @keywords internal
.aug_to_Ltable <- function(df, max_t, brts, L_extant) {
  if (is.null(L_extant)) return(NULL)

  t_ext_tip_val <- 1e11

  sp  <- df[df$t_ext != 0.0, , drop = FALSE]
  aug <- sp[sp$parent_id != -1L & sp$t_ext < t_ext_tip_val, , drop = FALSE]

  if (nrow(aug) == 0L) return(L_extant)

  id_to_label <- function(pid) {
    if (pid < 0L || (pid + 2L) > length(brts)) return(NA_integer_)
    bwd   <- brts[pid + 2L]
    row_k <- which.min(abs(L_extant[, 1L] - bwd))
    as.integer(L_extant[row_k, 3L])
  }

  next_lbl   <- as.integer(max(abs(L_extant[, 3L]))) + 1L
  aug        <- aug[order(aug$brts), , drop = FALSE]
  aug_labels <- setNames(integer(0), character(0))
  new_rows   <- matrix(0.0, nrow = nrow(aug), ncol = 4L)

  for (k in seq_len(nrow(aug))) {
    pid     <- aug$parent_id[k]
    pid_key <- as.character(pid)

    p_lbl <- if (pid_key %in% names(aug_labels)) {
      aug_labels[[pid_key]]
    } else {
      id_to_label(pid)
    }
    if (is.na(p_lbl) || length(p_lbl) == 0L)
      p_lbl <- as.integer(L_extant[1L, 3L])

    lbl      <- if (p_lbl < 0L) -next_lbl else next_lbl
    next_lbl <- next_lbl + 1L

    new_rows[k, 1L] <- max_t - aug$brts[k]
    new_rows[k, 2L] <- as.double(p_lbl)
    new_rows[k, 3L] <- as.double(lbl)
    new_rows[k, 4L] <- max_t - aug$t_ext[k]

    aug_labels[as.character(aug$id[k])] <- lbl
  }

  rbind(L_extant, new_rows)
}


#' @keywords internal
.augment_tree_internal <- function(tree,
                                   pars,
                                   sample_size = 1L,
                                   max_missing = 1e4,
                                   max_lambda  = 500,
                                   soc         = 2L,
                                   num_threads = 1L) {
  brts  <- .extract_brts(tree)
  max_n <- max(100L, as.integer(10L * sample_size))
  n     <- length(pars)
  mc_loglik(brts        = brts,
            pars        = as.numeric(pars),
            sample_size = as.integer(sample_size),
            maxN        = max_n,
            soc         = as.integer(soc),
            max_missing = as.integer(max_missing),
            max_lambda  = as.numeric(max_lambda),
            lower_bound = rep(-1e6, n),
            upper_bound = rep(1e6, n),
            xtol_rel    = 1e-3,
            num_threads = as.integer(num_threads))
}


#' @keywords internal
calc_p <- function(l_table, t) {
  l_table[, 1] <- t - l_table[, 1]
  return(sum(DDD::L2phylo(l_table, dropextinct = TRUE)$edge.length))
}
