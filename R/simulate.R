#' Simulate a phylogenetic tree under the general diversification model
#'
#' Unified entry point for forward and conditional tree simulation.
#' When \code{tree = NULL} (default) a complete tree is simulated forward
#' from the crown; when \code{tree} is supplied the observed extant tree is
#' augmented with stochastically drawn extinct lineages (importance-sampling
#' augmentation).
#'
#' @section Model specification:
#' The model is specified by a 3-element binary vector
#' \code{c(use_N, use_P, use_E)}. Named shortcuts:
#' \tabular{ll}{
#'   \code{"cr"} \tab \code{c(0,0,0)} constant rate\cr
#'   \code{"dd"} \tab \code{c(1,0,0)} diversity dependence\cr
#'   \code{"pd"} \tab \code{c(0,1,0)} phylogenetic diversity dependence\cr
#'   \code{"ep"} \tab \code{c(0,0,1)} evolutionary pendant\cr
#' }
#' Mixed models (e.g. \code{c(1,1,0)}) are fully supported.
#'
#' @section Parameter vector:
#' \preformatted{c(beta_0,  [beta_N],  [beta_P],  [beta_E],
#'   gamma_0, [gamma_N], [gamma_P], [gamma_E])}
#' Only active covariates are included; length = \code{2 + 2 * sum(model)}.
#'
#' @section Output — forward simulation:
#' A named list:
#' \describe{
#'   \item{\code{tes}}{Extant-only phylogeny (\code{phylo}).}
#'   \item{\code{tas}}{Full phylogeny with extinct lineages (\code{phylo}).}
#'   \item{\code{L}}{L-table (DDD format).}
#'   \item{\code{brts}}{Branching times of the extant tree.}
#'   \item{\code{status}}{\code{"done"}, \code{"extinct"}, or
#'     \code{"too_large"}.}
#'   \item{\code{model}}{Resolved binary vector.}
#'   \item{\code{pars}}{Parameter vector as supplied.}
#' }
#' When \code{pars} is a matrix, a list of such lists (one per row).
#'
#' @section Output — conditional simulation (augmentation):
#' When \code{n_trees = 1} the above fields are returned, extended with
#' importance-sampling statistics:
#' \describe{
#'   \item{\code{logf}}{Log-likelihood of the augmented tree under the model.}
#'   \item{\code{logg}}{Log sampling probability of the augmentation.}
#'   \item{\code{log_w}}{\code{logf - logg} — importance weight.}
#'   \item{\code{fhat}}{Monte Carlo log-likelihood estimate of the observed
#'     extant tree.}
#' }
#' When \code{n_trees > 1}, a list with:
#' \describe{
#'   \item{\code{trees}}{List of \code{n_trees} augmented-tree results
#'     (each with \code{tes}, \code{tas}, \code{L}, \code{brts},
#'     \code{status}).}
#'   \item{\code{logf}, \code{logg}, \code{log_w}}{Numeric vectors of length
#'     \code{n_trees}.}
#'   \item{\code{fhat}}{Scalar MC log-likelihood estimate.}
#'   \item{\code{tes}, \code{brts}, \code{model}, \code{pars}}{From the
#'     observed input tree.}
#' }
#'
#' @param tree Optional extant tree for conditional simulation. Accepts
#'   \code{NULL} (forward sim), a \code{simulate_tree} result list, a
#'   \code{phylo} object, or a numeric branching-time vector.
#' @param pars Numeric parameter vector or matrix. When a matrix, each row
#'   is one simulation and the function returns a list of results.
#' @param max_t Crown age. Required for forward simulation.
#' @param model String (\code{"cr"}, \code{"dd"}, \code{"pd"}, \code{"ep"})
#'   or length-3 binary integer vector. Default \code{"cr"}.
#' @param max_lin Maximum lineages before declaring the tree too large.
#'   Default \code{1e6}.
#' @param max_tries Maximum retries after extinction or overflow.
#'   Default \code{100}.
#' @param useDDD Convert L-table to \code{phylo} via \pkg{DDD}.
#'   Default \code{TRUE}.
#' @param n_trees Number of augmented trees to draw (conditional simulation
#'   only). Default \code{1L}. When \code{> 1} the output structure changes
#'   to include importance-sampling statistics across all draws.
#' @examples
#' \dontrun{
#' # --- Forward simulation ---
#' tr <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")
#' tr$brts
#'
#' # Batch: 50 CR simulations varying lambda and mu
#' pm <- cbind(beta_0  = runif(50, 0.3, 0.8),
#'             gamma_0 = runif(50, 0.0, 0.2))
#' sims <- simulate_tree(pars = pm, max_t = 5, model = "cr")
#' length(sims)          # 50
#' sims[[1]]$brts
#'
#' # --- Conditional simulation (augmentation + statistics) ---
#' # Single augmented tree with IS statistics
#' aug1 <- simulate_tree(tree = tr, pars = c(0.5, 0.1),
#'                       model = "cr", n_trees = 1L)
#' aug1$tas              # augmented phylo
#' aug1$logf             # log-likelihood of the augmented tree
#' aug1$fhat             # MC log-likelihood of the extant tree
#'
#' # Many augmented trees — E-step
#' augN <- simulate_tree(tree = tr, pars = c(0.5, 0.1),
#'                       model = "cr", n_trees = 200L)
#' augN$fhat             # MC log-likelihood
#' augN$log_w            # importance weights (length 200)
#' augN$trees[[1]]$tas   # first augmented phylo
#' }
#' @export
simulate_tree <- function(tree      = NULL,
                          pars,
                          max_t     = NULL,
                          model     = "cr",
                          max_lin   = 1e6,
                          max_tries = 100,
                          useDDD    = TRUE,
                          n_trees   = 1L) {

  model_bin <- .resolve_model(model)

  # Batch: pars is a matrix — one simulation per row
  if (is.matrix(pars)) {
    if (!is.numeric(pars) || nrow(pars) == 0L)
      stop("'pars' must be a non-empty numeric matrix.")
    return(lapply(seq_len(nrow(pars)), function(i)
      simulate_tree(tree, pars[i, ], max_t, model_bin,
                    max_lin, max_tries, useDDD, n_trees)))
  }

  if (!is.numeric(pars) || length(pars) == 0L)
    stop("'pars' must be a non-empty numeric vector.")

  # Conditional simulation (augmentation)
  if (!is.null(tree)) {
    return(.sim_tree_conditional(tree, pars, model_bin,
                                 as.integer(n_trees), useDDD))
  }

  # Forward simulation
  if (is.null(max_t)) stop("'max_t' is required for forward simulation.")
  if (!is.numeric(max_t) || length(max_t) != 1L || max_t <= 0)
    stop("'max_t' must be a positive number.")
  if (!is.numeric(max_lin) || length(max_lin) != 1L || max_lin <= 0)
    stop("'max_lin' must be a positive number.")

  expected_n <- 2L + 2L * sum(model_bin)
  if (length(pars) != expected_n) stop(.pars_error_msg(model_bin, expected_n))

  raw <- simulate_div_tree_cpp(
    .expand_pars(pars, model_bin), model_bin,
    max_t, as.integer(max_lin), as.integer(max_tries)
  )

  tes <- tas <- brts <- NULL
  if (raw$status == "done" && useDDD) {
    tes  <- tryCatch(DDD::L2phylo(raw$Ltable, dropextinct = TRUE),
                     error = function(e) NULL)
    tas  <- tryCatch(DDD::L2phylo(raw$Ltable, dropextinct = FALSE),
                     error = function(e) NULL)
    brts <- tryCatch(DDD::L2brts(raw$Ltable,  dropextinct = TRUE),
                     error = function(e) NULL)
  }

  list(tes = tes, tas = tas, L = raw$Ltable, brts = brts,
       status = raw$status, model = model_bin, pars = pars)
}


# --------------------------------------------------------------------------- #
#  Internal helpers                                                            #
# --------------------------------------------------------------------------- #

#' @keywords internal
.resolve_model <- function(model) {
  shortcuts <- list(cr = c(0L, 0L, 0L), dd = c(1L, 0L, 0L),
                    pd = c(0L, 1L, 0L), ep = c(0L, 0L, 1L))
  if (is.character(model)) {
    return(shortcuts[[match.arg(model, names(shortcuts))]])
  }
  model <- as.integer(model)
  if (length(model) != 3L || !all(model %in% 0:1)) {
    stop(paste0("'model' must be \"cr\", \"dd\", \"pd\", \"ep\", ",
                "or a length-3 binary integer vector."))
  }
  model
}

# Expand compact pars to full 8-element vector for C++.
# Layout: c(beta_0, beta_N, beta_P, beta_E, gamma_0, gamma_N, gamma_P, gamma_E)
#' @keywords internal
.expand_pars <- function(pars, model_bin) {
  expected_n <- 2L + 2L * sum(model_bin)
  if (length(pars) != expected_n) stop(.pars_error_msg(model_bin, expected_n))
  active <- which(model_bin == 1L)
  n_lam  <- 1L + length(active)
  beta   <- c(pars[1L],         rep(0.0, 3L))
  gamma  <- c(pars[n_lam + 1L], rep(0.0, 3L))
  if (length(active) > 0L) {
    beta [active + 1L] <- pars[2L:n_lam]
    gamma[active + 1L] <- pars[(n_lam + 2L):length(pars)]
  }
  c(beta, gamma)
}

#' @keywords internal
.pars_error_msg <- function(model_bin, expected_n) {
  lam <- paste(c("beta_0",
                 if (model_bin[1]) "beta_N",
                 if (model_bin[2]) "beta_P",
                 if (model_bin[3]) "beta_E"), collapse = ", ")
  mu  <- paste(c("gamma_0",
                 if (model_bin[1]) "gamma_N",
                 if (model_bin[2]) "gamma_P",
                 if (model_bin[3]) "gamma_E"), collapse = ", ")
  sprintf("model = c(%s) requires %d parameters: c(%s, %s)",
          paste(model_bin, collapse = ", "), expected_n, lam, mu)
}


# --------------------------------------------------------------------------- #
#  Conditional simulation (augmentation)                                      #
# --------------------------------------------------------------------------- #

#' @keywords internal
.sim_tree_conditional <- function(tree, pars, model_bin,
                                  n_trees = 1L, useDDD = TRUE) {
  brts  <- .extract_brts(tree)
  max_t <- brts[1L]

  expected_n <- 2L + 2L * sum(model_bin)
  if (length(pars) != expected_n) stop(.pars_error_msg(model_bin, expected_n))

  tes      <- .extract_tes(tree)
  L_extant <- .extract_Ltable(tree)

  aug <- tryCatch(
    .augment_tree_internal(tree, pars = pars, sample_size = n_trees,
                           max_missing = 1e4, max_lambda = 500,
                           num_threads = 1L),
    error = function(e) NULL
  )

  # Failure path — same structure regardless of n_trees
  if (is.null(aug) || length(aug$trees) == 0L) {
    base <- list(tes = tes, brts = brts, status = "failed",
                 model = model_bin, pars = pars,
                 logf = NA_real_, logg = NA_real_,
                 log_w = NA_real_, fhat = NA_real_)
    if (n_trees == 1L) {
      return(c(base, list(tas = NULL, L = NULL)))
    }
    return(c(base, list(trees = list())))
  }

  # Build augmented-tree results (shared for both n_trees = 1 and > 1)
  aug_trees <- lapply(aug$trees, function(df) {
    L   <- .aug_to_Ltable(df, max_t, brts, L_extant)
    tas <- if (!is.null(L) && useDDD) {
      tryCatch(DDD::L2phylo(L, dropextinct = FALSE), error = function(e) NULL)
    }
    list(tes    = tes,
         tas    = tas,
         L      = L,
         brts   = brts,
         status = if (!is.null(tas)) "done" else "failed")
  })

  stats <- list(logf  = aug$logf,
                logg  = aug$logg,
                log_w = aug$weights,
                fhat  = aug$fhat)

  # Single augmented tree: flat structure (backwards-compatible) + stats
  if (n_trees == 1L) {
    tr <- aug_trees[[1L]]
    return(c(tr, list(model = model_bin, pars = pars), stats))
  }

  # Multiple augmented trees: trees list + stats vectors
  c(list(trees = aug_trees,
         tes   = tes,
         brts  = brts,
         model = model_bin,
         pars  = pars),
    stats)
}


# --------------------------------------------------------------------------- #
#  Tree extraction helpers                                                     #
# --------------------------------------------------------------------------- #

#' @keywords internal
.extract_tes <- function(tree) {
  if (inherits(tree, "phylo"))  return(tree)
  if (!is.null(tree$tes))       return(tree$tes)
  if (!is.null(tree$tas))       return(prune_to_extant(tree$tas))
  NULL
}

#' @keywords internal
.extract_Ltable <- function(tree) {
  if (!is.null(tree$L)) return(tree$L)
  phy <- .extract_tes(tree)
  if (inherits(phy, "phylo")) {
    return(tryCatch(DDD::phylo2L(phy), error = function(e) NULL))
  }
  NULL
}

# Merge augmented extinct branches into the extant L-table.
#' @keywords internal
.aug_to_Ltable <- function(df, max_t, brts, L_extant) {
  if (is.null(L_extant)) return(NULL)

  t_ext_tip <- 1e11
  aug <- df[df$t_ext != 0.0 &
              df$parent_id != -1L &
              df$t_ext < t_ext_tip, , drop = FALSE]

  if (nrow(aug) == 0L) return(L_extant)

  id_to_label <- function(pid) {
    if (pid < 0L || (pid + 2L) > length(brts)) return(NA_integer_)
    as.integer(
      L_extant[which.min(abs(L_extant[, 1L] - brts[pid + 2L])), 3L]
    )
  }

  next_lbl   <- as.integer(max(abs(L_extant[, 3L]))) + 1L
  aug        <- aug[order(aug$brts), , drop = FALSE]
  aug_labels <- setNames(integer(0), character(0))
  new_rows   <- matrix(0.0, nrow = nrow(aug), ncol = 4L)

  for (k in seq_len(nrow(aug))) {
    pid   <- aug$parent_id[k]
    p_lbl <- aug_labels[[as.character(pid)]] %||% id_to_label(pid)
    if (is.na(p_lbl) || length(p_lbl) == 0L) p_lbl <- L_extant[1L, 3L]

    lbl      <- if (p_lbl < 0L) -next_lbl else next_lbl
    next_lbl <- next_lbl + 1L

    new_rows[k, ] <- c(max_t - aug$brts[k], p_lbl, lbl,
                       max_t - aug$t_ext[k])
    aug_labels[as.character(aug$id[k])] <- lbl
  }

  rbind(L_extant, new_rows)
}

`%||%` <- function(x, y) if (!is.null(x)) x else y


# --------------------------------------------------------------------------- #
#  Augmentation via mc_loglik                                                  #
# --------------------------------------------------------------------------- #

#' @keywords internal
.augment_tree_internal <- function(tree,
                                   pars,
                                   sample_size = 1L,
                                   max_missing = 1e4,
                                   max_lambda  = 500,
                                   num_threads = 1L) {
  brts <- .extract_brts(tree)
  mc_loglik(
    brts        = brts,
    pars        = as.numeric(pars),
    sample_size = as.integer(sample_size),
    maxN        = max(2000L, 200L * as.integer(sample_size)),
    max_missing = as.integer(max_missing),
    max_lambda  = as.numeric(max_lambda),
    lower_bound = rep(-1e6, length(pars)),
    upper_bound = rep( 1e6, length(pars)),
    xtol_rel    = 1e-3,
    num_threads = as.integer(num_threads)
  )
}
