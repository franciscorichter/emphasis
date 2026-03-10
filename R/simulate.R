#' Simulate a phylogenetic tree under the general diversification model
#'
#' Unified entry point for forward and conditional tree simulation.
#' When \code{tree = NULL} (default) a complete tree is simulated forward
#' from the crown; when \code{tree} is supplied the observed extant tree is
#' augmented with stochastically drawn extinct lineages.
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
#' @section Output -- forward simulation (single):
#' \describe{
#'   \item{\code{tes}}{Extant-only phylogeny (\code{phylo}), or \code{NULL}.}
#'   \item{\code{tas}}{Full phylogeny with extinct lineages (\code{phylo}),
#'     or \code{NULL}.}
#'   \item{\code{L}}{L-table (DDD format), or \code{NULL} if simulation
#'     failed.}
#'   \item{\code{status}}{\code{"done"}, \code{"extinct"}, or
#'     \code{"too_large"}.}
#'   \item{\code{survival_prob}}{Empirical survival probability: \code{1 /
#'     n_attempts} when the simulation succeeded, \code{0} otherwise.}
#' }
#'
#' @section Output -- forward simulation (batch, matrix \code{pars}):
#' \describe{
#'   \item{\code{simulations}}{List of individual results (one per row of
#'     \code{pars}), each with the fields above.}
#'   \item{\code{survival_prob}}{Fraction of simulations with
#'     \code{status == "done"}.}
#' }
#'
#' @section Output -- conditional simulation / augmentation:
#' When \code{n_trees = 1}:
#' \describe{
#'   \item{\code{tas}}{Augmented phylogeny (\code{phylo}) including extinct
#'     lineages, or \code{NULL} on failure or when \code{useDDD = FALSE}.}
#'   \item{\code{log_q}}{Log sampling probability \eqn{\log q(z \mid
#'     \text{obs}, \theta)} of the augmentation.  \code{NA} on failure.}
#' }
#' When \code{n_trees > 1}:
#' \describe{
#'   \item{\code{trees}}{List of \code{tas} phylogenies (length up to
#'     \code{n_trees}; failed draws are \code{NULL}).}
#'   \item{\code{log_q}}{Numeric vector of log sampling probabilities,
#'     one per successful draw.}
#' }
#'
#' @param tree Optional observed extant tree for augmentation. Accepts
#'   \code{NULL} (forward sim), a \code{simulate_tree} result, a
#'   \code{phylo} object, or a numeric branching-time vector.
#' @param pars Numeric parameter vector or matrix. When a matrix each row
#'   produces one forward simulation.
#' @param max_t Crown age (forward simulation only). Default \code{1}.
#' @param model Model specification (string, formula, or binary vector).
#'   Default \code{"cr"}.
#' @param max_lin Maximum lineages before declaring the tree too large.
#'   Default \code{1e6}.
#' @param max_tries Maximum additional attempts after extinction or overflow
#'   (forward simulation only). Retries are tracked at the R level to
#'   compute \code{survival_prob}. Default \code{100}.
#' @param useDDD Convert L-table to \code{phylo} via \pkg{DDD}.
#'   Default \code{TRUE}.
#' @param n_trees Number of augmented trees to draw (augmentation only).
#'   Default \code{1L}.
#' @param link Link function: \code{"linear"} (default) or
#'   \code{"exponential"}.
#' @param max_missing Maximum extinct lineages per augmented tree.
#'   Default \code{1e4}.
#' @param max_lambda Maximum speciation rate (thinning bound) for augmentation.
#'   Default \code{500}.
#' @param maxN Maximum total augmentation attempts. \code{NULL} (default)
#'   sets this to \code{max(2000, 200 * n_trees)}.
#' @param num_threads Threads for parallel augmentation. Default \code{1L}.
#' @examples
#' \dontrun{
#' # --- Forward simulation ---
#' tr <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")
#' tr$status          # "done"
#' tr$survival_prob   # e.g. 0.25 if took 4 attempts
#'
#' # Batch: 200 CR simulations
#' pm <- cbind(beta_0 = runif(200, 0.3, 0.8), gamma_0 = runif(200, 0.05, 0.3))
#' sims <- simulate_tree(pars = pm, max_t = 8, model = "cr")
#' sims$survival_prob              # fraction that succeeded
#' sims$simulations[[1]]$tes      # first extant phylo
#'
#' # --- Augmentation ---
#' aug <- simulate_tree(tree = tr, pars = c(0.5, 0.1), model = "cr")
#' aug$tas     # augmented phylo
#' aug$log_q   # log q(z | obs, theta)
#'
#' # Many augmented trees
#' augN <- simulate_tree(tree = tr, pars = c(0.5, 0.1), model = "cr",
#'                       n_trees = 100L)
#' augN$log_q          # vector of length 100
#' augN$trees[[1]]     # first augmented phylo
#' }
#' @export
simulate_tree <- function(tree        = NULL,
                          pars,
                          max_t       = 1,
                          model       = "cr",
                          max_lin     = 1e6,
                          max_tries   = 1,
                          useDDD      = TRUE,
                          n_trees     = 1L,
                          link        = "linear",
                          max_missing = 1e4,
                          max_lambda  = 500,
                          maxN        = NULL,
                          num_threads = 1L) {

  model_bin <- .resolve_model(model)
  link_int  <- .resolve_link(link)

  if (!is.numeric(pars) || length(pars) == 0L)
    stop("'pars' must be a non-empty numeric vector or matrix.")

  # ---------------------------------------------------------------------- #
  #  Batch forward simulation (matrix pars)                                 #
  # ---------------------------------------------------------------------- #
  if (is.matrix(pars)) {
    if (nrow(pars) == 0L) stop("'pars' matrix must have at least one row.")
    sims <- lapply(seq_len(nrow(pars)), function(i)
      simulate_tree(tree, pars[i, ], max_t, model_bin,
                    max_lin, max_tries, useDDD, n_trees, link_int,
                    max_missing, max_lambda, maxN, num_threads))
    if (is.null(tree)) {
      # Forward batch: aggregate survival_prob
      surv <- mean(vapply(sims, `[[`, 0.0, "survival_prob"))
      return(list(simulations = sims, survival_prob = surv))
    }
    # Conditional batch: just return the list
    return(sims)
  }

  # ---------------------------------------------------------------------- #
  #  Conditional simulation (augmentation)                                  #
  # ---------------------------------------------------------------------- #
  if (!is.null(tree)) {
    return(.sim_tree_conditional(tree, pars, model_bin,
                                 as.integer(n_trees), useDDD, link_int,
                                 max_missing, max_lambda, maxN,
                                 as.integer(num_threads)))
  }

  # ---------------------------------------------------------------------- #
  #  Forward simulation (single)                                            #
  # ---------------------------------------------------------------------- #
  if (!is.numeric(max_t) || length(max_t) != 1L || max_t <= 0)
    stop("'max_t' must be a positive number.")
  if (!is.numeric(max_lin) || length(max_lin) != 1L || max_lin <= 0)
    stop("'max_lin' must be a positive number.")

  expected_n <- 2L + 2L * sum(model_bin)
  if (length(pars) != expected_n) stop(.pars_error_msg(model_bin, expected_n))

  pars8       <- .expand_pars(pars, model_bin)
  max_lin_i   <- as.integer(max_lin)
  max_tries_i <- as.integer(max_tries)

  # Retry loop at R level to count attempts -> survival_prob
  n_attempts <- 0L
  raw        <- list(status = "extinct")
  while (raw$status != "done" && n_attempts <= max_tries_i) {
    raw        <- simulate_div_tree_cpp(pars8, model_bin, max_t, max_lin_i, 0L, link_int)
    n_attempts <- n_attempts + 1L
  }

  survival_prob <- if (raw$status == "done") 1.0 / n_attempts else 0.0

  tes <- tas <- NULL
  if (raw$status == "done" && useDDD) {
    tes <- tryCatch(DDD::L2phylo(raw$Ltable, dropextinct = TRUE),
                    error = function(e) NULL)
    tas <- tryCatch(DDD::L2phylo(raw$Ltable, dropextinct = FALSE),
                    error = function(e) NULL)
  }

  list(tes           = tes,
       tas           = tas,
       L             = if (raw$status == "done") raw$Ltable else NULL,
       status        = raw$status,
       survival_prob = survival_prob)
}


# --------------------------------------------------------------------------- #
#  Internal helpers                                                            #
# --------------------------------------------------------------------------- #

#' @keywords internal
.resolve_link <- function(link) {
  if (is.numeric(link)) return(as.integer(link))
  link <- match.arg(link, c("linear", "exponential"))
  switch(link, linear = 0L, exponential = 1L)
}

#' @keywords internal
.resolve_model <- function(model) {
  shortcuts <- list(cr = c(0L, 0L, 0L), dd = c(1L, 0L, 0L),
                    pd = c(0L, 1L, 0L), ep = c(0L, 0L, 1L))
  if (is.character(model)) {
    return(shortcuts[[match.arg(model, names(shortcuts))]])
  }
  if (inherits(model, "formula")) {
    return(.parse_model_formula(model))
  }
  model <- as.integer(model)
  if (length(model) != 3L || !all(model %in% 0:1)) {
    stop(paste0("'model' must be a formula (e.g. ~ N + PD), a string ",
                "(\"cr\", \"dd\", \"pd\", \"ep\"), or a length-3 binary ",
                "integer vector."))
  }
  model
}

#' @keywords internal
.parse_model_formula <- function(formula) {
  terms <- attr(stats::terms(formula), "term.labels")
  known <- c(N = 1L, PD = 2L, EP = 3L, E = 3L)
  terms_upper <- toupper(terms)
  model_bin <- c(0L, 0L, 0L)
  for (tm in terms_upper) {
    idx <- known[tm]
    if (is.na(idx)) {
      stop(sprintf("Unknown covariate '%s' in model formula. Use N, PD, and/or EP.", tm))
    }
    model_bin[idx] <- 1L
  }
  model_bin
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
#  Conditional simulation (augmentation) -- public output                      #
# --------------------------------------------------------------------------- #

#' @keywords internal
.sim_tree_conditional <- function(tree, pars, model_bin,
                                  n_trees = 1L, useDDD = TRUE, link = 0L,
                                  max_missing = 1e4, max_lambda = 500,
                                  maxN = NULL, num_threads = 1L) {
  brts  <- .extract_brts(tree)
  max_t <- brts[1L]

  expected_n <- 2L + 2L * sum(model_bin)
  if (length(pars) != expected_n) stop(.pars_error_msg(model_bin, expected_n))

  L_extant <- .extract_Ltable(tree)

  aug <- tryCatch(
    .augment_tree_internal(tree, pars = pars, model_bin = model_bin,
                           sample_size = n_trees,
                           max_missing = max_missing, max_lambda = max_lambda,
                           maxN = maxN, num_threads = num_threads, link = link),
    error = function(e) NULL
  )

  # Failure path
  if (is.null(aug) || length(aug$trees) == 0L) {
    if (n_trees == 1L) return(list(tas = NULL, log_q = NA_real_))
    return(list(trees = list(), log_q = numeric(0L)))
  }

  # Build augmented phylo for each valid draw
  tas_list <- lapply(aug$trees, function(df) {
    if (!useDDD) return(NULL)
    L <- .aug_to_Ltable(df, max_t, brts, L_extant)
    if (is.null(L)) return(NULL)
    tryCatch(DDD::L2phylo(L, dropextinct = FALSE), error = function(e) NULL)
  })

  log_q <- aug$logg   # log q(z | obs, theta) for each draw

  if (n_trees == 1L) {
    return(list(tas = tas_list[[1L]], log_q = log_q[1L]))
  }

  list(trees = tas_list, log_q = log_q)
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
  aug_labels <- list()
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
#  Internal augmentation (raw C++ output -- for inference use only)            #
# --------------------------------------------------------------------------- #

#' @keywords internal
.augment_tree_internal <- function(tree,
                                   pars,
                                   model_bin   = c(0L, 0L, 0L),
                                   sample_size = 1L,
                                   max_missing = 1e4,
                                   max_lambda  = 500,
                                   maxN        = NULL,
                                   num_threads = 1L,
                                   link        = 0L) {
  brts  <- .extract_brts(tree)
  pars8 <- .expand_pars(pars, model_bin)
  if (is.null(maxN)) maxN <- max(2000L, 200L * as.integer(sample_size))
  augment_trees(
    brts        = brts,
    pars        = as.numeric(pars8),
    sample_size = as.integer(sample_size),
    maxN        = as.integer(maxN),
    max_missing = as.integer(max_missing),
    max_lambda  = as.numeric(max_lambda),
    num_threads = as.integer(num_threads),
    model       = as.integer(model_bin),
    link        = as.integer(link)
  )
}
