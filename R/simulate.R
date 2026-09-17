#' Simulate a phylogenetic tree under the general diversification model
#'
#' Unified entry point for forward and conditional tree simulation.
#' When \code{tree = NULL} (default) a complete tree is simulated forward
#' from the crown; when \code{tree} is supplied the observed extant tree is
#' augmented with stochastically drawn extinct lineages.
#'
#' @section Model specification:
#' The model's covariates are diversity \code{N}, age-imbalance
#' \code{D = E - M} (a lineage's isolation relative to the clade mean
#' \code{M}, which is used only internally to centre \code{D}) and
#' evolutionary distinctiveness \code{ED} (the fair-proportion score on the
#' complete tree: the branches from the crown to the lineage, each divided by
#' the alive lineages below it; its last term is the pendant age and its sum
#' over lineages is Faith's phylogenetic diversity).  Selected by named
#' shortcuts or a formula:
#' \tabular{ll}{
#'   \code{"cr"}  \tab constant rate\cr
#'   \code{"dd"}  \tab diversity dependence (N)\cr
#'   \code{"d"}   \tab age-imbalance (D)\cr
#'   \code{"nd"}  \tab N + D\cr
#'   \code{"ed"}  \tab evolutionary distinctiveness (ED)\cr
#'   \code{"ned"} \tab N + ED\cr
#' }
#' Formulas \code{~ N}, \code{~ N + D} and \code{~ N + ED} are also accepted
#' (\code{ep} is a legacy alias for \code{d}).  An \code{ED} model needs the
#' tree's topology (a \code{phylo}, not a branching-time vector) and is not
#' available under the gaussian link.
#'
#' @section Parameter vector:
#' Compact layout; \code{beta_D} is the age-imbalance coefficient and
#' \code{beta_ED} the distinctiveness coefficient (\code{M} is internal only,
#' never part of the user vector):
#' \preformatted{c(beta_0,  [beta_N],  [beta_D],  [beta_ED],
#'   gamma_0, [gamma_N], [gamma_D], [gamma_ED])}
#' Only active covariates are included; length = \code{2 + 2 * sum(model)}.
#'
#' @section Reproducibility:
#' Every call draws its seed for the C++ sampler from R's generator, so
#' \code{set.seed()} fixes both the forward simulator and the thinning
#' augmenter: two runs under one seed return the same tree, two runs under
#' different seeds do not. The retries a forward simulation makes draw a seed
#' each, so they are distinct trees rather than one tree redrawn. Under
#' \code{parallel::mclapply} each child has its own R stream and therefore its
#' own seeds. Augmentation with \code{num_threads > 1} is not reproducible:
#' each worker draws from its own substream of the seed, but the sample keeps
#' the first \code{n_trees} augmentations to finish.
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
#'   \code{phylo} object, or a numeric branching-time vector.  A \code{phylo}
#'   also supplies the topology the pendant-age covariates \code{M} and
#'   \code{D} are computed from; with a bare branching-time vector every
#'   observed lineage is recorded as dating from the crown and \code{D = 0} at
#'   every observed branching event.
#' @param pars Numeric parameter vector or matrix. When a matrix each row
#'   produces one forward simulation.
#' @param max_t Crown age (forward simulation only). Default \code{1}.
#' @param model Model specification (string, formula, or binary vector).
#'   Default \code{"cr"}.
#' @param max_lin Maximum lineages before declaring the tree too large.
#'   Default \code{1e6}.
#' @param max_tries Maximum additional attempts after extinction or overflow
#'   (forward simulation only). Retries are tracked at the R level to
#'   compute \code{survival_prob}. Default \code{1}.
#' @param useDDD Convert L-table to \code{phylo} via \pkg{DDD}.
#'   Default \code{TRUE}.
#' @param n_trees Number of augmented trees to draw (augmentation only).
#'   Default \code{1L}.
#' @param link Link function: \code{"linear"} (default),
#'   \code{"exponential"}, or \code{"gaussian"}.
#' @param max_missing Maximum extinct lineages per augmented tree.
#'   Default \code{1e4}.
#' @param max_lambda Maximum speciation rate (thinning bound) for augmentation.
#'   Default \code{500}.
#' @param maxN Maximum total augmentation attempts. \code{NULL} (default)
#'   sets this to \code{max(2000, 200 * n_trees)}.
#' @param num_threads Threads for parallel augmentation. Default \code{1L}.
#' @param rho Sampling fraction in \code{(0, 1]}. Default \code{1} (complete
#'   sampling). With \code{rho < 1} a random fraction of extant tips is dropped
#'   (forward simulation) or unsampled extant lineages are inserted
#'   (augmentation), matching incomplete taxon sampling.
#' @param method Augmentation method: \code{"bdi"} (default, exact BDI
#'   sampler) or \code{"thinning"} (C++ thinning proposal). The BDI sampler
#'   draws from the exact conditional distribution under constant rates
#'   (ESS = S, zero variance) and uses a self-consistent backward-forward
#'   iteration under diversity dependence. It covers N-only models
#'   (\code{"cr"}, \code{"dd"}) on the linear and exponential links;
#'   any other model or link falls back to thinning.
#'
#'   \strong{\code{method = "bdi"} returns \code{NULL} for part of its draws.}
#'   The BDI sampler has no notion of the two crown lineages, so a lineage born
#'   before the first observed branching is recorded as a daughter of observed
#'   node 0, which is not yet born then.  Such a draw cannot be turned into a
#'   \code{tas}: it is refused, and that element of \code{trees} (or
#'   \code{tas}, with \code{n_trees = 1}) is \code{NULL}.  On seven trees
#'   measured for this release, 0 to 168 of 200 draws built; a tree whose first
#'   observed branching is late loses most of them.  \code{log_q} is returned
#'   for every draw either way.  It was previously silent and wrong -- every
#'   such draw produced a \code{phylo} with an edge of negative length.  The
#'   \code{"thinning"} method draws the parent from every lineage alive, the
#'   two crown lineages included, and is not affected.
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
#' @importFrom stats rbinom
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
                          num_threads = 1L,
                          rho         = 1.0,
                          method      = "bdi") {

  model_bin <- .resolve_model(model)
  link_int  <- .resolve_link(link)
  .check_rho(rho, "rho")

  if (!is.numeric(pars) || length(pars) == 0L)
    stop("'pars' must be a non-empty numeric vector or matrix.")

  # ---------------------------------------------------------------------- #
  #  Batch forward simulation (matrix pars)                                 #
  # ---------------------------------------------------------------------- #
  if (is.matrix(pars)) {
    if (nrow(pars) == 0L) stop("'pars' matrix must have at least one row.")
    do_sim <- function(i)
      simulate_tree(tree, pars[i, ], max_t, model_bin,
                    max_lin, max_tries, useDDD, n_trees, link_int,
                    max_missing, max_lambda, maxN, 1L, rho, method)
    if (num_threads > 1L && .Platform$OS.type == "unix") {
      sims <- parallel::mclapply(seq_len(nrow(pars)), do_sim,
                                 mc.cores = num_threads)
    } else {
      sims <- lapply(seq_len(nrow(pars)), do_sim)
    }
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
  method <- match.arg(method, c("thinning", "bdi"))

  if (!is.null(tree)) {
    return(.sim_tree_conditional(tree, pars, model_bin,
                                 as.integer(n_trees), useDDD, link_int,
                                 max_missing, max_lambda, maxN,
                                 as.integer(num_threads), rho = rho,
                                 method = method))
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
    raw        <- simulate_div_tree_cpp(pars8, model_bin, max_t, max_lin_i, 0L,
                                        link_int, seed = .draw_seed())
    n_attempts <- n_attempts + 1L
  }

  survival_prob <- if (raw$status == "done") 1.0 / n_attempts else 0.0

  # Apply incomplete sampling: randomly drop (1-rho) of extant tips
  L <- if (raw$status == "done") raw$Ltable else NULL
  if (!is.null(L) && rho < 1.0) {
    extant_idx <- which(L[, 4] == -1)
    # Must keep at least 2 tips for a valid tree
    n_drop <- rbinom(1, length(extant_idx), 1.0 - rho)
    n_drop <- min(n_drop, length(extant_idx) - 2L)
    if (n_drop > 0L) {
      drop <- sample(extant_idx, n_drop)
      L[drop, 4] <- max_t  # mark as "extinct at present" = removed
    }
  }

  tes <- tas <- NULL
  if (!is.null(L) && useDDD) {
    tes <- tryCatch(DDD::L2phylo(L, dropextinct = TRUE),
                    error = function(e) NULL)
    tas <- tryCatch(DDD::L2phylo(L, dropextinct = FALSE),
                    error = function(e) NULL)
  }

  list(tes           = tes,
       tas           = tas,
       L             = L,
       status        = raw$status,
       survival_prob = survival_prob)
}


# --------------------------------------------------------------------------- #
#  Internal helpers                                                            #
# --------------------------------------------------------------------------- #

# One seed for one C++ sampler call, drawn from R's generator.
#
# The C++ samplers used to seed themselves from the wall clock XOR the thread
# id, so set.seed() did not reach them and a forked child inherited its
# parent's engine state byte for byte (H47).  Drawing the seed here puts them
# on R's stream: set.seed() fixes them, and two forked children, whose R
# streams differ, draw different trees.  Each call gets its own draw, so the
# retries of a forward simulation are distinct trees rather than one tree
# repeated.
#' @keywords internal
.draw_seed <- function() sample.int(.Machine$integer.max, 1L)

#' @keywords internal
.resolve_link <- function(link) {
  if (is.numeric(link)) return(as.integer(link))
  link <- match.arg(link, c("linear", "exponential", "gaussian"))
  switch(link, linear = 0L, exponential = 1L, gaussian = 2L)
}

#' @keywords internal
.resolve_model <- function(model) {
  # Covariate slots c(use_N, use_M, use_D, use_ED):
  #   N = diversity, D = E - M age-imbalance, ED = evolutionary distinctiveness
  #   (fair proportion on the complete tree; see inst/include/ed_covariate.hpp).
  #   M (slot 2) is retained only as the internal centering reference for D and
  #   is not user-selectable.
  # Shortcuts: "dd" -> N, "d" -> D, "nd" -> N + D, "ed" -> ED, "ned" -> N + ED.
  # "ep"/"rd" are legacy D aliases.  A length-3 vector is the pre-ED layout
  # and is padded.
  shortcuts <- list(cr = c(0L, 0L, 0L, 0L), dd = c(1L, 0L, 0L, 0L),
                    d  = c(0L, 0L, 1L, 0L), nd = c(1L, 0L, 1L, 0L),
                    ed = c(0L, 0L, 0L, 1L), ned = c(1L, 0L, 0L, 1L),
                    rd = c(0L, 0L, 1L, 0L), ep = c(0L, 0L, 1L, 0L))
  if (is.character(model)) {
    return(shortcuts[[match.arg(model, names(shortcuts))]])
  }
  if (inherits(model, "formula")) {
    return(.parse_model_formula(model))
  }
  model <- as.integer(model)
  if (!(length(model) %in% c(3L, 4L)) || !all(model %in% 0:1)) {
    stop(paste0("'model' must be a formula (e.g. ~ N + D, ~ N + ED), a string ",
                "(\"cr\", \"dd\", \"d\", \"nd\", \"ed\", \"ned\"), or a binary ",
                "integer vector of length 3 or 4."))
  }
  .pad_model_bin(model)
}

# The canonical 4-slot model vector from a 3- or 4-slot one.
#' @keywords internal
.pad_model_bin <- function(model_bin) {
  model_bin <- as.integer(model_bin)
  if (length(model_bin) == 3L) model_bin <- c(model_bin, 0L)
  model_bin
}

# Number of covariate slots in the canonical layout, and the full parameter
# vector's length: c(beta_0, beta_N, beta_M, beta_D, gamma_0, gamma_N,
# gamma_M, gamma_D, beta_ED, gamma_ED).  The ED coefficients are appended so
# that an 8-element vector is exactly "ED absent".
.n_slots <- 4L
.n_full  <- 10L

#' @keywords internal
.parse_model_formula <- function(formula) {
  terms <- attr(stats::terms(formula), "term.labels")
  # User covariates N, D (slot 3) and ED (slot 4); legacy aliases EP/E for D.
  # M (slot 2) is internal only and not user-selectable.
  known <- c(N = 1L, D = 3L, EP = 3L, E = 3L, ED = 4L)
  terms_upper <- toupper(terms)
  model_bin <- c(0L, 0L, 0L, 0L)
  for (tm in terms_upper) {
    idx <- known[tm]
    if (is.na(idx)) {
      stop(sprintf("Unknown covariate '%s' in model formula. Use N, D and/or ED.", tm))
    }
    model_bin[idx] <- 1L
  }
  model_bin
}

# Where each covariate slot's coefficient sits in the full vector: beta at
# 2, 3, 4 for N, M, D and 9 for ED; gamma at 6, 7, 8 and 10.  The intercepts
# are 1 and 5.
.slot_beta  <- c(2L, 3L, 4L, 9L)
.slot_gamma <- c(6L, 7L, 8L, 10L)

# Expand compact pars to the full vector for C++.
# Layout: c(beta_0, beta_N, beta_M, beta_D, gamma_0, gamma_N, gamma_M, gamma_D,
#           beta_ED, gamma_ED).  With ED inactive the result is the first 8
# slots, the layout every pre-ED caller expects; with ED active all 10.
#' @keywords internal
.expand_pars <- function(pars, model_bin) {
  model_bin  <- .pad_model_bin(model_bin)
  expected_n <- 2L + 2L * sum(model_bin)
  if (length(pars) != expected_n) stop(.pars_error_msg(model_bin, expected_n))
  active <- which(model_bin == 1L)
  n_lam  <- 1L + length(active)
  full   <- numeric(.n_full)
  full[1L] <- pars[1L]
  full[5L] <- pars[n_lam + 1L]
  if (length(active) > 0L) {
    full[.slot_beta[active]]  <- pars[2L:n_lam]
    full[.slot_gamma[active]] <- pars[(n_lam + 2L):length(pars)]
  }
  if (model_bin[4L] == 1L) full else full[1:8]
}

#' @keywords internal
.pars_error_msg <- function(model_bin, expected_n) {
  model_bin <- .pad_model_bin(model_bin)
  lam <- paste(c("beta_0",
                 if (model_bin[1]) "beta_N",
                 if (model_bin[2]) "beta_M",
                 if (model_bin[3]) "beta_D",
                 if (model_bin[4]) "beta_ED"), collapse = ", ")
  mu  <- paste(c("gamma_0",
                 if (model_bin[1]) "gamma_N",
                 if (model_bin[2]) "gamma_M",
                 if (model_bin[3]) "gamma_D",
                 if (model_bin[4]) "gamma_ED"), collapse = ", ")
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
                                  maxN = NULL, num_threads = 1L,
                                  rho = 1.0, method = "thinning") {
  brts  <- .extract_brts(tree)
  max_t <- brts[1L]

  expected_n <- 2L + 2L * sum(model_bin)
  if (length(pars) != expected_n) stop(.pars_error_msg(model_bin, expected_n))

  L_extant <- .extract_Ltable(tree)

  if (method == "bdi" && !.bdi_supported(model_bin, link, rho)) method <- "thinning"

  if (method == "bdi") {
    aug <- tryCatch(
      .augment_tree_bdi(tree, pars = pars, model_bin = model_bin,
                        sample_size = n_trees,
                        max_missing = max_missing, link = link,
                        rho = rho),
      error = function(e) NULL
    )
  } else {
    aug <- tryCatch(
      .augment_tree_internal(tree, pars = pars, model_bin = model_bin,
                             sample_size = n_trees,
                             max_missing = max_missing, max_lambda = max_lambda,
                             maxN = maxN, num_threads = num_threads, link = link,
                             rho = rho),
      error = function(e) NULL
    )
  }

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

# Reserved lineage ids of the two crown lineages
#
# The two crown lineages carry no node, so the augmentation cannot name them by
# a node id; \code{inst/include/model_helpers.hpp} gives them ids of their own.
# \code{crown_id_a} is the one the forward sweep splits at the first observed
# branching event, which is what \code{.aug_to_Ltable} maps it onto.
#' @keywords internal
.crown_id_a <- -2L
#' @keywords internal
.crown_id_b <- -3L

# Merge augmented branches into the extant L-table.
#
# Every augmented lineage names the lineage it was drawn from, so every one of
# them becomes a row: an id of a lineage born at an observed branching event, or
# one of the two crown ids.  Nothing is dropped, and an attachment to a lineage
# born after the lineage attached to it is refused rather than turned into a
# phylo with a negative edge.
#' @keywords internal
.aug_to_Ltable <- function(df, max_t, brts, L_extant) {
  if (is.null(L_extant) || nrow(L_extant) < 2L) return(NULL)

  t_ext_tip      <- 1e11
  t_ext_unsampled <- 5e10
  aug <- df[df$t_ext != 0.0 & df$t_ext < t_ext_tip, , drop = FALSE]

  if (nrow(aug) == 0L) return(L_extant)
  if (any(aug$parent_id == -1L)) return(NULL)

  # Label and birth time of every observed lineage, by the id the C++ layer
  # gives it.  The crown lineages are the two oldest rows of the L-table; the
  # one that speciates first is the parent of the oldest row below them, and is
  # the one the sweep splits at the first observed event.
  crown <- which(L_extant[, 1L] == max(L_extant[, 1L]))
  if (length(crown) != 2L) return(NULL)
  crown_lbl <- as.integer(L_extant[crown, 3L])
  rest <- setdiff(seq_len(nrow(L_extant)), crown)
  if (length(rest)) {
    first <- as.integer(L_extant[rest[which.max(L_extant[rest, 1L])], 2L])
    if (first %in% crown_lbl) crown_lbl <- c(first, setdiff(crown_lbl, first))
  }

  lbl <- list(); birth <- list()
  put <- function(id, label, b) {
    lbl[[as.character(id)]]   <<- label
    birth[[as.character(id)]] <<- b
  }
  put(.crown_id_a, crown_lbl[1L], max_t)
  put(.crown_id_b, crown_lbl[2L], max_t)
  # Observed node id k is the lineage born at the k-th observed branching event
  # in forward time, which is entry k + 2 of the crown-first branching times.
  for (k in seq_len(length(brts) - 1L) - 1L) {
    r <- which.min(abs(L_extant[, 1L] - brts[k + 2L]))
    put(k, as.integer(L_extant[r, 3L]), L_extant[r, 1L])
  }

  next_lbl   <- as.integer(max(abs(L_extant[, 3L]))) + 1L
  aug        <- aug[order(aug$brts), , drop = FALSE]
  new_rows   <- matrix(0.0, nrow = nrow(aug), ncol = 4L)

  for (k in seq_len(nrow(aug))) {
    pk    <- as.character(aug$parent_id[k])
    p_lbl <- lbl[[pk]]
    if (is.null(p_lbl)) return(NULL)
    b <- max_t - aug$brts[k]
    if (b > birth[[pk]] + 1e-9) return(NULL)   # parent born after its daughter

    l        <- if (p_lbl < 0L) -next_lbl else next_lbl
    next_lbl <- next_lbl + 1L

    # Unsampled extant species (t_ext == t_ext_unsampled) are alive at present:
    # set L-table extinction time to -1 (DDD convention for extant tips)
    ext_k <- if (aug$t_ext[k] == t_ext_unsampled) -1 else max_t - aug$t_ext[k]
    new_rows[k, ] <- c(b, p_lbl, l, ext_k)
    put(aug$id[k], l, b)
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
                                   link        = 0L,
                                   rho         = 1.0,
                                   seed        = .draw_seed()) {
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
    link        = as.integer(link),
    rho         = as.numeric(rho),
    parent_tip_start = .pts(brts),
    seed        = as.integer(seed),
    parent_id   = .pid(brts)
  )
}
