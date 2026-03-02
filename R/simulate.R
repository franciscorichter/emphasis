#' Simulate a phylogenetic tree under a diversification model
#'
#' Unified wrapper for tree simulation. When \code{tree = NULL} (default),
#' simulates a complete tree from scratch (forward simulation). When
#' \code{tree} is supplied, performs conditional simulation: augments the
#' provided extant tree with stochastically drawn extinct lineages.
#'
#' @param tree Optional. An observed extant tree for conditional simulation.
#'   Accepts:
#'   \itemize{
#'     \item \code{NULL} (default): forward simulation from scratch.
#'     \item the list returned by \code{\link{simulate_tree}} (uses
#'       \code{$brts} and \code{$tes}),
#'     \item a \code{phylo} object (extant tips only), or
#'     \item a numeric vector of branching times (crown age first,
#'       sorted decreasing).
#'   }
#' @param pars Numeric parameter vector. Length and meaning depend on
#'   \code{model}:
#'   \describe{
#'     \item{\code{"pd"}}{4 values: \code{c(mu, lambda0, betaN, betaP)}.
#'       Speciation rate = \code{max(0, lambda0 + betaN*N + betaP*P/N)};
#'       extinction rate = \code{mu}.}
#'     \item{\code{"dd"}}{6 values: \code{c(A0, An, Ap, B0, Bn, Bp)}.
#'       Speciation rate = \code{exp(A0 + An*N + Ap*P)};
#'       extinction rate = \code{exp(B0 + Bn*N + Bp*P)}.}
#'     \item{\code{"cr"}}{2 values: \code{c(mu, lambda0)}.
#'       Constant-rate birth-death (special case of \code{"pd"} with
#'       \code{betaN = betaP = 0}).}
#'   }
#' @param max_t Crown age (simulation time). Required for forward simulation
#'   (\code{tree = NULL}); ignored when \code{tree} is provided (crown age is
#'   extracted from the tree).
#' @param model Diversification model: \code{"pd"} (phylogenetic diversity
#'   dependence), \code{"dd"} (diversity dependence, exponential rates), or
#'   \code{"cr"} (constant rate). Default \code{"pd"}. Only \code{"pd"} is
#'   supported for conditional simulation.
#' @param fast Logical. Use the fast C++ implementation (\code{TRUE}, default)
#'   or the slower pure-R implementation (\code{FALSE}). The \code{"dd"} model
#'   has no R implementation and always uses C++.
#' @param max_lin Maximum number of lineages before the simulation is
#'   considered non-extinct. Default \code{1e6}.
#' @param max_tries Maximum attempts to obtain a non-extinct, non-oversized
#'   tree. Default \code{100}.
#' @param useDDD Logical. Convert the L-table to \code{phylo} objects via
#'   \pkg{DDD}. Only applies when \code{fast = TRUE}. Default \code{TRUE}.
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{tes}}{Reconstructed phylogeny (extant tips only).}
#'     \item{\code{tas}}{Full phylogeny including extinct lineages;
#'       \code{NULL} when topology reconstruction fails.}
#'     \item{\code{L}}{L-table matrix (forward simulation only).}
#'     \item{\code{brts}}{Branching times of the extant tree.}
#'     \item{\code{status}}{\code{"done"} or \code{"failed"}.}
#'     \item{\code{model}}{The model used.}
#'     \item{\code{pars}}{The parameter vector supplied.}
#'   }
#' @examples
#' \dontrun{
#' # Forward simulation (PD model, fast C++)
#' tr <- simulate_tree(pars = c(0.1, 0.5, -0.01, 0.01), max_t = 5)
#'
#' # Conditional simulation — augment with extinct lineages
#' aug <- simulate_tree(tree = tr, pars = c(0.1, 0.5, -0.01, 0.01))
#' aug$tas  # phylo with extinct branches
#'
#' # Conditional from a plain phylo object
#' aug2 <- simulate_tree(tree = tr$tes, pars = c(0.1, 0.5, -0.01, 0.01))
#' }
#' @export
simulate_tree <- function(tree      = NULL,
                          pars,
                          max_t     = NULL,
                          model     = c("pd", "dd", "cr", "ep"),
                          fast      = TRUE,
                          max_lin   = 1e6,
                          max_tries = 100,
                          useDDD    = TRUE) {
  model <- match.arg(model)
  if (!is.null(tree)) return(.sim_tree_conditional(tree, pars, model))

  if (is.null(max_t)) stop("'max_t' is required for forward simulation (tree = NULL).")

  expected_n  <- c(pd = 4L, dd = 6L, cr = 2L, ep = 5L)
  param_names <- c(pd = "c(mu, lambda0, betaN, betaP)",
                   dd = "c(A0, An, Ap, B0, Bn, Bp)",
                   cr = "c(mu, lambda0)",
                   ep = "c(mu, lambda0, betaN, betaP, betaE)")
  if (length(pars) != expected_n[[model]]) {
    stop(sprintf("'%s' model requires %d parameters: %s",
                 model, expected_n[[model]], param_names[[model]]))
  }

  if (!fast && model == "dd") {
    message("No R implementation for 'dd' model; using fast C++ instead.")
    fast <- TRUE
  }

  if (!fast && model == "ep") {
    message("No R implementation for 'ep' model; using fast C++ instead.")
    fast <- TRUE
  }

  result <- if (fast) {
    .sim_tree_fast(pars, max_t, model, max_lin, max_tries, useDDD)
  } else {
    .sim_tree_slow(pars, max_t, model, max_tries)
  }

  result$model <- model
  result$pars  <- pars
  result
}

.sim_tree_fast <- function(pars, max_t, model, max_lin, max_tries, useDDD) {
  switch(model,
    pd = sim_tree_pd_cpp(pars, max_t,
                         max_lin   = max_lin,
                         max_tries = max_tries,
                         useDDD    = useDDD),
    dd = {
      raw  <- simulate_div_tree_cpp(pars, max_t, max_lin, max_tries)
      L    <- raw$tree
      tes  <- tas <- brts <- NULL
      if (raw$status == "done" && useDDD) {
        tes  <- DDD::L2phylo(L, dropextinct = TRUE)
        tas  <- DDD::L2phylo(L, dropextinct = FALSE)
        brts <- DDD::L2brts(L, dropextinct = TRUE)
      }
      list(tes = tes, tas = tas, L = L, brts = brts, status = raw$status)
    },
    cr = sim_tree_pd_cpp(c(pars[1], pars[2], 0, 0), max_t,
                         max_lin   = max_lin,
                         max_tries = max_tries,
                         useDDD    = useDDD),
    ep = sim_tree_ep_cpp(pars, max_t,
                         max_lin   = max_lin,
                         max_tries = max_tries,
                         useDDD    = useDDD)
  )
}

.sim_tree_slow <- function(pars, max_t, model, max_tries) {
  pd_pars <- switch(model,
    pd = pars,
    cr = c(pars[1], pars[2], 0, 0)
  )
  raw    <- NULL
  status <- "extinct"
  for (attempt in seq_len(max_tries)) {
    raw <- sim_tree_pd_R(pd_pars, max_t)
    # result = c(N, t, P); non-extinct iff N >= 2 and t reached max_t
    if (raw$result[1] >= 2 && abs(raw$result[2] - max_t) < sqrt(.Machine$double.eps)) {
      status <- "done"
      break
    }
  }
  list(tes = NULL, tas = raw$phy, L = NULL, brts = NULL, status = status)
}


# --------------------------------------------------------------------------- #
#  Conditional simulation helpers                                              #
# --------------------------------------------------------------------------- #

#' @keywords internal
.sim_tree_conditional <- function(tree, pars, model) {
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

  # Extant L-table needed for topology reconstruction of the augmented tree.
  # Required when tree is a simulate_tree() result; plain phylo inputs are
  # unsupported (tas will be NULL).
  L_extant <- if (is.list(tree) && !inherits(tree, "phylo") && !is.null(tree$L)) {
    tree$L
  } else if (inherits(tree, "phylo")) {
    tryCatch(DDD::phylo2L(tree), error = function(e) NULL)
  } else {
    NULL
  }

  # One augmented draw (sample_size = 1)
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
                status = "failed", model = model, pars = pars))
  }

  df <- aug$trees[[1L]]   # data frame: brts, n, t_ext, pd, id, parent_id

  # Build extended L-table and convert to phylo
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
       model  = model,
       pars   = pars)
}


# Convert an augmented-tree data frame (from e_cpp) into a DDD-compatible
# L-table that includes both the extant tree and the simulated extinct lineages.
#
# Design notes:
#   - `L_extant` (the forward-sim L-table stored in tree$L) provides the base
#     topology including the two mandatory root rows at backward time = max_t.
#   - emphasis init-tree node with id=k (C++ 0-indexed, sorted by ascending
#     forward brts) was born at backward time brts[k+2] (R 1-indexed; brts[1]
#     = max_t crown age).  We match to L_extant by backward birth time to
#     retrieve the DDD label for each init-tree node.
#   - Augmented extinct lineages (parent_id != -1, t_ext < 1e11) are appended
#     as new rows, inheriting the clade sign of their parent.
#' @keywords internal
.aug_to_Ltable <- function(df, max_t, brts, L_extant) {
  if (is.null(L_extant)) return(NULL)

  t_ext_tip_val <- 1e11

  # Keep only speciation rows (drop t_ext == 0 extinction-event rows)
  sp <- df[df$t_ext != 0.0, , drop = FALSE]

  # Augmented extinct lineages: have a parent AND go extinct before present
  aug <- sp[sp$parent_id != -1L & sp$t_ext < t_ext_tip_val, , drop = FALSE]

  if (nrow(aug) == 0L) return(L_extant)  # nothing to add

  # Map emphasis init-tree id (0-indexed) → DDD label in L_extant.
  # id=k → backward birth = brts[k+2] (brts[1]=max_t, brts[2..n+1] = internal)
  id_to_label <- function(pid) {
    if (pid < 0L || (pid + 2L) > length(brts)) return(NA_integer_)
    bwd   <- brts[pid + 2L]
    row_k <- which.min(abs(L_extant[, 1L] - bwd))
    as.integer(L_extant[row_k, 3L])
  }

  # Next unused label magnitude
  next_lbl <- as.integer(max(abs(L_extant[, 3L]))) + 1L

  # Process in chronological (ascending forward brts) order so that a node's
  # parent is always assigned a label before the node itself.
  aug <- aug[order(aug$brts), , drop = FALSE]

  # Labels assigned to augmented nodes (which may themselves parent later ones)
  aug_labels <- setNames(integer(0), character(0))

  new_rows <- matrix(0.0, nrow = nrow(aug), ncol = 4L)

  for (k in seq_len(nrow(aug))) {
    pid     <- aug$parent_id[k]
    pid_key <- as.character(pid)

    p_lbl <- if (pid_key %in% names(aug_labels)) {
      aug_labels[[pid_key]]
    } else {
      id_to_label(pid)
    }
    if (is.na(p_lbl) || length(p_lbl) == 0L) {
      p_lbl <- as.integer(L_extant[1L, 3L])  # fallback: first root label
    }

    lbl      <- if (p_lbl < 0L) -next_lbl else next_lbl
    next_lbl <- next_lbl + 1L

    new_rows[k, 1L] <- max_t - aug$brts[k]    # backward birth time
    new_rows[k, 2L] <- as.double(p_lbl)        # parent label
    new_rows[k, 3L] <- as.double(lbl)          # own label
    new_rows[k, 4L] <- max_t - aug$t_ext[k]   # backward death time

    # Store this node's label so later nodes can look it up as a parent
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
  brts   <- .extract_brts(tree)
  max_n  <- max(100L, as.integer(10L * sample_size))
  n      <- length(pars)
  e_cpp(brts        = brts,
        init_pars   = as.numeric(pars),
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
  return(sum(DDD::L2phylo(l_table, dropextinct = T)$edge.length))
}

#' simulation function to simulate a tree under the pd model, returns the full
#' tree.
#' @description super fast function to simulate the process of diversification
#' with diversity dependence and phylogenetic diversity dependence
#' @param pars parameter vector with c(mu, lambda_0, beta_N, beta_P)
#' @param max_t crown age
#' @return a list with the phylogeny, and a vector with population size,
#' time of extinction (equal to the crown #' age in the absence of extinction),
#' and the phylogenetic diversity at the time of extinction (or crown age).
#' @export
sim_tree_pd_R <- function(pars, max_t) {
  N1 <- 1
  N2 <- 1
  N <- N1 + N2
  t <- 0
  tree <- matrix(nrow = 2, ncol = 4)

  tree[1, ] <- c(0, 0, -1, -1)
  tree[2, ] <- c(0, -1, 2, -1)
  tree_ID <- 3

  mu <- pars[1]
  P <- 0
  while (t < max_t &&
         N1 >= 1 && N2 >= 1) {

    N <- N1 + N2
    spec_rate <- max(0, pars[2] +
                       pars[3] * N  +
                       ((P + N * (max_t - t) - t) / N) * pars[4])
    total_rate <- (spec_rate + mu) * N

    if (total_rate == 0) {
      t <- max_t
      break
    }

    next_event_time <- t + stats::rexp(n = 1, rate = total_rate)

    P <- calc_p(tree, t)

    if (next_event_time < max_t) {
      focal_spec <- max(0, pars[2] +
                          pars[3]*N  +
                          ((P + N * (next_event_time - t) - t) / N) * pars[4])

      pt = ((focal_spec + mu) * N ) / total_rate

      if (stats::runif(1) < pt) {
        if (stats::runif(1) < focal_spec / (focal_spec + mu)) {
          parent <- sample(which(tree[, 4] == -1), 1)

          new_ID <- tree_ID
          tree_ID <- tree_ID + 1
          if (tree[parent, 3] < 0) new_ID <- new_ID * -1

          new_spec <- c(next_event_time, tree[parent, 3], new_ID, -1)

          tree <- rbind(tree, new_spec)

          if (new_ID < 0) {
            N2 <- N2 + 1
          } else {
            N1 <- N1 + 1
          }
        } else {
          to_remove <- sample(which(tree[, 4] == -1), 1)

          tree[to_remove, 4] <- next_event_time

          if (tree[to_remove, 3] < 0) {
            N2 <- N2 - 1
          } else {
            N1 <- N1 - 1
          }
        }
      }
    }
    t <- next_event_time
  }

  if (N1 >= 1 && N2 >= 1) P <- calc_p(tree, t)
  N <- N1 + N2

  tree[, 1] <- max_t - tree[, 1]
  extinct <- which(tree[, 4] != -1)
  tree[extinct, 4] <- max_t - tree[extinct, 4]

  for_ddd <- as.matrix(tree[, 1:4])

  phy <- DDD::L2phylo(as.matrix(tree[, 1:4]), dropextinct = FALSE)

  t <- min(t, max_t)

  return(list("phy" = phy,
              "result" = c(N, t, P)))
}

#' simulation function to simulate a tree under the pd model, returning whether
#' the tree went extinct before max_t or not.
#' This function does not return a phylogenetic tree to improve computation
#' speed.
#' @description super fast function to simulate the process of diversification
#' with diversity dependence and phylogenetic diversity dependence
#' @param pars parameter vector with c(mu, lambda_0, beta_N, beta_P)
#' @param max_t crown age
#' @param num_repl number of replicates
#' @param max_lin number of lineages past which non-extinction is assumed.
#' @return a tibble with five columns: 1) whether or not the tree went extinct,
#' 2) the time of extinction (equal to the crown #' age in the absence of
#' extinction), 3) the number of tips,4) the phylogenetic diversity at the
# time of extinction (or crown age), and 5) the break condition indicating why
# the simulation was stopped (options: no break (none), exceeded max_t,
# extinction, or exceeded max_lin).
#' @export
sim_tree_is_extinct_pd <- function(pars, max_t, num_repl = 1, max_lin) {
    result <- simulate_pd_trees_cpp(pars, max_t, num_repl, max_lin)
    colnames(result) <- c("is_extinct", "t", "N", "P", "break_condition")
    result <- tibble::as_tibble(result)
    result$break_condition[result$break_condition == 3] <- "max_lin"
    result$break_condition[result$break_condition == 2] <- "extinction"
    result$break_condition[result$break_condition == 1] <- "max_t"
    result$break_condition[result$break_condition == 0] <- "none"

    return(result)
}

#' simulation function to simulate a tree under the pd model
#' @description super fast function to simulate the process of diversification
#' with diversity dependence and phylogenetic diversity dependence
#' @param pars parameter vector with c(mu, lambda_0, beta_N, beta_P)
#' @param max_t crown age
#' @param max_lin number of lineages past which non-extinction is assumed.
#' @param max_tries maximum number of tries to get a non-extinct tree.
#' @param useDDD whether to use DDD for tree conversion (default TRUE)
#' @return list with: 1) tes - reconstructed tree, 2) tas - tree with extinct
#' lineages, 3) L = Ltable and 4) brts - branching times of the reconstructed
#' tree and 5) status of simulation, options 1) "extinct", 2) "too_large" or 3)
#' "done".
#' @examples
#' set.seed(123)
#' pars <- c(mu = 0.1, lambda0 = 0.4, betaN = -0.05, betaP = 0.02)
#' tree_res <- try(sim_tree_pd_cpp(
#'   pars = pars,
#'   max_t = 5,
#'   max_lin = 1e4,
#'   max_tries = 5,
#'   useDDD = TRUE
#' ))
#' if (!inherits(tree_res, "try-error") && !is.null(tree_res$tes)) {
#'   head(tree_res$tes$tip.label)
#' }
#' @export
sim_tree_pd_cpp <- function(pars,
                            max_t,
                            max_lin = 1e6,
                            max_tries = 100,
                            useDDD=TRUE) {

  result <- simulate_single_pd_tree_cpp(pars,
                                        max_t,
                                        max_lin,
                                        max_tries)
  tes <- NULL
  tas <- NULL
  brts <- NULL

  if (result$status == "extinct") {
    warning("could not simulate tree, all trees went extinct, try increasing max_tries")
  }
  if (result$status == "too_large") {
    warning("could not simulate tree, all trees were too large, try increasing max_lin")
  }
  if (result$status == "done" & useDDD) {
    tes <- DDD::L2phylo(result$Ltable, dropextinct = TRUE)
    tas <- DDD::L2phylo(result$Ltable, dropextinct = FALSE)
    brts = DDD::L2brts(result$Ltable, dropextinct = TRUE)
  }

  out = list(tes = tes, tas = tas, L = result$Ltable, brts = brts,
             status = result$status)

  return(out)
}


#' Simulate a single phylogenetic tree under the EP (evolutionary pendant) model
#' @description Fast C++ Gillespie direct simulation where each lineage's
#' speciation rate depends on its own pendant edge length (time since last
#' divergence).
#' @param pars parameter vector with c(mu, lambda0, betaN, betaP, betaE)
#' @param max_t crown age
#' @param max_lin number of lineages past which non-extinction is assumed.
#' @param max_tries maximum number of tries to get a non-extinct tree.
#' @param useDDD whether to use DDD for tree conversion (default TRUE)
#' @return list with: 1) tes - reconstructed tree, 2) tas - tree with extinct
#' lineages, 3) L = Ltable and 4) brts - branching times of the reconstructed
#' tree and 5) status of simulation, options 1) "extinct", 2) "too_large" or 3)
#' "done".
#' @export
sim_tree_ep_cpp <- function(pars, max_t, max_lin = 1e6, max_tries = 100, useDDD = TRUE) {
  result <- simulate_single_ep_tree_cpp(pars, max_t, max_lin, max_tries)
  tes <- tas <- brts <- NULL
  if (result$status == "extinct")
    warning("could not simulate tree, all trees went extinct")
  if (result$status == "too_large")
    warning("could not simulate tree, all trees were too large")
  if (result$status == "done" && useDDD) {
    tes  <- DDD::L2phylo(result$Ltable, dropextinct = TRUE)
    tas  <- DDD::L2phylo(result$Ltable, dropextinct = FALSE)
    brts <- DDD::L2brts(result$Ltable, dropextinct = TRUE)
  }
  list(tes = tes, tas = tas, L = result$Ltable, brts = brts, status = result$status)
}


#' simulation function to simulate a tree under the pd model
#' @description super fast function to simulate the process of diversification
#' with diversity dependence and phylogenetic diversity dependence
#' @param mu_vec vector of extinction rates to explore
#' @param lambda_vec vector of lambda rates to explore
#' @param b_n_vec vector of B_n rates to explore
#' @param b_p_vec vector of B_p rates to explore
#' @param max_t crown age
#' @param max_N maximum number of lineages past which non-extinction is assumed.
#' @param num_repl number of replicates
#' @return a list with the used parameter combinations,
#' extinction of the phylogeny, and a vector with population size,
#' time of extinction (equal to the crown #' age in the absence of extinction),
#' and the phylogenetic diversity at the time of extinction (or crown age).
#' @export
sim_tree_pd_grid <- function(mu_vec,
                             lambda_vec,
                             b_n_vec,
                             b_p_vec,
                             max_t,
                             num_repl,
                             max_N) {
  result <- explore_grid_cpp(mu_vec,
                             lambda_vec,
                             b_n_vec,
                             b_p_vec,
                             max_t,
                             num_repl,
                             max_N)
  colnames(result) <- c("is_extinct", "t", "N", "P")
  result <- tibble::as_tibble(result)
  return(result)
}

