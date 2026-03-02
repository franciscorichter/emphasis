#' Simulate a phylogenetic tree under a diversification model
#'
#' Unified wrapper for tree simulation, dispatching to a C++ (fast) or pure-R
#' (slow) implementation and to the appropriate model-specific backend.
#'
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
#' @param max_t Crown age (simulation time).
#' @param model Diversification model: \code{"pd"} (phylogenetic diversity
#'   dependence), \code{"dd"} (diversity dependence, exponential rates), or
#'   \code{"cr"} (constant rate). Default \code{"pd"}.
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
#'     \item{\code{tes}}{Reconstructed phylogeny (extant tips only);
#'       \code{NULL} when \code{useDDD = FALSE}, \code{fast = FALSE}, or
#'       the simulation failed.}
#'     \item{\code{tas}}{Full phylogeny including extinct lineages;
#'       \code{NULL} when \code{useDDD = FALSE} or the simulation failed.
#'       For \code{fast = FALSE} this is the tree returned directly by the
#'       R simulator.}
#'     \item{\code{L}}{L-table matrix; \code{NULL} for pure-R runs.}
#'     \item{\code{brts}}{Branching times of the reconstructed tree;
#'       \code{NULL} when \code{useDDD = FALSE} or the simulation failed.}
#'     \item{\code{status}}{\code{"done"}, \code{"extinct"}, or
#'       \code{"too_large"}.}
#'     \item{\code{model}}{The model used.}
#'     \item{\code{pars}}{The parameter vector supplied.}
#'   }
#' @examples
#' \dontrun{
#' # PD model — fast C++
#' tree <- simulate_tree(c(0.1, 0.5, -0.01, 0.01), max_t = 5)
#'
#' # DD model — fast C++ only
#' tree <- simulate_tree(c(0.1, -0.01, 0, -2.5, 0, 0), max_t = 5, model = "dd")
#'
#' # CR model — slow pure R
#' tree <- simulate_tree(c(0.1, 0.5), max_t = 5, model = "cr", fast = FALSE)
#' }
#' @export
simulate_tree <- function(pars, max_t,
                          model    = c("pd", "dd", "cr", "ep"),
                          fast     = TRUE,
                          max_lin  = 1e6,
                          max_tries = 100,
                          useDDD   = TRUE) {
  model <- match.arg(model)

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

#' Augment an extant tree with simulated missing extinct lineages
#'
#' Given an observed extant phylogenetic tree, draws one or more augmented
#' copies in which the missing (extinct) lineages are filled in by stochastic
#' simulation under the model. This is the "conditional simulation" counterpart
#' of \code{\link{simulate_tree}}: instead of simulating a complete tree from
#' scratch, it conditions on the observed extant topology and adds back the
#' unobserved extinct branches.
#'
#' @param tree Observed extant tree. Accepts:
#'   \itemize{
#'     \item the list returned by \code{\link{simulate_tree}}
#'       (uses \code{$brts}),
#'     \item a \code{phylo} object (extant tips only), or
#'     \item a numeric vector of branching times (crown age first,
#'       sorted decreasing).
#'   }
#' @param pars Numeric parameter vector \code{c(mu, lambda0, betaN, betaP)}.
#' @param sample_size Number of independent augmented trees to draw.
#'   Default \code{1}.
#' @param max_missing Maximum number of missing (extinct) lineages added per
#'   augmented tree. Default \code{1e4}.
#' @param max_lambda Maximum per-lineage speciation rate during augmentation.
#'   Default \code{500}.
#' @param soc Stem (1) or crown (2) age condition. Default \code{2}.
#' @param num_threads Number of threads for parallel augmentation.
#'   Default \code{1}.
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{trees}}{List of \code{sample_size} augmented trees, each a
#'       data frame with columns \code{brts} (time), \code{n} (richness),
#'       \code{t_ext} (extinction time; \code{1e11} = extant tip), \code{pd}.}
#'     \item{\code{weights}}{Log importance weights for each augmented tree.}
#'     \item{\code{fhat}}{Monte Carlo log-likelihood estimate.}
#'     \item{\code{logf}}{Log model likelihood of each augmented tree.}
#'     \item{\code{logg}}{Log sampling probability of each augmented tree.}
#'   }
#' @seealso \code{\link{simulate_tree}} for full forward simulation;
#'   \code{\link{estimate_rates}} which uses augmentation internally.
#' @examples
#' \dontrun{
#' set.seed(42)
#' sim <- simulate_tree(c(0.1, 0.5, -0.02, 0.01), max_t = 5)
#'
#' # Draw 10 augmented copies of the extant tree
#' aug <- augment_tree(sim, pars = c(0.1, 0.5, -0.02, 0.01), sample_size = 10)
#' length(aug$trees)     # 10 augmented trees
#' aug$fhat              # Monte Carlo log-likelihood estimate
#' head(aug$trees[[1]])  # first augmented tree as data frame
#' }
#' @export
augment_tree <- function(tree,
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
