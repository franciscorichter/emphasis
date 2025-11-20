#' Safe C++ Tree Simulation Wrapper (for debugging)
#'
#' This function is a safe, minimal wrapper around the C++ tree simulation to avoid segfaults.
#' It checks input, output, and catches errors. It does not modify your original C++ or R functions.
#'
#' @param pars Numeric vector of parameters (length 6)
#' @param max_t Maximum simulation time
#' @param max_N Maximum number of taxa
#' @param max_tries Maximum number of tries
#' @return Number of species (tips) or NA on error
#' @export
#' Safe C++ Tree Simulation Wrapper (returns tree matrix and status)
#' Safe C++ Tree Simulation Wrapper (returns tree matrix and status)
safe_simulate_div_tree_cpp <- function(pars = c(0.1, 0.05, 0, 0, 0, 0), max_t = 10, max_N = 1000, max_tries = 10) {
  out <- tryCatch({
    res <- emphasis::safe_simulate_div_tree_cpp(pars, max_t, max_N, max_tries)
    # Defensive checks
    if (!is.list(res) || is.null(res$tree)) stop("No tree returned")
    tree_mat <- res$tree
    if (!is.matrix(tree_mat) || ncol(tree_mat) < 4) stop("Tree matrix malformed")
    return(list(tree = tree_mat, status = res$status))
  }, error = function(e) {
    message("Error: ", e$message)
    return(NULL)
  })
  return(out)
}

#' Convert simulation matrix to phylo object (requires ape)
sim_matrix_to_phylo <- function(tree_mat) {
  # Expecting columns: parent, time, id
  requireNamespace("ape")
  # Remove root row if present (parent==0)
  tree_mat <- tree_mat[tree_mat[,1]!=0,,drop=FALSE]
  n_tips <- length(unique(tree_mat[,3]))
  edges <- cbind(match(tree_mat[,1], tree_mat[,3]), seq_len(nrow(tree_mat)))
  phy <- list(
    edge = edges,
    Nnode = nrow(tree_mat) - n_tips + 1,
    tip.label = as.character(tree_mat[,3]),
    edge.length = tree_mat[,2]
  )
  class(phy) <- "phylo"
  phy <- ape::reorder.phylo(phy)
  phy
}

#' Plot a simulated tree from the C++ function
plot_safe_cpp_tree <- function(pars = c(0.1, 0.05, 0, 0, 0, 0), max_t = 10, max_N = 1000, max_tries = 10) {
  sim <- safe_simulate_div_tree_cpp(pars, max_t, max_N, max_tries)
  if (is.null(sim) || is.null(sim$tree)) {
    message("Simulation failed.")
    return(NULL)
  }
  phy <- sim_matrix_to_phylo(sim$tree)
  plot(phy, main = paste("Simulated tree (status:", sim$status, ")"))
  return(phy)
}

# Minimal test runner
run_safe_cpp_test <- function(n = 5) {
  cat("Running safe C++ simulation test...\n")
  sapply(seq_len(n), function(i) safe_simulate_div_tree_cpp())
}
