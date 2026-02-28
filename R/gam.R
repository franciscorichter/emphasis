#' Train a Generalized Additive Model on Simulation Data
#'
#' This function fits a Generalized Additive Model (GAM) to the provided simulation
#' results. It models the log-likelihood of simulation success (`loglik`) as a function
#' of the simulation parameters, using smooth terms for each parameter. Only simulations
#' that completed successfully are used for the model fitting.
#'
#' @param results A list containing simulation results with the following components:
#'   - `loglik_estimation`: A list of log-likelihood values or `try-error` objects.
#'   - `trees`: A list of data frames, each containing tree data.
#'   - `param`: A list of named vectors containing the simulation parameters `mu`,
#'     `lambda`, `betaN`, and `betaP`.
#'
#' @return A `gam` object representing the fitted model.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'results' is your list of simulation results:
#' gam_model <- train_GAM(results)
#' summary(gam_model)
#' }
#'
train_GAM <- function(results){

  completed_indices <- which(!sapply(results$loglik_estimation, inherits, "try-error"))
  trees = results$trees[completed_indices]
  num_extinct_per_tree <- sapply(trees, count_extinct_species)

  sim_data <- data.frame(mu = results$param$mu,
                         lambda = results$param$lambda,
                         betaN = results$param$betaN,
                         betaP = results$param$betaP,
                         loglik = unlist(results$loglik_estimation[completed_indices]),
                         nspecies =  num_extinct_per_tree
  )

  cat("Training GAM...")
  gam_loglik = mgcv::gam(loglik~s(mu)+s(lambda)+s(betaN)+s(betaP),data=sim_data,family= mgcv::scat(link="identity"))

  pvals <- summary(gam_loglik)$p.values
  significant <- pvals < 0.05
  if(all(significant)) {
    cat("All parameters are significant.\n")
  } else {
    cat("Not all parameters are significant.\n")
  }

  return(gam_loglik)
}

#' @keywords internal
count_extinct_species <- function(tree) {
  extinction_times = as.data.frame(tree)["t_ext"]
  sum(extinction_times < 1e+10)
}
