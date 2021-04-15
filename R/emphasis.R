#' Emphasis main function
#' @description Main function of emphasis, that uses an E-M approach to fit a
#' diversification model to a phylogenetic tree.
#' @param brts vector of branching times of the tree for which the model has to
#' be fitted
#' @param init_par initial parameter values of the model
#' @param soc number of species at the root (1) or crown (2). Default is 2.
#' @param model model to be used
#' @param lower_bound vector of the lower limit of parameter values used by the
#' model. Set to -Infinity if left empty.
#' @param upper_bound vector of the upper limit of parameter values used by the
#' model. Set to +Infinity if left empty
#' @param max_lambda maximum speciation rate, default is 500. Should not be set 
#' too high to avoid extremely long run times
#' @param xtol tolerance of step size in the M step
#' @param em_tol tolerance of step size in cycling through EM
#' @param sample_size_tol tolerance in determining the sample size
#' @param verbose if TRUE, provides textual output of intermediate steps 
#' @param return_trees boolean, if TRUE the simulated trees are returned as well
#' @param max_missing maximum number of tips a tree can be augmented with.
#' @param burnin_sample_size sample size during burn-in
#' @param pilot_sample_size vector of sample sizes used to determine the true 
#' sampling size
#' @param burnin_iterations number of iterations of the EM algorithm to discard
#' as burn-in
#' @param num_threads number of threads to be used. If set to 0, the maximum 
#' number of threads available is chosen. 
#' @param conditional a function that takes a parameter set as argument and returns
#' conditional probability. 
#' @export
#' @return a list with two components: 1) \code{pars} contains the average parameter
#' estimate and 2) \code{MCEM} matrix of parameter estimates and likelihoods.
emphasis <- function(brts,
                     init_par,
                     soc = 2,
                     model,
                     lower_bound = numeric(0),
                     upper_bound = numeric(0),
                     max_lambda = 500,
                     xtol = 0.001,
                     em_tol = 0.25,
                     sample_size_tol = 0.005,
                     verbose = FALSE,
                     return_trees = FALSE,
                     max_missing = 10000,  # maximum tree size
                     burnin_sample_size = 200,
                     pilot_sample_size = seq(100, 1000, by = 100),
                     burnin_iterations = 20,
                     num_threads = 0,
                     conditional = NULL) {
  
  if (NULL != conditional) {
    stopifnot(class(conditional) == "function")
  }
  if (class(brts) == "phylo") {
    cat("You have provided the full phylogeny instead of the branching times\n")
    cat("Emphasis will extract the branching times for your convenience\n")
    brts <- ape::branching.times(brts)
  }
  
  msg1 <- paste("Initializing emphasis...")
  msg2 <- paste("Age of the tree: ", max(brts))
  msg3 <- paste("Number of speciations: ", length(brts))
  msg4 <- paste("Diversification model to fit:", model)
  msg5 <- "######################################"
  cat(msg1, msg2, msg3, msg4, msg5, sep = "\n")
  
  cat("Performing Phase 1: burn-in", sep = "\n")
  mc <- mcEM_step(brts = brts,
                  pars = init_par,
                  sample_size = burnin_sample_size,
                  model = model,
                  soc = soc,
                  max_missing = max_missing,
                  max_lambda = max_lambda,
                  lower_bound = lower_bound,
                  upper_bound = upper_bound,
                  xtol = xtol,
                  num_threads = num_threads,
                  return_trees = FALSE,
                  verbose = FALSE,
                  tol = em_tol,
                  burnin = burnin_iterations,
                  conditional)
  
  M <- mc$mcem
  pars <- c(mean(utils::tail(M$par1, n = nrow(M) / 2)),
            mean(utils::tail(M$par2, n = nrow(M) / 2)),
            mean(utils::tail(M$par3, n = nrow(M) / 2)),
            mean(utils::tail(M$par4, n = nrow(M) / 2)))
  
  cat("\n", msg5, sep = "\n")
  cat("Phase 2: Assesing required MC sampling size \n")
  
  for (s in pilot_sample_size) {
    cat(paste("\n Sampling size: ", as.character(s), "\n"))
    mc <- mcEM_step(brts = brts,
                    pars = pars,
                    sample_size = s,
                    model = model,
                    soc = soc,
                    max_missing = max_missing,
                    max_lambda = max_lambda,
                    lower_bound = lower_bound,
                    upper_bound = upper_bound,
                    xtol = xtol,
                    num_threads = num_threads,
                    return_trees = FALSE,
                    verbose = FALSE,
                    tol = em_tol,
                    burnin = 10,
                    conditional)
    
    ta <- tail(mc$mcem, n = nrow(M) / 2)
    pars <- c(mean(ta$par1),
              mean(ta$par2),
              mean(ta$par3),
              mean(ta$par4))
    M <- rbind(M, mc$mcem)
  }
  n.r <- get_required_sampling_size(M[ -(1:burnin_iterations), ],
                                    tol = sample_size_tol)
  sample_size <- max(pilot_sample_size + 2, n.r)
  n.r_old <- -1
  j <- 1
  while (n.r_old < n.r) {
    msg6 <- paste0("Required sampling size: ", n.r)
    msg7 <- paste0("Phase 3: Performing metaiteration: ", j)
    cat("\n", msg5, msg7, msg6, sep = "\n")
    mc <- mcEM_step(brts = brts,
                    pars = pars,
                    sample_size = sample_size,
                    model = model,
                    soc = soc,
                    max_missing = max_missing,
                    max_lambda = max_lambda,
                    lower_bound = lower_bound,
                    upper_bound = upper_bound,
                    xtol = xtol,
                    num_threads = num_threads,
                    return_trees = FALSE,
                    verbose = FALSE,
                    tol = em_tol,
                    burnin = 2,
                    conditional)
    
    M <- rbind(M, mc$mcem)
    n.r_old <- n.r
    j <- j + 1
    n.r <- get_required_sampling_size(M[ -(1:burnin_iterations), ],
                                      tol = sample_size_tol)
    pars <- as.numeric(colMeans(mc$mcem)[1:4])
    sample_size <- n.r
  }
  
  cat(pars)
  return(list(pars = pars, MCEM = M))
}

#' @keywords internal
#' this is an internal function
mcEM_step <- function(brts,
                      pars,
                      sample_size,
                      model,
                      soc,
                      max_missing,
                      max_lambda,
                      lower_bound,
                      upper_bound,
                      xtol, # tolerance in M step
                      tol,   #tolerance of mcEM step
                      burnin,
                      num_threads,
                      return_trees,
                      verbose,
                      conditional) {
  mcem <- NULL
  sde <- 10
  i <- 0
  times <- NULL
  while (sde > tol) {
    i <- i + 1
    results <- remphasis::em_cpp(brts,
                                 pars,
                                 sample_size,
                                 maxN = 10 * sample_size,                   
                                 locate_plugin(model),           
                                 soc,
                                 max_missing,           
                                 max_lambda,           
                                 lower_bound,
                                 upper_bound,
                                 xtol_rel = xtol,                   
                                 num_threads,
                                 return_trees,
                                 conditional)
    pars <- results$estimates
    mcem <- rbind(mcem, data.frame(par1 = pars[1],
                                   par2 = pars[2],
                                   par3 = pars[3],
                                   par4 = pars[4],
                                   fhat = results$fhat,
                                   sample_size = sample_size))
    
    if (verbose) {
      print(paste("(mean of) loglikelihood estimation: ", mean(mcem$fhat)))
    }
    
    times <- c(times, mcem$time)
    time_p_it <- mean(times)
    
    if (i > burnin) {
      mcem_est <- mcem[floor(nrow(mcem) / 2):nrow(mcem), ]
      mcem_est <- mcem_est[is.finite(mcem_est$fhat), ]
      #  sde0 = sde
      sde <- stats::sd(mcem_est$fhat) / sqrt(nrow(mcem_est))
      #  mde = mean(mcem_est$fhat)
      msg <- paste("Iteration:", i, " SE of the loglikelihood: ", sde)
      cat("\r", msg)
    } else {
      msg <- paste("Remaining time (burn-in): ",
                   round(time_p_it * (burnin - i), digits = 0), "sec")
      cat("\r", msg)
    }
  }
  return(list(mcem = mcem))
}

#' @keywords internal
#' this is an internal function
get_required_sampling_size <- function(M, tol = .05) {
  hlp <- MASS::rlm(M$fhat ~ I(1 / M$sample_size), weights = M$sample_size)
  ab <- stats::coef(hlp)
  f_r <- ab[1] - tol
  n_r <- ceiling(ab[2] / (f_r - ab[1]))
  return(n_r)
}
