## GAM estimation for model 


## leanrning 

# optim with GAM min. (Maybe DE)
# install.packages("svMisc")
library(svMisc)

# Parameters for generatePhyloPD
n_trees_desired <- 100


mu_interval <- c(0, 0.5)
lambda_interval <- c(0.5, 3)
betaN_interval <- c(-0.01, 0.01) 
betaP_interval <- c(-0.001, 0.001) 

load("~/Library/CloudStorage/Dropbox/Pancho/51 - emphasisLD/data/FamilyAllTrees.RData")

phylo = FamilyAllTrees$Nectariniidae$tree
length(ape:::branching.times(phylo))

# Call the function
results <- AugmentMultiplePhyloPD(phylo = phylo,
                                  n_trees = 100,
                                  mu_interval = mu_interval,
                                  lambda_interval = lambda_interval,
                                  betaN_interval = betaN_interval,
                                  betaP_interval =betaP_interval,
                                  max_lin = 100000,
                                  max_tries = 10)

trees_from_results <- results$trees
param_from_results <- results$param
rejected_overruns_from_results <- results$rejected_overruns
rejected_lambda_from_results <- results$rejected_lambda
rejected_zero_weights_from_results <- results$rejected_zero_weights
times_from_results <- results$times
loglik_estimation_from_results <- as.numeric(unlist(results$loglik_estimation))
logf_from_results <- results$logf
logg_from_results <- results$logg


# Further analysis can be done using the 'results' list.
dat3 = data.frame(srv=loglik_estimation_from_results[!is.na(loglik_estimation_from_results)],
                  p1=param_from_results$mu,
                  p2=param_from_results$lambda,
                  p3=param_from_results$betaN,
                  p4=param_from_results$betaP)

gam_loglik_pd = mgcv::gam(srv~s(p1)+s(p3)+s(p2)+s(p4),data=dat3)


