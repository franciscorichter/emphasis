## GAM estimation for model 

# optim with GAM min. (Maybe DE)
# install.packages("svMisc")
library(svMisc)

# Parameters for generatePhyloPD
n_trees_desired <- 10000

mu_interval <- c(0, 0.5)
lambda_interval <- c(0.5, 3)
betaN_interval <- c(-0.01, 0.01) 
betaP_interval <- c(-0.01, 0.01) 

#load("~/Library/CloudStorage/Dropbox/Pancho/51 - emphasisLD/data/FamilyAllTrees.RData")

#phylo = FamilyAllTrees$Nectariniidae$tree

phylo = FamilyBirdTrees$Nectariniidae$tree
length(ape:::branching.times(phylo))

# Call the function
results <- AugmentMultiplePhyloPD(phylo = phylo,
                                  n_trees = n_trees_desired,
                                  mu_interval = mu_interval,
                                  lambda_interval = lambda_interval,
                                  betaN_interval = betaN_interval,
                                  betaP_interval =betaP_interval,
                                  max_lin = 10000,
                                  max_tries = 10)

trees_from_results <- results$trees
param_from_results <- results$param
rejected_overruns_from_results <- results$rejected_overruns
rejected_lambda_from_results <- results$rejected_lambda
rejected_zero_weights_from_results <- results$rejected_zero_weights
times_from_results <- results$times
loglik_estimation_from_results <- as.numeric(unlist(results$loglik_estimation))
logf_from_results <- as.numeric(unlist(results$logf))
logg_from_results <- as.numeric(unlist(results$logg))


# Further analysis can be done using the 'results' list.
dat3 = data.frame(srv=loglik_estimation_from_results[!is.na(loglik_estimation_from_results)],
                  p1=param_from_results$mu,
                  p2=param_from_results$lambda,
                  p3=param_from_results$betaN,
                  p4=param_from_results$betaP)

gam_loglik_pd = mgcv::gam(srv~s(p1)+s(p3)+s(p2)+s(p4),data=dat3)


summary(gam_loglik_pd)


vals = data.frame(p1 = 0.1,p2 = 0.5,p3 = -0.01,p4 = 0)
predict(gam_loglik_pd,vals)

#pars2 = c(1000,1,0,1,1,2)
#DDD:::dd_loglik(pars1 = c(0.5,0.1,0.4/0.01),
#                pars2 = pars2,
#                brts = ape:::branching.times(phylo),
#                missnumspec = 0)

lambda_sample <- runif(1000, lambda_interval[1], lambda_interval[2])
mu_sample <- runif(1000, mu_interval[1], mu_interval[2])
betaN_sample <- runif(1000, betaN_interval[1], betaN_interval[2])
betaP_sample <- runif(1000, betaP_interval[1], betaP_interval[2])
sim.param <- data.frame(p1 = mu_sample, 
                        p2 = lambda_sample,
                        p3 = betaN_sample,
                        p4 = betaP_sample)

pred =predict(gam_loglik_pd,sim.param)

