## GAM estimation for model 

# optim with GAM min. (Maybe DE)
# install.packages("svMisc")
library(svMisc)
library(emphasis)
# Parameters for generatePhyloPD
time = proc.time()

n_trees <- 200000

mu_interval <- c(0, 0.5)
lambda_interval <- c(0.5, 3)
betaN_interval <- c(-0.02, 0) 
betaP_interval <- c(-0.02, 0.02) 
max_lin = 10000
max_tries = 10

load("~/Library/CloudStorage/Dropbox/Pancho/51 - emphasisLD/data/FamilyBirdTrees.Rdata")
j=101
phylo = FamilyBirdTrees[[j]]$tree

#DDD:::dd_ML(brts=ape::branching.times(phylo))
#plot(phylo)

# Call the function
results10 <- AugmentMultiplePhyloPD(phylo = phylo,
                                  n_trees = n_trees,
                                  mu_interval = mu_interval,
                                  lambda_interval = lambda_interval,
                                  betaN_interval = betaN_interval,
                                  betaP_interval =betaP_interval,
                                  max_lin = max_lin,
                                  max_tries = max_tries)

sim_time10 = get.time(time,mode = "min")

trees_from_results <- results$trees
param_from_results <- results$param
rejected_overruns_from_results <- unlist(results$rejected_overruns)


rejected_lambda_from_results <- unlist(results$rejected_lambda)

rejected_zero_weights_from_results <- unlist(results$rejected_zero_weights)


times_from_results <- results$times
loglik_estimation_from_results <- as.numeric(unlist(results$loglik_estimation))
logf_from_results <- as.numeric(unlist(results$logf))
logg_from_results <- as.numeric(unlist(results$logg))


# Further analysis can be done using the 'results' list.
wi = !is.na(loglik_estimation_from_results)
dat8 = data.frame(srv=loglik_estimation_from_results[wi],
                  p1=param_from_results$mu[wi],
                  p2=param_from_results$lambda[wi],
                  p3=param_from_results$betaN[wi],
                  p4=param_from_results$betaP[wi],
                  rzw = rejected_zero_weights_from_results[wi],
                  rl = rejected_lambda_from_results[wi],
                  ro = rejected_overruns_from_results[wi])

dat_all <- rbind(dat3,data.frame(srv="NaN",p1=param_from_results$mu[!wi],
                                 p2=param_from_results$lambda[!wi],
                                 p3=param_from_results$betaN[!wi],
                                 p4=param_from_results$betaP[!wi]))

library(ggplot2)
ggplot(dat5) + geom_point(aes(x=p1,y=p2,colour=srv)) + theme_minimal()
ggplot(dat5) + geom_point(aes(x=p3,y=p4,colour=srv)) + theme_minimal()

ggplot(dat3) + geom_point(aes(x=p1,y=p2,colour=rzw)) + theme_minimal()
ggplot(dat3) + geom_point(aes(x=p3,y=p4,colour=rzw)) + theme_minimal()

ggplot(dat3) + geom_point(aes(x=p1,y=p2,colour=rl)) + theme_minimal()
ggplot(dat3) + geom_point(aes(x=p3,y=p4,colour=rl)) + theme_minimal()

ggplot(dat3) + geom_point(aes(x=p1,y=p2,colour=ro)) + theme_minimal()
ggplot(dat3) + geom_point(aes(x=p3,y=p4,colour=ro)) + theme_minimal()


gam_loglik_pd = mgcv::gam(srv~s(p1)+s(p3)+s(p2)+s(p4),data=dat5)

summary(gam_loglik_pd)

gam_loglik_second_order = mgcv::gam(srv~s(p1)+s(p3)+s(p2)+s(p4)+
                            s(p1,p2)+s(p1,p3)+s(p1,p4)+
                            s(p2,p3)+s(p2,p4)+s(p3,p4),data=dat8)


summary(gam_loglik_second_order)


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

