 ### DE estimation for PD 
library(emphasis)

load("100simPD.RData")
trees = results$trees
BRTS = results$brts


num_iterations = 1000
num_points = 100
max_missing = 20000
lower_bound = c(0,0.5,-0.05,-0.05)
upper_bound = c(0.5,2,0.01,0.03)
maxN = 100
max_lambda = 1000

disc_prop = 0.75
sd_vec = c(0.1,0.5,0.01,0.01)

ESTIMATIONS = list()
for(i in 1:length(trees)){
  brts = BRTS[[i]]
  est_de = emphasis_de(brts = brts,
                       num_iterations = num_iterations,
                       num_points = num_points,
                       max_missing = max_missing,
                       sd_vec = sd_vec,
                       lower_bound = lower_bound,
                       upper_bound = upper_bound,
                       maxN = maxN,
                       max_lambda = max_lambda,
                       disc_prop = disc_prop,
                       verbose = TRUE)
  ESTIMATIONS[[i]] = est_de
  save.image(file=paste0("tree",i,".RData"))
}