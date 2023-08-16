# Set fixed parameters
# Example: mu = 0.5, lambda = 0.6, betaN = 0.1, betaP = 0.2
rm(list=ls())
set.seed(123)
rm(list=ls())
sim.param <- c(0.1, 0.8, -0.0125, 0.005)

# Generate 10 trees using the fixed parameters
results <- generatePhyloPD_fixedParam(n_trees = 100, sim.param = sim.param,max_t = 20)


mat = NULL
for(i in 1:length(results$trees)){
  coords = ape::ltt.plot.coords(results$trees[[i]])
  df = cbind(coords,data.frame(tree=i))
  mat = rbind(mat,df)
}

ggplot(mat) + geom_line(aes(x=time,y=N,as.factor(tree)),alpha=0.5) +theme_minimal()


