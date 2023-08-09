#   Pipeline for full emphasis

## Input 

rm(list=ls())
#source("~/Dropbox/Pancho/05 - SystematicBiology (paper 2)/DEphasis/scripts/functions.R")
load("~/Dropbox/Pancho/52 - EmphasisComplete/06-data/Phylogenies/FamilyAllTrees.RData")
#theme_set(theme_bw())
#l = load_brts()

clade = "Megapodiidae"
brts = ape::branching.times(FamilyAllTrees$Megapodiidae$tree)

num_iterations = 100
num_points = 100
max_missing = 2000
lower_bound = c(0,0.5,-0.05,-0.05)
upper_bound = c(0.5,1,0.05,0.05)
maxN = 100
max_lambda = 1000

disc_prop = 0.75
sd_vec = c(0.5,1,0.1,0.1)
#alpha = sd_vec/n_it



### Prove how to use EM algorithm 


## How to use DE algorithm 

est_de = emphasis_de(brts = brts,
            num_iterations = num_iterations,
            num_points = num_points,
            max_missing = max_missing,
            sd_vec = sd_vec,
            lower_bound = lower_bound,
            upper_bound = upper_bound,maxN = maxN,
            max_lambda = max_lambda,
            disc_prop = disc_prop,
            verbose = TRUE)

#### How to vizualize results 

viz = viz_de(est_de)

### Prove how to calculate the GAM version 


### Use Generalized 