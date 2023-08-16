#   Pipeline for full emphasis
library(emphasis)
## Input 

rm(list=ls())
#source("~/Dropbox/Pancho/05 - SystematicBiology (paper 2)/DEphasis/scripts/functions.R")
load("~/Dropbox/Pancho/52 - EmphasisComplete/06-data/Phylogenies/FamilyAllTrees.RData")
#theme_set(theme_bw())
#l = load_brts()

clade = "Megapodiidae"
brts = ape::branching.times(FamilyAllTrees$Megapodiidae$tree)

num_iterations = 50000
num_points = 100
max_missing = 20000
lower_bound = c(0,0.5,-0.05,-0.05)
upper_bound = c(0.5,2,0.01,0.03)
maxN = 100
max_lambda = 1000

disc_prop = 0.75
sd_vec = c(0.1,0.1,0.01,0.01)
#alpha = sd_vec/n_it



### Prove how to use EM algorithm 


## How to use DE algorithm 



#save.image(file="results9Aug.RData")
#### How to vizualize results 
library(tidyr)
library(purrr)
library(dplyr)
library(ggplot2)
viz = viz_de(est_de)

### Prove how to calculate the GAM version 


### Use Generalized 