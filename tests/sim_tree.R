library(emphasis)
pars = c(0.1,0.8,-0.0175,0)
max_t = 20

pars = pars*max_t
max_lin = 1000000
max_tries = 1

#sim_single_tree_pd_cpp
tree=simulate_single_pd_tree_cpp(pars = pars,
                       max_t = 1,
                       max_lin = max_lin,
                       max_tries = max_tries)
