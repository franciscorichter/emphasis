library(emphasis)
pars = c(0.1,0.8,-0.0175,0)
max_t = 20

pars = pars*max_t
max_lin = 1000000
max_tries = 1

#sim_single_tree_pd_cpp
tree = sim_tree_pd_cpp(pars = pars,
                       max_t = 5,
                       max_lin = max_lin,
                       max_tries = max_tries)

# this is not for official testing:
if (requireNamespace("treestats")){
  treestats::crown_age(tree$tes) # 5
  resc_tree <- emphasis::rescale_tree(tree$tes, new_crown_age = 2)
  treestats::crown_age(resc_tree) # 2
  
  resc_tree <- emphasis::rescale_tree(tree$tes, new_crown_age = 1)
  treestats::crown_age(resc_tree) # 1
}