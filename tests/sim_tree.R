library(emphasis)
#pars = c(0.1,0.8,-0.0175,0)
pars = c(0.24624818,  0.81232484, -0.01140782,  0.00000000)
max_t = 20

max_lin = 100000
max_tries = 1

#sim_single_tree_pd_cpp
tree=sim.tree(pars = pars*max_t,
                       max_t = 1,
                       max_lin = max_lin,
                       max_tries = max_tries)

tree

sim.tree <- function(pars,
                     max_t = 1,
                     max_lin,
                     max_tries){
  it <- try(sim_tree_pd_cpp(pars = pars*max_t,
                            max_t = 1,
                            max_lin = max_lin,
                            max_tries = max_tries),silent = TRUE)
  if(class(it)!="try-error"){
    result = it
  }else{
    result = NaN
  }
  return(result)
}
