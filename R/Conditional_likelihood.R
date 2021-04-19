# load package to fit gam
#library(mgcv)
# Simulate survival
create_grid <- function(llim,rlim,n.grid){
  p = length(llim)
  theta.range = cbind(llim,rlim)
  pars = NULL
  for (i in 1:p){
    v = rep(rep(seq(theta.range[i,1],
                    theta.range[i,2],
                    length.out = n.grid),
                each=n.grid^(p-i)),n.grid^(i-1))
    pars = cbind(pars,v)
  }
  return(pars)
}


Simulation_step <- function(grid,model,ct,timeLimit){
  srv = vector(mode="numeric",length=nrow(grid))
  Trees = vector(mode="list",length=nrow(grid))
  for (i in 1:nrow(grid)){
    svMisc:::progress(i,max.value = nrow(grid))
    tau = try(emphasis::sim_survival(diversification_model = list(pars=grid[i,],model=model),
                           ct=ct,
                           timeLimit = timeLimit),silent = TRUE)
    if(class(tau)=="try-error"){
      srv[i] = -1
      Trees[[i]] = "error"
    }else{
      srv[i] = tau$srv
      Trees[[i]] = tau$tree
    }
    
  }
  return(list(srv = srv, Trees = Trees))
}

sim_survival <- function (diversification_model, ct, timeLimit=10){
  pars = diversification_model$pars
  model = diversification_model$model
  if(model == "rpd1"){
    tree = simTree_dd(pars,ct,timeLimit = timeLimit)
  }
  if(model == "rpd5"){
    tree = simTree_pd(pars,ct,timeLimit = timeLimit)
  }
  if (length(tree) == 1) {
    srv = 0
  }else{
    srv =1
  }
  return(list(tree=tree,srv=srv))
}


fit_gam_survival <- function(simulations,splines="bivariate"){
  
  if(splines=="bivariate"){
    srv.gam = mgcv::gam(srv~s(p1,p2)+s(p1,p3)+s(p2,p3),family = binomial,data=simulations)
  }
  if(splines=="univariate"){
    srv.gam = mgcv::gam(srv~s(p1)+s(p3)+s(p2),family = binomial,data=simulations)
  }

  return(srv.gam)
}


  