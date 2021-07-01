####  Sample fulll trees 

sample_tree_full <- function(diversification_model,ct,nhpp_step=1){  
  # test1Ok
  cbt  = 1e-10
  tree = list(extant=data.frame(brts=c(0,0),
                                parent=c(1,2),
                                child=c(2,3),
                                clade=c(0,1)),
              extinct=data.frame(brts=NULL,
                                 parent=NULL,
                                 child=NULL,
                                 t_ext=NULL),
              ct=ct)
  next_bt = 0 
  while(cbt < ct & sum(tree$extant$clade==1)>0 & sum(tree$extant$clade==0)>0){
    
    next_bt = ifelse((ct-next_bt)>2,(next_bt+ct)/2,ct)
    
    ### Draw speciation 
    event_time = draw_event_time(cbt,
                                 next_bt,
                                 diversification_model,
                                 tree=tree)
    if(event_time<next_bt){
      
      ## resolve allocation and update tree
      allocation = draw_allocation_full(event_time,
                                        diversification_model,
                                        tree)
      row_extant = which(tree$extant$child==allocation$species)
      if(allocation$event=="e"){
        to_add = tree$extant[tree$extant$child==allocation$species,]
        tree$extinct = rbind(tree$extinct,data.frame(brts=to_add$brts,
                                                     parent=to_add$parent,
                                                     child=to_add$child,
                                                     t_ext=event_time))
        tree$extant = tree$extant[-row_extant,]
      }
      
      if(allocation$event=="s"){
        next_child = max(tree$extant$child,tree$extinct$child)+1
        to_add = data.frame(brts=event_time,parent=allocation$species,child=next_child,clade=tree$extant$clade[row_extant])
        tree$extant = rbind(tree$extant,to_add)
      }
      
    }
    cbt = min(event_time, next_bt)
  }
  if(sum(tree$extant$clade==1)>0 & sum(tree$extant$clade==0)>0){
    tree$extant$clade=NULL
    class(tree)="etree"
  }else{
    tree = NULL
  }
  return(tree)
}

draw_allocation_full <- function(event_time,
                                 diversification_model,
                                 tree){
  
  # calculate rates
  l = speciation_rate(tm = event_time,
                      tree = tree,
                      diversification_model = diversification_model)
  m = extinction_rate(tm = event_time,
                      tree = tree,
                      diversification_model = diversification_model)
  mu = sum(m)
  lambda = sum(l)
  ## choose speciation or extinction 
  to = sample(c("e","s"),1,prob = c(mu/(mu+lambda),lambda/(mu+lambda))) 
  
  ## choose species
  current_species <- get_current_species(tm = event_time,
                                         tree = tree)
  if(to=="e") probs = m
  if(to=="s") probs = l

  species = sample(current_species,prob = probs,size=1)
  
  return(list(event=to,species=species))
}


draw_event_time <- function(cbt,
                            next_bt,
                            diversification_model,
                            tree){
  
  
  nsr = sum_of_rates
  
  key = 0 
  while(key == 0 & cbt < next_bt){
    lambda_max = optim(cbt,
                       fn = nsr,
                       tree = tree,
                       diversification_model = diversification_model,
                       lower = cbt,
                       upper = next_bt,
                       method ="L-BFGS-B",
                       control=list(fnscale=-1,maxit=1))$value
    
    u1 = runif(1)
    if(lambda_max==0){
      cbt = Inf
    }else{
      cbt = cbt - log(x = u1)/lambda_max
    }
    
    if(cbt < next_bt){
      u2 = runif(1)
      
      pt = nsr(cbt,
               tree,
               diversification_model)/lambda_max
      
      if(u2<pt){
        key = 1
      }
    }
  }
  
  return(cbt)
  
}
