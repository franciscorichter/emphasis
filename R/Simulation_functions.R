
### simulation of trees 

simTree_dd <- function(pars,ct,timeLimit){
  setTimeLimit(timeLimit)
  tree =  data.frame(brts=c(0,0),
                     to=c(1,1),
                     t_ext=c(Inf,Inf),
                     lambda=c(0,0), 
                     clade = c("a","b"))
  cbt = 0
  N1 = 1 
  N2 = 1
  mu = max(0,pars[1])
  while((cbt < ct)  &  (N1 >= 1) &  (N2 >= 1)){
    N = N1 + N2 
    lambda_ct = max(0,pars[2] + pars[3]*N)  # diversity dependence only 
    rate_max = (lambda_ct+mu)*N
    u1 = runif(1)
    next_event_time = cbt-log(x = u1)/rate_max
    
    if(next_event_time < ct){

        to = sample(c(1,0),size=1,prob=c(lambda_ct,mu)/(lambda_ct+mu))
        clade = sample(c("a","b"),size=1,prob=c(N1,N2)/(N1+N2))
        if(to == 1){
          if(clade=="a"){
            N1 = N1+1
          }else{
            N2 = N2+1
          }
        }else{
          if(clade=="a"){
            N1 = N1-1
          }else{
            N2 = N2-1
          } 
          if((N1 >= 1) &  (N2 >= 1)){
            ext_spec = sample(which(tree$to == 1 & tree$t_ext == Inf & tree$clade == clade),1)
            tree$t_ext[ext_spec] = next_event_time
          }
        }
        tree = rbind(tree,data.frame(brts=next_event_time,to=to,t_ext=Inf,lambda=lambda_ct, clade = clade))
    }
    cbt = next_event_time
  }

  if((N1 == 0) |  (N2 == 0)){
    tree = cbt
  }else{
    tree = rbind(tree,data.frame(brts=ct,to=1,t_ext=Inf,lambda=lambda_ct, clade=0))
  }
  return(tree)
}


simTree_pd <- function(pars,ct,timeLimit){
  setTimeLimit(timeLimit)
  tree =  data.frame(brts=c(0,0),
                     to=c(1,1),
                     t_ext=c(Inf,Inf),
                     clade = c("a","b"))
  cbt = 0
  N1 = 1 
  N2 = 1
  mu = max(0,pars[1])
  P=0
  while((cbt < ct)  &  (N1 >= 1) &  (N2 >= 1)){
    next_bt = min(c(tree$brts[tree$brts>cbt],ct))
    N = N1+N2
    lambda_mx = max(0,pars[2] + pars[3]*N  +  ((P + N*(next_bt-cbt)-cbt)/N )*pars[4])
    rate_max = (lambda_mx+mu)*N 
    u1 = runif(1)
    next_event_time = cbt-log(x = u1)/rate_max
    P = P + N*(next_event_time-cbt)
    
    if(next_event_time < next_bt){
      u2 = runif(1)
      lambda_ct = max(0,pars[2] + pars[3]*N  +  ((P+N*(next_event_time-cbt) - next_event_time)/N)*pars[4])
      pt =( (lambda_ct+mu)*N )/rate_max  
      
      if(u2<pt){
        to = sample(c(1,0),size=1,prob=c(lambda_ct,mu)/(lambda_ct+mu))
        clade = sample(c("a","b"),size=1,prob=c(N1,N2)/(N1+N2))
        if(to == 1){
          if(clade=="a"){
            N1 = N1+1
          }else{
            N2 = N2+1
          }
        }else{
          if(clade=="a"){
            N1 = N1-1
          }else{
            N2 = N2-1
          }
          if((N1 >= 1) &  (N2 >= 1)){
            ext_spec = sample(which(tree$to == 1 & tree$t_ext == Inf & tree$clade == clade),1)
            tree$t_ext[ext_spec] = next_event_time
            removed_branch_length = tree$t_ext[ext_spec] - tree$brts[ext_spec]
            P = P - removed_branch_length
          }
        }
            
        tree = rbind(tree,data.frame(brts=next_event_time,to=to,t_ext=Inf,clade=clade))
        
      }
    }
    
    cbt = next_event_time
  }
  if((N1 == 0) |  (N2 == 0)){
    tree = cbt
  }else{
    tree = rbind(tree,data.frame(brts=ct,to=1,t_ext=Inf, clade=0))
  }
  return(tree)
}
