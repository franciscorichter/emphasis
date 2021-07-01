get_extant <- function(tm,tree){
  origin = setdiff(c(tree$extant$parent,tree$extinct$parent),c(tree$extant$child,tree$extinct$child))
  # obtain extant tree from full tree
  if (nrow(tree$extinct)>0){
    if (tm==0){tm<-1e-10}
    tree$extant$clade<-NULL
    extinct = tree$extinct[tree$extinct$brts<tm & tree$extinct$t_ext<tm,]
    extinct = extinct[order(extinct$brts),]
    extant = rbind(tree$extant[tree$extant$brts<tm,],
                   tree$extinct[tree$extinct$brts<tm & tree$extinct$t_ext>=tm,-4])
    extant = extant[order(extant$brts),]
    extinct$t_ext = NULL
    if (nrow(extinct)>0){
      i<-1
      while (i <= nrow(extant)){
        if ((sum(extant$parent[i]==extant$child)==0) & (extant$parent[i]!=origin)){
          ind = which(extant$parent[i]==extinct$child)
          ind2 = which(extinct$parent==extant$parent[i] & extinct$brts < extant$brts[i])
          child = extant$child[i]
          new.extant.row<-extinct[ind,]
          new.extant.row[3]<-child
          new.extinct.row <- extant[i,c(1,3,2)]
          extant[i,]<-new.extant.row
          extinct[ind,]<-new.extinct.row
          extinct$parent[ind2]<- child
        } else {
          i <- i+1
        }
      }
      extant$parent[1]<-origin
      extant$parent[2]<-extant$child[1]
    } 
  } else {
    extant = tree$extant
  }
  extant = extant[order(extant$brts),]
  if (extant$parent[2]==origin){
    extant[c(1,2),]=extant[c(2,1),]
  }
  extant = extant[extant$brts<=tm,]
  return(extant)
}




transf <- function(name_spe,vec){
  which(vec==name_spe)
}

newick<- function(tree,CT){
  n<-nrow(tree)
  child.nms<-as.character(tree$child)
  parent.nms<-as.character(tree$parent)
  species.nms<-unique(child.nms,parent.nms)
  n.species<-length(species.nms)
  CT<-rep(CT,n.species)
  for (i in seq(n,1)){
    nw<-paste("(",parent.nms[i],":",as.character(CT[which(species.nms==parent.nms[i])]-tree$brts[i]),",",child.nms[i],":",as.character(CT[which(species.nms==child.nms[i])]-tree$brts[i]),")", sep = "")
    j<-which(parent.nms[i]==child.nms)
    rp<-which(parent.nms==child.nms[j])
    if (length(rp)>0){
      parent.nms[rp]<-nw
    }
    species.nms[which(species.nms==child.nms[j])]<-nw
    child.nms[j]<-nw
    CT[j]<-CT[j]-(CT[which(species.nms==parent.nms[i])]-tree$brts[i])
  }
  return(paste(child.nms[1],";",sep=""))
  #return(child.nms)
}

PDTT_plot <- function(tree){
  ct = max(tree$brts)
  times = seq(0,ct,length.out = 1000)
  times = sort(times,decreasing = T)
  PD=NULL
  for(i in 1:length(times)){
    G=GPD(times[i],tree[-1,])
    
    PD=rbind(PD,data.frame(time=rep(times[i],nrow(G)),
                         P=colSums(G)/(nrow(G)-1),
                         lineage=as.character(1:nrow(G))))
  }
  g1 = ggplot(PD)+
    geom_line(aes(x=time,y=P,colour=lineage,alpha=0.5))+
    theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
return(g1)
}

etree2phylo <- function(etree){
  ext = get_extant(tm = etree$ct, tree = etree)
  nw = newick(ext,CT = etree$ct)
  tr = ape:::read.tree(text=nw)
  return(tr)
}

phylo2etree <- function(phylo){
  #transformation of ultrametric trees into data frame
  tree = DDD::phylo2L(phylo)
  brts_dd = tree[,1]
  brts = cumsum(-diff(c(brts_dd,0)))
  
  tree = list(extant = data.frame(brts = c(0,brts[-length(brts)]),
                                  parent=c(1,abs(tree[,2][-1])),
                                  child=abs(tree[,3])),
              extinct = data.frame(brts = numeric(),
                                   parent = numeric(),
                                   child = numeric(),
                                   t_ext = numeric()),
              
              ct=brts_dd[1])
  tree$extant[3:nrow(tree$extant),"parent"]=tree$extant[3:nrow(tree$extant),"parent"]+1
  tree$extant$child = tree$extant$child + 1 
  tree$extant$parent[1:2]=1
  tree$extinct$parent=tree$extinct$parent+1
  tree$extinct$child=tree$extinct$child+1
  #NOTE: This next line is a hack
  tree$extant[2,2]<-2
  class(tree)="etree"
  return(tree)
}



GPD<- function(tree,tm){
  # input: an ultramedric tree defined by a data.frame
  # with columns brts, parent, child
  # the first two rows are 
  n<-nrow(tree)
  child.nms<-as.character(tree$child)
  parent.nms<-as.character(tree$parent)
  species.nms<-child.nms
  gpd<-matrix(0,ncol=n,nrow=n)
  dimnames(gpd)<-list(species.nms,species.nms)
  species<-as.list(1:n)
  for (i in seq(n,2)){
    p.set<- species[[which(species.nms==parent.nms[i])]]
    c.set<- species[[which(species.nms==child.nms[i])]]
    gpd[p.set,c.set]<- tm-tree$brts[i]
    species[[which(species.nms==parent.nms[i])]]<-c(p.set,c.set)
  }
  gpd<-gpd+t(gpd)
  return(gpd)
}



# more utilities  (emphasis)

n_from_time <- function(tm,tree){
  # return N at tm.
  if(tm==0){
    N=2
  }else{
    extended_tree = extend_tree(tree)

    n = cumsum(extended_tree$event)+cumsum(extended_tree$event-1)+1
    brts = extended_tree$brts
    if(tm==0) tm = 0.000000000000001
    N = n[max(which(brts < tm))]
  }
  return(N)
} 

n_for_all_bt <- function(tree){
  brts = c(tree$extant$brts,
           tree$extinct$brts,
           tree$extinct$t_ext)
  brts = c(0,sort(brts[brts!=0]))
  n = sapply(brts,n_from_time,tree)
  return(n)
}

extend_tree <- function(tree){
  if(is.null(tree$extinct)){
    tree = list(extant=tree,extinct=data.frame(brts=NULL,t_ext=NULL))
  } 
  extended_tree = data.frame(brts = c(tree$extant$brts,
                                      tree$extinct$brts,
                                      tree$extinct$t_ext),
                             event = c(rep(1,nrow(tree$extant)),
                                       rep(1,nrow(tree$extinct)),
                                       rep(0,nrow(tree$extinct))))
  extended_tree = extended_tree[order(extended_tree$brts),]
  extended_tree = rbind(extended_tree,data.frame(brts=tree$ct,event=2))
  if(extended_tree$brts[2]==0){
    extended_tree = extended_tree[-1,]
  }
  return(extended_tree)
}


foo <- function(phylo, metric = "colless") {
  if(!is.null(phylo$tree)) phylo = phylo$tree
  if (metric == "colless") {
    xx <- apTreeshape:::as.treeshape(phylo)  # convert to apTreeshape format
    apTreeshape:::colless(xx, "yule")  # calculate colless' metric
  } else if (metric == "gamma") {
    ape:::gammaStat(phylo)
  } else stop("metric should be one of colless or gamma")
}

phylodiversity <- function(tree,tm){
  # input: an ultrametric tree
  if(!is.null(tree$extant)){
    tree = get_extant(tree,tm)
  }else{
    tree = tree[tree$brts<tm,]
  }
  dt<-diff(c(tree$brts[-1],tm))
  return(sum(dt*(2:(length(dt)+2-1))))
}



get_current_species <- function(tm,tree){
  species = c(tree$extant$child[tree$extant$brts<tm],
              tree$extinct$child[tree$extinct$t_ext>tm])
  return(species)
}


AIC_llik <- function(LogLik,k){
  aic <- (2*k)-(2*LogLik)
  return(aic)
}

AICw <- function(l1,l2,k1,k2){
  IC <- AIC_llik(c(l1,l2),c(k1,k2))
  bestmodelIC <- min(IC)
  weights <- exp(-0.5*(IC-bestmodelIC))
  weights <- weights/sum(weights)
  return(weights[1])
}



# time calculation
get.time <- function(time,mode='sec'){
  dif = proc.time()-time
  ti = as.numeric(dif[3])
  if(mode == 'min')  ti = ti/60
  if(mode == 'hou') ti = ti/3600
  return(ti)
}


