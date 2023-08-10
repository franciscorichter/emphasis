## Visualizing results 

minpars = as.data.frame(est_de$min_pars)
names(minpars) = c("par1","par2","par3","par4")
meanpars = as.data.frame(est_de$mean_pars)
names(meanpars) = c("par1","par2","par3","par4")

learning.curves(minpars,meanpars)

viz_de_vectors(est_de$minloglik,est_de$meanloglik)
