.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
suppressMessages(library(emphasis))
aug <- emphasis:::.augment_tree_bdi
b0 <- c(5,4.745201,4.530461,4.067871,3.622029,2.929002,1.468698,1.386875,1.302139,0.044729)
for (s in c(1, 20, 200, 1e4)) {
  brts <- b0*s
  for (p in list(c(0.5,0.3)/s, c(0.3,0.5)/s, c(0.4,0.4)/s)) {
    set.seed(11)
    a <- tryCatch(aug(brts,pars=p,model_bin=c(0L,0L,0L),sample_size=20L,link=0L,rho=1),
                  error=function(e) paste("ERROR:",conditionMessage(e)))
    if (is.character(a)) { cat(sprintf("scale=%-6g lam=%.3g mu=%.3g %s\n", s,p[1],p[2],a)); next }
    bl <- DDD::bd_loglik(pars1=c(p,0,0),pars2=c(0,0,1,0,2),brts=brts,missnumspec=0)
    cat(sprintf("scale=%-6g lam=%-9.3g mu=%-9.3g ntree=%2d rng(lw)=%.2e fhat-bd=%+.2e\n",
        s,p[1],p[2],length(a$trees),diff(range(a$weights)), a$fhat-bl))
  }
}
