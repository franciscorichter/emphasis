.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
f <- function(p8, lab) {
  r <- replicate(200, { x <- emphasis:::simulate_div_tree_cpp(p8, c(0L,0L,0L), 5, 1000L, 0L, 2L); c(x$status, nrow(x$Ltable)) })
  cat(lab, ":\n"); print(table(status = r[1,], rows = r[2,]))
}
f(c(-0.5,0,0,0,0.1,0,0,0), "lam=-0.5*e^-.5, mu=0.1*e^-.5 (total<0)")
f(c(-0.5,0,0,0,1.0,0,0,0), "lam=-0.5*e^-.5, mu=1.0*e^-.5 (total>0, p=-1)")
f(c(0.5,0,0,0,-0.1,0,0,0), "lam>0, mu<0 (p>1)")
