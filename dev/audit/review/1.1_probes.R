lib <- commandArgs(TRUE)[1]; .libPaths(c(lib, .libPaths())); library(emphasis)
mk <- function(k, tp = k+1) data.frame(brts=c(seq_len(k),tp), n=c(2:(k+1),k+2), t_ext=1e11,
  pd=0, tip_start=0, id=c(seq_len(k)-1L,-1L), parent_id=-1L)
ev <- function(p,tr,m,l) emphasis:::eval_logf(p, list(tr), model=as.integer(m), link=as.integer(l), rho=1)$logf
p8 <- function(...) {v <- c(...); c(v, rep(0, 8-length(v)))}
cat(sprintf("%-46s %s\n", "probe", "logf"))
f <- function(lab, x) cat(sprintf("%-46s %s\n", lab, format(x)))
# exp link, eta underflows to exactly 0 (true logf finite ~ -800*k)
f("exp link eta=-800 (exp underflows to 0), k=3", ev(p8(-800,0,0,0,-800), mk(3), c(0,0,0), 1))
f("exp link eta=-400 (no underflow), k=3",        ev(p8(-400,0,0,0,-400), mk(3), c(0,0,0), 1))
f("exp link eta=-745.2 boundary, k=3",            ev(p8(-745.2,0,0,0,-745.2), mk(3), c(0,0,0), 1))
# gaussian negative intercept, EVEN number of scored speciation nodes
f("gauss beta0=-0.4, k=4 (even # negative rates)", ev(p8(-0.4,0,0,0,0.2), mk(4), c(0,0,0), 2))
f("gauss beta0=-0.4, k=3 (odd)",                   ev(p8(-0.4,0,0,0,0.2), mk(3), c(0,0,0), 2))
f("gauss beta0=0 (rate exactly 0), k=3",           ev(p8(0,0,0,0,0.2), mk(3), c(0,0,0), 2))
# mu floor asymmetry: extinction node with mu = 0 (linear, gamma0 = 0)
tr <- data.frame(brts=c(1,2,3), n=c(2,3,2), t_ext=c(1e11,2,1e11), pd=0, tip_start=0,
                 id=c(0L,1L,-1L), parent_id=-1L)
f("linear, mu=0 at an extinction node",            ev(p8(0.5,0,0,0,0), tr, c(0,0,0), 0))
f("linear, mu=1e-300 at an extinction node",       ev(p8(0.5,0,0,0,1e-300), tr, c(0,0,0), 0))
# lambda exactly 0 at the SAME tree (speciation) for contrast
f("linear, lambda=0, mu=0.1",                      ev(p8(0,0,0,0,0.1), tr, c(0,0,0), 0))
# huge rates: prod_ overflow with negative sum_
f("linear lambda=1e300 k=3 (prod overflow)",       ev(p8(1e300,0,0,0,1), mk(3), c(0,0,0), 0))
f("NaN parameter",                                 ev(p8(NaN,0,0,0,0.2), mk(3), c(0,0,0), 0))
