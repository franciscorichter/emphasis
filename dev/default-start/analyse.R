suppressMessages(devtools::load_all("~/pancho/emphasis-conv", quiet = TRUE)); library(ape)
o <- "~/pancho/conv-out"; tr <- readRDS(file.path(o, "trees.rds"))
rd <- function(set, i, v) readRDS(file.path(o, sprintf("%s_%02d_%s.rds", set, i, v)))
ex <- t(sapply(tr$cr, function(t) .cr_mle(branching.times(t), cond = TRUE)))
tab <- do.call(rbind, lapply(1:40, function(i) { a <- rd("cr", i, "old"); b <- rd("cr", i, "new")
  data.frame(i = i, n = Ntip(tr$cr[[i]]), lam_ex = ex[i,1], mu_ex = ex[i,2],
             lam_old = exp(a$pars[1]), mu_old = exp(a$pars[2]), it_old = a$iter, s_old = a$secs,
             lam_new = exp(b$pars[1]), mu_new = exp(b$pars[2]), it_new = b$iter, s_new = b$secs, stop_new = b$stop) }))
cat("== cr, 40 trees: |error| vs exact conditioned MLE ==\n")
f <- function(x) sprintf("median %.3f  max %.3f", median(x), max(x))
cat("old  lam", f(abs(tab$lam_old - tab$lam_ex)), " | mu", f(abs(tab$mu_old - tab$mu_ex)), " | iters median", median(tab$it_old), " secs median", median(tab$s_old), "\n")
cat("new  lam", f(abs(tab$lam_new - tab$lam_ex)), " | mu", f(abs(tab$mu_new - tab$mu_ex)), " | iters median", median(tab$it_new), " secs median", median(tab$s_new), "\n")
cat("medians: exact mu", round(median(tab$mu_ex),3), " old", round(median(tab$mu_old),3), " new", round(median(tab$mu_new),3), "\n")
print(table(tab$stop_new)); print(summary(tab$it_new))
cat("tree 4:\n"); print(round(tab[4, c("lam_ex","mu_ex","lam_old","mu_old","lam_new","mu_new")], 3))
cat("\n== cr, rule only (midpoint start, projected rule), 10 trees ==\n")
for (i in 1:10) { r <- rd("cr", i, "rule_only"); a <- rd("cr", i, "old")
  cat(sprintf("%2d exact %.3f %.3f | old %.3f %.3f (%3d it) | rule_only %.3f %.3f (%3d it, %s, %4.0fs)\n", i, ex[i,1], ex[i,2],
      exp(a$pars[1]), exp(a$pars[2]), a$iter, exp(r$pars[1]), exp(r$pars[2]), r$iter, r$stop, r$secs)) }
cat("\n== dd, 10 trees (pars: beta0 betaN gamma0 gammaN; loglik) ==\n")
for (i in 1:10) { s <- lapply(c("old","new","ref"), function(v) rd("dd", i, v))
  cat(sprintf("%2d n=%3d", i, Ntip(tr$dd[[i]])))
  for (k in 1:3) cat(sprintf(" | %s %s ll %.2f %3d it %4.0fs %s", c("old","new","ref")[k],
      paste(sprintf("%.3f", s[[k]]$pars), collapse=" "), s[[k]]$loglik, s[[k]]$iter, s[[k]]$secs, substr(s[[k]]$stop,1,4)))
  cat("\n") }
saveRDS(tab, file.path(o, "cr_table.rds"))
