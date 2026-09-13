.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape); library(DDD) })
pars <- c(0.6, 0.25); max_t <- 6
sim <- NULL
for (i in 1:300) { s <- simulate_tree(pars = pars, max_t = max_t, model = "cr", max_tries = 50)
  if (s$status == "done" && !is.null(s$tes) && Ntip(s$tes) >= 12 && Ntip(s$tes) <= 30) { sim <- s; break } }
L <- sim$L; brts <- emphasis:::.extract_brts(sim)
colnames(L) <- c("birth_age","parent","label","ext_age")
cat("sim$L (raw simulator L-table):\n"); print(round(L, 5))
cat("\nobserved brts:\n"); print(round(brts, 5))
ext <- which(L[,4] != -1)
for (r in ext) {
  m <- which(abs(brts - L[r,1]) < 1e-6)
  if (length(m)) {
    same <- which(abs(L[,1] - L[r,1]) < 1e-6 & seq_len(nrow(L)) != r)
    cat(sprintf("extinct row %d (label %d, parent %d, born %.5f, died %.5f) coincides with brts[%d]; other rows born at the same age: %s\n",
        r, L[r,3], L[r,2], L[r,1], L[r,4], m, paste(sprintf("row %d (label %d, parent %d, ext %.3f)", same, L[same,3], L[same,2], L[same,4]), collapse="; ")))
  }
}
# DDD reference: does phylo2L(tas) show the same coincidences?
L2 <- DDD::phylo2L(sim$tas)
cat(sprintf("\nphylo2L(tas): %d rows, %d extinct; duplicated birth ages among all rows: %d; in sim$L: %d\n",
    nrow(L2), sum(L2[,4] != -1), sum(duplicated(round(L2[,1],6))), sum(duplicated(round(L[,1],6)))))
# what does .aug_to_Ltable's lookup return for each observed branching, on each base?
Le <- DDD::phylo2L(sim$tes)
for (j in 2:length(brts)) {
  rf <- which.min(abs(L[,1] - brts[j])); re <- which.min(abs(Le[,1] - brts[j]))
  if (L[rf,4] != -1) cat(sprintf("brts[%d]=%.5f -> sim$L row %d is EXTINCT (died %.4f); phylo2L(tes) row %d ext %.1f\n", j, brts[j], rf, L[rf,4], re, Le[re,4]))
}
