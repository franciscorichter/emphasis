# Old against new forward simulator on raw-ED generators: tip count, gamma,
# Colless, median branching time, survival; 400 surviving clades each.
#   Rscript dev/simcheck_ed.R <lib> <tag> <out.csv>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
suppressMessages({ library(emphasis); library(ape) })
cat("emphasis from", find.package("emphasis"), "\n")
colless <- function(phy) {
  n <- Ntip(phy); if (n < 3) return(NA_real_)
  desc <- function(node) { if (node <= n) return(1L); kids <- phy$edge[phy$edge[, 1] == node, 2]; sum(vapply(kids, desc, 1L)) }
  s <- 0; for (nd in (n + 1):(n + phy$Nnode)) { kids <- phy$edge[phy$edge[, 1] == nd, 2]; s <- s + abs(desc(kids[1]) - desc(kids[2])) }
  s / ((n - 1) * (n - 2) / 2)
}
gammastat <- function(phy) tryCatch(gammaStat(phy), error = function(e) NA_real_)
cases <- list(
  ned_exp_neg = list(pars = c(log(0.5), -0.0088, -0.0899, log(0.05), 0, 0), model = "ned", link = "exponential"),
  ned_exp_pos = list(pars = c(log(0.5), -0.0200, +0.0500, log(0.05), 0, 0), model = "ned", link = "exponential"),
  ed_exp_neg  = list(pars = c(log(0.5), -0.3380, log(0.05), 0), model = "ed", link = "exponential"),
  ned_lin_neg = list(pars = c(0.5, -0.0020, -0.0300, 0.05, 0, 0), model = "ned", link = "linear"))
out <- NULL
for (nm in names(cases)) {
  cs <- cases[[nm]]; got <- 0; att <- 0; t0 <- proc.time()[["elapsed"]]
  while (got < 400 && att < 4000) {
    att <- att + 1; set.seed(att)
    s <- simulate_tree(pars = cs$pars, max_t = 15, model = cs$model, link = cs$link, rho = 1, max_lin = 5000L, max_tries = 0, num_threads = 1L)
    if (!identical(s$status, "done") || is.null(s$tes)) next
    got <- got + 1; phy <- s$tes
    out <- rbind(out, data.frame(case = nm, tag = args[2], tips = Ntip(phy), tas = Ntip(s$tas), gamma = gammastat(phy),
                                 colless = colless(phy), brt = median(branching.times(phy)) / max(branching.times(phy)),
                                 survival = NA_real_))
  }
  cat(sprintf("%-12s %d of %d attempts survived, %.1f s\n", nm, got, att, proc.time()[["elapsed"]] - t0))
  out$survival[out$case == nm] <- got / att
}
write.csv(out, args[3], row.names = FALSE)
