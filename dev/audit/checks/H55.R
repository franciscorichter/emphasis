.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(5)
brts <- sort(c(5, runif(14, 0, 5)), decreasing = TRUE)
# (1) 8-element pars with non-zero slot 3/4 accepted by .augment_tree_bdi?
pars8 <- c(0.5, 0, 0.1, 0.2, 0.2, 0, 0, 0)
r <- tryCatch(emphasis:::.augment_tree_bdi(brts, pars = pars8, model_bin = c(0L,0L,0L), sample_size = 5L, link = 0L, rho = 1), error = function(e) e)
cat("8-elem pars, CR: ", if (inherits(r,"error")) conditionMessage(r) else sprintf("accepted, %d trees", length(r$trees)), "\n")
r2 <- tryCatch(emphasis:::.augment_tree_bdi(brts, pars = pars8, model_bin = c(1L,0L,1L), sample_size = 5L, link = 0L, rho = 1), error = function(e) e)
cat("8-elem pars, model (1,0,1): ", if (inherits(r2,"error")) conditionMessage(r2) else sprintf("accepted, %d trees", length(r2$trees)), "\n")
# (2) rate convention: R sampler vs C++ (M = P/N, D = E - M)
N <- 4; P <- 6; E <- 2.5
cat("R .bdi_lam(slot3=0.1 on P, slot4=0.2 on E):", emphasis:::.bdi_lam(pars8, N, P, E, c(0,0,0), 0L),
    "  C++ convention beta0 + 0.1*M + 0.2*D:", 0.5 + 0.1*(P/N) + 0.2*(E - P/N), "\n")
# (3) inertness of the Gaussian-closure branch under reachable configs (compact pars, DD, slots 3/4 = 0)
pars8_dd <- emphasis:::.expand_pars(c(0.6, -0.01, 0.2, 0.0), c(1L,0L,0L))
cat("expanded DD pars8:", pars8_dd, "\n")
bt <- sort(brts[1] - brts[-1]); tp <- brts[1]
a <- emphasis:::.bdi_iterate(pars8_dd, c(1L,0L,0L), 0L, bt, tp, use_gaussian_closure = TRUE)
b <- emphasis:::.bdi_iterate(pars8_dd, c(1L,0L,0L), 0L, bt, tp, use_gaussian_closure = FALSE)
g <- a$t_grid
cat("max|Nhat diff| =", max(abs(a$Nhat_fun(g) - b$Nhat_fun(g))), " max|p diff| =", max(abs(a$p_fun(g) - b$p_fun(g))),
    " max|Ehat diff| =", max(abs(a$Ehat_fun(g) - b$Ehat_fun(g))), " (Ehat differs but has slope 0)\n")
cat("vN nonzero?", any(a$vN_vals != 0), " cNP nonzero?", any(a$cNP_vals != 0), "\n")
# public gate
cat(".bdi_supported: ", emphasis:::.bdi_supported(c(1,0,0),0), emphasis:::.bdi_supported(c(0,1,0),0), emphasis:::.bdi_supported(c(1,0,1),1), emphasis:::.bdi_supported(c(1,0,0),2), "\n")
