suppressWarnings(suppressMessages(devtools::load_all(".", quiet = TRUE)))
ns <- asNamespace("emphasis")

set.seed(11)
phy  <- ape::rphylo(14, 0.7, 0.15)
phy  <- ns$prune_to_extant(phy)
brts <- ns$.extract_brts(phy)
pts  <- ns$.pts(brts)
cat("tips", ape::Ntip(phy), " brts", length(brts), " pts", length(pts), "\n")

# --- augmentation with real extinct lineages under the D model --------------
a <- augment_trees(as.numeric(brts), c(0.9, 0, 0, 0.10, 0.35, 0, 0, 0),
                   40L, 20000L, 300L, 1e6, 1L,
                   model = c(0L, 0L, 1L), link = 0L, rho = 1,
                   parent_tip_start = pts)
cat("trees:", length(a$trees), " with augmented lineages:",
    sum(vapply(a$trees, function(d) nrow(d) > length(brts), TRUE)), "\n")
allp <- unlist(lapply(a$trees, `[[`, "pd"))
alln <- unlist(lapply(a$trees, `[[`, "n"))
allt <- unlist(lapply(a$trees, `[[`, "brts"))
cat("pd  : min", min(allp), " max", max(allp), " any NaN", anyNA(allp), "\n")
cat("M   : range", paste(round(range(allp / alln), 3), collapse = " .. "), "\n")
cat("M/t : max  ", max(allp / alln / allt), " (must be <= 1)\n")
cat("P<0 : ", sum(allp < 0), "\n")
cat("logf finite:", sum(is.finite(a$logf)), "/", length(a$logf), "\n")

# D at the events of one augmented tree must still centre on the alive set.
d <- a$trees[[which.max(vapply(a$trees, nrow, 0L))]]
E <- ifelse(d$focal_tip_start >= 0, d$brts - d$focal_tip_start, d$pd / d$n)
cat("D range:", paste(round(range(E - d$pd / d$n), 3), collapse = " .. "), "\n")

# --- end-to-end: the thinning driver reaches C++ with the topology ----------
lf <- function(bd) {
  set.seed(3)
  eval_logf(c(0.9, 0, 0, bd, 0.35, 0, 0, 0), a$trees,
            model = c(0L, 0L, 1L), link = 0L, rho = 1)$logf
}
cat("logf sensitive to beta_D:", !isTRUE(all.equal(lf(0), lf(0.3))), "\n")

set.seed(5)
fit_d <- estimate_rates(phy, method = "mcem", model = "d",
                        init_pars = c(0.9, 0.0, 0.3, 0.0),
                        link = "exponential",
                        control = list(lower_bound = c(-2, -0.5, -4, -0.5),
                                       upper_bound = c(2, 0.5, 1, 0.5),
                                       num_trees = 40L, max_iter = 3L,
                                       sampling = "dynamic_fresh"))
print(fit_d)

set.seed(6)
fit_cem <- estimate_rates(phy, method = "cem", model = "d",
                          link = "exponential",
                          control = list(lower_bound = c(-2, -0.5, -4, -0.5),
                                         upper_bound = c(2, 0.5, 1, 0.5),
                                         max_iter = 3L, num_particles = 12L,
                                         num_trees = 5L))
cat("cem loglik:", fit_cem$loglik, "\n")
