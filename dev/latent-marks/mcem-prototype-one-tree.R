suppressMessages(devtools::load_all(quiet = TRUE))
mb <- .resolve_model("cr"); sg <- 0.5
mko <- function(x) { s <- exp(x[1]); .latent_rate_marks(s, log_alpha = -s^2 / 2, kmax = 15) }
set.seed(31)
repeat { r <- .simulate_marks(c(log(0.35), log(0.07)), "cr", 14, mko(log(sg)), max_N = 20000, max_tries = 1L)
  if (r$status == "done") { nt <- sum(r$Ltable[, 4] == -1); if (nt >= 80 && nt <= 150) break } }
phy <- DDD::L2phylo(r$Ltable, dropextinct = TRUE); cat("tips", ape::Ntip(phy), " lineages ever", nrow(r$Ltable), "\n")
brts <- .extract_brts(phy); th <- c(log(0.35), log(0.07))
a <- augment_trees(brts, .expand_pars(th, mb), 500L, 25000L, 5000L, 1e6, 1L, model = mb, link = 1L, rho = 1,
                   parent_tip_start = .pts(brts), seed = 3L, parent_id = .pid(brts))
ess <- function(lw) { lw <- lw[is.finite(lw)]; w <- exp(lw - max(lw)); sum(w)^2 / sum(w^2) }
lw0 <- a$logf - a$logg; lw1 <- .eval_logf_marks(.expand_pars(th, mb), a$trees, mb, mko(log(sg))) - a$logg
cat(sprintf("ESS of 500 at the truth: no marks %.1f, latent rate sigma 0.5 %.1f\n", ess(lw0), ess(lw1)))
t0 <- Sys.time()
f <- .mcem_marks(phy, "cr", c(log(0.3), log(0.1)), log(0.3), mko, n_iter = 15L, N = 200L, verbose = TRUE)
cat(sprintf("estimate: lambda0 %.3f mu0 %.3f sigma %.3f  (%.0f s)\n", exp(f$estimate[1]), exp(f$estimate[2]), exp(f$estimate[3]),
    as.numeric(Sys.time() - t0, units = "secs")))
