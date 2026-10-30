suppressMessages(devtools::load_all(quiet = TRUE))
mb <- .resolve_model("cr"); sg <- 0.5
mko <- function(x) { s <- exp(x[1]); .latent_rate_marks(s, log_alpha = -s^2 / 2, kmax = 15) }
set.seed(31)
repeat { r <- .simulate_marks(c(log(0.35), log(0.07)), "cr", 14, mko(log(sg)), max_N = 20000, max_tries = 1L)
  if (r$status == "done") { nt <- sum(r$Ltable[, 4] == -1); if (nt >= 80 && nt <= 150) break } }
phy <- DDD::L2phylo(r$Ltable, dropextinct = TRUE); brts <- .extract_brts(phy)
pf <- .expand_pars(c(log(0.35), log(0.07)), mb); mk <- mko(log(sg))
o <- .mark_offsets(pf, phy, "cr", mk)
cat(sprintf("tips %d; offsets on %d stretches, range %.2f to %.2f\n", ape::Ntip(phy), length(o$off), min(o$off), max(o$off)))
summ <- function(lw) { lw <- lw[is.finite(lw)]; m <- max(lw); w <- exp(lw - m)
  c(logL = m + log(mean(w)), se = sd(w) / (sqrt(length(w)) * mean(w)), ess = sum(w)^2 / sum(w^2)) }
for (aw in c(FALSE, TRUE)) {
  t0 <- Sys.time()
  a <- if (aw) augment_trees(brts, pf, 2000L, 100000L, 5000L, 1e6, 1L, model = mb, link = 1L, rho = 1,
                             parent_tip_start = .pts(brts), seed = 5L, parent_id = .pid(brts),
                             offset_id = o$id, offset_from = o$from, offset_to = o$to, offset = o$off)
       else augment_trees(brts, pf, 2000L, 100000L, 5000L, 1e6, 1L, model = mb, link = 1L, rho = 1,
                          parent_tip_start = .pts(brts), seed = 5L, parent_id = .pid(brts))
  s <- summ(.eval_logf_marks(pf, a$trees, mb, mk) - a$logg)
  cat(sprintf("%-12s logL %.3f (se %.3f)  ESS %.1f of %d  envelope violations %d  %.0f s\n", if (aw) "mark-aware" else "blind",
      s["logL"], s["se"], s["ess"], length(a$trees), a$envelope_violations, as.numeric(Sys.time() - t0, units = "secs")))
}
