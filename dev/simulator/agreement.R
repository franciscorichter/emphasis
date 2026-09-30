# Does the browser simulator (assets/simulator/emphasis-sim.js, the engine of
# simulator.html) draw from the same distribution as simulate_tree()?
#
#   Rscript dev/simulator/agreement.R [K] [scenarios]   # from the package root; K = 2000
#
# scenarios is an optional comma-separated subset of the names below; with a
# subset the check still prints its table but does not rewrite agreement.csv.
#
# Two checks.
#   ed    The JS fair-proportion routine against emphasis:::ed_fair_proportion()
#         on the same trees: equal to rounding, tip by tip.
#   draws K surviving trees per scenario from each simulator, seven summaries
#         per tree.  For every summary the difference of the two means in units
#         of its standard error (z), and a two-sample KS test on the extant-tip
#         count.  Agreement: every |z| < 4 and every KS p > 0.001.
#
# The two simulators do not share a random stream, so agreement is distributional.
# Writes dev/simulator/agreement.csv.

suppressMessages(devtools::load_all(".", quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
K    <- if (length(args)) as.integer(args[1]) else 2000L
node <- Sys.which("node")
if (!nzchar(node)) node <- path.expand("~/.local/bin/node")
tmp  <- tempfile("simagree"); dir.create(tmp)

# --------------------------------------------------------------------------- #
#  ed: the JS routine against the package's
# --------------------------------------------------------------------------- #
ed_csv <- file.path(tmp, "ed.csv")
system2(node, c("dev/simulator/js_ed.js", "40", "11", ed_csv))
ed <- read.csv(ed_csv)
ed_err <- do.call(rbind, lapply(split(ed, ed$tree), function(d) {
  r <- emphasis:::ed_fair_proportion(as.integer(d$parent), d$birth, d$alive == 1, d$t[1])
  keep <- d$alive == 1
  data.frame(n_tips = sum(keep), max_abs_err = max(abs(r$ed[keep] - d$ed_js[keep])))
}))
cat(sprintf("ed: %d trees, %d tips, largest |ED_js - ED_cpp| = %.2e\n",
            nrow(ed_err), sum(ed_err$n_tips), max(ed_err$max_abs_err)))
ed_ok <- max(ed_err$max_abs_err) < 1e-9

# --------------------------------------------------------------------------- #
#  draws
# --------------------------------------------------------------------------- #
scen <- list(
  list(name = "cr_linear",    model = "cr",  link = "linear",      pars = c(0.8, 0.3),                      maxT = 7,  rho = 1),
  list(name = "cr_linear_rho",model = "cr",  link = "linear",      pars = c(0.8, 0.3),                      maxT = 7,  rho = 0.5),
  list(name = "dd_linear",    model = "dd",  link = "linear",      pars = c(1.0, -0.015, 0.1, 0),           maxT = 10, rho = 1),
  list(name = "dd_gaussian",  model = "dd",  link = "gaussian",    pars = c(1.2, 0.02, 0.2, 0),             maxT = 8,  rho = 1),
  list(name = "d_linear",     model = "d",   link = "linear",      pars = c(0.6, -0.3, 0.15, 0),            maxT = 7,  rho = 1),
  list(name = "nd_exp",       model = "nd",  link = "exponential", pars = c(0, -0.02, -0.3, log(0.1), 0, 0.2), maxT = 8, rho = 1),
  list(name = "ed_linear",    model = "ed",  link = "linear",      pars = c(0.9, -0.1, 0.2, 0.02),          maxT = 6,  rho = 1),
  list(name = "ned_exp",      model = "ned", link = "exponential", pars = c(0, -0.015, -0.15, log(0.15), 0, 0.05), maxT = 7, rho = 1),
  list(name = "edc_linear",   model = "edc", link = "linear",      pars = c(0.9, -0.2, 0.2, 0.02),          maxT = 6,  rho = 1),
  list(name = "nedc_exp",     model = "nedc",link = "exponential", pars = c(0, -0.015, -0.2, log(0.15), 0, 0.05), maxT = 7, rho = 1)
)
only <- if (length(args) > 1) strsplit(args[2], ",")[[1]] else NULL
if (!is.null(only)) scen <- Filter(function(s) s$name %in% only, scen)
scen_json <- file.path(tmp, "scen.json")
writeLines(jsonlite::toJSON(scen, auto_unbox = TRUE, digits = NA), scen_json)

t0 <- Sys.time()
js_csv <- file.path(tmp, "js.csv")
system2(node, c("dev/simulator/js_draws.js", scen_json, K, "7", js_csv))
js <- read.csv(js_csv)
cat(sprintf("js draws: %.0f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# The same seven summaries from an R L-table.
summarise_L <- function(L, T) {
  dropped <- L[, 4] == T            # simulate_tree() marks a dropped tip with death age max_t
  extant  <- L[, 4] == -1 | dropped
  birth   <- T - L[, 1]
  death   <- ifelse(extant, Inf, T - L[, 4])
  lab     <- L[, 3]
  row_of  <- function(l) match(l, lab)
  par_row <- row_of(L[, 2]); par_row[1:2] <- NA     # crown rows have no parent row
  kids    <- split(seq_len(nrow(L))[-(1:2)], factor(par_row[-(1:2)], levels = seq_len(nrow(L))))
  tip_start <- vapply(seq_len(nrow(L)), function(i) max(c(birth[i], birth[kids[[i]]])), 0)
  surv <- extant
  for (i in rev(seq_len(nrow(L)))) if (!surv[i] && length(kids[[i]])) surv[i] <- any(surv[kids[[i]]])
  end_recon <- vapply(seq_len(nrow(L)), function(i) {
    if (!surv[i]) return(-Inf)
    if (extant[i]) return(Inf)
    k <- kids[[i]][surv[kids[[i]]]]
    max(birth[k])
  }, 0)
  recon <- function(tt) sum(surv & birth <= tt & end_recon > tt)
  c(n_extant = sum(extant), n_sampled = sum(L[, 4] == -1), n_rows = nrow(L),
    n_extinct = sum(!extant), pendant_mean = mean(T - tip_start[extant]),
    n_complete_half = sum(birth <= T / 2 & death > T / 2),
    n_recon_half = recon(T / 2), n_recon_80 = recon(0.8 * T))
}

set.seed(7)
t0 <- Sys.time()
rr <- do.call(rbind, lapply(scen, function(s) {
  out <- t(vapply(seq_len(K), function(k) {
    x <- simulate_tree(pars = s$pars, max_t = s$maxT, model = s$model, link = s$link,
                       rho = s$rho, max_tries = 10000, useDDD = FALSE)
    stopifnot(x$status == "done")
    c(attempts = 1 / x$survival_prob, summarise_L(x$L, s$maxT))
  }, numeric(9)))
  data.frame(scenario = s$name, out)
}))
cat(sprintf("R draws: %.0f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

stats <- c("attempts", "n_extant", "n_sampled", "n_rows", "n_extinct", "pendant_mean",
           "n_complete_half", "n_recon_half", "n_recon_80")
res <- do.call(rbind, lapply(scen, function(s) {
  a <- rr[rr$scenario == s$name, ]; b <- js[js$scenario == s$name, ]
  do.call(rbind, lapply(stats, function(v) {
    se <- sqrt(var(a[[v]]) / nrow(a) + var(b[[v]]) / nrow(b))
    z  <- if (se > 0) (mean(b[[v]]) - mean(a[[v]])) / se else 0
    data.frame(scenario = s$name, stat = v, mean_R = mean(a[[v]]), mean_js = mean(b[[v]]), z = z)
  }))
}))
ks <- vapply(scen, function(s) suppressWarnings(
  ks.test(rr$n_extant[rr$scenario == s$name], js$n_extant[js$scenario == s$name])$p.value), 0)
names(ks) <- vapply(scen, `[[`, "", "name")

res$ks_p_n_extant <- ifelse(res$stat == "n_extant", ks[res$scenario], NA)
if (is.null(only)) write.csv(transform(res, mean_R = signif(mean_R, 5), mean_js = signif(mean_js, 5), z = round(z, 2),
                    ks_p_n_extant = signif(ks_p_n_extant, 3)),
          "dev/simulator/agreement.csv", row.names = FALSE)

wide <- reshape(res[, c("scenario", "stat", "z")], idvar = "scenario", timevar = "stat",
                direction = "wide")
names(wide) <- sub("^z\\.", "", names(wide))
wide$ks_p <- signif(ks[wide$scenario], 3)
cat(sprintf("\nz = (mean_js - mean_R) / se, K = %d surviving trees per simulator and scenario\n", K))
print(format(wide, digits = 2), row.names = FALSE)

ok <- ed_ok && all(abs(res$z) < 4) && all(ks > 0.001)
cat(sprintf("\n%s: max |z| = %.2f, min KS p = %.3g, ED %s\n",
            if (ok) "AGREE" else "DISAGREE", max(abs(res$z)), min(ks),
            if (ed_ok) "exact" else "MISMATCH"))
quit(status = if (ok) 0 else 1)
