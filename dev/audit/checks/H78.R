## H78 — stale `(use_N, use_P, use_E)` / link `0/1` documentation on the C++ wrappers.
## Self-contained. Run: Rscript dev/audit/checks/H78.R
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
options(width = 120)
root <- "/Users/pancho/Code/emphasis"

cat("== Part A: static evidence — where the stale strings live ==\n")
files <- c(list.files(file.path(root, "R"), full.names = TRUE, pattern = "\\.R$"),
           list.files(file.path(root, "src"), full.names = TRUE, pattern = "\\.cpp$"),
           list.files(file.path(root, "inst/include"), full.names = TRUE, pattern = "\\.hpp$"),
           list.files(file.path(root, "man"), full.names = TRUE, pattern = "\\.Rd$"),
           list.files(file.path(root, "tests/testthat"), full.names = TRUE, pattern = "\\.R$"),
           file.path(root, "README.md"))
pat_stale <- "use_P|use_E|beta_E|beta_P|gamma_E|gamma_P"
pat_link  <- "0 = linear.*1 = exponential|\\\\code\\{0\\} = linear, \\\\code\\{1\\} = exponential"
pat_new   <- "use_M|use_D"
hits <- function(pat) {
  out <- character()
  for (f in files) {
    l <- readLines(f, warn = FALSE)
    i <- grep(pat, l)
    if (length(i)) out <- c(out, sprintf("%s:%d: %s", sub(paste0(root, "/"), "", f), i, trimws(l[i])))
  }
  out
}
cat("\n-- stale covariate names (use_P/use_E/beta_E/...):\n"); writeLines(hits(pat_stale))
cat("\n-- link docs listing only 0/1:\n"); writeLines(hits(pat_link))
cat("\n-- current names (use_M/use_D) for contrast:\n"); writeLines(hits(pat_new))

cat("\n-- installed help (Rd_db of the scratch build): does ?augment_trees / ?eval_logf say use_P?\n")
db <- tools::Rd_db("emphasis", lib.loc = "/Users/pancho/.claude/jobs/867af780/tmp/rlib")
for (nm in c("augment_trees.Rd", "eval_logf.Rd", "em_cpp.Rd", "m_cpp.Rd", "simulate_div_tree_cpp.Rd")) {
  txt <- paste(capture.output(tools::Rd2txt(db[[nm]], options = list(underline_titles = FALSE))), collapse = "\n")
  cat(sprintf("  %-26s use_P: %-5s  'use_D': %-5s  gaussian mentioned: %s\n", nm,
              grepl("use_P", txt), grepl("use_D", txt), grepl("gaussian", txt)))
}

cat("\n== Part B: behavioural evidence that the wrapper docs describe the wrong contract ==\n")
# Hand-built tree in the C++ convention (see H11.R): forward times, crown at 0,
# observed speciation nodes at brts 1..k with n = 2..k+1, closing sentinel at tp.
# `pd` is the pendant PD P at the node; C++ turns it into M = P/n.
mk_tree <- function(k, tp = k + 1, pd = 0) {
  data.frame(brts = c(seq_len(k), tp), n = c(2:(k + 1), k + 2), t_ext = 1e11,
             pd = pd, tip_start = 0, id = c(seq_len(k) - 1L, -1L), parent_id = -1L)
}
ev <- function(p8, tr, model, link) emphasis:::eval_logf(p8, list(tr), model = as.integer(model), link = as.integer(link), rho = 1)$logf
hand_loglik <- function(tr, lam_fun, mu_fun) {
  dt <- diff(c(0, tr$brts)); N <- tr$n; P <- tr$pd
  sum(log(lam_fun(N[-nrow(tr)], P[-nrow(tr)]))) - sum(dt * N * (lam_fun(N, P) + mu_fun(N, P)))
}
tr0 <- mk_tree(5)
tr1 <- mk_tree(5, pd = c(0.4, 1.1, 1.9, 2.5, 3.3, 4.0))

cat("\nB1. link = 2 (gaussian) is accepted by every wrapper whose roxygen says '0 = linear, 1 = exponential':\n")
p_g <- c(0.8, 0, 0, 0, 0.2, 0, 0, 0)            # CR under gaussian: lambda = b0 e^{-1/2}
p_l <- c(0.8 * exp(-0.5), 0, 0, 0, 0.2 * exp(-0.5), 0, 0, 0)  # same rates under linear
lf_g <- ev(p_g, tr0, c(0, 0, 0), 2); lf_l <- ev(p_l, tr0, c(0, 0, 0), 0); lf_e <- ev(p_g, tr0, c(0, 0, 0), 1)
cat(sprintf("   eval_logf link=2: %.6f | link=0 with b0*e^-1/2: %.6f | link=1 same pars: %.6f  -> link=2 is a third, distinct link (diff to linear-equiv = %.1e)\n",
            lf_g, lf_l, lf_e, abs(lf_g - lf_l)))
aug <- emphasis:::augment_trees(brts = c(4, 3, 2, 1), pars = p_g, sample_size = 20, maxN = 2000,
                                max_missing = 50, max_lambda = 100, num_threads = 1,
                                model = c(0L, 0L, 0L), link = 2L, rho = 1)
cat(sprintf("   augment_trees link=2: %d trees returned, all logf finite: %s\n", length(aug$trees), all(is.finite(aug$logf))))
sim <- emphasis:::simulate_div_tree_cpp(p_g, c(0L, 0L, 0L), max_t = 3, max_N = 500L, max_tries = 20L, link = 2L)
cat(sprintf("   simulate_div_tree_cpp link=2: status = %s, Ltable rows = %d\n", sim$status, nrow(sim$Ltable)))

cat("\nB2. model slot 2 is documented as use_P; C++ never reads it, and pars[3] multiplies M = P/N, not P:\n")
pM <- c(0.8, 0, 0.3, 0, 0.2, 0, 0, 0)            # beta_M = 0.3
v000 <- ev(pM, tr1, c(0, 0, 0), 0); v010 <- ev(pM, tr1, c(0, 1, 0), 0); v110 <- ev(pM, tr1, c(1, 1, 0), 0)
cat(sprintf("   eval_logf(beta_M=0.3): model=c(0,0,0) %.6f | c(0,1,0) %.6f | c(1,1,0) %.6f  -> slot 2 gates nothing (max diff %.1e)\n",
            v000, v010, v110, max(abs(c(v000 - v010, v000 - v110)))))
h_M <- hand_loglik(tr1, function(N, P) pmax(0, 0.8 + 0.3 * P / N), function(N, P) rep(0.2, length(N)))
h_P <- hand_loglik(tr1, function(N, P) pmax(0, 0.8 + 0.3 * P),     function(N, P) rep(0.2, length(N)))
cat(sprintf("   hand formula with M=P/N: %.6f (|diff| = %.1e) ; with raw P: %.6f (|diff| = %.1e)\n",
            h_M, abs(h_M - v000), h_P, abs(h_P - v000)))

cat("\nB3. model slot 3 is documented as use_E; C++ uses D = E - M (for parent_id = -1 nodes E := M, so D = 0):\n")
pD <- c(0.8, 0, 0, 0.5, 0.2, 0, 0, 0)            # beta_D = 0.5, pd != 0 so raw E = M != 0
v_D  <- ev(pD, tr1, c(0, 0, 1), 0); v_0 <- ev(c(0.8, 0, 0, 0, 0.2, 0, 0, 0), tr1, c(0, 0, 1), 0)
h_E  <- hand_loglik(tr1, function(N, P) pmax(0, 0.8 + 0.5 * P / N), function(N, P) rep(0.2, length(N)))
cat(sprintf("   eval_logf(beta_D=0.5, model=c(0,0,1)) = %.6f ; beta_D=0 = %.6f (|diff| = %.1e -> D = 0 on these nodes)\n",
            v_D, v_0, abs(v_D - v_0)))
cat(sprintf("   if slot 3 meant raw E (= pendant age = M here) the hand value would be %.6f (|diff| = %.1e)\n", h_E, abs(h_E - v_D)))

cat("\nB4. The public layer already uses the current names — the stale text is confined to the wrappers:\n")
cat(sprintf("   .resolve_model('dd') = %s ; .resolve_model('nd') = %s ; .resolve_link('gaussian') = %d\n",
            paste(emphasis:::.resolve_model("dd"), collapse = ","), paste(emphasis:::.resolve_model("nd"), collapse = ","),
            emphasis:::.resolve_link("gaussian")))

cat("\nB5. Effect on the validation-study configurations (cr / dd, link 0, rho = 1): none — the wrappers receive\n")
cat("    model = c(0,0,0) / c(1,0,0) and link = 0L regardless of what their roxygen says; the text is comments only.\n")
cat(sprintf("   eval_logf cr: %.6f (link 0 doc'd and used) ; dd c(1,0,0): %.6f\n",
            ev(c(0.8, 0, 0, 0, 0.2, 0, 0, 0), tr0, c(0, 0, 0), 0), ev(c(0.8, -0.02, 0, 0, 0.2, 0, 0, 0), tr0, c(1, 0, 0), 0)))
