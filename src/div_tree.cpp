#include "general_tree.hpp"

#include <Rcpp.h>
#include <cstdint>

// src/rcpp_mce.cpp — the `seed` argument into the integer the samplers run on.
namespace emphasis { uint64_t resolve_seed(int seed); }

//' Simulate a single phylogenetic tree under the general diversification model
//'
//' @param pars Numeric vector of length 8:
//'   \code{c(beta_0, beta_N, beta_M, beta_D, gamma_0, gamma_N, gamma_M, gamma_D)}.
//'   Inactive parameters (those not selected by \code{model}) should be set to 0.
//' @param model Integer vector of length 3: \code{c(use_N, use_P, use_E)} (each 0 or 1).
//' @param max_t Crown age (forward simulation end time).
//' @param max_N Maximum number of lineages before simulation is declared too large.
//' @param max_tries Maximum retries after extinction; an overflow is not retried.
//' @param link Link function: 0 = linear (max(0,...)), 1 = exponential.
//' @param seed Positive integer seeding the simulator. \code{0} (the default)
//'   draws one from R's generator, so \code{set.seed()} reaches the simulator
//'   whether or not a seed is passed. Note that \code{max_tries} retries
//'   continue the same stream: a retry is a fresh tree, not the same one again.
//' @param marks_S,log_clam,log_cmu,marks_C,marks_P,marks_pi Latent marks
//'   (dev/latent-marks/design.md; the arguments of \code{eval_logf_marks}).
//'   The default, one mark, simulates without marks.
//' @return A named list with \code{Ltable} (4-column numeric matrix in DDD format),
//'   \code{status} (one of \code{"done"}, \code{"extinct"}, \code{"too_large"})
//'   and \code{marks}, each L-table row's mark at the end (0-based).
// [[Rcpp::export]]
Rcpp::List simulate_div_tree_cpp(Rcpp::NumericVector  pars,
                                  Rcpp::IntegerVector  model,
                                  double               max_t,
                                  int                  max_N,
                                  int                  max_tries,
                                  int                  link = 0,
                                  int                  seed = 0,
                                  int                  marks_S = 1,
                                  Rcpp::NumericVector  log_clam = Rcpp::NumericVector::create(0.0),
                                  Rcpp::NumericVector  log_cmu = Rcpp::NumericVector::create(0.0),
                                  Rcpp::NumericVector  marks_C = Rcpp::NumericVector::create(1.0),
                                  Rcpp::NumericVector  marks_P = Rcpp::NumericVector::create(),
                                  Rcpp::NumericVector  marks_pi = Rcpp::NumericVector::create(1.0)) {
  const uint64_t rseed = emphasis::resolve_seed(seed);
  // 8, 10 or 12 parameters, 3 to 5 model flags: the ED and K slots are
  // appended and default to absent (model.hpp).
  if (pars.size() != 8 && pars.size() != 10 && pars.size() != 12) {
    Rcpp::stop("simulate_div_tree_cpp: pars must have length 8, 10 or 12 (got %d)", (int) pars.size());
  }
  if (model.size() < 3 || model.size() > 5) {
    Rcpp::stop("simulate_div_tree_cpp: model must have length 3, 4 or 5 (got %d)", (int) model.size());
  }
  if (link != 0 && link != 1) {
    Rcpp::stop("simulate_div_tree_cpp: link must be 0 (linear) or 1 (exponential) (got %d)", link);
  }
  std::array<double, 12> p{};
  for (int i = 0; i < pars.size(); ++i) p[static_cast<size_t>(i)] = pars[i];
  std::array<int, 5> m{};
  for (int i = 0; i < model.size(); ++i) m[static_cast<size_t>(i)] = model[i];

  sim_tree::general_div sim(max_t, p, m, static_cast<size_t>(max_N), link, rseed);
  if (marks_S > 1) {
    const size_t S = static_cast<size_t>(marks_S);
    const std::vector<double> lcl(log_clam.begin(), log_clam.end()), lcm(log_cmu.begin(), log_cmu.end()),
      mC(marks_C.begin(), marks_C.end()), mP(marks_P.begin(), marks_P.end()), mpi(marks_pi.begin(), marks_pi.end());
    if (link != 1) Rcpp::stop("simulate_div_tree_cpp: latent marks need the exponential link");
    if (lcl.size() != S || lcm.size() != S || mpi.size() != S ||
        (mP.empty() && mC.size() != S * S * S) || (!mP.empty() && mP.size() != S * S))
      Rcpp::stop("simulate_div_tree_cpp: marks do not match marks_S");
    sim.set_marks(marks_S, lcl, lcm, mC, mP, mpi);
  }
  sim.simulate_tree_ltable();

  int tries = 0;
  while (sim.break_type == sim_tree::extinction && tries < max_tries) {
    sim.simulate_tree_ltable();
    ++tries;
  }

  const auto& tree = sim.ltable;
  Rcpp::NumericMatrix out(static_cast<int>(tree.size()), 4);
  const double crown_age = sim.max_t;

  for (size_t i = 0; i < tree.size(); ++i) {
    out(static_cast<int>(i), 0) = crown_age - static_cast<double>(tree[i].start_date);
    out(static_cast<int>(i), 1) = tree[i].parent_label;
    out(static_cast<int>(i), 2) = tree[i].label;
    const double ed = static_cast<double>(tree[i].end_date);
    out(static_cast<int>(i), 3) = (ed >= 0.0) ? (crown_age - ed) : -1.0;
  }

  std::string status = "done";
  if (sim.break_type == sim_tree::extinction)    status = "extinct";
  if (sim.break_type == sim_tree::maxN_exceeded) status = "too_large";

  Rcpp::IntegerVector marks(static_cast<int>(tree.size()));
  for (size_t i = 0; i < tree.size(); ++i) marks[static_cast<int>(i)] = tree[i].mark;

  return Rcpp::List::create(Rcpp::Named("Ltable") = out,
                             Rcpp::Named("status") = status,
                             Rcpp::Named("marks")  = marks);
}
