#include "general_tree.hpp"

#include <Rcpp.h>

//' Simulate a single phylogenetic tree under the general diversification model
//'
//' @param pars Numeric vector of length 8:
//'   \code{c(beta_0, beta_N, beta_P, beta_E, gamma_0, gamma_N, gamma_P, gamma_E)}.
//'   Inactive parameters (those not selected by \code{model}) should be set to 0.
//' @param model Integer vector of length 3: \code{c(use_N, use_P, use_E)} (each 0 or 1).
//' @param max_t Crown age (forward simulation end time).
//' @param max_N Maximum number of lineages before simulation is declared too large.
//' @param max_tries Maximum retries after extinction or overflow.
//' @param link Link function: 0 = linear (max(0,...)), 1 = exponential.
//' @return A named list with \code{Ltable} (4-column numeric matrix in DDD format)
//'   and \code{status} (one of \code{"done"}, \code{"extinct"}, \code{"too_large"}).
// [[Rcpp::export]]
Rcpp::List simulate_div_tree_cpp(Rcpp::NumericVector  pars,
                                  Rcpp::IntegerVector  model,
                                  double               max_t,
                                  int                  max_N,
                                  int                  max_tries,
                                  int                  link = 0) {
  std::array<double, 8> p = {
    pars[0], pars[1], pars[2], pars[3],
    pars[4], pars[5], pars[6], pars[7]
  };
  std::array<int, 3> m = { model[0], model[1], model[2] };

  sim_tree::general_div sim(max_t, p, m, static_cast<size_t>(max_N), link);
  sim.simulate_tree_ltable();

  int tries = 0;
  while ((sim.break_type == sim_tree::extinction ||
          sim.break_type == sim_tree::maxN_exceeded) &&
         tries < max_tries) {
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

  return Rcpp::List::create(Rcpp::Named("Ltable") = out,
                             Rcpp::Named("status") = status);
}
