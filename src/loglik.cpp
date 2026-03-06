// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>
#include "emphasis.hpp"
#include "model.hpp"
#include "rinit.h"
using namespace Rcpp;

namespace loglik {

emphasis::tree_t pack(const Rcpp::DataFrame& r_tree) {
  emphasis::tree_t new_tree;
  Rcpp::NumericVector brts      = r_tree["brts"];
  Rcpp::NumericVector n         = r_tree["n"];
  Rcpp::NumericVector t_ext     = r_tree["t_ext"];
  Rcpp::NumericVector pd        = r_tree["pd"];
  Rcpp::NumericVector tip_start_v;
  bool has_tip_start = r_tree.containsElementNamed("tip_start");
  if (has_tip_start) tip_start_v = r_tree["tip_start"];

  for (int i = 0; i < r_tree.nrow(); ++i) {
    emphasis::node_t entry;
    entry.brts      = brts[i];
    entry.n         = n[i];
    entry.t_ext     = t_ext[i];
    entry.pd        = pd[i];
    entry.tip_start = has_tip_start ? tip_start_v[i] : 0.0;
    entry.clade     = 0;
    entry.id        = -1;
    entry.parent_id = -1;
    new_tree.push_back(entry);
  }
  return new_tree;
}

}  // namespace loglik


//' Evaluate log p(obs, z | theta) for a set of pre-computed augmented trees
//'
//' Given a list of augmented trees (produced by \code{\link{augment_trees}})
//' and a parameter vector, computes \code{logf[i] = log p(obs, z_i | theta)}
//' for each tree.  IS aggregation (\code{fhat}) is left to the R layer via
//' \code{.is_fhat}.
//'
//' @param pars Numeric vector of 8 model parameters.
//' @param trees List of augmented-tree data frames (output of
//'   \code{\link{augment_trees}}).
//' @param model Integer vector \code{c(use_N, use_P, use_E)}.
//' @param link Link function: \code{0} = linear, \code{1} = exponential.
//' @return A named list with one element:
//' \describe{
//'   \item{logf}{Numeric vector of \code{log p(obs, z_i | theta)}, length =
//'     \code{length(trees)}.}
//' }
//' @export
// [[Rcpp::export(name = "eval_logf")]]
Rcpp::List eval_logf_cpp(const std::vector<double>& pars,
                         const Rcpp::List& trees,
                         Rcpp::IntegerVector model = Rcpp::IntegerVector::create(0, 0, 0),
                         int link = 0) {
  std::vector<int> model_bin = {model[0], model[1], model[2]};
  emphasis::param_t lb8(8, -1e6), ub8(8, 1e6);
  auto mdl = emphasis::Model(lb8, ub8, model_bin, link);

  std::vector<double> logf(trees.size());
  for (int i = 0; i < trees.size(); ++i) {
    auto local_tree = loglik::pack(Rcpp::as<Rcpp::DataFrame>(trees[i]));
    logf[i] = mdl.loglik(pars, local_tree);
  }

  return Rcpp::List::create(Rcpp::Named("logf") = logf);
}