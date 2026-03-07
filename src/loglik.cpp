
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


//' Evaluate log p(obs, z | theta) and log q(z | obs, theta) for augmented trees
//'
//' Given a list of augmented trees (produced by \code{\link{augment_trees}})
//' and a parameter vector, computes for each tree:
//' \itemize{
//'   \item \code{logf[i] = log p(obs, z_i | theta)} (model log-likelihood)
//'   \item \code{logg[i] = log q(z_i | obs, theta)} (proposal log-probability)
//' }
//' Both are evaluated at the supplied \code{pars}, regardless of which
//' parameters were used to simulate the trees.  This is essential for
//' shared-tree (Mode 2) evaluation where trees simulated at one particle
//' are scored at another.
//'
//' @param pars Numeric vector of 8 model parameters.
//' @param trees List of augmented-tree data frames (output of
//'   \code{\link{augment_trees}}).
//' @param model Integer vector \code{c(use_N, use_P, use_E)}.
//' @param link Link function: \code{0} = linear, \code{1} = exponential.
//' @return A named list:
//' \describe{
//'   \item{logf}{Numeric vector of \code{log p(obs, z_i | theta)}.}
//'   \item{logg}{Numeric vector of \code{log q(z_i | obs, theta)}.}
//' }
//' @keywords internal
// [[Rcpp::export(name = "eval_logf")]]
Rcpp::List eval_logf_cpp(const std::vector<double>& pars,
                         const Rcpp::List& trees,
                         Rcpp::IntegerVector model = Rcpp::IntegerVector::create(0, 0, 0),
                         int link = 0) {
  std::vector<int> model_bin = {model[0], model[1], model[2]};
  emphasis::param_t lb8(8, -1e6), ub8(8, 1e6);
  auto mdl = emphasis::Model(lb8, ub8, model_bin, link);

  std::vector<double> logf(trees.size());
  std::vector<double> logg(trees.size());
  for (int i = 0; i < trees.size(); ++i) {
    auto local_tree = loglik::pack(Rcpp::as<Rcpp::DataFrame>(trees[i]));
    logf[i] = mdl.loglik(pars, local_tree);
    logg[i] = mdl.sampling_prob(pars, local_tree);
  }

  return Rcpp::List::create(
    Rcpp::Named("logf") = logf,
    Rcpp::Named("logg") = logg
  );
}