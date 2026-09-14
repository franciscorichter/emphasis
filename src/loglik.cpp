
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
  Rcpp::NumericVector tip_start_v, focal_tip_start_v;
  Rcpp::IntegerVector clade_v, id_v, parent_id_v;
  bool has_tip_start = r_tree.containsElementNamed("tip_start");
  bool has_focal_ts  = r_tree.containsElementNamed("focal_tip_start");
  bool has_clade     = r_tree.containsElementNamed("clade");
  bool has_id        = r_tree.containsElementNamed("id");
  bool has_parent_id = r_tree.containsElementNamed("parent_id");
  if (has_tip_start) tip_start_v = r_tree["tip_start"];
  if (has_focal_ts)  focal_tip_start_v = r_tree["focal_tip_start"];
  if (has_clade)     clade_v = r_tree["clade"];
  if (has_id)        id_v = r_tree["id"];
  if (has_parent_id) parent_id_v = r_tree["parent_id"];

  for (int i = 0; i < r_tree.nrow(); ++i) {
    emphasis::node_t entry;
    entry.brts      = brts[i];
    entry.n         = n[i];
    entry.t_ext     = t_ext[i];
    entry.pd        = pd[i];
    entry.tip_start = has_tip_start ? tip_start_v[i] : 0.0;
    entry.id        = has_id ? id_v[i] : -1;
    entry.parent_id = has_parent_id ? parent_id_v[i] : -1;
    // A tree built without the column carries no splitting lineage, so the
    // event is mean-field (D = 0) exactly where it was before: at a node whose
    // parent is not on record.
    entry.focal_tip_start = has_focal_ts ? focal_tip_start_v[i]
                          : (emphasis::has_parent(entry.parent_id) ? 0.0
                                                                  : emphasis::ts_unknown);
    // Absent, the tree is read under the legacy convention, which is what a
    // data frame built without the column has always meant.
    entry.clade     = has_clade ? clade_v[i] : 0;
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
                         int link = 0,
                         double rho = 1.0) {
  std::vector<int> model_bin = {model[0], model[1], model[2]};
  emphasis::param_t lb8(8, -1e6), ub8(8, 1e6);
  auto mdl = emphasis::Model(lb8, ub8, model_bin, link, rho);

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


//' The pendant PD and thinning rate the sampler sees at arbitrary times
//'
//' \code{Model::pendant_pd_at} reads P off the node that governs the segment a
//' candidate time falls in; \code{Model::nh_rate} is the intensity of the
//' sampler's birth process there.  This exposes both, so a test can hold P
//' against an independent recomputation from the augmented tree and against
//' \code{pd + n * (t - brts)} of the governing node, and can scan the intensity
//' over a segment to check that its value at the segment start dominates it.
//'
//' @param pars Numeric vector of 8 model parameters.
//' @param tree One augmented-tree data frame.
//' @param times Numeric vector of times at which to evaluate.
//' @param model Integer vector \code{c(use_N, use_M, use_D)}.
//' @param link Link function: 0 = linear, 1 = exponential, 2 = gaussian.
//' @param rho Sampling fraction.
//' @return A named list with \code{pd} (the pendant PD used) and \code{nh}
//'   (the non-homogeneous thinning rate) at each time.
//' @keywords internal
// [[Rcpp::export(name = "eval_nh_rate")]]
Rcpp::List eval_nh_rate_cpp(const std::vector<double>& pars,
                            const Rcpp::DataFrame& tree,
                            const std::vector<double>& times,
                            Rcpp::IntegerVector model = Rcpp::IntegerVector::create(0, 0, 0),
                            int link = 0,
                            double rho = 1.0) {
  std::vector<int> model_bin = {model[0], model[1], model[2]};
  emphasis::param_t lb8(8, -1e6), ub8(8, 1e6);
  auto mdl = emphasis::Model(lb8, ub8, model_bin, link, rho);
  auto local_tree = loglik::pack(tree);

  std::vector<double> pd(times.size()), nh(times.size());
  for (size_t i = 0; i < times.size(); ++i) {
    pd[i] = mdl.pendant_pd_at(times[i], local_tree);
    nh[i] = mdl.nh_rate(times[i], pars, local_tree);
  }
  return Rcpp::List::create(Rcpp::Named("pd") = pd, Rcpp::Named("nh") = nh);
}
