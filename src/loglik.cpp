// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>
#include "emphasis.hpp"
#include "plugin.hpp"
#include "rinit.h"
using namespace Rcpp;

namespace loglik {

emphasis::tree_t pack(const Rcpp::DataFrame& r_tree) {
  emphasis::tree_t new_tree;
  
  Rcpp::NumericVector brts = r_tree["brts"];
  Rcpp::NumericVector n = r_tree["n"];
  Rcpp::NumericVector t_ext = r_tree["t_ext"];
  
  for (size_t i = 0; i < r_tree.nrow(); ++i) {
    emphasis::node_t entry;
    entry.brts = brts[i];
    entry.n    = n[i];
    entry.t_ext = t_ext[i];
    new_tree.push_back(entry);
  }
  return new_tree;
}

}

//' function to calculate log likelihood of pars for a tree set,
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector loglikelihood(const Rcpp::NumericVector& pars_r,
                                  const Rcpp::List& trees,
                                  const std::string& plugin) {
  
  auto model = emphasis::create_plugin_model(plugin);  
  
  Rcpp::NumericVector logf(trees.size());
  
  std::vector<double> pars(pars_r.begin(), pars_r.end()); // not sure if this can't be sugared directly on input.
  
  // unpack list and convert to tree_t
  for (size_t i = 0; i < trees.size(); ++i) {
    auto local_tree = loglik::pack(Rcpp::as<Rcpp::DataFrame>(trees[i]));
    logf[i] = model->loglik(pars, local_tree);
  }
  
  return logf;
}