// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>
#include "emphasis.hpp"
#include "model.hpp"
#include "rinit.h"
#include "precision_weights.hpp"
using namespace Rcpp;

namespace loglik {

emphasis::tree_t pack(const Rcpp::DataFrame& r_tree) {
  emphasis::tree_t new_tree;
  
  // NOTE 10/8/22: this can probably be micro-optimized, as it would make more 
  // sense to avoid copying everything, and instead traverse the DataFrame as
  // a matrix in a row wise fashion.
  
  Rcpp::NumericVector brts  = r_tree["brts"];
  Rcpp::NumericVector n     = r_tree["n"];
  Rcpp::NumericVector t_ext = r_tree["t_ext"];
  Rcpp::NumericVector pd    = r_tree["pd"];
  
  for (size_t i = 0; i < r_tree.nrow(); ++i) {
    emphasis::node_t entry;
    entry.brts = brts[i];
    entry.n    = n[i];
    entry.t_ext = t_ext[i];
    entry.pd    = pd[i];
    new_tree.push_back(entry);
  }
  return new_tree;
}

}

//' function to calculate log likelihood of pars for a tree set,
//' @param pars vector of parameter values
//' @param trees list of trees, e.g. a multiPhylo object
//' @param logg vector of logg values
//' @param plugin name of used plugin
//' @param num_rejected number of rejected trees
//' @return list with the following entries: 
//' \itemize{
//'  \item{logf}{logf values}
//'  \item{log_w}{log of weight}
//'  \item{fhat}{fhat values}
//'  \item{N}{number of trees}
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List loglikelihood(const std::vector<double>& pars,
                         const Rcpp::List& trees,
                         const Rcpp::NumericVector& logg,
                         const std::string& plugin,
                         const int num_rejected) {
  
  auto model = emphasis::create_model();  
  
  std::vector<double> logf(trees.size());
  std::vector<double> log_w(trees.size());
  
  // unpack list and convert to tree_t
  for (size_t i = 0; i < trees.size(); ++i) {
    auto local_tree = loglik::pack(Rcpp::as<Rcpp::DataFrame>(trees[i]));
    logf[i] = model->loglik(pars, local_tree);
    log_w[i] = logf[i] - logg[i];
  }

  double fhat = logf.front();
  
  if (trees.size() > 1) {
  
    const double max_log_w = *std::max_element(log_w.cbegin(), log_w.cend());
    double sum_w = calc_sum_w(log_w.begin(), log_w.end(), max_log_w);
    
    fhat = std::log(sum_w / (trees.size() + num_rejected)) + max_log_w;
  }
  
  Rcpp::List ret;
  ret["logf"]    = logf;
  ret["weights"] = log_w;
  ret["fhat"]    = fhat;
  ret["N"]       = trees.size() + num_rejected;
  
  return ret;
}