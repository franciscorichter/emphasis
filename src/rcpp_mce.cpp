// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>
#include <array>
#include <vector>
#include <tbb/tbb.h>
#include "emphasis.hpp"
#include "model.hpp"
#include "rinit.h"
#include "unpack.h"
using namespace Rcpp;


//' Monte Carlo log-likelihood via importance sampling with augmented trees
//'
//' Augments an observed tree (represented by its branching times) with
//' stochastically drawn extinct lineages, computes importance-sampling
//' weights, and returns a Monte Carlo estimate of the log-likelihood.
//'
//' @param brts Numeric vector of branching times (crown age first,
//'   sorted decreasing).
//' @param pars Numeric vector of model parameters.
//' @param sample_size Number of valid augmented trees to collect.
//' @param maxN Maximum total augmentation attempts (including failures).
//' @param max_missing Maximum number of extinct lineages per tree.
//' @param max_lambda Upper bound on the speciation rate (thinning bound).
//' @param lower_bound Lower bounds for parameter optimisation (same length
//'   as \code{pars}).
//' @param upper_bound Upper bounds for parameter optimisation.
//' @param xtol_rel Relative tolerance for internal optimisation.
//' @param num_threads Number of threads for parallel augmentation.
//' @param model Integer vector of length 3: c(use_N, use_P, use_E).
//' @param link Link function: 0 = linear (max(0,...)), 1 = exponential.
//' @return A named list:
//' \describe{
//'   \item{trees}{List of augmented-tree data frames (brts, n, t_ext, pd,
//'     id, parent_id).}
//'   \item{weights}{Numeric vector of log importance weights
//'     (\code{logf - logg}).}
//'   \item{fhat}{Monte Carlo log-likelihood estimate.}
//'   \item{logf}{Per-tree log-likelihood values.}
//'   \item{logg}{Per-tree log sampling probabilities.}
//'   \item{rejected}{Number of trees rejected (unhandled).}
//'   \item{rejected_overruns}{Rejected due to too many missing lineages.}
//'   \item{rejected_lambda}{Rejected due to lambda bound exceeded.}
//'   \item{rejected_zero_weights}{Rejected due to zero weight.}
//'   \item{time}{Elapsed time in milliseconds.}
//' }
//' @export
// [[Rcpp::export(name = "mc_loglik")]]
List rcpp_mce(const std::vector<double>& brts,
              const std::vector<double>& pars,
              int sample_size,
              int maxN,
              int max_missing,
              double max_lambda,
              const std::vector<double>& lower_bound,
              const std::vector<double>& upper_bound,
              double xtol_rel,
              int num_threads,
              Rcpp::IntegerVector model = Rcpp::IntegerVector::create(0, 0, 0),
              int link = 0)
{
  std::vector<int> model_bin = {model[0], model[1], model[2]};
  auto mdl = emphasis::Model(lower_bound, upper_bound, model_bin, link);

  auto E = emphasis::E_step(sample_size,
                            maxN,
                            pars,
                            brts,
                            mdl,
                            max_missing,
                            max_lambda,
                            num_threads);
  List ret;
  List trees;
  for (const emphasis::tree_t& tree : E.trees) {
    trees.push_back(unpack(tree));
  }
  ret["trees"] = trees;
  ret["rejected"] = E.info.rejected;
  ret["rejected_overruns"] = E.info.rejected_overruns;
  ret["rejected_lambda"] = E.info.rejected_lambda;
  ret["rejected_zero_weights"] = E.info.rejected_zero_weights;
  ret["time"] = E.info.elapsed;
  ret["weights"] = E.weights;
  ret["fhat"] = E.info.fhat;
  ret["logf"] = E.logf_;
  ret["logg"] = E.logg_;
  return ret;
}
