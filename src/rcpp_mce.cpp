
#include <Rcpp.h>
#include <array>
#include <vector>
#include <tbb/tbb.h>
#include "emphasis.hpp"
#include "model.hpp"
#include "rinit.h"
#include "unpack.h"
using namespace Rcpp;


//' Draw augmented trees via importance sampling
//'
//' Augments an observed extant tree (given by its branching times) with
//' stochastically drawn extinct lineages, and returns the simulated trees
//' together with their log proposal probabilities.  Scoring (\code{logf})
//' and IS aggregation (\code{fhat}) are intentionally left to the R layer
//' via \code{\link{eval_logf}} and \code{.is_fhat}.
//'
//' @param brts Numeric vector of branching times (crown age first,
//'   sorted decreasing).
//' @param pars Numeric vector of 8 model parameters
//'   \code{c(beta_0, beta_N, beta_P, beta_E, gamma_0, gamma_N, gamma_P, gamma_E)}.
//' @param sample_size Number of valid augmented trees to collect.
//' @param maxN Maximum total augmentation attempts (including failures).
//' @param max_missing Maximum extinct lineages per augmented tree.
//' @param max_lambda Upper bound on the speciation rate (thinning bound).
//' @param num_threads Threads for parallel augmentation.
//' @param model Integer vector \code{c(use_N, use_P, use_E)}.
//' @param link Link function: \code{0} = linear, \code{1} = exponential.
//' @return A named list:
//' \describe{
//'   \item{trees}{List of augmented-tree data frames.}
//'   \item{logf}{Per-tree log p(obs, z | theta).}
//'   \item{logg}{Per-tree log q(z | obs, theta).}
//'   \item{rejected}{Unhandled rejections.}
//'   \item{rejected_overruns}{Rejected: too many extinct lineages.}
//'   \item{rejected_lambda}{Rejected: lambda bound exceeded.}
//'   \item{rejected_zero_weights}{Rejected: zero IS weight.}
//'   \item{time}{Elapsed time (ms).}
//' }
//' @keywords internal
// [[Rcpp::export(name = "augment_trees")]]
List rcpp_mce(const std::vector<double>& brts,
              const std::vector<double>& pars,
              int sample_size,
              int maxN,
              int max_missing,
              double max_lambda,
              int num_threads,
              Rcpp::IntegerVector model = Rcpp::IntegerVector::create(0, 0, 0),
              int link = 0,
              double rho = 1.0)
{
  std::vector<int> model_bin = {model[0], model[1], model[2]};
  // Bounds are only needed for M-step (nlopt); E-step does not use them.
  std::vector<double> lb8(8, -1e6), ub8(8, 1e6);
  auto mdl = emphasis::Model(lb8, ub8, model_bin, link, rho);

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
  ret["trees"]                 = trees;
  ret["logf"]                  = E.logf_;
  ret["logg"]                  = E.logg_;
  ret["rejected"]              = E.info.rejected;
  ret["rejected_overruns"]     = E.info.rejected_overruns;
  ret["rejected_lambda"]       = E.info.rejected_lambda;
  ret["rejected_zero_weights"] = E.info.rejected_zero_weights;
  ret["time"]                  = E.info.elapsed;
  return ret;
}
