
#include <stdexcept>
#include <string>
#include <Rcpp.h>
#include "emphasis.hpp"
#include "model.hpp"
#include "rinit.h"
using namespace Rcpp;

// src/rcpp_mce.cpp — the `seed` argument into the integer the samplers run on.
namespace emphasis { uint64_t resolve_seed(int seed); }


namespace {

  DataFrame unpack(const emphasis::tree_t& tree)
  {
    NumericVector brts, n, t_ext;
    for (const emphasis::node_t& node : tree) {
      brts.push_back(node.brts);
      n.push_back(node.n);
      t_ext.push_back(node.t_ext);
    }
    return DataFrame::create(Named("brts") = brts, Named("n") = n, Named("t_ext") = t_ext);
  }

}

//' function to perform one step of the E-M algorithm
//' @param brts vector of branching times
//' @param init_pars vector of initial parameter files
//' @param sample_size number of samples
//' @param maxN maximum number of failed trees
//' @param max_missing maximum number of species missing
//' @param max_lambda maximum speciation rate
//' @param lower_bound vector of lower bound values for optimization, should
//' be equal in length to the vector of init_pars
//' @param upper_bound vector of upper bound values for optimization, should
//' be equal in length to the vector of init_pars
//' @param xtol_rel relative tolerance for optimization
//' @param num_threads number of threads used.
//' @param copy_trees if set to true, the trees generated are returned as well
//' @param model integer vector of length 3: c(use_N, use_P, use_E)
//' @param link link function: 0 = linear (max(0,...)), 1 = exponential
//' @param rho sampling fraction in (0, 1]
//' @param parent_tip_start tip start of the lineage that splits at each
//' observed branching event, in the same forward-time order as \code{brts};
//' empty (the default) when only branching times are available, in which case
//' every observed lineage is recorded as dating from the crown and D = 0 at
//' every observed event.
//' @param seed positive integer seeding the C++ sampler; 0 (the default) draws
//' one from R's generator, so set.seed() reaches the sampler whether or not a
//' seed is passed.
//' @param rconditional R function that evaluates the GAM function.
//' @return a list with the following components:
//' \describe{
//'  \item{trees}{list of trees}
//'  \item{rejected}{number of rejected trees}
//'  \item{rejected_overruns}{number of trees rejected due to too large size}
//'  \item{rejected_lambda}{number of trees rejected due to lambda errors}
//'  \item{rejected_zero_weights}{number of trees rejected due to zero weight (log weight -Inf)}
//'  \item{rejected_nonfinite}{number of trees rejected due to a +Inf or NaN log weight}
//'  \item{num_trees}{number of trees the E-step returned}
//'  \item{estimates}{vector of estimates}
//'  \item{nlopt}{nlopt status}
//'  \item{fhat}{vector of fhat values}
//'  \item{time}{time elapsed}
//'  \item{weights}{vector of weights}
//'  \item{logf}{vector of log p(obs, z_i | theta) for each valid tree}
//'  \item{logg}{vector of log q(z_i | obs, theta) for each valid tree}
//' }
//' @keywords internal
// [[Rcpp::export(name = "em_cpp")]]
List rcpp_mcem(const std::vector<double>& brts,
               const std::vector<double>& init_pars,
               int sample_size,
               int maxN,
               int max_missing,
               double max_lambda,
               const std::vector<double>& lower_bound,
               const std::vector<double>& upper_bound,
               double xtol_rel,
               int num_threads,
               bool copy_trees,
               Rcpp::IntegerVector model = Rcpp::IntegerVector::create(0, 0, 0),
               int link = 0,
               double rho = 1.0,
               Nullable<Function> rconditional = R_NilValue,
               Rcpp::NumericVector parent_tip_start = Rcpp::NumericVector::create(),
               int seed = 0,
               Rcpp::IntegerVector parent_id = Rcpp::IntegerVector::create())
{
  emphasis::rng::set_seed(emphasis::resolve_seed(seed));
  const std::vector<double> pts(parent_tip_start.begin(), parent_tip_start.end());
  const std::vector<int> pid(parent_id.begin(), parent_id.end());
  // 8 slots means ED absent; the internal layout is 10 (model.hpp).
  const bool ok_len = (init_pars.size() == 8 || init_pars.size() == emphasis::n_params) &&
                      init_pars.size() == lower_bound.size() && init_pars.size() == upper_bound.size();
  if (!ok_len) {
    throw std::invalid_argument("em_cpp: init_pars, lower_bound and upper_bound must have length 8 or 10, equal (got " +
      std::to_string(init_pars.size()) + ", " + std::to_string(lower_bound.size()) + ", " +
      std::to_string(upper_bound.size()) + ")");
  }
  const std::vector<double> init_pars10 = emphasis::pad_params(init_pars);
  const std::vector<double> lower10 = emphasis::pad_params(lower_bound);
  const std::vector<double> upper10 = emphasis::pad_params(upper_bound);
  std::vector<int> model_bin(model.begin(), model.end());
  // The proposal reads M = P/N at the segment start. Without the observed
  // topology P is the legacy quantity (every observed lineage dated from the
  // crown), which the scorer cannot reproduce from the finished tree, so
  // log q would not be the density the sampler drew from. cr, dd, d and nd
  // never reach this; an M-active model passed as a binary vector does.
  if (model_bin[1] == 1 && pts.empty()) {
    throw std::invalid_argument(
      "em_cpp: an M-dependent model (model[2] == 1) needs the observed topology. "
      "Pass parent_tip_start; without it log q is not the density the sampler "
      "draws from.");
  }
  if (model_bin.size() > 3 && model_bin[3] == 1 && (pid.empty() || pts.empty())) {
    throw std::invalid_argument(
      "em_cpp: the ED covariate needs the tree's topology: pass a phylo object "
      "(or a simulate_tree() result), not a bare branching-time vector.");
  }
  auto mdl = emphasis::Model(lower10, upper10, model_bin, link, rho);

  emphasis::conditional_fun_t conditional{};
  if (rconditional.isNotNull()) {
    // The R-side conditional receives the layout the caller used.
    const std::size_t n_out = init_pars.size();
    conditional = [cond= Function(rconditional), n_out](const emphasis::param_t& pars) {
      return as<double>( cond(NumericVector(pars.cbegin(), pars.cbegin() + static_cast<long>(n_out))) );
    };
  }
  auto mcem = emphasis::mcem(sample_size,
                             maxN,
                             init_pars10,
                             brts,
                             mdl,
                             max_missing,
                             max_lambda,
                             lower10,
                             upper10,
                             xtol_rel,
                             num_threads,
                             conditional ? &conditional : nullptr,
                             pts,
                             pid);
  if (mcem.e.trees.empty()) {
    throw std::runtime_error("no trees, no optimization");
  }
  List ret;
  if (copy_trees) {
    List trees;
    for (const emphasis::tree_t& tree : mcem.e.trees) {
      trees.push_back(unpack(tree));
    }
    ret["trees"] = trees;
  } else {
    ret["trees"] = static_cast<int>(mcem.e.trees.size());
  }
  ret["rejected"] = mcem.e.info.rejected;
  ret["rejected_overruns"] = mcem.e.info.rejected_overruns;
  ret["rejected_lambda"] = mcem.e.info.rejected_lambda;
  ret["rejected_zero_weights"] = mcem.e.info.rejected_zero_weights;
  ret["rejected_nonfinite"] = mcem.e.info.rejected_nonfinite;
  ret["num_trees"] = mcem.e.info.num_trees;
  {
    // Returned in the layout the caller used: 8 slots when it passed 8.
    const long n_out = static_cast<long>(std::min(init_pars.size(), mcem.m.estimates.size()));
    ret["estimates"] = NumericVector(mcem.m.estimates.begin(), mcem.m.estimates.begin() + n_out);
  }
  ret["nlopt"] = mcem.m.opt;
  ret["fhat"]  = mcem.e.info.fhat;
  ret["time"]  = mcem.e.info.elapsed + mcem.m.elapsed;
  ret["weights"] = mcem.e.weights;
  ret["logf"] = NumericVector(mcem.e.logf_.begin(), mcem.e.logf_.end());
  ret["logg"] = NumericVector(mcem.e.logg_.begin(), mcem.e.logg_.end());
  return ret;
}
