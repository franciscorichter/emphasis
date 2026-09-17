
#include <stdexcept>
#include <string>
#include <Rcpp.h>
#include "emphasis.hpp"
#include "model.hpp"
#include "rinit.h"
using namespace Rcpp;


namespace {

  // The covariate columns are carried over when the data frame has them: the
  // M-step scores the same trees the E-step drew, and pd / tip_start /
  // focal_tip_start are what M and D are read from.  Dropping them here made
  // every tree reaching the M-step look like one with P = 0 and no parent.
  std::vector<emphasis::tree_t> pack(List rtrees)
  {
    std::vector<emphasis::tree_t> trees;
    for (auto it = rtrees.cbegin(); it != rtrees.cend(); ++it) {
      emphasis::tree_t tree;
      auto df = DataFrame(*it);
      auto brts = as<NumericVector>(df["brts"]);
      auto n = as<NumericVector>(df["n"]);
      auto t_ext = as<NumericVector>(df["t_ext"]);
      NumericVector pd, tip_start, focal_tip_start;
      IntegerVector clade, id, parent_id;
      const bool has_pd        = df.containsElementNamed("pd");
      const bool has_ts        = df.containsElementNamed("tip_start");
      const bool has_focal_ts  = df.containsElementNamed("focal_tip_start");
      const bool has_clade     = df.containsElementNamed("clade");
      const bool has_id        = df.containsElementNamed("id");
      const bool has_parent_id = df.containsElementNamed("parent_id");
      if (has_pd)        pd = as<NumericVector>(df["pd"]);
      if (has_ts)        tip_start = as<NumericVector>(df["tip_start"]);
      if (has_focal_ts)  focal_tip_start = as<NumericVector>(df["focal_tip_start"]);
      if (has_clade)     clade = as<IntegerVector>(df["clade"]);
      if (has_id)        id = as<IntegerVector>(df["id"]);
      if (has_parent_id) parent_id = as<IntegerVector>(df["parent_id"]);
      for (auto i = 0; i < brts.size(); ++i) {
        const int pid = has_parent_id ? parent_id[i] : -1;
        emphasis::node_t node{};
        node.brts = brts[i];
        node.n = n[i];
        node.t_ext = t_ext[i];
        node.pd = has_pd ? pd[i] : 0.0;
        node.tip_start = has_ts ? tip_start[i] : 0.0;
        node.focal_tip_start = has_focal_ts ? focal_tip_start[i]
                             : (emphasis::has_parent(pid) ? 0.0
                                                          : emphasis::ts_unknown);
        node.clade = has_clade ? clade[i] : 0;
        node.id = has_id ? id[i] : -1;
        node.parent_id = pid;
        tree.push_back(node);
      }
      trees.emplace_back(std::move(tree));
    }
    return trees;
  }

}

//' function to perform one step of the E-M algorithm
//' @param e_step result of e_step function, as a list
//' @param init_pars vector of initial parameter values
//' @param plugin string indicating plugin used, currently available: 'rpd1' and
//' 'rpd5c'
//' @param lower_bound vector of lower bound values for optimization, should
//' be equal in length to the vector of init_pars
//' @param upper_bound vector of upper bound values for optimization, should
//' be equal in length to the vector of init_pars
//' @param xtol_rel relative tolerance for optimization
//' @param num_threads number of threads used.
//' @param model integer vector of length 3: c(use_N, use_P, use_E)
//' @param link link function: 0 = linear (max(0,...)), 1 = exponential
//' @param rconditional R function that evaluates the GAM function.
//' @return list with the following entries:
//' \describe{
//'  \item{estimates}{vector of estimates}
//'  \item{nlopt}{nlopt status}
//'  \item{time}{used computation time}
//' }
//' @keywords internal
// [[Rcpp::export(name = "m_cpp")]]
List rcpp_mcm(List e_step,
              const std::vector<double>& init_pars,
              const std::string& plugin,
              const std::vector<double>& lower_bound,
              const std::vector<double>& upper_bound,
              double xtol_rel,
              int num_threads,
              Rcpp::IntegerVector model = Rcpp::IntegerVector::create(0, 0, 0),
              int link = 0,
              double rho = 1.0,
              Nullable<Function> rconditional = R_NilValue)
{
  // The internal parameter layout has 10 slots (the last two the ED
  // coefficients); an 8-element vector means ED absent and is padded.
  const bool ok_len = (init_pars.size() == 8 || init_pars.size() == emphasis::n_params) &&
                      init_pars.size() == lower_bound.size() && init_pars.size() == upper_bound.size();
  if (!ok_len) {
    throw std::invalid_argument("m_cpp: init_pars, lower_bound and upper_bound must have length 8 or 10, equal (got " +
      std::to_string(init_pars.size()) + ", " + std::to_string(lower_bound.size()) + ", " +
      std::to_string(upper_bound.size()) + ")");
  }
  const std::vector<double> init_pars10 = emphasis::pad_params(init_pars);
  const std::vector<double> lower10 = emphasis::pad_params(lower_bound);
  const std::vector<double> upper10 = emphasis::pad_params(upper_bound);
  auto E = emphasis::E_step_t{};
  E.trees = pack(as<List>(e_step["trees"]));
  E.weights = as<std::vector<double>>(e_step["weights"]);
  if (E.weights.size() != E.trees.size()) {
    throw std::invalid_argument("m_cpp: weights (" + std::to_string(E.weights.size()) +
      ") and trees (" + std::to_string(E.trees.size()) + ") differ in length");
  }
  E.info.num_trees = static_cast<int>(E.trees.size());
  E.info.rejected = as<int>(e_step["rejected"]);
  E.info.rejected_overruns = as<int>(e_step["rejected_overruns"]);
  E.info.rejected_lambda = as<int>(e_step["rejected_lambda"]);
  E.info.rejected_zero_weights = as<int>(e_step["rejected_zero_weights"]);
  E.info.elapsed = as<double>(e_step["time"]);
  E.info.fhat = as<double>(e_step["fhat"]);

  if (E.trees.empty()) {
    throw std::runtime_error("no trees, no optimization");
  }
  std::vector<int> model_bin(model.begin(), model.end());
  auto mdl = emphasis::Model(lower10, upper10, model_bin, link, rho);
  emphasis::conditional_fun_t conditional{};
  if (rconditional.isNotNull()) {
    // The R-side conditional receives the layout it was written for: 8
    // slots when the caller passed 8, the full 10 otherwise.
    const std::size_t n_out = init_pars.size();
    conditional = [cond= Function(rconditional), n_out](const emphasis::param_t& pars) {
      return as<double>( cond(NumericVector(pars.cbegin(), pars.cbegin() + static_cast<long>(n_out))) );
    };
  }
  auto M = emphasis::M_step(init_pars10,
                            E.trees,
                            E.weights,
                            mdl,
                            lower10,
                            upper10,
                            xtol_rel,
                            num_threads,
                            conditional ? &conditional : nullptr);
  List ret;
  // Returned in the layout the caller used: 8 slots when it passed 8.
  const long n_out = static_cast<long>(std::min(init_pars.size(), M.estimates.size()));
  ret["estimates"] = NumericVector(M.estimates.begin(), M.estimates.begin() + n_out);
  ret["nlopt"] = M.opt;
  ret["time"]  = M.elapsed;
  return ret;
}
