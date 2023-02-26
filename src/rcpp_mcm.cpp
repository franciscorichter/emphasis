// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>
#include "emphasis.hpp"
#include "model.hpp"
#include "rinit.h"
using namespace Rcpp;


namespace {

  std::vector<emphasis::tree_t> pack(List rtrees)
  {
    std::vector<emphasis::tree_t> trees;
    for (auto it = rtrees.cbegin(); it != rtrees.cend(); ++it) {
      emphasis::tree_t tree;
      auto df = DataFrame(*it);
      auto brts = as<NumericVector>(df["brts"]);
      auto n = as<NumericVector>(df["n"]);
      auto t_ext = as<NumericVector>(df["t_ext"]);
      for (auto i = 0; i < brts.size(); ++i) {
        tree.push_back(emphasis::node_t{brts[i], n[i], t_ext[i], 0.0});
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
//' @param rconditional R function that evaluates the GAM function.
//' @return list with the following entries:
//' \itemize{
//'  \item{estimates}{vector of estimates}
//'  \item{nlopt}{nlopt status}
//'  \item{time}{used computation time} 
//' }
//' @export
// [[Rcpp::export(name = "m_cpp")]]
List rcpp_mcm(List e_step,       
              const std::vector<double>& init_pars,      
              const std::string& plugin,             
              const std::vector<double>& lower_bound,  
              const std::vector<double>& upper_bound,  
              double xtol_rel,                     
              int num_threads,
              Nullable<Function> rconditional = R_NilValue)
{
  auto E = emphasis::E_step_t{};
  E.trees = pack(as<List>(e_step["trees"]));
  E.weights = as<std::vector<double>>(e_step["weights"]);
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
  auto model = emphasis::Model(lower_bound, upper_bound);
  emphasis::conditional_fun_t conditional{};
  if (rconditional.isNotNull()) {
    conditional = [cond= Function(rconditional)](const emphasis::param_t& pars) {
      return as<double>( cond(NumericVector(pars.cbegin(), pars.cend())) );
    };
  }
  auto M = emphasis::M_step(init_pars, 
                            E.trees, 
                            E.weights, 
                            model, 
                            lower_bound, 
                            upper_bound, 
                            xtol_rel, 
                            num_threads, 
                            conditional ? &conditional : nullptr);
  List ret;
  ret["estimates"] = NumericVector(M.estimates.begin(), M.estimates.end());
  ret["nlopt"] = M.opt;
  ret["time"]  = M.elapsed;
  return ret;
}
