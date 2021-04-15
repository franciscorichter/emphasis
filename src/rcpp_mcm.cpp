// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>
#include "emphasis.hpp"
#include "plugin.hpp"
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
  E.rejected = as<int>(e_step["rejected"]);
  E.rejected_overruns = as<int>(e_step["rejected_overruns"]);
  E.rejected_lambda = as<int>(e_step["rejected_lambda"]);
  E.rejected_zero_weights = as<int>(e_step["rejected_zero_weights"]);
  E.elapsed = as<double>(e_step["time"]);
  E.fhat = as<double>(e_step["fhat"]);

  if (E.trees.empty()) {
    throw std::runtime_error("no trees, no optimization");
  }
  auto model = emphasis::create_plugin_model(plugin);
  emphasis::conditional_fun_t conditional{};
  if (rconditional.isNotNull()) {
    conditional = [cond= Function(rconditional)](const emphasis::param_t& pars) {
      return as<double>( cond(NumericVector(pars.cbegin(), pars.cend())) );
    };
  }
  auto M = emphasis::M_step(init_pars, 
                            E.trees, 
                            E.weights, 
                            model.get(), 
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
