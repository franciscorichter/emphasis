// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>
#include "emphasis.hpp"
#include "plugin.hpp"
#include "rinit.h"
using namespace Rcpp;


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


// [[Rcpp::export(name = "e_cpp")]]
List rcpp_mce(const std::vector<double>& brts,       
              const std::vector<double>& init_pars,      
              int sample_size,
              int maxN,
              const std::string& plugin,             
              int soc,
              int max_missing,               
              double max_lambda,             
              const std::vector<double>& lower_bound,  
              const std::vector<double>& upper_bound,  
              double xtol_rel,                     
              int num_threads)
{
  auto model = emphasis::create_plugin_model(plugin);
  auto E = emphasis::E_step(sample_size,
                            maxN,
                            init_pars,
                            brts,
                            model.get(),
                            soc,
                            max_missing,
                            max_lambda,
                            num_threads);
  List ret;
  List trees;
  for (const emphasis::tree_t& tree : E.trees) {
    trees.push_back(unpack(tree));
  }
  ret["trees"] = trees;
  ret["rejected"] = E.rejected;
  ret["rejected_overruns"] = E.rejected_overruns;
  ret["rejected_lambda"] = E.rejected_lambda;
  ret["rejected_zero_weights"] = E.rejected_zero_weights;
  ret["time"] = E.elapsed;
  ret["weights"] = E.weights;
  ret["fhat"] = E.fhat;
  return ret;
}
