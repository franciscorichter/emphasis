// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>
#include "emphasis.hpp"
#include "plugin.hpp"
#include "rinit.h"
using namespace Rcpp;


namespace {

  DataFrame unpack(const emphasis::tree_t& tree)
  {
    NumericVector brts, n, t_ext, pd;
    for (const emphasis::node_t& node : tree) {
      brts.push_back(node.brts);
      n.push_back(node.n);
      t_ext.push_back(node.t_ext);
      pd.push_back(node.pd);
    }
    return DataFrame::create(Named("brts") = brts, 
                             Named("n") = n, 
                             Named("t_ext") = t_ext,
                             Named("pd") = pd);
  }

}

//' function to perform one step of the E-M algorithm
//' @param brts vector of branching times
//' @param init_pars vector of initial parameter files
//' @param sample_size number of samples
//' @param maxN maximum number of failed trees
//' @param plugin string indicating plugin used, currently available: 'rpd1' and
//' 'rpd5c'
//' @param soc number of lineages at the root/crown (1/2)
//' @param max_missing maximum number of species missing
//' @param max_lambda maximum speciation rate
//' @param lower_bound vector of lower bound values for optimization, should
//' be equal in length to the vector of init_pars
//' @param upper_bound vector of upper bound values for optimization, should
//' be equal in length to the vector of init_pars
//' @param xtol_rel relative tolerance for optimization
//' @param num_threads number of threads used.
//' @return a list with the following components: 
//' \itemize{
//'  \item{trees}{list of trees}
//'  \item{rejected}{number of rejected trees}
//'  \item{rejected_overruns}{number of trees rejected due to too large size}
//'  \item{rejected_lambda}{number of trees rejected due to lambda errors}
//'  \item{rejected_zero_weights}{number of trees rejected to zero weight}
//'  \item{time_elapsed}{time used}
//'  \item{weights}{vector of weights}
//'  \item{fhat}{vector of fhat values}
//'  \item{logf}{vector of logf values}
//'  \item{logg}{vector of logg values}
//' }
//' @export
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
  try {
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
    ret["logf"] = E.logf_;
    ret["logg"] = E.logg_;
    return ret;
  } catch (emphasis::emphasis_error e) {
    std::cerr << e.what() << "\n";
    return NA_REAL;
  } catch (...) {
    std::cerr << "unknown exception\n";
    return NA_REAL;
  }
}
