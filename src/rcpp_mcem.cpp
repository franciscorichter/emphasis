// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>
#include "emphasis.hpp"
#include "model.hpp"
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

//' function to perform one step of the E-M algorithm
//' @param brts vector of branching times
//' @param init_pars vector of initial parameter files
//' @param sample_size number of samples
//' @param maxN maximum number of failed trees
//' @param soc number of lineages at the root/crown (1/2)
//' @param max_missing maximum number of species missing
//' @param max_lambda maximum speciation rate
//' @param lower_bound vector of lower bound values for optimization, should
//' be equal in length to the vector of init_pars
//' @param upper_bound vector of upper bound values for optimization, should
//' be equal in length to the vector of init_pars
//' @param xtol_rel relative tolerance for optimization
//' @param num_threads number of threads used.
//' @param copy_trees if set to true, the trees generated are returned as well
//' @param rconditional R function that evaluates the GAM function.
//' @return a list with the following components: 
//' \itemize{
//'  \item{trees}{list of trees}
//'  \item{rejected}{number of rejected trees}
//'  \item{rejected_overruns}{number of trees rejected due to too large size}
//'  \item{rejected_lambda}{number of trees rejected due to lambda errors}
//'  \item{rejected_zero_weights}{number of trees rejected to zero weight}
//'  \item{estimates}{vector of estimates}
//'  \item{nlopt}{nlopt status}
//'  \item{fhat}{vector of fhat values}
//'  \item{time}{time elapsed}
//'  \item{weights}{vector of weights}
//' }
//' @export
// [[Rcpp::export(name = "em_cpp")]]
List rcpp_mcem(const std::vector<double>& brts,       
               const std::vector<double>& init_pars,      
               int sample_size,
               int maxN,
               int soc,
               int max_missing,               
               double max_lambda,             
               const std::vector<double>& lower_bound,  
               const std::vector<double>& upper_bound,  
               double xtol_rel,                     
               int num_threads,
               bool copy_trees,
               Nullable<Function> rconditional = R_NilValue) 
{
  auto model = emphasis::Model(lower_bound, upper_bound);
  
  emphasis::conditional_fun_t conditional{};
  if (rconditional.isNotNull()) {
    conditional = [cond= Function(rconditional)](const emphasis::param_t& pars) {
      return as<double>( cond(NumericVector(pars.cbegin(), pars.cend())) );
    };
  }
  auto mcem = emphasis::mcem(sample_size,
                             maxN,
                             init_pars,
                             brts,
                             model,
                             soc,
                             max_missing,
                             max_lambda,
                             lower_bound,
                             upper_bound,
                             xtol_rel,
                             num_threads,
                             conditional ? &conditional : nullptr);
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
  ret["rejected"] = mcem.e.rejected;
  ret["rejected_overruns"] = mcem.e.rejected_overruns;
  ret["rejected_lambda"] = mcem.e.rejected_lambda;
  ret["rejected_zero_weights"] = mcem.e.rejected_zero_weights;
  ret["estimates"] = NumericVector(mcem.m.estimates.begin(), mcem.m.estimates.end());
  ret["nlopt"] = mcem.m.opt;
  ret["fhat"]  = mcem.e.fhat;
  ret["time"]  = mcem.e.elapsed + mcem.m.elapsed;
  ret["weights"] = mcem.e.weights;
  return ret;
}