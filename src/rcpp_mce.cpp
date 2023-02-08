// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>
#include <vector>
#include <tbb/tbb.h>
#include "emphasis.hpp"
#include "model.hpp"
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
              int soc,
              int max_missing,               
              double max_lambda,             
              const std::vector<double>& lower_bound,  
              const std::vector<double>& upper_bound,  
              double xtol_rel,                     
              int num_threads)
{
  auto model = emphasis::Model(lower_bound, upper_bound);

  auto E = emphasis::E_step(sample_size,
                            maxN,
                            init_pars,
                            brts,
                            model,
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
}


std::vector<double> rcpp_mce_nothrow(const std::vector<double>& brts,       
                                     const std::vector<double>& init_pars,      
                                     int sample_size,
                                     int maxN,
                                     int soc,
                                     int max_missing,               
                                     double max_lambda,             
                                     const std::vector<double>& lower_bound,  
                                     const std::vector<double>& upper_bound,  
                                     double xtol_rel,                     
                                     int num_threads)
{
  auto model = emphasis::Model(lower_bound, upper_bound);
  
  try {
  auto E = emphasis::E_step(sample_size,
                            maxN,
                            init_pars,
                            brts,
                            model,
                            soc,
                            max_missing,
                            max_lambda,
                            num_threads);
  std::vector<double> out = {E.fhat, 
                             double(E.rejected_lambda),
                             double(E.rejected_overruns), 
                             double(E.rejected_zero_weights)};
  return out;
  } catch (const emphasis::emphasis_error_E& E) {
    return {-1.0, //E.E_.fhat, 
            double(E.E_.rejected_lambda),
            double(E.E_.rejected_overruns), 
            double(E.E_.rejected_zero_weights)};
  }
}


// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_mce_grid(const Rcpp::NumericMatrix pars_R,
                                  const std::vector<double>& brts,       
                                  int sample_size,
                                  int maxN,
                                  int soc,
                                  int max_missing,               
                                  double max_lambda,             
                                  const std::vector<double>& lower_bound,  
                                  const std::vector<double>& upper_bound,  
                                  double xtol_rel,                     
                                  int num_threads) {
  
  std::vector< std::vector<double> > results(pars_R.nrow(), std::vector<double>(4, 0.0));
  
  // copy parameters to C++ object, for multithreading
  std::vector< std::vector<double> > pars(pars_R.nrow());
  std::vector<double> row_entry(pars_R.ncol());
  for (int i = 0; i < pars_R.nrow(); ++i) {
    for(int j = 0; j < pars_R.ncol(); ++j) {
      row_entry[j] = pars_R(i, j);
    }
    pars[i] = row_entry;
  }
  
  const int grainsize = pars.size() / std::max<unsigned>(1, std::min<unsigned>(std::thread::hardware_concurrency(), num_threads));
  tbb::parallel_for(tbb::blocked_range<unsigned>(0, pars.size(), grainsize), [&](const tbb::blocked_range<unsigned>& r) {
    for (unsigned i = r.begin(); i < r.end(); ++i) {
      std::vector<double> local_pars = pars[i];
      results[i] = rcpp_mce_nothrow(brts,       
                                     local_pars,      
                                     sample_size,
                                     maxN,
                                     soc,
                                     max_missing,               
                                     max_lambda,             
                                     lower_bound,  
                                     upper_bound,  
                                     xtol_rel,                     
                                     num_threads);
     }
    });
  
  // and now back to Rcpp::NumericMatrix...
  Rcpp::NumericMatrix out(pars.size(), pars[0].size());
  for (int i = 0; i < pars.size(); ++i) {
    for (int j = 0; j < pars[i].size(); ++j) {
      out(i, j) = pars[i][j];
    }
  }
  return out;
}



