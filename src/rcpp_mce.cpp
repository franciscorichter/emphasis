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
                                  int num_threads) 
{                     
  std::vector<std::array<double, 6>> results(pars_R.nrow(), {0});
  
  // copy parameters to C++ object, for multithreading
  std::vector<std::vector<double>> pars(pars_R.nrow());
  std::vector<double> row_entry(pars_R.ncol());
  for (int i = 0; i < pars_R.nrow(); ++i) {
    for(int j = 0; j < pars_R.ncol(); ++j) {
      row_entry[j] = pars_R(i, j);
    }
    pars[i] = row_entry;
  }
  
  const int grainsize = pars.size() / std::max<unsigned>(1, num_threads);
  
  tbb::task_arena arena(num_threads);
  
  tbb::parallel_for(tbb::blocked_range<unsigned>(0, pars.size(), grainsize), [&](const tbb::blocked_range<unsigned>& r) {
    for (unsigned i = r.begin(); i < r.end(); ++i) {
      auto model = emphasis::Model(lower_bound, upper_bound);
      auto EI = emphasis::E_step_info(sample_size,
                                      maxN,
                                      pars[i],
                                      brts,
                                      model,
                                      soc,
                                      max_missing,
                                      max_lambda);
      results[i] = { EI.fhat, 
                     double(EI.rejected_lambda),
                     double(EI.rejected_overruns), 
                     double(EI.rejected_zero_weights),
                     double(EI.logf),
                     double(EI.logg)};
    }
  });
  
  // and now back to Rcpp::NumericMatrix...
  Rcpp::NumericMatrix out(results.size(), results[0].size());
  for (size_t i = 0; i < results.size(); ++i) {
    for (size_t j = 0; j < results[i].size(); ++j) {
      if (std::isnan(results[i][j])) {
        out(i, j) = NA_REAL;
      } else {
        out(i, j) = results[i][j];
      }
    }
  }
  return out;
}

using mat2d = std::vector< std::vector<double>>;

template<typename T>
Rcpp::NumericMatrix two_d_vec_to_mat(const std::vector<T> from_cpp) {
  // and now back to Rcpp::NumericMatrix...
  Rcpp::NumericMatrix out(from_cpp.size(), from_cpp[0].size());
  for (size_t i = 0; i < from_cpp.size(); ++i) {
    for (size_t j = 0; j < from_cpp[i].size(); ++j) {
      if (std::isnan(from_cpp[i][j])) {
        out(i, j) = NA_REAL;
      } else {
        out(i, j) = from_cpp[i][j];
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List rcpp_mce_grid_factorial(const Rcpp::NumericMatrix pars_R,
                                  const std::vector<double>& brts,       
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
  std::vector<std::array<double, 6>> results(pars_R.nrow(), {0});
  
  // copy parameters to C++ object, for multithreading
  std::vector<std::vector<double>> pars(pars_R.nrow());
  std::vector<double> row_entry(pars_R.ncol());
  for (int i = 0; i < pars_R.nrow(); ++i) {
    for(int j = 0; j < pars_R.ncol(); ++j) {
      row_entry[j] = pars_R(i, j);
    }
    pars[i] = row_entry;
  }
  
  mat2d logf_grid = std::vector< std::vector< double >>(pars.size(), std::vector<double>(pars.size(), 0.0));
  mat2d logg_grid = logf_grid;
  
  const int grainsize = pars.size() / std::max<unsigned>(1, num_threads);
  
  tbb::task_arena arena(num_threads);
  
  tbb::parallel_for(tbb::blocked_range<unsigned>(0, pars.size(), grainsize), [&](const tbb::blocked_range<unsigned>& r) {
    for (unsigned i = r.begin(); i < r.end(); ++i) {
      auto model = emphasis::Model(lower_bound, upper_bound);
      auto EI = emphasis::E_step_info_grid(maxN,
                                           pars[i],
                                           pars,
                                           brts,
                                           model,
                                           soc,
                                           max_missing,
                                           max_lambda);
      results[i] = { EI.fhat, 
                     double(EI.rejected_lambda),
                     double(EI.rejected_overruns), 
                     double(EI.rejected_zero_weights),
                     double(EI.logf),
                     double(EI.logg)};
      logf_grid[i] = EI.logf_grid;
      logg_grid[i] = EI.logg_grid;
    }
  });
  
  
  Rcpp::List out;
  out["results"]   = two_d_vec_to_mat(results);
  out["logf_grid"] = two_d_vec_to_mat(logf_grid);
  out["logg_grid"] = two_d_vec_to_mat(logg_grid);
  return out;
}
