#include "phylodiv_tree.h"

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix simulate_pd_trees_cpp(Rcpp::NumericVector pars,
                                 float max_t,
                                 size_t repl) {
  
  phylodiv sim_tree(max_t,  {pars[0], pars[1], pars[2], pars[3]});
  
  Rcpp::NumericMatrix results(repl, 4);
  for (size_t r = 0; r < repl; ++r) {
    bool is_extant = sim_tree.simulate_tree();
    results(repl, 0) = is_extant;
    results(repl, 1) = sim_tree.t;
    results(repl, 2) = sim_tree.N;
    results(repl, 3) = sim_tree.P;
  }
  return results;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix explore_grid_cpp(Rcpp::NumericVector par1,
                                     Rcpp::NumericVector par2,
                                     Rcpp::NumericVector par3,
                                     Rcpp::NumericVector par4,
                                     float max_t,
                                     int num_repl) {
  
  size_t total_number_simulations = par1.size() * par2.size() * 
                                    par3.size() * par4.size() * num_repl;
  
  Rcpp::NumericMatrix output(total_number_simulations, 8);
  int row = 0;
  for (auto a : par1) {
    for (auto b : par2) {
      for (auto c : par3) {
        for (auto d : par4) {
          
          phylodiv sim_tree(max_t,  {static_cast<double>(a), 
                                     static_cast<double>(b), 
                                     static_cast<double>(c), 
                                     static_cast<double>(d)});
          
          for (size_t r = 0; r < num_repl; ++r) {
            bool is_extant = sim_tree.simulate_tree();
            
            output(row, 0) = a;
            output(row, 1) = b;
            output(row, 2) = c;
            output(row, 3) = d;
            
            output(row, 4) = is_extant;
            output(row, 5) = sim_tree.t;
            output(row, 6) = sim_tree.N;
            output(row, 7) = sim_tree.P;
            row++;
          }
        }
      }
    }
  }
  return output;
}