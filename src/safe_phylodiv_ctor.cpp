#include <Rcpp.h>
#include "div_tree.hpp"
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List safe_phylodiv_ctor_cpp(NumericVector pars, float max_t, size_t max_N) {
  // Only construct the object, do not call simulate_tree
  sim_tree::phylodiv simulation(max_t, {pars[0], pars[1], pars[2], pars[3], pars[4], pars[5]}, max_N);
  NumericMatrix out(1, 4);
  out(0, 0) = 0;
  out(0, 1) = 0;
  out(0, 2) = 1;
  out(0, 3) = 0;
  return List::create(Named("tree") = out, Named("status") = "constructed");
}
