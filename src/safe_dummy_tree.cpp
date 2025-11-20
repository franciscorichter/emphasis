#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List safe_dummy_tree_cpp(NumericVector pars, double max_t, int max_N, int max_tries) {
  NumericMatrix tree(2, 3);
  tree(0, 0) = 1; tree(0, 1) = 0; tree(0, 2) = 1;
  tree(1, 0) = 2; tree(1, 1) = 0.5; tree(1, 2) = 2;
  return List::create(Named("tree") = tree, Named("status") = "dummy");
}
