#include <Rcpp.h>
#include "div_tree.hpp"
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List safe_simulate_div_tree_cpp(NumericVector pars, float max_t, size_t max_N, int max_tries) {
  sim_tree::phylodiv simulation(max_t, {pars[0], pars[1], pars[2], pars[3], pars[4], pars[5]}, max_N);

  bool is_extinct = simulation.simulate_tree();
  size_t tries = 0;
  while ((simulation.break_type == sim_tree::breaks::extinction ||
          simulation.break_type == sim_tree::breaks::maxN_exceeded) && tries < (size_t)max_tries) {
    is_extinct = simulation.simulate_tree();
    tries++;
  }

  auto tree = simulation.tree;
  size_t nrow = tree.size();
  Rcpp::NumericMatrix out(nrow ? nrow : 1, 4);

  auto crown_age = simulation.t;

  if (nrow) {
    size_t cnt = 0;
    for (const auto& i : tree) {
      out(cnt, 0) = i.parent_label;
      out(cnt, 1) = crown_age - i.start_date;
      out(cnt, 2) = i.label;
      auto end_date = (i.end_date == -1) ? crown_age : (crown_age - i.end_date);
      out(cnt, 3) = end_date;
      cnt++;
    }
  } else {
    // Dummy row if simulation failed
    out(0, 0) = 0;
    out(0, 1) = 0;
    out(0, 2) = 1;
    out(0, 3) = 0;
  }

  std::string status = "done";
  if (simulation.break_type == sim_tree::breaks::extinction) status = "extinct";
  if (simulation.break_type == sim_tree::breaks::maxN_exceeded) status = "too_large";

  Rcpp::List res;
  res["tree"] = out;
  res["status"] = status;
  return res;
}
