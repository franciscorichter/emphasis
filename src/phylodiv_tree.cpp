#include "phylodiv_tree.hpp"

#include <Rcpp.h>

Rcpp::List simulate_single_pd_tree_cpp(Rcpp::NumericVector pars,
                            float max_t,
                            float max_N,
                            int max_tries) {
  
  sim_tree::phylodiv simulation(max_t,  {pars[0], pars[1], pars[2], pars[3]}, max_N);
  
  simulation.simulate_tree_ltable();
  size_t tries = 0;
  while (simulation.break_type == sim_tree::breaks::extinction ||
         simulation.break_type == sim_tree::breaks::maxN_exceeded) {
    simulation.simulate_tree_ltable();
    tries++;
    if (tries > max_tries) break;
  }
  
  
  auto tree = simulation.ltable;
  Rcpp::NumericMatrix out(tree.size(), 4);
  size_t cnt = 0;
  
  auto crown_age = simulation.t;
  
  for (const auto& i : tree) {
    out(cnt, 0) = crown_age - i.start_date;
    out(cnt, 1) = i.parent_label;
    out(cnt, 2) = i.label;
    
    auto end_date = i.end_date;
    if (end_date >= 0.0) end_date = crown_age - end_date;
    
    out(cnt, 3) = end_date;
    cnt++;
  }
  
  std::string status = "done";
  if (simulation.break_type == sim_tree::breaks::extinction) status = "extinct";
  if (simulation.break_type == sim_tree::breaks::maxN_exceeded) status = "too_large";
  
  Rcpp::List res;
  res["Ltable"] = out;
  res["status"] = status;
  
  return res;
}


Rcpp::List simulate_single_ep_tree_cpp(Rcpp::NumericVector pars,
                                        float max_t,
                                        float max_N,
                                        int max_tries) {
  sim_tree::pendant_div simulation(max_t,
    {pars[0], pars[1], pars[2], pars[3], pars[4]}, (size_t)max_N);

  simulation.simulate_tree_ltable();
  size_t tries = 0;
  while (simulation.break_type == sim_tree::breaks::extinction ||
         simulation.break_type == sim_tree::breaks::maxN_exceeded) {
    simulation.simulate_tree_ltable();
    if (++tries > (size_t)max_tries) break;
  }

  auto& tree = simulation.ltable;
  Rcpp::NumericMatrix out(tree.size(), 4);
  float crown_age = simulation.t;
  for (size_t cnt = 0; cnt < tree.size(); ++cnt) {
    out(cnt, 0) = crown_age - tree[cnt].start_date;
    out(cnt, 1) = tree[cnt].parent_label;
    out(cnt, 2) = tree[cnt].label;
    auto end_date = tree[cnt].end_date;
    if (end_date >= 0.f) end_date = crown_age - end_date;
    out(cnt, 3) = end_date;
  }

  std::string status = "done";
  if (simulation.break_type == sim_tree::breaks::extinction) status = "extinct";
  if (simulation.break_type == sim_tree::breaks::maxN_exceeded) status = "too_large";

  return Rcpp::List::create(Rcpp::Named("Ltable") = out,
                             Rcpp::Named("status") = status);
}


Rcpp::NumericMatrix simulate_pd_trees_cpp(Rcpp::NumericVector pars,
                                 float max_t,
                                 size_t repl,
                                 float max_N) {
  
  sim_tree::phylodiv sim_tree(max_t,  {pars[0], pars[1], pars[2], pars[3]}, max_N);
  
  Rcpp::NumericMatrix results(repl, 5);
  for (size_t r = 0; r < repl; ++r) {
    bool is_extinct = sim_tree.simulate_tree();
    results(r, 0) = is_extinct;
    results(r, 1) = sim_tree.t;
    results(r, 2) = sim_tree.N;
    results(r, 3) = sim_tree.P;
    results(r, 4) = sim_tree.break_type;
  }
  return results;
}

Rcpp::NumericMatrix explore_grid_cpp(Rcpp::NumericVector par1,
                                     Rcpp::NumericVector par2,
                                     Rcpp::NumericVector par3,
                                     Rcpp::NumericVector par4,
                                     float max_t,
                                     int num_repl,
                                     int max_N) {
  
  size_t total_number_simulations = par1.size() * par2.size() * 
                                    par3.size() * par4.size() * num_repl;
  

  Rcpp::NumericMatrix output(total_number_simulations, 9);
  int row = 0;
  for (auto a : par1) {
    for (auto b : par2) {
      for (auto c : par3) {
        for (auto d : par4) {
          
          sim_tree::phylodiv sim_tree(max_t,  {static_cast<double>(a), 
                                     static_cast<double>(b), 
                                     static_cast<double>(c), 
                                     static_cast<double>(d)},
                            max_N);
          
          for (int r = 0; r < num_repl; ++r) {
            bool is_extinct = sim_tree.simulate_tree();
            
            output(row, 0) = a;
            output(row, 1) = b;
            output(row, 2) = c;
            output(row, 3) = d;
            
            output(row, 4) = is_extinct;
            output(row, 5) = sim_tree.t;
            output(row, 6) = sim_tree.N;
            output(row, 7) = sim_tree.P;
            output(row, 8) = sim_tree.break_type;
            row++;
          }
        }
      }
    }
  }
  return output;
}