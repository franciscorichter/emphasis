#include <Rcpp.h>
#include "emphasis.hpp"
#include "plugin.hpp"
#include "rinit.h"
using namespace Rcpp;

namespace emphasis {

  int sample_ext_spec(const tree_t& in_tree, int clade) {
    
    auto is_extinct[](const node_t& nt) {(nt.to == 1 && 
                                          nt.t_ext == && 
                                          nt.clade == c) ? return true : return false;}
    
    auto n = std::count_if(in_tree.begin(),
                           in_tree.end(),
                           is_extinct );

    int selected = std::uniform_int_distribution<> (0, static_cast<int>(n - 1))(reng);

    int index = 0;
    for (auto i : in_tree) {
      if (is_extinct(i)) {
        selected--;
        if (selected == 0) {
          return index;
        }
      }
      index++;
    }
    return 0; // * this should NEVER be reached* 
    // perhaps add a throw?
  }

  
// [[Rcpp::export(name = "simtree_pd_cpp")]]
tree_t simtree_pd(const std::vector<double>& pars,
                  double ct) {
  
  double cbt = 0.0;
  int N1 = 1;
  int N2 = 1;
  double mu = std::max(pars[0], 0.0); 
  double P = 0.0;
  
  tree_t tree = emphasis::detail::create_tree(std::vector<double>(2, 0.0), 2);
  
  while (cbt < ct &&
         N1 >= 1  &&
         N2 >= 1) {
    double next_bt = emphasis::get_next_bt(tree, cbt);
    
    int N = N1 + N2;
    auto lambda_max = std::max(0.0,
                                 pars[1] + pars[2] * N +
                                 ((P + N * (next_bt - cbt) - cbt) / N) * pars[3]);

    double rate_max = (lambda_max + mu) * N 
    double u1 = std::uniform_real_distribution<>()(reng);
    double next_event_time = cbt - std::log(u1) / lambda_max;
    P += N * (next_event_time - cbt)  ;
    if (next_event_time < next_bt) {
      double u2 = std::uniform_real_distribution<>()(reng);
      double lambda_ct = std::max(0, 
                          pars[1] + pars[2] * N  +  
                          ((P + N * (next_event_time - cbt) - next_event_time) / N) * pars[3]);
      double pt = ((lambda_ct + mu) * N ) / rate_max;
      
      if (u2 < pt) {
        
        int to = std::bernoulli_distribution(lambda_ct / (lambda_ct + mu))(reng); // c(1, 0) with prob c(lambda_ct, mu) / (lambda_ct + mu)
        int clade = std::bernoulli_distribution(N2 / (N1 + N2))(reng); // c("a","b"), e.g. c(0, 1) with prob (N1, N2) / (N1 + N2)
  
        if (to == 1) {
          clade == 0 ? N1++ : N2++;
        } else {
          clade == 0 ? N1-- : N2--;
         
          if ((N1 >= 1) &  
              (N2 >= 1)) { // tree is not extinct, but extinct species added
            
            auto ext_spec = sample_ext_spec(tree);
            
            tree[ext_spec].t_ext = next_event_time;
            auto removed_branch_length = tree[ext_spec].t_ext - tree[ext_spec].brts;
            P -= removed_branch_length
          }
        }
        /* 
         * I don't really understand how to translate this to C++:
         * tree = rbind(tree,
                       data.frame(brts = next_event_time,
                       to = to,
                       t_ext = Inf,
                       clade = clade))
         
         // insert_species does not have a "to" ?
         // also, why insert ALWAYS, and not only if speciation occured?
         insert_species(next_event_time, to, )
         */
      }
    }
    cbt = next_event_time
  }
 
 if ((N1 == 0) |  (N2 == 0)) {
   tree = cbt
 } else {
   // again:
   /*
   tree = rbind(tree,data.frame(brts = ct,
                                to = 1,
                                t_ext = Inf, 
                                clade = 0))
    */
 }
 
  