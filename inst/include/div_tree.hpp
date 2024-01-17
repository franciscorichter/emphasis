//
//  div_tree.hpp
//  div_tree
//

#ifndef div_tree_h
#define div_tree_h

#include <array>
#include <vector>
#include <random>
#include <thread>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace sim_tree {

struct rnd_t {
  std::mt19937 rndgen_;
  
  rnd_t() : rndgen_(std::random_device()()) {}
  
  float expon(float lambda) {
    if (lambda == 0.f) return 1e20f;
    return std::exponential_distribution<float>(lambda)(rndgen_);
  }
  
  bool bernouilli(double p) {
    std::bernoulli_distribution d(p);
    return d(rndgen_);
  }
  
  size_t random_number(size_t n) {
    if (n <= 1) return 0;
    return std::uniform_int_distribution<size_t>(0, n - 1)(rndgen_);
  }
};

struct branch {
  
  branch(float bd, int pl, int lab, float ext) :
  start_date(bd),
  parent_label(pl),
  label(lab),
  end_date(ext)
  {}
  
  float start_date;
  int parent_label;
  int label;
  float end_date;
  std::vector< int > daughters;
  
  void remove_daughter(int to_remove) {
    // code below can be made explicit to speed up.
    
    if (daughters.empty()) { // exception, should not happen
      return;
    }
    
    if (daughters.size() == 1) {
      daughters.clear();
      return;
    }
    
    if (daughters.size() == 2) {
      if (daughters[0] == to_remove) {
        daughters[0] = daughters[1];
      }
      daughters.pop_back();
      return;
    }
    
    if (daughters.size() > 2) {
      // this should not happen normally.
      for (size_t i = 0; i < daughters.size(); ++i) {
        if (daughters[i] == to_remove) {
          daughters[i] = daughters.back();
          daughters.pop_back();
          return;
        }
      }
    }
    return;
  }
  
};

enum breaks { none, finished, extinction, maxN_exceeded };

struct phylodiv {
  float max_t, P, t;
  size_t max_N, N;
  std::array<double, 6> pars; // A_0, A_n, A_p, B_0, B_n, B_p
  
  std::vector<branch> tree, ltable;
  breaks break_type;
  rnd_t rndgen;
  
  phylodiv(float total_time, const std::array<double, 6>& p, size_t maxN)
    : max_t(total_time), max_N(maxN), pars(p), P(0.f), t(0.f), rndgen(), break_type(none) {}
  
  float calculate_full_phylodiv() {
    float ph = 0.f;
    for (const auto& i : tree) {
      float bl = (i.end_date == -1) ? (t - i.start_date) : (i.end_date - i.start_date);
      ph += bl;
    }
    return ph;
  }
  
  float calculateRate(const std::array<double, 3>& rate_params, size_t N, float P) {
    return exp(rate_params[0] + rate_params[1] * N + rate_params[2] * P);
  }
  
  bool simulate_tree() {
    size_t N1 = 1, N2 = 1;
    N = 2;
    tree.clear();
    tree.emplace_back(0.f, 0, -1, -1);
    tree.emplace_back(0.f, -1, 2, -1);
    int tree_id = 3;
    
    while (true) {
      N = N1 + N2;
      if (N >= max_N) {
        break_type = maxN_exceeded;
        break;
      }
      
      float spec_rate = calculateRate({pars[0], pars[1], pars[2]}, N, P);
      float ext_rate = calculateRate({pars[3], pars[4], pars[5]}, N, P);
      float total_rate = spec_rate + ext_rate;
      float next_event_time = t + rndgen.expon(total_rate);
      
      P += (t - next_event_time) * N;
      if (next_event_time < max_t) {
        if (rndgen.bernouilli(spec_rate / total_rate)) {
          size_t parent = sample_tip(N);
          int new_id_1 = tree_id++;
          int new_id_2 = tree_id++;
          
          if (tree[parent].label < 0) {
            new_id_1 *= -1;
            new_id_2 *= -1;
            N2 += 2;
          } else {
            N1 += 2;
          }
          
          tree[parent].end_date = next_event_time;
          ltable[parent].end_date = next_event_time;
          tree.emplace_back(next_event_time, tree[parent].label, new_id_1, -1);
          tree.emplace_back(next_event_time, tree[parent].label, new_id_2, -1);
        } else {
          size_t to_remove = sample_tip(N);
          tree[to_remove].end_date = next_event_time;
          if (tree[to_remove].label < 0) N2--;
          else N1--;
          
          P -= purge_tree_record(to_remove);
        }
      }
      
      t = next_event_time;
      if (N1 < 1 || N2 < 1) {
        break_type = extinction;
        break;
      }
      if (t >= max_t) {
        break_type = finished;
        break;
      }
    }
    
    N = N1 + N2;
    P = calculate_full_phylodiv(); // Final check
    return break_type == extinction;
  }

    double purge_tree_record(size_t focal_index) {
      double bt_removed = 0.f;
      auto label_removed = tree[focal_index].label;
      auto parent = tree[focal_index].parent_label;
      // first, we remove the focal tree:
      bt_removed += tree[focal_index].end_date - tree[focal_index].start_date;
      tree[focal_index] = tree.back();
      tree.pop_back();
      // now, we have to remove it from the parent:
      auto parent_loc = std::find_if(tree.begin(), tree.end(),
                                     [parent](const branch& other){return other.label == parent;});
  
      if (parent_loc != tree.end()) { // root has no parents.
        parent_loc->remove_daughter(label_removed);
        if (parent_loc->daughters.empty()) {
          auto parent_index = std::distance(tree.begin(), parent_loc);
          bt_removed += purge_tree_record(parent_index);
        }
      }
      return bt_removed;
    }
  
    size_t sample_tip(size_t N) {
      size_t index = 0; 
      if (N * 10 < tree.size()) {
        // there are many extinct tips!
        std::vector< size_t > alive;
        for (size_t i = 0; i < tree.size(); ++i) {
          if (tree[i].end_date == -1) alive.push_back(i);
        }
        index = alive[ rndgen.random_number(alive.size())];
      } else {
        index = rndgen.random_number(tree.size());
        while(tree[index].end_date != -1) {
          index = rndgen.random_number(tree.size());
        }
      }
      
      return index;
    }
    
    
  };
}


#endif /* div_tree_h */
