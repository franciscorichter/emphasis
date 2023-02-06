//
//  phylodiv_tree.h
//  phylodiv_tree
//
//  Created by Thijs Janzen on 14/09/2021.
//  Copyright © 2021 Thijs Janzen. All rights reserved.
//

#ifndef phylodiv_tree_h
#define phylodiv_tree_h

#include <array>
#include <vector>
#include <random>
#include <thread>
#include <algorithm>

namespace sim_tree {

  
  struct rnd_t {
    std::mt19937 rndgen_;
  
    rnd_t() {
      std::mt19937 rndgen_t(get_seed());
      rndgen_ = rndgen_t;
    }
  
    int get_seed() {
      const auto tt = static_cast<int64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
      auto tid = std::this_thread::get_id();
      const uint64_t e3{ std::hash<std::remove_const_t<decltype(tid)>>()(tid) };
      auto output = static_cast<int>(tt + e3);
      if (output < 0) output *= -1;
      return output;
    }
  
    size_t random_number(size_t n)    {
      if(n <= 1) return 0;
      return static_cast<size_t>(std::uniform_int_distribution<> (0, static_cast<int>(n - 1))(rndgen_));
    }
  
    bool bernouilli(double p) {
      std::bernoulli_distribution d(p);
      return(d(rndgen_));
    }
  
    float expon(float lambda) {
      if (lambda == 0.f) return 1e20f;
      return std::exponential_distribution<float>(lambda)(rndgen_);
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
        for (int i = 0; i < daughters.size(); ++i) {
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
  
  enum breaks {none, finished, extinction, maxN_exceeded};
  
  struct phylodiv {
    float max_t;
    float P; // phylogenetic diversity
    float t; // time
    float max_N; // number of lineages after which we assume non-extinction
    size_t N; // number of extant species
    std::array< double, 4> pars; // mu, lambda, B_n, B_p
  
    std::vector< branch > tree;
  
    breaks break_type;
  
    rnd_t rndgen;
  
    phylodiv(float total_time,
             const std::array<double, 4>& p,
             int maxN) :
    max_t(total_time),
    max_N(maxN),
    pars(p) {
      P = 0.f;
      t = 0.f;
      rndgen = rnd_t();
    }
  
    float calculate_full_phylodiv(float t) {
      float ph = 0.f;
      for (const auto& i : tree) {
        if (i.end_date == -1) {
          float bl = t - i.start_date;
          ph += bl;
        } else {
          float bl = i.end_date - i.start_date;
          ph += bl;
        }
      }
      return ph;
    }
  
    bool simulate_tree() {
     size_t N1 = 1;
     size_t N2 = 1;
     N = 2;
     
     float prev_t = 0.f;
     t = 0.f;
  
     tree.clear();
     tree.push_back( branch(t, 0, -1, -1));
     tree.push_back( branch(t, -1, 2, -1));
  
     int tree_id = 3;
     P = 0.f;
     
     float mu = pars[0];

     break_type = breaks::none;
     
     while(true) {
  
       N = N1 + N2;
       if (N >= max_N) {
         break_type = breaks::maxN_exceeded;
         break;
       }

       
       float spec_rate = pars[1] + pars[2] * N  +
                         ((P + N * (max_t - t) - t) / N) * pars[3];
       if (spec_rate < 0.f) spec_rate = 0.f;
  
       float total_rate = (spec_rate + mu ) * N;
  
       float next_event_time = t + rndgen.expon(total_rate);

       P += (t - prev_t) * N;
       
       //P = calculate_full_phylodiv(t);
       
       if (next_event_time < max_t) {
         float focal_spec = pars[1] +
                            pars[2] * N  +
                            ((P + N * (next_event_time - t) - t) / N) * pars[3];
         float pt = ((focal_spec + mu) * N ) / total_rate;
         
         if (rndgen.bernouilli(pt)) {
           // event is accepted
           if (rndgen.bernouilli(focal_spec / (focal_spec + mu))) {
             // speciation
             size_t parent = sample_tip(N);
             int new_id_1 = tree_id;
             tree_id++;
             int new_id_2 = tree_id;
             tree_id++;
  
             if (tree[parent].label < 0) {
               new_id_1 *= -1;
               new_id_2 *= -1;
               N2++;
             } else {
               N1++;
             }
  
             tree[parent].end_date = next_event_time;
             tree[parent].daughters.push_back(new_id_1);
             tree[parent].daughters.push_back(new_id_2);
             tree.push_back( branch(next_event_time, tree[parent].label, new_id_1, -1));
             tree.push_back( branch(next_event_time, tree[parent].label, new_id_2, -1));
           } else {
             // extinction
             size_t to_remove = sample_tip(N);
             tree[to_remove].end_date = next_event_time;
             if (tree[to_remove].label < 0) {
               N2--;
             } else {
               N1--;
             }
  
             P -= purge_tree_record(to_remove);
           }
         }
         prev_t = t;
         t = next_event_time;
       } else {
         t = max_t;
       }
     
       if (N1 < 1 || N2 < 1) {
         break_type = breaks::extinction;
         break;
       }
       if (t >= max_t) {
         break_type = breaks::finished;
         break;
       }
     }
     
     N = N1 + N2;
     P = calculate_full_phylodiv(t); // final check.
     
     if (break_type == breaks::extinction) {
       return true;
     }
     
     return false;
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


#endif /* phylodiv_tree_h */
