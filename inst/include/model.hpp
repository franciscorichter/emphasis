/* Copyright (c) 2007-2014 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

//
// C++-API for plug-in diversification model
// Hanno Hildenbrandt 2020
//

#ifndef EMPHASIS_PLUGIN_HPP_INCLUDED
#define EMPHASIS_PLUGIN_HPP_INCLUDED

#include <memory>
#include <vector>
#include <functional>
#include "model_helpers.hpp"

using namespace emphasis::detail;



namespace emphasis {
  using param_t = std::vector<double>;                  // unspecific parameters
  using tree_t = std::vector<node_t>;                   // tree, sorted by note_t::brts

  namespace {  
    using reng_t = std::mt19937_64;   // we need doubles
    static thread_local reng_t reng_ = make_random_engine<reng_t>();
  }  


// diversification model
class Model
{
public:
  Model() = default;
  ~Model() = default;
  
  const char* description() const { return "rpd5c model"; }   // textual description of the model
  bool is_threadsafe() const { return true; }                  // is this implementation thread-safe?
  bool numerical_max_lambda() const { return false; }
  int nparams() const {return 4;}                               // number of parameters
  
  // diversification model
  double extinction_time(double t_speciation, const param_t& pars, const tree_t& tree) const {
    return t_speciation + emphasis::detail::trunc_exp(tree[tree.size() - 1].brts - t_speciation, pars[0], reng_);
  }
  
  double speciation_rate(const param_t& pars, const node_t& node) const {
    const double lambda = pars[1] + pars[2] * node.n + pars[3] * node.pd / node.n;
    return std::max(0.0, lambda);
  }
  

  double nh_rate(double t, const param_t& pars, const tree_t& tree) const {
    auto it = lower_bound_node(t, tree.size(), &tree[0]);
    const double pd = calculate_pd(t, tree.size(), &tree[0]);
    const double lambda = std::max(0.0, pars[1] + pars[2] * it->n + pars[3] * pd / it->n);
    return lambda * it->n * (1.0 - std::exp(-pars[0] * (tree[tree.size() - 1].brts - t)));
  }
  
  
  double sampling_prob(const param_t& pars, const tree_t& tree) const {
    mu_integral muint(pars[0], tree[tree.size() - 1].brts);
    double inte = 0;
    double logg = 0;
    double prev_brts = 0;
    double tips = tree[0].n;
    double Ne = 0.0;
    for (unsigned i = 0; i < tree.size(); ++i) {
      const auto& node = tree[i];
      const double lambda = speciation_rate(pars, node);
      inte += node.n * lambda * muint(prev_brts, node.brts);
      tips += is_tip(node);
      Ne -= is_extinction(node);
      if (is_missing(node)) {
        const double lifespan = node.t_ext - node.brts;
        logg += std::log(node.n * pars[0] * lambda) - pars[0] * lifespan - std::log(2.0 * tips + Ne++);
      }
      prev_brts = node.brts;
    }
    return logg - inte;
  }
  
  double loglik(const param_t& pars, const tree_t& tree) const {
    log_sum log_lambda{};
    int cex = 0;
    double inte = 0.0;
    double prev_brts = 0.0;
    for (unsigned i = 0; i < tree.size(); ++i) {
      const auto& node = tree[i];
      const double lambda = speciation_rate(pars, node);
      if (is_extinction(node)) {
        ++cex;
      }
      else if (i != tree.size() - 1) {
        log_lambda += lambda;
      }
      inte += (node.brts - prev_brts) * node.n * (lambda + pars[0]);
      prev_brts = node.brts;
    }
    const double loglik = std::log(pars[0]) * cex + log_lambda.result() - inte;
    return loglik;
  }
  
  // optional hints for optimization step
  param_t lower_bound() const { return param_t(); }
  param_t upper_bound() const { return param_t(); }
};

}

#endif
