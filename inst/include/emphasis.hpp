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

// public C++ API
// Hanno 2020

#ifndef EMPHASIS_EMPHASIS_HPP_INCLUDED
#define EMPHASIS_EMPHASIS_HPP_INCLUDED

#include <stdexcept>
#include <limits>
#include <memory>
#include <vector>
#include <functional>
#include "plugin.hpp"


namespace emphasis {

  static constexpr int default_max_missing_branches = 10000;
  static constexpr double default_max_aug_lambda = 500.0;


  class emphasis_error : public std::runtime_error
  {
  public:
    explicit emphasis_error(const std::string& what) : std::runtime_error(what) {}
    explicit emphasis_error(const char* what) : std::runtime_error(what) {}
  };


  using brts_t = std::vector<double>;     // input tree


  // results from e
  struct E_step_t
  {
    E_step_t(E_step_t&&) = default;
    E_step_t& operator=(E_step_t&&) = default;
    E_step_t(const E_step_t&) = delete;
    E_step_t& operator=(const E_step_t&) = delete;
    E_step_t() {};
    ~E_step_t() {};

    std::vector<tree_t> trees;          // augmented trees
    std::vector<double> weights;
    std::vector<double> logg_;
    std::vector<double> logf_;
    
    double fhat;                        // mean, unscaled, weight
    int rejected_overruns = 0;          // # trees rejected because overrun of missing branches
    int rejected_lambda = 0;            // # trees rejected because of lambda overrun
    int rejected_zero_weights = 0;      // # trees rejected because of zero-weight
    int rejected = 0;                   // # trees rejected because of unhandled exception
    double elapsed = 0;                 // elapsed runtime [ms]
  };


  E_step_t E_step(int N,      // sample size
                  int maxN,   // max number of augmented trees (incl. invalid)
                  const param_t& pars,
                  const brts_t& brts,
                  class Model* model,
                  int soc = 2,
                  int max_missing = default_max_missing_branches,
                  double max_lambda = default_max_aug_lambda,
                  int num_threads = 0);


  // results from m
  struct M_step_t
  {
    param_t estimates;
    int opt = -1;                       // nlopt result
    double minf = 0.0;
    double elapsed = 0.0;               // elapsed runtime [ms]
  };


  using conditional_fun_t = std::function<double(const param_t&)>;
  
  
  M_step_t M_step(const param_t& pars,
                  const std::vector<tree_t>& trees,          // augmented trees
                  const std::vector<double>& weights,
                  class Model* model,
                  const param_t& lower_bound = {}, // overrides model.lower_bound
                  const param_t& upper_bound = {}, // overrides model.upper.bound
                  double xtol_rel = 0.001,
                  int num_threads = 0,
                  conditional_fun_t* conditional = nullptr);


  // results from mcem
  struct mcem_t
  {
    mcem_t(mcem_t&&) = default;
    mcem_t& operator=(mcem_t&&) = default;
    mcem_t(const mcem_t&) = delete;
    mcem_t& operator=(const mcem_t&) = delete;
    mcem_t() {};
    ~mcem_t() {};

    E_step_t e;
    M_step_t m;
  };

  
  mcem_t mcem(int N,      // sample size
              int maxN,   // max. number of augmented trees (incl. invalid)
              const param_t& pars,
              const brts_t& brts,
              class Model* model,
              int soc = 2,
              int max_missing = default_max_missing_branches,
              double max_lambda = default_max_aug_lambda,
              const param_t& lower_bound = {}, // overrides model.lower_bound
              const param_t& upper_bound = {}, // overrides model.upper.bound
              double xtol_rel = 0.001,
              int num_threads = 0,
              conditional_fun_t* conditional = nullptr);


  std::unique_ptr<class Model> create_plugin_model(const std::string& model_dll);

}

#endif
