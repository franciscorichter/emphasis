#include <thread>
#include <mutex>
#include <chrono>
#include <numeric>
#include <cmath>  
#include <algorithm>
#include <functional>
#include <memory>
#include <tbb/tbb.h>
#include "plugin.hpp"
#include "emphasis.hpp"
#include "sbplx.hpp"


namespace emphasis {

  namespace {

    struct nlopt_f_data
    {
      nlopt_f_data(const Model* M, 
                   const std::vector<tree_t>& Trees, 
                   const std::vector<double>& W,
                   conditional_fun_t* Conditional)
        : model(M), trees(Trees), w(W), conditional(Conditional)
      {
      }

      ~nlopt_f_data()
      {
        auto empty = tree_t{};
      }

      const Model* model;
      const std::vector<tree_t>& trees;
      const std::vector<double>& w;
      conditional_fun_t* conditional;
    };


    double objective(unsigned int n, const double* x, double*, void* func_data)
    {
      auto psd = reinterpret_cast<nlopt_f_data*>(func_data);
      param_t pars(x, x + n);
      const double Q = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, psd->trees.size()), 0.0, 
        [&](const tbb::blocked_range<size_t>& r, double q) -> double {
          for (size_t i = r.begin(); i < r.end(); ++i) {
            const double loglik = psd->model->loglik(pars, psd->trees[i]);
            q += loglik * psd->w[i];
          }
          return q;
        },
        std::plus<double>{}
      );
      if (nullptr == psd->conditional) {
        return -Q;
      }
      return -Q * psd->conditional->operator()(pars);
    }
    
  }


  M_step_t M_step(const param_t& pars,
                  const std::vector<tree_t>& trees,          // augmented trees
                  const std::vector<double>& weights,
                  class Model* model,
                  const param_t& lower_bound, // overrides model.lower_bound
                  const param_t& upper_bound, // overrides model.upper.bound
                  double xtol_rel,
                  int num_threads,
                  conditional_fun_t* conditional)
  {
    if (!model->is_threadsafe()) num_threads = 1;
    tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);
    auto T0 = std::chrono::high_resolution_clock::now();
    nlopt_f_data sd{ model, trees, weights, conditional };
    auto M = M_step_t{};
    sbplx nlopt(pars.size());
    M.estimates = pars;
    nlopt.set_xtol_rel(xtol_rel);
    auto lower = lower_bound.empty() ? model->lower_bound() : lower_bound;
    if (!lower.empty()) nlopt.set_lower_bounds(lower);
    auto upper = upper_bound.empty() ? model->upper_bound() : upper_bound;
    if (!upper.empty()) nlopt.set_upper_bounds(upper);
    nlopt.set_min_objective(objective, &sd);
    M.minf = nlopt.optimize(M.estimates);
    M.opt = static_cast<int>(nlopt.result());
    auto T1 = std::chrono::high_resolution_clock::now();
    M.elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count());
    return M;
  }

}
