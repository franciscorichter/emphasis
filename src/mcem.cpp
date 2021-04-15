#include <mutex>
#include <atomic>
#include <algorithm>
#include <chrono>
#include <functional>
#include "emphasis.hpp"
#include "plugin.hpp"
#include "augment_tree.hpp"
#include "model_helpers.hpp"


namespace emphasis {


  mcem_t mce(int N,      // sample size
            int maxN,   // max. number augmented trees (incl. invalid)
            const param_t& pars,
            const brts_t& brts,
            class Model* model,
            int soc,
            int max_missing,
            double max_lambda,
            const param_t& lower_bound, // overrides model.lower_bound
            const param_t& upper_bound, // overrides model.upper.bound
            double xtol,
            int num_threads)
  {
    auto EM = mcem_t();
    EM.e = E_step(N, maxN, pars, brts, model, soc, max_missing, max_lambda, num_threads);
    return EM;
  }
  

  mcem_t mcem(int N,      // sample size
              int maxN,   // max. number augmented trees (incl. invalid)
              const param_t& pars,
              const brts_t& brts,
              class Model* model,
              int soc,
              int max_missing,
              double max_lambda,
              const param_t& lower_bound, // overrides model.lower_bound
              const param_t& upper_bound, // overrides model.upper.bound
              double xtol,
              int num_threads,
              conditional_fun_t* conditional)
  {
    auto EM = mcem_t();
    EM.e = E_step(N, maxN, pars, brts, model, soc, max_missing, max_lambda, num_threads);
    // optimize
    if (!EM.e.trees.empty()) {
      EM.m = M_step(pars, EM.e.trees, EM.e.weights, model, lower_bound, upper_bound, xtol, num_threads, conditional);
    }
    return EM;
  }

}
