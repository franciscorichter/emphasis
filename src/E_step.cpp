#include <mutex>
#include <atomic>
#include <algorithm>
#include <numeric>
#include <string>
#include <tbb/tbb.h>
#include "model.hpp"
#include "emphasis.hpp"
#include "augment_tree.hpp"
#include "model_helpers.hpp"
#include "precision_weights.hpp"


namespace emphasis {
  namespace detail {
    // this little addition reduces the load to memory allocator massively.
    tree_t thread_local pooled_tree;

    void inplace_cumsum_of_diff(brts_t& input)
    {
      double sum = 0.0;
      for (size_t i = 1; i < input.size(); ++i) {
        input[i - 1] = sum += (input[i - 1] - input[i]);
      }
      input.back() += sum;
    }


    tree_t create_tree(brts_t brts, double soc)
    {
      inplace_cumsum_of_diff(brts);
      tree_t tree;
      for (size_t i = 0; i < brts.size(); ++i) {
        tree.push_back({ brts[i], soc + i, t_ext_tip });
      }
      std::sort(tree.begin(), tree.end(), detail::node_less{});
      return(tree);
    }

  }


E_step_t E_step(int N,               
                int maxN,
                const param_t& pars,
                const brts_t& brts,
                const Model& model,
                int soc,
                int max_missing,
                double max_lambda,
                int num_threads)
{
  if (!model.is_threadsafe()) num_threads = 1;
  num_threads = std::max(1, std::min(num_threads, static_cast<int>(std::thread::hardware_concurrency())));
  std::mutex mutex;
  std::atomic<bool> stop{ false };
  
  tree_t init_tree = detail::create_tree(brts, static_cast<double>(soc));
  auto E = E_step_t{};
  auto T0 = std::chrono::high_resolution_clock::now();
  
  const int grainsize = maxN / num_threads;
  tbb::task_arena arena(num_threads);
  
  tbb::parallel_for(tbb::blocked_range<unsigned>(0, maxN, grainsize), [&](const tbb::blocked_range<unsigned>& r) {
    for (unsigned i = r.begin(); i < r.end(); ++i) {
      try {
        if (!stop) {
          // reuse tree from pool
          auto& pool_tree = detail::pooled_tree;
          
          emphasis::augment_tree(pars, init_tree, model, max_missing, max_lambda, pool_tree);
          
          double log_w = 0.0;
          double logf = 0.0;
          double logg = 0.0;
          {
            logf = model.loglik(pars, pool_tree);
            logg = model.sampling_prob(pars, pool_tree);
            log_w = logf - logg;
          }
          if (std::isfinite(log_w) && (0.0 < std::exp(log_w))) {
            if (!stop) {
              std::lock_guard<std::mutex> _(mutex);
              E.trees.emplace_back(pool_tree.cbegin(), pool_tree.cend());
              E.weights.push_back(log_w);
              E.logf_.push_back(logf);
              E.logg_.push_back(logg);
              stop = (static_cast<int>(E.trees.size()) == N);
            }
          }
          else {
            std::lock_guard<std::mutex> _(mutex);
            ++E.info.rejected_zero_weights;
          }
        }
      }
      catch (const augmentation_overrun&) {
        std::lock_guard<std::mutex> _(mutex);
        ++E.info.rejected_overruns;
      }
      catch (const augmentation_lambda&) {
        std::lock_guard<std::mutex> _(mutex);
        ++E.info.rejected_lambda;
      }
      catch (...) {
        std::lock_guard<std::mutex> _(mutex);
        ++E.info.rejected;
      }
    }
  });
  E.info.num_trees = static_cast<int>(E.trees.size());
  auto T1 = std::chrono::high_resolution_clock::now();
  E.info.elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count());
  
  if (E.info.num_trees < N) {
    throw emphasis_error_E(E);
  }
  
  const double max_log_w = *std::max_element(E.weights.cbegin(), E.weights.cend());
  double sum_w = calc_sum_w(E.weights.begin(), E.weights.end(), max_log_w);
  E.info.fhat = std::log(sum_w / N) + max_log_w;
  return E;
}

  E_step_info_t E_step_info_grid(int maxN,
                                const param_t& focal_pars,
                                const std::vector<param_t>& all_pars,
                                const brts_t& brts,
                                const Model& model,
                                int soc,
                                int max_missing,
                                double max_lambda)
  {
    int N = 1; // perhaps we can recode this a bit nicer.
    tree_t init_tree = detail::create_tree(brts, static_cast<double>(soc));
    
    auto E = E_step_t{};
    auto T0 = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < maxN; ++i) {
      try {
        // reuse tree from pool
        auto& pool_tree = detail::pooled_tree;
        emphasis::augment_tree(focal_pars, init_tree, model, max_missing, max_lambda, pool_tree);
        
        double log_w = 0.0;
        double logf = 0.0;
        double logg = 0.0;
        {
          logf = model.loglik(focal_pars, pool_tree);
          logg = model.sampling_prob(focal_pars, pool_tree);
          log_w = logf - logg;
        }
        if (std::isfinite(log_w) && (0.0 < std::exp(log_w))) {
         // E.trees.emplace_back(pool_tree.cbegin(), pool_tree.cend());
          E.trees.push_back(pool_tree);
          E.weights.push_back(log_w);
          E.logf_.push_back(logf);
          E.logg_.push_back(logg);
          if (static_cast<int>(E.trees.size()) == N) {
            // now, we need to traverse the grid and update *all* logf/logg values
            
            for (const auto& par : all_pars) {
              E.info.logf_grid.push_back(model.loglik(par, pool_tree));
              E.info.logg_grid.push_back(model.sampling_prob(par, pool_tree));                             
            }

            break;
          }
        }
        else {
          ++E.info.rejected_zero_weights;
        }
      }
      catch (const augmentation_overrun&) {
        ++E.info.rejected_overruns;
      }
      catch (const augmentation_lambda&) {
        ++E.info.rejected_lambda;
      }
      catch (...) {
        ++E.info.rejected;
      }
    }
    E.info.num_trees = static_cast<int>(E.trees.size());
    auto T1 = std::chrono::high_resolution_clock::now();
    E.info.elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count());
    
    if (!E.logf_.empty()) {
      if (E.logf_.size() > 1) {
        E.info.logf = std::accumulate(E.logf_.begin(), E.logf_.end(), 0.0) * 1.0 / E.logf_.size();
        E.info.logg = std::accumulate(E.logg_.begin(), E.logg_.end(), 0.0) * 1.0 / E.logg_.size();
      } else {
        E.info.logf = E.logf_.front();
        E.info.logg = E.logg_.front();
      }
    }
    
    if (E.info.num_trees < N) {
      E.info.fhat = std::nan("");    // ?????  --> TJ: this is NA. Perhaps Inf is better?
    } else {  
      const double max_log_w = *std::max_element(E.weights.cbegin(), E.weights.cend());
      double sum_w = calc_sum_w(E.weights.begin(), E.weights.end(), max_log_w);
      E.info.fhat = std::log(sum_w / N) + max_log_w;
    }
    return E.info;
  }
}
