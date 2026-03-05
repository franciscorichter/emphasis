#include <cassert>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <atomic>
#include <tuple>
#include <memory>
#include <thread>
#include <unordered_map>
#include <tbb/tbb.h>
#include "model.hpp"
#include "augment_tree.hpp"
#include "model_helpers.hpp"
#include "sbplx.hpp"


namespace emphasis {

  namespace {


    template <typename IT>
    inline IT make_node(IT it, double t, double n, double t_ext, int id = -1, int parent_id = -1)
    {
      it->brts = t; it->n = n; it->t_ext = t_ext; it->pd = 0.0;
      it->tip_start = t; it->clade = 0; it->id = id; it->parent_id = parent_id;
      return it;
    }


    template <typename IT>
    inline IT make_extinct_node(IT it, double t, double n, int id = -1, int parent_id = -1)
    {
      it->brts = t; it->n = n; it->t_ext = t_ext_extinct; it->pd = 0.0;
      it->tip_start = 0.0; it->clade = 0; it->id = id; it->parent_id = parent_id;
      return it;
    }


    class maximize_lambda
    {
      using ml_state = std::tuple<const param_t&, const tree_t&, const Model&>;

      static double dx(unsigned int n, const double* x, double*, void* func_data)
      {
        auto ml = *reinterpret_cast<ml_state*>(func_data);
        return std::max(0.0, std::get<2>(ml).nh_rate(*x, std::get<0>(ml), std::get<1>(ml)));
      }

    public:
      explicit maximize_lambda() : nlopt_() {}

      double operator()(double t0, double t1, const param_t& pars, tree_t& tree, const Model& model)
      {
        auto x0 = std::min(t0, t1);
        auto x1 = std::max(t0, t1);
        nlopt_.set_lower_bounds(x0);
        nlopt_.set_upper_bounds(x1);
        nlopt_.set_xtol_rel(0.0001);
        ml_state mls(pars, tree, model);
        nlopt_.set_max_objective(dx, &mls);
        return nlopt_.optimize(x0);
      }

    private:
      sbplx_1 nlopt_;
    };


    auto thread_local reng = detail::make_random_engine<std::default_random_engine>();
    maximize_lambda thread_local tlml;


    // After augmentation, assign tip_start for each speciation node and
    // compute pendant PD per node.
    // tip_start of a lineage = its birth time (brts of its speciation event).
    // For initial tree nodes (no parent), tip_start = 0 (crown age start).
    void compute_pendant_pd(tree_t& tree)
    {
      // First pass: set tip_start for speciation nodes.
      // Speciation nodes are those that are NOT extinction nodes.
      // tip_start = brts (birth time of the lineage).
      // For initial tree nodes (parent_id == -1), tip_start = 0.
      for (auto& node : tree) {
        if (!detail::is_extinction(node)) {
          node.tip_start = (node.parent_id == -1) ? 0.0 : node.brts;
        } else {
          node.tip_start = 0.0;
        }
      }

      // Second pass: compute pendant PD for each node using
      // calculate_pendant_pd which sums (brts - tip_start_i) for alive lineages.
      for (auto& node : tree) {
        node.pd = detail::calculate_pendant_pd(node.brts, tree);
      }
    }


    double get_next_bt(const tree_t& tree, double cbt)
    {
      auto it = std::upper_bound(tree.cbegin(), tree.cend(), cbt, detail::node_less{});
      return (it != tree.cend()) ? it->brts : tree.back().brts;
    }


    // insert speciation node_t before t_spec,
    // inserts extinction node_t before t_ext
    // and tracks n
    void insert_species(double t_spec, double t_ext, tree_t& tree, int id = -1, int parent_id = -1)
    {
      auto n_after = [](tree_t::iterator it) {
        const auto to = detail::is_extinction(*it) ? -1.0 : 1.0;
        return it->n + to;
      };
      tree.reserve(tree.size() + 2);   // keep iterators valid
      auto first = std::lower_bound(tree.begin(), tree.end(), t_spec, detail::node_less{});
      auto n = (first != tree.begin()) ? n_after(first - 1) : tree.front().n;
      first = make_node(tree.emplace(first), t_spec, n, t_ext, id, parent_id);
      // recalculate dirty range
      for (++first; first->brts < t_ext; ++first) {
        first->n = n_after(first - 1);
      }
      make_extinct_node(tree.emplace(first), t_ext, n_after(first - 1), id, parent_id);
    }


    void do_augment_tree(const param_t& pars, tree_t& tree, const Model& model, int max_missing, double max_lambda, int& next_id)
    {
      double cbt = 0;
      tree.reserve(5 * tree.size());    // just a guess, should cover most 'normal' cases
      int num_missing_branches = 0;
      const double b = tree.back().brts;
      auto& ml = tlml;
      while (cbt < b) {
        double next_bt = get_next_bt(tree, cbt);
        double lambda_max = ml(cbt, next_bt, pars, tree, model);
        if (lambda_max > max_lambda) throw augmentation_lambda{};
        double next_speciation_time = next_bt;
        if (0.0 != lambda_max) {
          const double u1 = std::uniform_real_distribution<>()(reng);
          next_speciation_time = cbt - std::log(u1) / lambda_max;
        }
        if (next_speciation_time < next_bt) {
          double u2 = std::uniform_real_distribution<>()(reng);
          // calc pd(next_speciation_time)
          double pt = std::max(0.0, model.nh_rate(next_speciation_time, pars, tree)) / lambda_max;
          if (u2 < pt) {
            double extinction_time = model.extinction_time(next_speciation_time, pars, tree);
            // find lineages alive at next_speciation_time and pick one as parent
            std::vector<int> alive_ids;
            for (const auto& node : tree) {
              if (!detail::is_extinction(node) &&
                  node.brts < next_speciation_time &&
                  node.t_ext > next_speciation_time) {
                alive_ids.push_back(node.id);
              }
            }
            int chosen_parent_id = -1;
            if (!alive_ids.empty()) {
              std::uniform_int_distribution<size_t> uid(0, alive_ids.size() - 1);
              chosen_parent_id = alive_ids[uid(reng)];
            }
            int new_id = next_id++;
            insert_species(next_speciation_time, extinction_time, tree, new_id, chosen_parent_id);
            num_missing_branches++;
            if (num_missing_branches > max_missing) {
              throw augmentation_overrun{};
            }
          }
        }
        cbt = std::min(next_speciation_time, next_bt);
      }
      compute_pendant_pd(tree);
    }


    void do_augment_tree_cont(const param_t& pars, tree_t& tree, const Model& model, int max_missing, double max_lambda, int& next_id)
    {
      double cbt = 0;
      tree.reserve(5 * tree.size());    // just a guess, should cover most 'normal' cases
      int num_missing_branches = 0;
      const double b = tree.back().brts;
      double lambda2 = 0.0;
      bool dirty = true;
      while (cbt < b) {
        double next_bt = get_next_bt(tree, cbt);
        double lambda1 = (!dirty) ? lambda2 : std::max(0.0, model.nh_rate(cbt, pars, tree));
        lambda2 = std::max(0.0, model.nh_rate(next_bt, pars, tree));
        double lambda_max = std::max<double>(lambda1, lambda2);
        if (lambda_max > max_lambda) throw augmentation_lambda{};
        double next_speciation_time = next_bt;
        if (0.0 != lambda_max) {
          const double u1 = std::uniform_real_distribution<>()(reng);
          next_speciation_time = cbt - std::log(u1) / lambda_max;
        }
        dirty = false;
        if (next_speciation_time < next_bt) {
          double u2 = std::uniform_real_distribution<>()(reng);
          double pt = std::max(0.0, model.nh_rate(next_speciation_time, pars, tree)) / lambda_max;
          if (u2 < pt) {
            double extinction_time = model.extinction_time(next_speciation_time, pars, tree);
            // find lineages alive at next_speciation_time and pick one as parent
            std::vector<int> alive_ids;
            for (const auto& node : tree) {
              if (!detail::is_extinction(node) &&
                  node.brts < next_speciation_time &&
                  node.t_ext > next_speciation_time) {
                alive_ids.push_back(node.id);
              }
            }
            int chosen_parent_id = -1;
            if (!alive_ids.empty()) {
              std::uniform_int_distribution<size_t> uid(0, alive_ids.size() - 1);
              chosen_parent_id = alive_ids[uid(reng)];
            }
            int new_id = next_id++;
            insert_species(next_speciation_time, extinction_time, tree, new_id, chosen_parent_id);
            num_missing_branches++;
            if (num_missing_branches > max_missing) {
              throw augmentation_overrun{};
            }
            dirty = true;   // tree changed
          }
        }
        cbt = std::min(next_speciation_time, next_bt);
      }
      compute_pendant_pd(tree);
    }

  } // namespace augment


  void augment_tree(const param_t& pars, const tree_t& input_tree, const Model& model, int max_missing, double max_lambda, tree_t& pooled)
  {
    pooled.resize(input_tree.size());
    std::copy(input_tree.cbegin(), input_tree.cend(), pooled.begin());
    // assign sequential IDs to initial tree nodes; augmented nodes get IDs starting after
    int next_id = 0;
    for (auto& node : pooled) {
      if (!detail::is_extinction(node)) {
        node.id = next_id++;
        node.parent_id = -1;
      } else {
        node.id = -1;
        node.parent_id = -1;
      }
    }
    if (model.numerical_max_lambda()) {
      do_augment_tree(pars, pooled, model, max_missing, max_lambda, next_id);
    }
    else {
      do_augment_tree_cont(pars, pooled, model, max_missing, max_lambda, next_id);
    }
  }


  // returns one augmented tree per vpars
  // failures results in empty tree
  std::vector<tree_t> augment_trees(const std::vector<param_t>& vpars, const tree_t& input_tree, const Model& model, int max_missing, double max_lambda, int num_threads)
  {
    if (!model.is_threadsafe()) num_threads = 1;
    num_threads = std::max(1, std::min(num_threads, static_cast<int>(std::thread::hardware_concurrency())));
    tbb::task_arena arena(num_threads);
    std::vector<tree_t> trees(vpars.size(), input_tree);
    const size_t grainsize = vpars.size() / num_threads;
    tbb::parallel_for(tbb::blocked_range<size_t>(0ull, vpars.size(), grainsize), [&](const tbb::blocked_range<size_t>& r) {
      for (size_t i = r.begin(); i < r.end(); ++i) {
        try {
          int next_id = 0;
          for (auto& node : trees[i]) {
            if (!detail::is_extinction(node)) {
              node.id = next_id++;
              node.parent_id = -1;
            } else {
              node.id = -1;
              node.parent_id = -1;
            }
          }
          if (model.numerical_max_lambda()) {
            do_augment_tree(vpars[i], trees[i], model, max_missing, max_lambda, next_id);
          }
          else {
            do_augment_tree_cont(vpars[i], trees[i], model, max_missing, max_lambda, next_id);
          }
        } catch (...) {
          trees[i].clear();
        }
      }
    });
    return trees;
  }


}
