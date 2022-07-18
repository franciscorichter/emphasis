#include <cassert>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <atomic>
#include <tuple>
#include <memory>
#include "plugin.hpp"
#include "augment_tree.hpp"
#include "model_helpers.hpp"
#include "sbplx.hpp"


namespace emphasis {

  namespace {


    template <typename IT>
    inline IT make_node(IT it, double t, double n, double t_ext)
    {
      it->brts = t; it->n = n; it->t_ext = t_ext; it->pd = 0.0;
      return it;
    }


    template <typename IT>
    inline IT make_extinct_node(IT it, double t, double n)
    {
      it->brts = t; it->n = n; it->t_ext = t_ext_extinct; it->pd = 0.0;
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


    double get_next_bt(const tree_t& tree, double cbt)
    {
      auto it = std::upper_bound(tree.cbegin(), tree.cend(), cbt, detail::node_less{});
      return (it != tree.cend()) ? it->brts : tree.back().brts;
    }


    // insert speciation node_t before t_spec,
    // inserts extinction node_t before t_ext
    // and tracks n
    void insert_species(double t_spec, double t_ext, tree_t& tree)
    {
      auto n_after = [](tree_t::iterator it) { 
        const auto to = detail::is_extinction(*it) ? -1.0 : 1.0;
        return it->n + to; 
      };
      tree.reserve(tree.size() + 2);   // keep iterators valid
      auto first = std::lower_bound(tree.begin(), tree.end(), t_spec, detail::node_less{});
      auto n = (first != tree.begin()) ? n_after(first - 1) : tree.front().n;
      first = make_node(tree.emplace(first), t_spec, n, t_ext);
      // recalculate dirty range
      for (++first; first->brts < t_ext; ++first) {
        first->n = n_after(first - 1);
      }
      make_extinct_node(tree.emplace(first), t_ext, n_after(first - 1));
    }


    void do_augment_tree(const param_t& pars, tree_t& tree, const Model& model, int max_missing, double max_lambda)
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
        double u1 = std::uniform_real_distribution<>()(reng);
        double next_speciation_time = cbt - std::log(u1) / lambda_max;
        if (next_speciation_time < next_bt) {
          double u2 = std::uniform_real_distribution<>()(reng);
          // calc pd(next_speciation_time)
          double pt = std::max(0.0, model.nh_rate(next_speciation_time, pars, tree)) / lambda_max;
          if (u2 < pt) {
            double extinction_time = model.extinction_time(next_speciation_time, pars, tree);
            insert_species(next_speciation_time, extinction_time, tree);
            num_missing_branches++;
            if (num_missing_branches > max_missing) {
              throw augmentation_overrun{};
            }
          }
        }
        cbt = std::min(next_speciation_time, next_bt);
      }
      for (auto& node : tree) {
        node.pd = detail::calculate_pd(node.brts, static_cast<unsigned>(tree.size()), tree.data());
      }
    }


    void do_augment_tree_cont(const param_t& pars, tree_t& tree, const Model& model, int max_missing, double max_lambda)
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
        const double u1 = std::uniform_real_distribution<>()(reng);
        double next_speciation_time = cbt - std::log(u1) / lambda_max;
        dirty = false;
        if (next_speciation_time < next_bt) {
          double u2 = std::uniform_real_distribution<>()(reng);
          double pt = std::max(0.0, model.nh_rate(next_speciation_time, pars, tree)) / lambda_max;
          if (u2 < pt) {
            double extinction_time = model.extinction_time(next_speciation_time, pars, tree);
            insert_species(next_speciation_time, extinction_time, tree);
            num_missing_branches++;
            if (num_missing_branches > max_missing) {
              throw augmentation_overrun{};
            }
            dirty = true;   // tree changed
          }
        }
        cbt = std::min(next_speciation_time, next_bt);
      }
      for (auto& node : tree) {
        node.pd = detail::calculate_pd(node.brts, static_cast<unsigned>(tree.size()), tree.data());
      }
    }

  } // namespace augment


  void augment_tree(const param_t& pars, const tree_t& input_tree, Model* model, int max_missing, double max_lambda, tree_t& pooled)
  {
    pooled.resize(input_tree.size());
    std::copy(input_tree.cbegin(), input_tree.cend(), pooled.begin());
    if (model->numerical_max_lambda()) {
      do_augment_tree(pars, pooled, *model, max_missing, max_lambda);
    }
    else {
      do_augment_tree_cont(pars, pooled, *model, max_missing, max_lambda);
    }
  }


}
