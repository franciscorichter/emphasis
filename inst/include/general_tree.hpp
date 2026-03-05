//
//  general_tree.hpp
//
//  General diversification model with optional N, P, E dependence.
//
//  Rates:
//    lambda(s,t) = max(0, beta_0  + beta_N  * N(t) + beta_P  * P(t) + beta_E  * E(s,t))
//    mu(s,t)     = max(0, gamma_0 + gamma_N * N(t) + gamma_P * P(t) + gamma_E * E(s,t))
//
//  where E(s,t) = t - tip_start(s) is the current pendant edge length of lineage s,
//  P(t) = sum of pendant lengths of alive lineages = N*t - sum(tip_start),
//  and the model binary vector selects which covariates are active.
//
//  Parameters (always 8, zero-padded for inactive terms):
//    pars = { beta_0, beta_N, beta_P, beta_E, gamma_0, gamma_N, gamma_P, gamma_E }
//  Model:
//    model = { use_N, use_P, use_E }  (each 0 or 1)
//

#ifndef general_tree_h
#define general_tree_h

#include <array>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <thread>

namespace sim_tree {

struct rnd_t {
  std::mt19937 rndgen_;

  rnd_t() {
    const auto tt = static_cast<int64_t>(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());
    auto tid = std::this_thread::get_id();
    const uint64_t e3{ std::hash<std::remove_const_t<decltype(tid)>>()(tid) };
    rndgen_.seed(static_cast<unsigned int>(tt ^ static_cast<int64_t>(e3)));
  }

  double expon(double lambda) {
    if (lambda <= 0.0) return 1e20;
    return std::exponential_distribution<double>(lambda)(rndgen_);
  }

  bool bernouilli(double p) {
    return std::bernoulli_distribution(p)(rndgen_);
  }

  size_t random_number(size_t n) {
    if (n <= 1) return 0;
    return std::uniform_int_distribution<size_t>(0, n - 1)(rndgen_);
  }

  double uniform() {
    return std::uniform_real_distribution<double>(0.0, 1.0)(rndgen_);
  }
};

// Branch record for the L-table.
// start_date : forward birth time of this lineage (used for L-table output)
// tip_start  : forward time when this lineage last became a pendant tip
//              (= start_date at birth; reset to speciation time when parent speciates)
// parent_label, label, end_date : standard L-table fields
struct gbranch {
  float start_date;
  float tip_start;
  int   parent_label;
  int   label;
  float end_date;   // -1 if still alive

  gbranch(float bd, float ts, int pl, int lab, float ext)
    : start_date(bd), tip_start(ts), parent_label(pl), label(lab), end_date(ext) {}
};

enum breaks { none, finished, extinction, maxN_exceeded };

// General diversification simulator.
struct general_div {
  double max_t;
  float  t;
  size_t max_N, N;
  std::array<double, 8> pars;   // beta_0..beta_E, gamma_0..gamma_E
  std::array<int, 3>    model;  // {use_N, use_P, use_E}
  int    link;                  // 0 = linear (max(0,...)), 1 = exponential

  std::vector<gbranch> ltable;
  breaks break_type;
  rnd_t  rndgen;

  general_div(double total_time,
              const std::array<double, 8>& p,
              const std::array<int, 3>&    m,
              size_t maxN,
              int link_type = 0)
    : max_t(total_time), t(0.f), max_N(maxN), N(0),
      pars(p), model(m), link(link_type), break_type(none), rndgen() {}

  // Apply link function to linear predictor
  double apply_link(double eta) const {
    if (link == 1) return std::exp(eta);
    return eta < 0.0 ? 0.0 : eta;
  }

  // Per-lineage speciation rate
  double compute_lambda(double Nval, double Pval, double Eval) const {
    double eta = pars[0] + pars[1] * Nval + pars[2] * Pval + pars[3] * Eval;
    return apply_link(eta);
  }

  // Per-lineage extinction rate
  double compute_mu(double Nval, double Pval, double Eval) const {
    double eta = pars[4] + pars[5] * Nval + pars[6] * Pval + pars[7] * Eval;
    return apply_link(eta);
  }

  // Uniform random alive lineage (rejection sampling)
  size_t sample_tip() {
    size_t idx = rndgen.random_number(ltable.size());
    while (ltable[idx].end_date != -1.f)
      idx = rndgen.random_number(ltable.size());
    return idx;
  }

  // Forward Gillespie simulation; fills ltable.
  // Returns true if the tree went extinct.
  bool simulate_tree_ltable() {
    size_t N1 = 1, N2 = 1;
    N = 2;
    t = 0.f;

    ltable.clear();
    ltable.push_back(gbranch(0.f, 0.f,  0, -1, -1.f));
    ltable.push_back(gbranch(0.f, 0.f, -1,  2, -1.f));
    int tree_id = 3;
    break_type  = none;

    // P(t) = N*t - sum_tip_start  (pendant PD of alive lineages)
    double sum_tip_start = 0.0;

    while (true) {
      N = N1 + N2;
      if (N >= max_N) { break_type = maxN_exceeded; break; }

      const double Nval = static_cast<double>(N);
      const double Pval = Nval * static_cast<double>(t) - sum_tip_start;

      // ------------------------------------------------------------------
      // Compute total event rate
      // ------------------------------------------------------------------
      double total_rate = 0.0;

      if (!model[2]) {
        // No E-dependence: all lineages share the same per-lineage rates
        const double lam = compute_lambda(Nval, Pval, 0.0);
        const double mu  = compute_mu   (Nval, Pval, 0.0);
        total_rate = Nval * (lam + mu);
      } else {
        // E-dependence: sum per-lineage rates
        for (const auto& br : ltable) {
          if (br.end_date != -1.f) continue;
          const double E = static_cast<double>(t) - static_cast<double>(br.tip_start);
          total_rate += compute_lambda(Nval, Pval, E) + compute_mu(Nval, Pval, E);
        }
      }

      const double next_time = static_cast<double>(t) + rndgen.expon(total_rate);

      if (next_time >= max_t) {
        t = static_cast<float>(max_t);
        break_type = finished;
        break;
      }
      t = static_cast<float>(next_time);

      // Recompute P at new t for event selection
      const double Pval2 = Nval * static_cast<double>(t) - sum_tip_start;

      // ------------------------------------------------------------------
      // Choose focal lineage and event type
      // ------------------------------------------------------------------
      size_t focal;
      bool   is_spec;

      if (!model[2]) {
        const double lam = compute_lambda(Nval, Pval2, 0.0);
        const double mu  = compute_mu   (Nval, Pval2, 0.0);
        focal   = sample_tip();
        is_spec = rndgen.bernouilli(lam / (lam + mu));
      } else {
        // Weighted selection: precompute per-lineage rates at new t, then sample
        std::vector<size_t> alive_idx;
        std::vector<double> lam_vec, mu_vec;
        double rate_sum = 0.0;
        for (size_t i = 0; i < ltable.size(); ++i) {
          if (ltable[i].end_date != -1.f) continue;
          const double E   = static_cast<double>(t) - static_cast<double>(ltable[i].tip_start);
          const double lam = compute_lambda(Nval, Pval2, E);
          const double mu  = compute_mu   (Nval, Pval2, E);
          alive_idx.push_back(i);
          lam_vec.push_back(lam);
          mu_vec.push_back(mu);
          rate_sum += lam + mu;
        }

        const double u = rndgen.uniform() * rate_sum;
        double cum = 0.0;
        size_t sel = alive_idx.size() - 1;  // fallback: last alive
        for (size_t k = 0; k < alive_idx.size(); ++k) {
          cum += lam_vec[k] + mu_vec[k];
          if (u <= cum) { sel = k; break; }
        }
        focal   = alive_idx[sel];
        is_spec = rndgen.bernouilli(lam_vec[sel] / (lam_vec[sel] + mu_vec[sel]));
      }

      // ------------------------------------------------------------------
      // Apply event
      // ------------------------------------------------------------------
      if (is_spec) {
        // Speciation: parent continues (tip_start reset), new daughter added
        const double old_ts = static_cast<double>(ltable[focal].tip_start);
        ltable[focal].tip_start = t;
        sum_tip_start += -old_ts + static_cast<double>(t) + static_cast<double>(t);

        int new_id = tree_id++;
        if (ltable[focal].label < 0) { new_id = -new_id; N2++; } else { N1++; }

        ltable.push_back(gbranch(t, t, ltable[focal].label, new_id, -1.f));
      } else {
        // Extinction
        sum_tip_start -= static_cast<double>(ltable[focal].tip_start);
        ltable[focal].end_date = t;
        if (ltable[focal].label < 0) N2--; else N1--;
      }

      if (N1 < 1 || N2 < 1) { break_type = extinction; break; }
    }

    N = N1 + N2;
    return break_type == extinction;
  }
};

}  // namespace sim_tree

#endif /* general_tree_h */
