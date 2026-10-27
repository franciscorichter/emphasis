//
//  general_tree.hpp
//
//  General diversification model with optional N, M, D dependence.
//
//  Orthogonal covariate basis:
//    N(t)     = number of alive lineages
//    M(t)=P/N = mean pendant age, P(t) = sum of pendant lengths = N*t - sum(tip_start)
//    D(s,t)   = E(s,t) - M(t),  E(s,t) = t - tip_start(s)   (Sum_s D = 0)
//
//  Rates:
//    lambda(s,t) = max(0, beta_0  + beta_N  * N + beta_M  * M + beta_D  * D)
//    mu(s,t)     = max(0, gamma_0 + gamma_N * N + gamma_M * M + gamma_D * D)
//
//  Parameters (always 8, zero-padded for inactive terms); slots 2/6 -> M, 3/7 -> D:
//    pars = { beta_0, beta_N, beta_M, beta_D, gamma_0, gamma_N, gamma_M, gamma_D }
//  Model:
//    model = { use_N, use_M, use_D }  (each 0 or 1)
//

#ifndef general_tree_h
#define general_tree_h

#include <array>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include "ed_covariate.hpp"

namespace sim_tree {

// The simulator's engine, seeded from the caller's integer rather than from
// the wall clock XOR the thread id, so that set.seed() in R reaches it
// (src/div_tree.cpp draws the integer; see the rng notes in model.hpp).  Same
// engine type and same seed_seq recipe as the augmenter's streams.
struct rnd_t {
  std::mt19937_64 rndgen_;

  explicit rnd_t(uint64_t seed) {
    const uint32_t w[5] = {
      static_cast<uint32_t>(seed), static_cast<uint32_t>(seed >> 32),
      0u, 0u, 0x9e3779b9u
    };
    std::seed_seq sseq(w, w + 5);
    rndgen_.seed(sseq);
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
  int   k;          // ancestral splits: speciation events on the path from the crown
  int   mark;       // latent mark (dev/latent-marks/design.md); 0 without marks

  gbranch(float bd, float ts, int pl, int lab, float ext, int kk = 0, int mk = 0)
    : start_date(bd), tip_start(ts), parent_label(pl), label(lab), end_date(ext), k(kk), mark(mk) {}
};

enum breaks { none, finished, extinction, maxN_exceeded };

// General diversification simulator.
struct general_div {
  double max_t;
  float  t;
  size_t max_N, N;
  // {beta_0, beta_N, beta_M, beta_D, gamma_0, gamma_N, gamma_M, gamma_D,
  //  beta_ED, gamma_ED, beta_K, gamma_K}: the ED and K coefficients are
  //  appended, as in model.hpp.
  std::array<double, 12> pars;
  std::array<int, 5>     model;  // {use_N, use_M, use_D, use_ED, use_K}
  int    link;                   // 0 = linear, 1 = exponential

  std::vector<gbranch> ltable;
  breaks break_type;
  rnd_t  rndgen;

  // Latent marks (model.hpp, marks_t): S marks multiplying a lineage's rates,
  // drawn at the crown from pi and at every split for both daughters, from
  // the joint kernel C (S^3, read symmetrised by a fair swap) or, when P
  // (S x S) is set, independently.  With S = 1 nothing is drawn, so the
  // random stream and every tree are what they were without marks.
  int mS = 1;
  std::vector<double> m_clam{1.0}, m_cmu{1.0}, m_C, m_P, m_pi{1.0};

  void set_marks(int S, const std::vector<double>& log_clam, const std::vector<double>& log_cmu,
                 const std::vector<double>& C, const std::vector<double>& P,
                 const std::vector<double>& pi) {
    mS = S; m_C = C; m_P = P; m_pi = pi;
    m_clam.resize(static_cast<size_t>(S)); m_cmu.resize(static_cast<size_t>(S));
    for (size_t m = 0; m < static_cast<size_t>(S); ++m) {
      m_clam[m] = std::exp(log_clam[m]); m_cmu[m] = std::exp(log_cmu[m]);
    }
  }
  bool has_marks() const { return mS > 1; }
  // An index drawn from n weights starting at w (not necessarily normalised).
  int draw_index(const double* w, size_t n) {
    double tot = 0.0;
    for (size_t i = 0; i < n; ++i) tot += w[i];
    const double u = rndgen.uniform() * tot;
    double cum = 0.0;
    for (size_t i = 0; i < n; ++i) { cum += w[i]; if (u < cum) return static_cast<int>(i); }
    return static_cast<int>(n - 1);
  }
  // The two daughters' marks from a parent's mark m.
  void draw_daughters(int m, int& a, int& b) {
    const size_t S = static_cast<size_t>(mS), mm = static_cast<size_t>(m);
    if (!m_P.empty()) {
      a = draw_index(m_P.data() + mm * S, S);
      b = draw_index(m_P.data() + mm * S, S);
    } else {
      const int ab = draw_index(m_C.data() + mm * S * S, S * S);
      a = ab / mS; b = ab % mS;
      if (rndgen.uniform() < 0.5) std::swap(a, b);
    }
  }

  general_div(double total_time,
              const std::array<double, 12>& p,
              const std::array<int, 5>&     m,
              size_t maxN,
              int link_type,
              uint64_t seed)
    : max_t(total_time), t(0.f), max_N(maxN), N(0),
      pars(p), model(m), link(link_type), break_type(none), rndgen(seed) {}

  // Evolutionary distinctiveness of every row of the L-table at time tt:
  // ED_i(tt) for an alive row, NaN otherwise.  The first two rows are the
  // crown lineages and have no parent; every other row's parent is the row
  // whose label its parent_label names.  Same routine as the likelihood's
  // (ed_covariate.hpp), so the simulator and the scorer agree on ED.
  std::vector<double> ed_now(double tt) const {
    const size_t n = ltable.size();
    std::vector<int>    parent(n, emphasis::ed::no_parent_idx);
    std::vector<double> birth(n);
    std::vector<char>   alive(n);
    std::vector<int>    row_of_label;   // label -> row, labels are small ints of either sign
    int lo = 0, hi = 0;
    for (const auto& br : ltable) { lo = std::min(lo, br.label); hi = std::max(hi, br.label); }
    row_of_label.assign(static_cast<size_t>(hi - lo + 1), -1);
    for (size_t i = 0; i < n; ++i) row_of_label[static_cast<size_t>(ltable[i].label - lo)] = static_cast<int>(i);
    for (size_t i = 0; i < n; ++i) {
      birth[i] = static_cast<double>(ltable[i].start_date);
      alive[i] = (ltable[i].end_date == -1.f) ? 1 : 0;
      if (i >= 2) {
        const int pl = ltable[i].parent_label;
        if (pl >= lo && pl <= hi) parent[i] = row_of_label[static_cast<size_t>(pl - lo)];
      }
    }
    emphasis::ed::result_t r;
    emphasis::ed::fair_proportion(parent, birth, alive, tt, r);
    std::vector<double> out(n, std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < n; ++i) if (alive[i]) out[i] = r.c[i] + (tt - r.ts[i]);
    // Centred covariate (model[3] >= 3, model.hpp Model::ed_centred): ED minus
    // its mean over the lineages alive now.
    if (model[3] >= 3) {
      double s = 0.0; int m = 0;
      for (size_t i = 0; i < n; ++i) if (alive[i]) { s += out[i]; ++m; }
      const double mean = m ? s / m : 0.0;
      for (size_t i = 0; i < n; ++i) if (alive[i]) out[i] -= mean;
    }
    return out;
  }

  // Apply link function to linear predictor (linear and exponential only)
  double apply_link(double eta) const {
    if (link == 1) return std::exp(eta);
    return eta < 0.0 ? 0.0 : eta;
  }

  // Orthogonal covariate basis {N, M=P/N, D=E-M}. Inputs are still (N, P, E);
  // the transform to (M, D) happens here so callers pass raw covariates.

  // Per-lineage speciation rate; EDval is the lineage's evolutionary
  // distinctiveness (0 when the ED covariate is inactive)
  double compute_lambda(double Nval, double Pval, double Eval, double EDval = 0.0,
                        double Kval = 0.0) const {
    const double M = (Nval > 0.0) ? Pval / Nval : 0.0;
    const double D = Eval - M;
    double eta = pars[0] + pars[1] * Nval + pars[2] * M + pars[3] * D + pars[8] * EDval
               + pars[10] * Kval;
    return apply_link(eta);
  }

  // Per-lineage extinction rate
  double compute_mu(double Nval, double Pval, double Eval, double EDval = 0.0,
                    double Kval = 0.0) const {
    const double M = (Nval > 0.0) ? Pval / Nval : 0.0;
    const double D = Eval - M;
    double eta = pars[4] + pars[5] * Nval + pars[6] * M + pars[7] * D + pars[9] * EDval
               + pars[11] * Kval;
    return apply_link(eta);
  }

  // The rates of row i of the L-table, its mark's multiplier included.
  double lambda_of(size_t i, double Nval, double Pval, double Eval, double EDval, double Kval) const {
    return compute_lambda(Nval, Pval, Eval, EDval, Kval) * m_clam[static_cast<size_t>(ltable[i].mark)];
  }
  double mu_of(size_t i, double Nval, double Pval, double Eval, double EDval, double Kval) const {
    return compute_mu(Nval, Pval, Eval, EDval, Kval) * m_cmu[static_cast<size_t>(ltable[i].mark)];
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
    if (has_marks()) {
      int a = 0, b = 0;
      draw_daughters(draw_index(m_pi.data(), static_cast<size_t>(mS)), a, b);
      ltable[0].mark = a; ltable[1].mark = b;
    }

    // P(t) = N*t - sum_tip_start  (pendant PD of alive lineages)
    double sum_tip_start = 0.0;
    // Running sums for O(1) total-rate under exponential EP:
    //   sum_exp_bE = Σ_{alive s} exp(-β_E · tip_start_s)
    //   sum_exp_gE = Σ_{alive s} exp(-γ_E · tip_start_s)
    // Allows: total_λ = exp(β₀+β_N·N+β_P·P+β_E·t) · sum_exp_bE
    // Per-lineage rates are needed for D and for ED; the factored running sums
    // serve D alone under the exponential link, since a lineage's ED constant
    // changes at every event in its clade and cannot be maintained
    // incrementally.
    const bool use_ed      = (model[3] != 0);
    const bool use_k       = (model[4] != 0);
    // K is per lineage like D and ED, and K alone takes the same per-lineage
    // loop the linear-link D takes.
    const bool per_lineage = (model[2] == 1) || use_ed || use_k || has_marks();
    const bool ep_exp = (model[2] == 1 && !use_ed && !use_k && !has_marks() && link == 1);
    double sum_exp_bE = ep_exp ? 2.0 * std::exp(-pars[3] * 0.0) : 0.0;  // 2 lineages at ts=0
    double sum_exp_gE = ep_exp ? 2.0 * std::exp(-pars[7] * 0.0) : 0.0;

    // The L-table keeps every lineage ever born, and the per-lineage rates and
    // ED read all of it at every event, so a clade that churns -- high
    // speciation and extinction over a long crown -- makes each event slower
    // than the last while N stays under max_N.  On a 30 Ma crown a corner of
    // auto_bounds()'s box ran one simulation for hours this way.  So the table
    // is capped too: ten times max_N lineages ever born is "too large" like N.
    const size_t max_rows = (max_N > std::numeric_limits<size_t>::max() / 10) ? max_N : 10 * max_N;

    while (true) {
      N = N1 + N2;
      if (N >= max_N || ltable.size() >= max_rows) { break_type = maxN_exceeded; break; }

      const double Nval = static_cast<double>(N);
      const double Pval = Nval * static_cast<double>(t) - sum_tip_start;

      // ------------------------------------------------------------------
      // Compute total event rate
      // ------------------------------------------------------------------
      double total_rate = 0.0;
      double ed_next = -1.0;   // set by the raw-ED branch, which draws its own waiting time

      if (!per_lineage) {
        // No E-dependence: all lineages share the same per-lineage rates — O(1)
        const double lam = compute_lambda(Nval, Pval, 0.0);
        const double mu  = compute_mu   (Nval, Pval, 0.0);
        total_rate = Nval * (lam + mu);
      } else if (use_ed) {
        // ED (with or without D and K): O(N) per-lineage loop over the alive
        // rows, each with its own pendant age and its own ED at the current time.
        const std::vector<double> edv = ed_now(static_cast<double>(t));
        for (size_t i = 0; i < ltable.size(); ++i) {
          const auto& br = ltable[i];
          if (br.end_date != -1.f) continue;
          const double E = static_cast<double>(t) - static_cast<double>(br.tip_start);
          total_rate += lambda_of(i, Nval, Pval, E, edv[i], br.k) + mu_of(i, Nval, Pval, E, edv[i], br.k);
        }
        // Raw ED grows at unit slope between events, so the hazard is not
        // constant until the next event: under the exponential link every
        // lineage's speciation rate carries the factor exp(beta_ED * dt).  The
        // waiting time is drawn by thinning (Ogata): propose from a bound on
        // the hazard over a horizon, accept with the ratio of the hazard at
        // the proposed time to the bound.  With the coefficients at or below
        // zero the hazard falls, so the bound is the current total and the
        // horizon is the rest of the crown age; with a positive coefficient
        // the horizon is the time over which the hazard can grow by half.
        // Holding the rate fixed instead (the code before this) put events
        // early by the rate's change over one waiting time.  The centred
        // covariate is constant between events and keeps the plain draw.
        if (model[3] < 3) {
          const double base = static_cast<double>(t);
          const double bpos = std::max(pars[8], 0.0), gpos = std::max(pars[9], 0.0);
          auto sums_at = [&](double dt, double& L, double& M) {
            const double tt = base + dt;
            const double Pv = Nval * tt - sum_tip_start;
            L = 0.0; M = 0.0;
            for (size_t i = 0; i < ltable.size(); ++i) {
              const auto& br = ltable[i];
              if (br.end_date != -1.f) continue;
              const double E = tt - static_cast<double>(br.tip_start);
              L += lambda_of(i, Nval, Pv, E, edv[i] + dt, br.k);
              M += mu_of    (i, Nval, Pv, E, edv[i] + dt, br.k);
            }
          };
          double off = 0.0; double found = -1.0;
          for (int guard = 0; guard < 1000000; ++guard) {
            const double remain = max_t - (base + off);
            if (remain <= 0.0) break;
            double L0, M0; sums_at(off, L0, M0);
            double h = remain;
            if (link == 1) {
              const double hp = std::max(bpos, gpos);
              if (hp > 0.0) h = std::min(h, std::log(1.5) / hp);
            } else if (link == 0) {
              const double g = Nval * (bpos + gpos);
              if (g > 0.0) h = std::min(h, std::max(0.5 * (L0 + M0), 1.0) / g);
            }
            double Lh = L0, Mh = M0;
            if (bpos > 0.0 || gpos > 0.0) sums_at(off + h, Lh, Mh);
            const double bound = (bpos > 0.0 ? Lh : L0) + (gpos > 0.0 ? Mh : M0);
            if (bound <= 0.0) { if (h >= remain) break; off += h; continue; }
            const double dt = rndgen.expon(bound);
            if (dt > h) { off += h; continue; }
            const double cand = off + dt;
            if (base + cand >= max_t) break;
            double Lc, Mc; sums_at(cand, Lc, Mc);
            if (rndgen.uniform() * bound <= Lc + Mc) { found = cand; break; }
            off = cand;
          }
          ed_next = (found >= 0.0) ? base + found : max_t;
        }
      } else if (ep_exp) {
        // EP + exponential link, orthogonal basis {N, M=P/N, D=E-M}: O(1) via
        // factored running sums.  With D = (t - ts_s) - M,
        //   λ(s,t) = exp(β₀+β_N·N+(β_M-β_D)·M+β_D·t) · exp(-β_D·ts_s)
        //   Σ_s λ(s,t) = exp(A'_λ + β_D·t) · sum_exp_bE,
        //   A'_λ = β₀+β_N·N+(β_M-β_D)·M,   M = P/N.
        const double td = static_cast<double>(t);
        const double M = (Nval > 0.0) ? Pval / Nval : 0.0;
        total_rate = std::exp(pars[0] + pars[1]*Nval + (pars[2]-pars[3])*M + pars[3]*td) * sum_exp_bE
                   + std::exp(pars[4] + pars[5]*Nval + (pars[6]-pars[7])*M + pars[7]*td) * sum_exp_gE;
      } else {
        // D under the linear link, or K: O(N) per-lineage loop (max(0,·)
        // clipping, or the per-lineage integer, prevents factoring)
        for (size_t i = 0; i < ltable.size(); ++i) {
          const auto& br = ltable[i];
          if (br.end_date != -1.f) continue;
          const double E = static_cast<double>(t) - static_cast<double>(br.tip_start);
          total_rate += lambda_of(i, Nval, Pval, E, 0.0, br.k) + mu_of(i, Nval, Pval, E, 0.0, br.k);
        }
      }

      const double next_time = (ed_next >= 0.0) ? ed_next
                             : static_cast<double>(t) + rndgen.expon(total_rate);

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

      if (!per_lineage) {
        const double lam = compute_lambda(Nval, Pval2, 0.0);
        const double mu  = compute_mu   (Nval, Pval2, 0.0);
        focal   = sample_tip();
        is_spec = rndgen.bernouilli(lam / (lam + mu));
      } else {
        // Weighted selection: precompute per-lineage rates at new t, then sample
        const std::vector<double> edv = use_ed ? ed_now(static_cast<double>(t))
                                               : std::vector<double>(ltable.size(), 0.0);
        std::vector<size_t> alive_idx;
        std::vector<double> lam_vec, mu_vec;
        double rate_sum = 0.0;
        for (size_t i = 0; i < ltable.size(); ++i) {
          if (ltable[i].end_date != -1.f) continue;
          const double E   = static_cast<double>(t) - static_cast<double>(ltable[i].tip_start);
          const double lam = lambda_of(i, Nval, Pval2, E, edv[i], ltable[i].k);
          const double mu  = mu_of    (i, Nval, Pval2, E, edv[i], ltable[i].k);
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
        const double td = static_cast<double>(t);
        ltable[focal].tip_start = t;
        sum_tip_start += -old_ts + td + td;
        // Maintain exponential running sums: remove old parent, add parent+daughter at new ts
        if (ep_exp) {
          sum_exp_bE += -std::exp(-pars[3] * old_ts) + 2.0 * std::exp(-pars[3] * td);
          sum_exp_gE += -std::exp(-pars[7] * old_ts) + 2.0 * std::exp(-pars[7] * td);
        }

        int new_id = tree_id++;
        if (ltable[focal].label < 0) { new_id = -new_id; N2++; } else { N1++; }

        // both daughters carry one more ancestral split, and draw their marks
        ltable[focal].k += 1;
        int ma = 0, mb = 0;
        if (has_marks()) draw_daughters(ltable[focal].mark, ma, mb);
        if (has_marks()) ltable[focal].mark = ma;
        ltable.push_back(gbranch(t, t, ltable[focal].label, new_id, -1.f, ltable[focal].k, mb));
      } else {
        // Extinction
        const double dead_ts = static_cast<double>(ltable[focal].tip_start);
        sum_tip_start -= dead_ts;
        if (ep_exp) {
          sum_exp_bE -= std::exp(-pars[3] * dead_ts);
          sum_exp_gE -= std::exp(-pars[7] * dead_ts);
        }
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
