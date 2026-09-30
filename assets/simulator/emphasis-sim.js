// Forward simulator of the endogenous diversification models in emphasis, for the browser.
//
// A port of sim_tree::general_div::simulate_tree_ltable()
// (inst/include/general_tree.hpp) and emphasis::ed::fair_proportion()
// (inst/include/ed_covariate.hpp).  It follows the C++ step for step: two
// crown lineages at t = 0, Gillespie waiting times from the total rate at the
// current time, the focal lineage and event type drawn from the rates at the
// event time, and a draw that ends as "extinct" as soon as either crown side
// has no lineage left.  Incomplete sampling follows simulate_tree() in
// R/simulate.R: Binomial(n, 1 - rho) extant tips are dropped, at least two kept.
//
// The random stream is not R's, so a seed here does not reproduce an R tree;
// dev/simulator/agreement.R checks that the two simulators draw from the same
// distribution.
//
// Works as a browser script (window.EmphasisSim) and as a Node module.

(function (root) {
  'use strict';

  // xoshiro128** seeded through splitmix32.
  function makeRng(seed) {
    let s = (seed >>> 0) || 1;
    const sm = () => {
      s = (s + 0x9e3779b9) >>> 0;
      let z = s;
      z = Math.imul(z ^ (z >>> 16), 0x85ebca6b) >>> 0;
      z = Math.imul(z ^ (z >>> 13), 0xc2b2ae35) >>> 0;
      return (z ^ (z >>> 16)) >>> 0;
    };
    let a = sm(), b = sm(), c = sm(), d = sm();
    const next = () => {
      const r = Math.imul(rotl(Math.imul(b, 5) >>> 0, 7), 9) >>> 0;
      const t = (b << 9) >>> 0;
      c ^= a; d ^= b; b ^= c; a ^= d; c ^= t; d = rotl(d, 11);
      return r;
    };
    // 53-bit uniform on [0, 1)
    const uniform = () => ((next() >>> 5) * 67108864 + (next() >>> 6)) / 9007199254740992;
    return {
      uniform,
      expon: rate => rate <= 0 ? 1e20 : -Math.log(1 - uniform()) / rate,
      bernoulli: p => uniform() < p,
      int: n => Math.floor(uniform() * n),
    };
  }
  function rotl(x, k) { return ((x << k) | (x >>> (32 - k))) >>> 0; }

  // Model shortcuts -> c(use_N, use_M, use_D, use_ED), as .resolve_model().
  const MODELS = {
    cr: [0, 0, 0, 0], dd: [1, 0, 0, 0], d: [0, 0, 1, 0],
    nd: [1, 0, 1, 0], ed: [0, 0, 0, 1], ned: [1, 0, 0, 1],
    edc: [0, 0, 0, 3], nedc: [1, 0, 0, 3],     // slot 4 = 3: ED centred on its mean over the living
  };
  const LINKS = { linear: 0, exponential: 1, gaussian: 2 };

  // Compact user vector -> the 10-slot layout (.expand_pars()).
  const SLOT_BETA = [1, 2, 3, 8], SLOT_GAMMA = [5, 6, 7, 9];
  function expandPars(pars, model) {
    const active = [0, 1, 2, 3].filter(k => model[k]);
    const nLam = 1 + active.length;
    if (pars.length !== 2 * nLam) throw new Error('pars must have length ' + 2 * nLam);
    const full = new Array(10).fill(0);
    full[0] = pars[0];
    full[4] = pars[nLam];
    active.forEach((k, j) => { full[SLOT_BETA[k]] = pars[1 + j]; full[SLOT_GAMMA[k]] = pars[nLam + 1 + j]; });
    return full;
  }

  // general_div::compute_lambda / compute_mu
  function rateFns(p, link) {
    const f = (i0, iN, iM, iD, iED) => (N, P, E, ED) => {
      const M = N > 0 ? P / N : 0;
      const D = E - M;
      if (link === 2) {
        const d = p[iN] * N + p[iM] * M + p[iD] * D - 1;
        return p[i0] * Math.exp(-0.5 * d * d);
      }
      const eta = p[i0] + p[iN] * N + p[iM] * M + p[iD] * D + p[iED] * ED;
      return link === 1 ? Math.exp(eta) : (eta < 0 ? 0 : eta);
    };
    return { lambda: f(0, 1, 2, 3, 8), mu: f(4, 5, 6, 7, 9) };
  }

  // ed::fair_proportion over the rows born by tt.  Rows are appended in birth
  // order and each row's daughters are appended in birth order, so no sort is
  // needed.  Returns ED per row (NaN for a row not alive).
  function edNow(rows, tt, out) {
    const n = rows.length;
    const nDesc = out.nDesc, inherit = out.inherit, ed = out.ed;
    for (let i = n - 1; i >= 0; i--) {
      const r = rows[i];
      let s = r.end < 0 ? 1 : 0;
      for (const c of r.kids) s += nDesc[c];
      nDesc[i] = s;
    }
    inherit[0] = 0; inherit[1] = 0;
    for (let i = 0; i < n; i++) {
      const r = rows[i], self = r.end < 0 ? 1 : 0;
      let suf = 0;
      for (const c of r.kids) suf += nDesc[c];
      let pre = 0, from = r.birth;
      const inh = inherit[i];
      for (const c of r.kids) {
        const to = rows[c].birth, desc = self + suf;
        if (desc > 0) pre += (to - from) / desc;
        from = to;
        inherit[c] = inh + pre;
        suf -= nDesc[c];
      }
      ed[i] = self ? inh + pre + (tt - from) : NaN;
    }
    return ed;
  }

  // One draw of simulate_tree_ltable().  Rows: {birth, tipStart, parent (row
  // index, -1 for a crown lineage), label, end (-1 while alive), kids}.
  function drawOnce(p, model, link, maxT, maxLin, rng, trace) {
    const R = rateFns(p, link);
    const useEd = model[3] !== 0;
    const centred = model[3] >= 3;             // general_div::ed_now(): EDc = ED - mean over the alive
    const perLineage = model[2] === 1 || useEd;
    const rows = [
      { birth: 0, tipStart: 0, parent: -1, label: -1, end: -1, kids: [] },
      { birth: 0, tipStart: 0, parent: -1, label: 2, end: -1, kids: [] },
    ];
    const alive = [0, 1];            // row indices of the alive lineages
    const pos = [0, 1];              // row -> position in `alive`
    const ws = { nDesc: [], inherit: [], ed: [] };
    let N1 = 1, N2 = 1, t = 0, sumTs = 0, treeId = 3, status = null;
    const rec = trace ? [] : null;

    const snapshot = (N, lamV, muV) => {
      let lo = Infinity, hi = -Infinity, sl = 0, sm = 0;
      for (let k = 0; k < lamV.length; k++) {
        lo = Math.min(lo, lamV[k]); hi = Math.max(hi, lamV[k]); sl += lamV[k]; sm += muV[k];
      }
      rec.push({ t, N, lam: sl / lamV.length, lamLo: lo, lamHi: hi, mu: sm / muV.length });
    };
    const lineageRates = (N, P) => {
      const edv = useEd ? edNow(rows, t, ws) : null;
      if (centred) {
        let m = 0;
        for (const i of alive) m += edv[i];
        m /= alive.length;
        for (const i of alive) edv[i] -= m;
      }
      const lamV = new Array(alive.length), muV = new Array(alive.length);
      for (let k = 0; k < alive.length; k++) {
        const i = alive[k], E = t - rows[i].tipStart, ED = useEd ? edv[i] : 0;
        lamV[k] = R.lambda(N, P, E, ED); muV[k] = R.mu(N, P, E, ED);
      }
      return [lamV, muV];
    };

    if (trace) {
      const [lv, mv] = perLineage ? lineageRates(2, 0) : [[R.lambda(2, 0, 0, 0)], [R.mu(2, 0, 0, 0)]];
      snapshot(2, lv, mv);
    }

    for (;;) {
      const N = N1 + N2;
      if (N >= maxLin) { status = 'too_large'; break; }
      const P = N * t - sumTs;

      let total = 0;
      if (!perLineage) {
        total = N * (R.lambda(N, P, 0, 0) + R.mu(N, P, 0, 0));
      } else {
        const [lv, mv] = lineageRates(N, P);
        for (let k = 0; k < lv.length; k++) total += lv[k] + mv[k];
      }

      const next = t + rng.expon(total);
      if (next >= maxT) { t = maxT; status = 'done'; break; }
      t = next;
      const P2 = N * t - sumTs;

      let focal, isSpec, lamV, muV;
      if (!perLineage) {
        const lam = R.lambda(N, P2, 0, 0), mu = R.mu(N, P2, 0, 0);
        focal = alive[rng.int(alive.length)];
        isSpec = rng.bernoulli(lam / (lam + mu));
        if (trace) { lamV = [lam]; muV = [mu]; }
      } else {
        [lamV, muV] = lineageRates(N, P2);
        let sum = 0;
        for (let k = 0; k < lamV.length; k++) sum += lamV[k] + muV[k];
        const u = rng.uniform() * sum;
        let cum = 0, sel = alive.length - 1;
        for (let k = 0; k < alive.length; k++) {
          cum += lamV[k] + muV[k];
          if (u <= cum) { sel = k; break; }
        }
        focal = alive[sel];
        isSpec = rng.bernoulli(lamV[sel] / (lamV[sel] + muV[sel]));
      }

      const f = rows[focal];
      if (isSpec) {
        sumTs += -f.tipStart + 2 * t;
        f.tipStart = t;
        let id = treeId++;
        if (f.label < 0) { id = -id; N2++; } else N1++;
        const r = rows.length;
        rows.push({ birth: t, tipStart: t, parent: focal, label: id, end: -1, kids: [] });
        f.kids.push(r);
        pos[r] = alive.length; alive.push(r);
      } else {
        sumTs -= f.tipStart;
        f.end = t;
        if (f.label < 0) N2--; else N1--;
        const k = pos[focal], last = alive.pop();
        if (last !== focal) { alive[k] = last; pos[last] = k; }
      }

      if (trace) {
        // rates after the event, at the event time, for the plot
        const Nn = N1 + N2;
        if (Nn > 0) {
          const Pn = Nn * t - sumTs;
          const [lv, mv] = perLineage ? lineageRates(Nn, Pn)
                                      : [[R.lambda(Nn, Pn, 0, 0)], [R.mu(Nn, Pn, 0, 0)]];
          snapshot(Nn, lv, mv);
        }
      }
      if (N1 < 1 || N2 < 1) { status = 'extinct'; break; }
    }
    if (trace && status === 'done') {
      const N = N1 + N2, P = N * t - sumTs;
      const [lv, mv] = perLineage ? lineageRates(N, P) : [[R.lambda(N, P, 0, 0)], [R.mu(N, P, 0, 0)]];
      snapshot(N, lv, mv);
    }
    return { rows, status, trace: rec, tEnd: t };
  }

  // simulate_tree(): retry until a draw survives or maxTries attempts are used,
  // then apply incomplete sampling.
  function simulate(opt) {
    const model = MODELS[opt.model] || opt.model;
    const link = typeof opt.link === 'number' ? opt.link : LINKS[opt.link || 'linear'];
    const p = opt.pars.length === 10 ? opt.pars.slice() : expandPars(opt.pars, model);
    const maxT = opt.maxT, maxLin = opt.maxLin || 1e6, maxTries = opt.maxTries || 1;
    const rho = opt.rho == null ? 1 : opt.rho;
    const rng = opt.rng || makeRng(opt.seed == null ? (Math.random() * 4294967296) >>> 0 : opt.seed);

    let res, attempts = 0;
    do {
      res = drawOnce(p, model, link, maxT, maxLin, rng, !!opt.trace);
      attempts++;
    } while (res.status !== 'done' && attempts < maxTries);

    for (const r of res.rows) r.sampled = r.end < 0;
    if (res.status === 'done' && rho < 1) {
      const ext = res.rows.map((r, i) => i).filter(i => res.rows[i].end < 0);
      let nDrop = 0;
      for (let k = 0; k < ext.length; k++) if (rng.uniform() < 1 - rho) nDrop++;
      nDrop = Math.min(nDrop, ext.length - 2);
      for (let k = 0; k < nDrop; k++) {        // partial Fisher-Yates: sample without replacement
        const j = k + rng.int(ext.length - k);
        [ext[k], ext[j]] = [ext[j], ext[k]];
        res.rows[ext[k]].sampled = false;
      }
    }
    return {
      rows: res.rows, status: res.status, attempts,
      survivalProb: res.status === 'done' ? 1 / attempts : 0,
      trace: res.trace, maxT, pars: p, model, link, rho,
    };
  }

  // Which rows carry a lineage that is alive at time tt (or, with
  // sampledOnly, a sampled extant tip) below them.
  function survivors(rows, tt, sampledOnly) {
    const s = new Array(rows.length);
    for (let i = rows.length - 1; i >= 0; i--) {
      const r = rows[i];
      if (r.birth > tt) { s[i] = false; continue; }
      let v = sampledOnly ? (r.end < 0 && r.sampled) : (r.end < 0 || r.end > tt);
      if (!v) for (const c of r.kids) if (s[c]) { v = true; break; }
      s[i] = v;
    }
    return s;
  }

  const api = { makeRng, simulate, survivors, expandPars, edNow, MODELS, LINKS };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.EmphasisSim = api;
})(typeof window !== 'undefined' ? window : this);
