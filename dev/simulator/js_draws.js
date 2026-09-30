// Per-tree summaries from the browser simulator, for dev/simulator/agreement.R.
//
//   node dev/simulator/js_draws.js scenarios.json K seed out.csv
//
// Each scenario is {name, model, link, pars, maxT, rho}.  Every draw is
// conditioned on survival the way simulate_tree(max_tries = ...) is, and the
// statistics are the ones agreement.R computes from the R L-table.
'use strict';
const fs = require('fs');
const path = require('path');
const sim = require(path.join(__dirname, '..', '..', 'assets', 'simulator', 'emphasis-sim.js'));

const [scenFile, Karg, seedArg, outFile] = process.argv.slice(2);
const scen = JSON.parse(fs.readFileSync(scenFile, 'utf8'));
const K = +Karg;
const rng = sim.makeRng(+seedArg);

const lines = ['scenario,attempts,n_extant,n_sampled,n_rows,n_extinct,pendant_mean,n_complete_half,n_recon_half,n_recon_80'];
for (const s of scen) {
  for (let k = 0; k < K; k++) {
    const r = sim.simulate({ model: s.model, link: s.link, pars: s.pars, maxT: s.maxT,
                             rho: s.rho, maxTries: 10000, maxLin: 1e6, rng });
    if (r.status !== 'done') throw new Error(s.name + ': draw did not survive in 10000 attempts');
    const rows = r.rows, T = s.maxT;
    const ext = rows.filter(x => x.end < 0);
    const at = tt => rows.filter(x => x.birth <= tt && (x.end < 0 || x.end > tt)).length;
    // reconstructed tree: rows with an extant descendant, alive at tt means the
    // row was born by tt and still carries such a descendant past tt
    const surv = sim.survivors(rows, T, false);
    const recon = tt => {
      let n = 0;
      rows.forEach((x, i) => {
        if (!surv[i] || x.birth > tt) return;
        // the row's reconstructed span ends at its own tip (extant) or at the
        // birth of its last surviving daughter
        let endR = x.end < 0 ? Infinity : -Infinity;
        for (const c of x.kids) if (surv[c]) endR = Math.max(endR, rows[c].birth);
        if (endR > tt) n++;
      });
      return n;
    };
    const pend = ext.reduce((a, x) => a + (T - x.tipStart), 0) / ext.length;
    lines.push([s.name, r.attempts, ext.length, ext.filter(x => x.sampled).length, rows.length,
                rows.length - ext.length, pend, at(T / 2), recon(T / 2), recon(0.8 * T)].join(','));
  }
}
fs.writeFileSync(outFile, lines.join('\n') + '\n');
