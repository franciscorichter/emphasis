// Trees from the browser simulator with the ED it assigns every extant tip at
// the present, for agreement.R to rescore with emphasis:::ed_fair_proportion().
//
//   node dev/simulator/js_ed.js n_trees seed out.csv
'use strict';
const fs = require('fs');
const path = require('path');
const sim = require(path.join(__dirname, '..', '..', 'assets', 'simulator', 'emphasis-sim.js'));

const [nArg, seedArg, outFile] = process.argv.slice(2);
const rng = sim.makeRng(+seedArg);
const lines = ['tree,row,parent,birth,alive,ed_js,t'];
for (let k = 0; k < +nArg; k++) {
  const T = 5 + 3 * rng.uniform();
  const r = sim.simulate({ model: 'ed', link: 'linear', pars: [0.9, -0.1, 0.35, 0.05],
                           maxT: T, maxTries: 10000, rng });
  const ws = { nDesc: [], inherit: [], ed: [] };
  const ed = sim.edNow(r.rows, T, ws);
  r.rows.forEach((x, i) => lines.push([k, i + 1, x.parent + 1, x.birth, x.end < 0 ? 1 : 0,
                                       Number.isNaN(ed[i]) ? 'NA' : ed[i], T].join(',')));
}
fs.writeFileSync(outFile, lines.join('\n') + '\n');
