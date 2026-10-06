'use strict';
/**
 * T14  Snapshot regression of the new defaults (unified model, signalFree training): fixed seeds and parameters, the six algorithms (FOURIER, MMSE-M, MMSE-P, SMI, DL,
 *      BEAMSPACE), SINR, EVM, SINR_opt and SER of 4 seeds, for the default training mode and for withSignal. The numbers are stored in tests/golden/default_snapshot.json.
 *      A change that moves them must update the file on purpose (node tests/make_golden.js) and say why in the commit.  Tolerance: 1e-9 (relative, 1e-12 absolute floor).
 */
const fs = require('fs'), path = require('path');
const { compute } = require('./golden_cases.js');
const U = require('./_util.js');
module.exports = {
    id: 'T14', title: 'golden snapshot of the default settings (unified, signalFree): SINR / EVM / SINR_opt / SER',
    async run() {
        const gold = JSON.parse(fs.readFileSync(path.join(__dirname, 'golden', 'default_snapshot.json'), 'utf8')), now = compute(), checks = [];
        const dflt = JSON.stringify(now.defaults) === JSON.stringify(gold.defaults);
        console.log('defaults:', JSON.stringify(now.defaults));
        checks.push(U.check('T14 default settings (model, trainMode, smiSingular, gammaRelDb, d_min, fc, df)', JSON.stringify(now.defaults), 'equal to the golden file', dflt));
        let worstAll = 0;
        for (const key of Object.keys(gold.cases)) {
            let worst = 0;
            gold.cases[key].forEach((g, i) => { const n = now.cases[key][i]; for (const f of ['sinrDb', 'evm', 'sinrOptDb', 'ser']) worst = Math.max(worst, Math.abs(n[f] - g[f]) / Math.max(Math.abs(g[f]), 1e-12)); });
            worstAll = Math.max(worstAll, worst);
            checks.push(U.check(`T14 ${key}`, U.e(worst), '< 1e-9 (relative)', worst < 1e-9));
        }
        console.log(`max relative deviation over all ${Object.keys(gold.cases).length} cases: ${U.e(worstAll)}`);
        return { id: this.id, title: this.title, checks };
    }
};
