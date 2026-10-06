'use strict';
// Shared by make_golden.js and t14_golden.js: the fixed cases of the snapshot regression (new defaults: unified model, signalFree training, 16-QAM).
const Core = require('../core.js');
const ALGOS = ['FOURIER', 'MMSE', 'MMSEP', 'SMI', 'DL', 'BEAMSPACE'];
const PARAMS = { N: 8, L: 100, aoaT: 10, aoaJ: 40, snr: 20, sir: -10, v: 100, kDb: 20, latMs: 1, calDeg: 5, taper: 'NONE', mod: 'QAM16' };   // model, trainMode, d_min, fc, df: the defaults
const SEEDS = [11, 22, 33, 44];
function compute(extra = {}) {
    const out = { params: PARAMS, defaults: {}, cases: {} };
    for (const tm of [null, 'withSignal']) {
        for (const algo of ALGOS) {
            const rows = [];
            for (const seed of SEEDS) {
                Core.setSeed(seed);
                const s = Core.createSys(); s.rollCal();
                Object.assign(s, PARAMS, { algo, freshRealization: true }, extra); if (tm) s.trainMode = tm;
                if (!out.defaults.model) out.defaults = { model: s.model, trainMode: Core.CONFIG.trainMode, smiSingular: s.smiSingular, gammaRelDb: s.gammaRelDb, d_min: s.d_min, fc: s.fc, scs: s.scs };
                s.snaps = []; s.computeMath();
                rows.push({ seed, sinrDb: s.sinrDb, evm: s.evm, sinrOptDb: s.sinrOptDb, ser: s.ser });
            }
            out.cases[`${tm || 'default(signalFree)'}/${algo}`] = rows;
        }
    }
    return out;
}
module.exports = { compute };
