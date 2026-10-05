'use strict';
/**
 * T6  (informational, no PASS/FAIL)  Diagonal-loading sweep of DL, unified model: K = 20 dB, v = 0, tau = 0, SNR = 20 dB, L = 100 (and L = 12 for a
 *     small-sample reference), gamma = gamma_rel * sigma_n^2, gamma_rel = -10 ... +30 dB in 2 dB steps, withSignal and signalFree.
 *     Trials are paired (same seed for every gamma_rel and both modes). The whole curve is printed and the best gamma_rel is marked.
 *     Expected for withSignal: rising then falling or saturating (loading fights the self-nulling of MPDR; too much loading approaches the matched filter).
 */
const Core = require('../core.js');
const U = require('./_util.js');

const BASE = { model: 'unified', freshRealization: true, algo: 'DL', N: 8, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, v: 0, latMs: 0, kDb: 20, calDeg: 0, taper: 'NONE', d_min: 30 };
const GRID = []; for (let g = -10; g <= 30; g += 2) GRID.push(g);

module.exports = {
    id: 'T6', title: 'gamma_rel sweep of DL (informational): K = 20 dB, v = 0',
    async run({ trials = 1000, seed = 616161 } = {}) {
        const info = [], sys = Core.createSys(); Core.setSeed(1); sys.rollCal();
        for (const L of [100, 12]) for (const tm of ['withSignal', 'signalFree']) {
            const means = [], ses = [];
            for (const g of GRID) {
                const a = [];
                for (let t = 0; t < trials; t++) { Core.setSeed(seed + t); Object.assign(sys, BASE, { L, trainMode: tm, gammaRelDb: g }); sys.snaps = []; sys.computeMath(); a.push(sys.sinrDb); }
                means.push(U.mean(a)); ses.push(U.se(a));
            }
            let best = 0; means.forEach((m, i) => { if (m > means[best]) best = i; });
            const lo = Math.min(...means), hi = Math.max(...means);
            console.log(`\nT6  DL, ${tm}, L = ${L}, ${trials} paired trials   (SNR = 20 dB: sigma_n^2 = 0.01, gamma = 10^(gamma_rel/10) * 0.01)`);
            console.log(U.pad('gamma_rel dB', 13), U.pad('gamma', 9), U.rpad('SINR dB', 9), U.rpad('SE', 7), ' curve');
            GRID.forEach((g, i) => console.log(U.pad(g, 13), U.pad(U.e(Math.pow(10, g / 10) * 0.01, 1), 9), U.rpad(U.f(means[i], 2), 9), U.rpad(U.f(ses[i], 2), 7), ' ' + '#'.repeat(Math.max(0, Math.round((means[i] - lo) / Math.max(hi - lo, 1e-9) * 40))) + (i === best ? '  <= best' : '')));
            const shape = best === 0 ? 'monotone decreasing' : best === GRID.length - 1 ? 'increasing up to the end of the range (saturating?)' : 'rises, peaks, then falls';
            console.log(`     best gamma_rel = ${GRID[best]} dB (SINR ${U.f(means[best], 2)} ± ${U.f(ses[best], 2)} dB); at 0 dB ${U.f(means[GRID.indexOf(0)], 2)}, at +10 dB ${U.f(means[GRID.indexOf(10)], 2)}, at ${GRID[0]} dB ${U.f(means[0], 2)}, at +30 dB ${U.f(means[GRID.length - 1], 2)}; shape: ${shape}`);
            info.push({ name: `T6 DL ${tm} L=${L}`, value: `best gamma_rel ${GRID[best]} dB (${U.f(means[best], 2)} dB)`,
                note: `${shape}; 0 dB: ${U.f(means[GRID.indexOf(0)], 2)}, +10 dB: ${U.f(means[GRID.indexOf(10)], 2)}, -10 dB: ${U.f(means[0], 2)}, +30 dB: ${U.f(means[GRID.length - 1], 2)}` });
        }
        return { id: this.id, title: this.title, checks: [], info };
    }
};
