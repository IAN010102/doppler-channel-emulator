'use strict';
/**
 * T3  Aging.
 *
 *   T3a  (check) signalFree training, K -> infinity, target ON the track (theta_1 = 0: no angle drift, a single path), jammer at theta_2 = -30 deg
 *        drifting with d_min. So the only aging source is the jammer's angle drift. Trials are paired (same seed for every tau).
 *        For SMI, DL, BEAMSPACE and d_min in {5, 30} m (small): mean SINR_inst must decrease monotonically as tau grows
 *        (tau = 0, 1, 2, 4, 6, 8, 10 ms). A step violates this only if the paired increase exceeds one standard error of that
 *        paired difference (mean(d) > SE(d), d = SINR(tau_k+1) - SINR(tau_k)). d_min = 500 m is listed for information.
 *   T3b  (check) signalFree, tau = 0, K = 20 dB, no angle drift (d_min -> infinity): v = 0 vs v = 300 km/h, paired trials.
 *        Pure Doppler must not change the (target-free) covariance statistics: |mean diff| <= 3 SE.
 *   T3c  (informational) withSignal, diffuse-only channel (K -> 0), no drift: tau sweep at v = 300 for SMI, DL, BEAMSPACE (signalFree next to it).
 *        Records whether the curve is non-monotone (expected: self-nulling of MPDR relaxes as the channel decorrelates from the training window).
 *   T3d  (check) FOURIER does not adapt: under the T3a scenario its SINR_inst must not depend on tau (1e-9 dB).
 */
const Core = require('../core.js');
const U = require('./_util.js');

const TAUS_MS = [0, 1, 2, 4, 6, 8, 10];
const BASE = { model: 'unified', freshRealization: true, N: 8, L: 100, snr: 20, sir: -10, calDeg: 0, taper: 'NONE', gammaDL: 0.01, gammaRelDb: 10, mod: 'QAM16' };
const A_SCEN = { kDb: 400, aoaT: 0, aoaJ: -30, v: 300, trainMode: 'signalFree' };       // jammer-only aging

function values(sys, over, trials, seed0) {
    const out = []; out.opt = [];
    for (let t = 0; t < trials; t++) {
        Core.setSeed(seed0 + t);
        Object.assign(sys, BASE, over); sys.snaps = [];
        sys.computeMath();
        out.push(sys.sinrDb); out.opt.push(sys.sinrOptDb);
    }
    return out;
}
const bar = (x, lo, hi, w = 30) => '#'.repeat(Math.max(0, Math.round((x - lo) / (hi - lo) * w)));

module.exports = {
    id: 'T3', title: 'aging: jammer-only drift (signalFree) monotone; Doppler-only invariance; MPDR tau sweep (info); FOURIER control',
    async run({ trials = 2000, seed = 424242 } = {}) {
        const checks = [], info = [];
        const sys = Core.createSys(); Core.setSeed(1); sys.rollCal();

        // ---------------------------------------------------------------- T3a
        console.log(`T3a  signalFree, K -> inf, theta1 = 0 (no drift), theta2 = -30, v = 300, ${trials} paired trials; mean SINR_inst [dB] vs tau [ms]`);
        console.log(U.pad('d_min', 7), U.pad('algorithm', 10), TAUS_MS.map(x => U.rpad(x, 8)).join(''), '   worst step (mean d / SE_d)');
        for (const dmin of [5, 30, 500]) for (const algo of ['SMI', 'DL', 'BEAMSPACE']) {
            const cols = TAUS_MS.map(tau => values(sys, Object.assign({}, A_SCEN, { d_min: dmin, algo, latMs: tau }), trials, seed));
            const means = cols.map(U.mean);
            let worst = -Infinity, ok = true, worstTxt = '';
            for (let k = 0; k + 1 < TAUS_MS.length; k++) {
                const d = cols[k + 1].map((x, i) => x - cols[k][i]), m = U.mean(d), s = U.se(d);
                const viol = m > s;                                    // paired increase larger than one standard error
                if (viol) ok = false;
                const r = s > 0 ? m / s : (m > 0 ? Infinity : -Infinity);
                if (r > worst) { worst = r; worstTxt = `${TAUS_MS[k]}->${TAUS_MS[k + 1]} ms: ${U.f(m, 4)} dB / ${U.f(s, 4)}`; }
            }
            console.log(U.pad('', 7), U.pad('  - SINR_opt', 10), cols.map(c => U.rpad(U.f(U.mean(c) - U.mean(c.opt), 2), 8)).join(''), '   (mean SINR minus mean SINR_opt [dB]; SINR_opt at tau = 0: ' + U.f(U.mean(cols[0].opt), 2) + ' dB)');
            const small = dmin < 100;
            console.log(U.pad(dmin, 7), U.pad(algo, 10), means.map(x => U.rpad(U.f(x, 2), 8)).join(''), '   ', worstTxt, small ? (ok ? 'PASS' : 'FAIL') : (ok ? '(info: monotone)' : '(info: not monotone)'));
            if (small) checks.push(U.check(`T3a d_min=${dmin} m ${algo}: mean SINR_inst non-increasing in tau`, `drop 0 -> 10 ms: ${U.f(means[0] - means[means.length - 1], 2)} dB`, 'each step: mean(d) <= SE(d)', ok, `worst step ${worstTxt}`));
            else info.push({ name: `T3a d_min=${dmin} m ${algo}`, value: `drop 0 -> 10 ms: ${U.f(means[0] - means[means.length - 1], 3)} dB`, note: ok ? 'monotone by the same rule' : `violates the rule: ${worstTxt}` });
        }

        // ---------------------------------------------------------------- T3b
        for (const algo of ['SMI', 'DL', 'BEAMSPACE']) {
            const o = { kDb: 20, aoaT: 60, aoaJ: -40, trainMode: 'signalFree', d_min: 1e30, algo, latMs: 0 };
            const a = values(sys, Object.assign({ v: 0 }, o), trials, seed), b = values(sys, Object.assign({ v: 300 }, o), trials, seed);
            const d = b.map((x, i) => x - a[i]), m = U.mean(d), s = U.se(d), pass = Math.abs(m) <= 3 * s;
            console.log(`T3b  ${U.pad(algo, 10)} signalFree, tau = 0, no drift: v = 0 -> ${U.f(U.mean(a), 3)} dB, v = 300 -> ${U.f(U.mean(b), 3)} dB, paired diff ${U.f(m, 4)} ± ${U.f(s, 4)} dB  ${pass ? 'PASS' : 'FAIL'}`);
            checks.push(U.check(`T3b ${algo}: signalFree, tau = 0: mean SINR_inst at v = 300 equals v = 0`, `${U.f(m, 4)} dB`, `3*SE = ${U.f(3 * s, 4)} dB`, pass));
        }

        // ---------------------------------------------------------------- T3c (informational)
        console.log(`\nT3c  (informational) diffuse-only (K -> 0), v = 300, theta1 = 60, no drift, ${trials} paired trials; mean SINR_inst [dB] vs tau [ms]`);
        const C_SCEN = { kDb: -300, aoaT: 60, aoaJ: -40, v: 300, d_min: 1e30 }, TC = [0, 0.25, 0.5, 1, 2, 4, 8];
        for (const algo of ['SMI', 'DL', 'BEAMSPACE']) for (const tm of ['withSignal', 'signalFree']) {
            const cols = TC.map(tau => values(sys, Object.assign({}, C_SCEN, { algo, trainMode: tm, latMs: tau }), trials, seed)), means = cols.map(U.mean);
            let rises = 0, maxRise = 0;
            for (let k = 0; k + 1 < TC.length; k++) { const d = cols[k + 1].map((x, i) => x - cols[k][i]), m = U.mean(d), s = U.se(d); if (m > s) { rises++; maxRise = Math.max(maxRise, m); } }
            const lo = Math.min(...means), hi = Math.max(...means);
            console.log(U.pad(algo, 10), U.pad(tm, 11), means.map(x => U.rpad(U.f(x, 2), 8)).join(''), '  ', rises ? `NON-MONOTONE (${rises} rising steps, max +${U.f(maxRise, 2)} dB)` : 'monotone');
            console.log(U.pad('', 22), means.map(x => U.rpad(bar(x, lo - 0.3, hi + 0.3, 7), 8)).join(''));
            info.push({ name: `T3c ${algo} ${tm}`, value: means.map(x => U.f(x, 2)).join(' / ') + ' dB', note: rises ? `non-monotone: ${rises} rising step(s), max +${U.f(maxRise, 2)} dB` : 'monotone' });
        }
        console.log('     (tau = ' + TC.join(' / ') + ' ms.)  Expected for withSignal: the MPDR self-nulling is strongest when the channel at t_app equals the training channel and relaxes as it decorrelates.');

        // ---------------------------------------------------------------- T3d
        let maxDev = 0;
        const ref = values(sys, Object.assign({}, A_SCEN, { d_min: 5, algo: 'FOURIER', latMs: 0 }), 200, seed);
        for (const tau of TAUS_MS) { const x = values(sys, Object.assign({}, A_SCEN, { d_min: 5, algo: 'FOURIER', latMs: tau }), 200, seed); x.forEach((val, i) => { maxDev = Math.max(maxDev, Math.abs(val - ref[i])); }); }
        console.log(`T3d  FOURIER, T3a scenario, d_min = 5: max |SINR_inst(tau) - SINR_inst(0)| = ${U.e(maxDev)} dB`);
        checks.push(U.check('T3d FOURIER: SINR_inst independent of tau (no adaptation)', U.e(maxDev) + ' dB', '1e-9 dB', maxDev < 1e-9));
        return { id: this.id, title: this.title, checks, info };
    }
};
