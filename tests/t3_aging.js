'use strict';
/**
 * T3  Aging monotonicity. Diffuse paths only (K -> 0), fixed L, no angle drift (R_min -> infinity). As tau grows, the trial
 *     mean of SINR_inst of the *adaptive* weights should decrease monotonically.  Trials are paired (same seed for every tau,
 *     so the realisation, symbols and noise are identical and only t_app moves).
 *
 *   T3a  monotone decrease: for consecutive taus,  mean(d) <= 2 * SE(d)  (d = paired difference, per trial). A strictly
 *        monotone sample mean is not demanded where the true curve has flattened (plateau after decorrelation).
 *   T3b  v = 0: SINR_inst must not depend on tau at all (static channel)  -> exact.
 *   T3c  literal "tau = 0 equals v = 0": mean SINR_inst at (v = 300, tau = 0) vs (v = 0, tau = 0), paired, |diff| <= 3 SE.
 *   T3d  control: FOURIER does not adapt, so no aging is expected (informational, not a check).
 */
const Core = require('../core.js');
const U = require('./_util.js');

const TAUS_MS = [0, 0.25, 0.5, 1, 2, 4, 8];
const BASE = { model: 'unified', freshRealization: true, kDb: -300, R_min: 1e30, N: 8, L: 100, aoaT: 60, aoaJ: -40, snr: 20, sir: -10, calDeg: 0,
    taper: 'NONE', gammaDL: 0.01, mod: 'QAM16' };

function trialValues(sys, algo, v, tauMs, trials, seed0) {
    const out = [];
    for (let t = 0; t < trials; t++) {
        Core.setSeed(seed0 + t);
        Object.assign(sys, BASE, { algo, v, latMs: tauMs }); sys.snaps = [];
        sys.computeMath();
        out.push(sys.sinrDb);
    }
    return out;
}

module.exports = {
    id: 'T3', title: 'aging monotonicity (diffuse-only), v = 0 invariance, tau = 0 vs v = 0',
    async run({ trials = 1000, seed = 424242 } = {}) {
        const checks = [];
        const sys = Core.createSys(); Core.setSeed(1); sys.rollCal();

        // ---------------------------------------------------------------- T3a: monotone decrease
        console.log(`T3a  v = 300 km/h, theta_1 = 60 deg, K -> 0, L = ${BASE.L}, ${trials} paired trials; mean SINR_inst [dB] vs tau [ms]`);
        console.log(U.pad('algorithm', 11), TAUS_MS.map(x => U.rpad(x, 8)).join(''), '   worst step (mean d / SE_d)');
        for (const algo of ['SMI', 'DL', 'BEAMSPACE']) {
            const cols = TAUS_MS.map(tau => trialValues(sys, algo, 300, tau, trials, seed));
            const means = cols.map(U.mean);
            let worst = -Infinity, ok = true, worstTxt = '';
            for (let k = 0; k + 1 < TAUS_MS.length; k++) {
                const d = cols[k + 1].map((x, i) => x - cols[k][i]);
                const m = U.mean(d), s = U.se(d), r = s > 0 ? m / s : (m > 0 ? Infinity : -Infinity);
                if (m > 2 * s) ok = false;
                if (r > worst) { worst = r; worstTxt = `${TAUS_MS[k]}->${TAUS_MS[k + 1]} ms: ${U.f(m, 3)} dB / ${U.f(s, 3)}`; }
            }
            console.log(U.pad(algo, 11), means.map(x => U.rpad(U.f(x, 2), 8)).join(''), '   ', worstTxt, ok ? 'PASS' : 'FAIL');
            checks.push(U.check(`T3a ${algo}: mean SINR_inst non-increasing in tau`, `drop 0 -> 8 ms: ${U.f(means[0] - means[means.length - 1], 2)} dB`, 'each step: mean(d) <= 2 SE(d)', ok, `worst step ${worstTxt}`));
        }

        // ---------------------------------------------------------------- T3b: v = 0 -> no tau dependence
        let maxDev = 0;
        for (const algo of ['SMI', 'DL', 'BEAMSPACE']) {
            const ref = trialValues(sys, algo, 0, 0, 50, seed);
            for (const tau of [1, 4, 8]) { const x = trialValues(sys, algo, 0, tau, 50, seed); x.forEach((val, i) => { maxDev = Math.max(maxDev, Math.abs(val - ref[i])); }); }
        }
        console.log(`T3b  v = 0: max |SINR_inst(tau) - SINR_inst(0)| over algorithms, tau = 1/4/8 ms, 50 trials = ${U.e(maxDev)} dB`);
        checks.push(U.check('T3b v = 0: SINR_inst independent of tau', U.e(maxDev) + ' dB', '1e-9 dB', maxDev < 1e-9));

        // ---------------------------------------------------------------- T3c: literal "tau = 0 equals v = 0"
        for (const algo of ['SMI', 'DL', 'BEAMSPACE']) {
            const a = trialValues(sys, algo, 0, 0, trials, seed), b = trialValues(sys, algo, 300, 0, trials, seed);
            const d = b.map((x, i) => x - a[i]), m = U.mean(d), s = U.se(d);
            const pass = Math.abs(m) <= 3 * s;
            console.log(`T3c  ${U.pad(algo, 10)} mean SINR_inst: v = 0 -> ${U.f(U.mean(a), 2)} dB,  v = 300 (tau = 0) -> ${U.f(U.mean(b), 2)} dB,  paired diff ${U.f(m, 3)} ± ${U.f(s, 3)} dB  ${pass ? 'PASS' : 'FAIL'}`);
            checks.push(U.check(`T3c ${algo}: mean SINR_inst at (v=300, tau=0) equals (v=0, tau=0)`, `${U.f(m, 3)} dB`, `3*SE = ${U.f(3 * s, 3)} dB`, pass,
                'Doppler makes the channel time-varying inside the training window, so R_hat (and the weights) differ from the static case even at tau = 0; the specification only guarantees equality of the aging term'));
        }

        // ---------------------------------------------------------------- T3d: control
        const f0 = trialValues(sys, 'FOURIER', 300, 0, trials, seed), f8 = trialValues(sys, 'FOURIER', 300, 8, trials, seed);
        console.log(`T3d  control, FOURIER (non-adaptive): mean SINR_inst ${U.f(U.mean(f0), 2)} dB at tau = 0, ${U.f(U.mean(f8), 2)} dB at tau = 8 ms  -> no aging expected`);
        return { id: this.id, title: this.title, checks };
    }
};
