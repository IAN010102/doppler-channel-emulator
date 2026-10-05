'use strict';
/**
 * T1  Regression: at v = 0 and tau = 0 the unified model should agree with the legacy model on the mean SINR,
 *     within Monte-Carlo error. Criterion (fixed in advance): |mean_unified - mean_legacy| <= 3 * SE_diff, evaluated on the
 *     dB scale (what the UI shows). The linear-scale means are printed as well.
 *
 *   T1   default model parameters (K = 20 dB, jammer with Rician spread in legacy, single-path jammer in unified)
 *   T1b  K -> infinity (kDb = 400): both models then reduce to the same LoS-only scenario - a diagnostic that isolates
 *        the differences between the models (it is NOT a substitute for T1).
 */
const Core = require('../core.js');
const U = require('./_util.js');

const ALGOS = ['FOURIER', 'MMSE', 'SMI', 'DL', 'BEAMSPACE'];
const BASE = { N: 8, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, v: 0, latMs: 0, calDeg: 0, sll: 35, gammaDL: 0.01, taper: 'NONE', mod: 'QAM16' };

function sample(model, algo, over, trials, seed) {
    Core.setSeed(seed);
    const sys = Core.createSys(); sys.rollCal();
    Object.assign(sys, BASE, over, { model, algo, freshRealization: true });
    const db = [];
    for (let t = 0; t < trials; t++) {
        sys.snaps = []; sys.snapKey = '';
        sys.computeMath();
        if (Number.isFinite(sys.sinrDb)) db.push(sys.sinrDb);
    }
    return db;
}

function runCase(label, over, trials, seed) {
    const checks = [];
    console.log(`\n${label}   (${trials} trials per cell, K = ${over.kDb} dB)`);
    console.log(U.pad('algorithm', 11), U.rpad('legacy dB', 15), U.rpad('unified dB', 15), U.rpad('diff dB', 9), U.rpad('SE_diff', 8), U.rpad('z', 7),
        U.rpad('legacy lin', 11), U.rpad('unified lin', 12), U.rpad('z_lin', 7), 'result');
    ALGOS.forEach((algo, k) => {
        const a = sample('legacy', algo, over, trials, seed + 101 * k);
        const b = sample('unified', algo, over, trials, seed + 101 * k + 50);
        const ma = U.mean(a), mb = U.mean(b), sa = U.se(a), sb = U.se(b);
        const diff = mb - ma, sed = Math.hypot(sa, sb), z = sed > 0 ? diff / sed : (diff === 0 ? 0 : Infinity);
        const la = a.map(U.lin10), lb = b.map(U.lin10);
        const zl = (U.mean(lb) - U.mean(la)) / Math.hypot(U.se(la), U.se(lb));
        const pass = Math.abs(diff) <= Math.max(3 * sed, 1e-6);       // 1e-6 dB floor: FOURIER is deterministic, SE = 0
        console.log(U.pad(algo, 11), U.rpad(`${U.f(ma, 2)} ± ${U.f(sa, 2)}`, 15), U.rpad(`${U.f(mb, 2)} ± ${U.f(sb, 2)}`, 15), U.rpad(U.f(diff, 2), 9), U.rpad(U.f(sed, 3), 8), U.rpad(U.f(z, 1), 7),
            U.rpad(U.f(U.mean(la), 2), 11), U.rpad(U.f(U.mean(lb), 2), 12), U.rpad(U.f(zl, 1), 7), pass ? 'PASS' : 'FAIL');
        checks.push(U.check(`${label.split(' ')[0]} ${algo}: mean SINR dB, unified - legacy`, `${U.f(diff, 3)} dB`, `3*SE = ${U.f(3 * sed, 3)} dB`, pass,
            `z = ${U.f(z, 1)}; legacy ${U.f(ma, 2)} ± ${U.f(sa, 2)}, unified ${U.f(mb, 2)} ± ${U.f(sb, 2)}`));
    });
    return checks;
}

module.exports = {
    id: 'T1', title: 'regression: unified vs legacy at v = 0, tau = 0',
    async run({ trials = 2000, seed = 20260101 } = {}) {
        const checks = [];
        checks.push(...runCase('T1  default K', { kDb: 20 }, trials, seed));
        checks.push(...runCase('T1b K -> inf (diagnostic)', { kDb: 400 }, trials, seed + 7));
        return { id: this.id, title: this.title, checks };
    }
};
