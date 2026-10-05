'use strict';
/**
 * T1  Regression of the unified model against the legacy model at v = 0, tau = 0.
 *
 *   T1a  (check)  K -> infinity: both models reduce to the same LoS-only scenario. Mean SINR of the five algorithms must agree:
 *        |mean_unified - mean_legacy| <= 3 * SE_diff on the dB scale (what the UI shows), with a 1e-6 dB floor because FOURIER is deterministic
 *        (SE = 0). Training data: the default (withSignal).
 *   T1b  (informational, no PASS/FAIL)  K = 20 dB. Five algorithms x {withSignal, signalFree} x {legacy, unified}. Decomposition of the gap
 *        gap_ws = U_ws - L_ws  and  gap_sf = U_sf - L_sf:
 *           gap_ws - gap_sf   = the part of the gap that disappears when the target is removed from the training data
 *                               (= the part explained by "training contains the signal"),
 *           gap_sf            = what remains (static fading vs per-snapshot fading, single-path jammer, window handling ...).
 *        The model-internal effect of the training data, U_sf - U_ws and L_sf - L_ws, is listed too.
 *        The DL loading is the same absolute value in both models (gamma = 0.01 = gamma_rel 0 dB at SNR = 20 dB).
 */
const Core = require('../core.js');
const U = require('./_util.js');

const ALGOS = ['FOURIER', 'MMSE', 'SMI', 'DL', 'BEAMSPACE'];
const BASE = { N: 8, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, v: 0, latMs: 0, calDeg: 0, sll: 35, gammaDL: 0.01, gammaRelDb: 0, taper: 'NONE', mod: 'QAM16' };

function sample(model, algo, over, trials, seed) {
    Core.setSeed(seed);
    const sys = Core.createSys(); sys.rollCal();
    Object.assign(sys, BASE, over, { model, algo, freshRealization: true });
    const db = []; db.opt = [];
    for (let t = 0; t < trials; t++) {
        sys.snaps = []; sys.snapKey = '';
        sys.computeMath();
        if (Number.isFinite(sys.sinrDb)) { db.push(sys.sinrDb); db.opt.push(sys.sinrOptDb); }
    }
    return db;
}
const fmt = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`;

module.exports = {
    id: 'T1', title: 'regression: unified vs legacy at v = 0, tau = 0 (T1a check; T1b informational decomposition)',
    async run({ trials = 2000, seed = 20260101 } = {}) {
        const checks = [], info = [];

        // ---------------------------------------------------------------- T1a
        const over = { kDb: 400 };
        console.log(`T1a  K -> inf (kDb = 400), v = 0, tau = 0, withSignal, ${trials} trials per cell`);
        console.log(U.pad('algorithm', 11), U.rpad('legacy dB', 15), U.rpad('unified dB', 15), U.rpad('diff dB', 9), U.rpad('SE_diff', 8), U.rpad('z', 7), 'result');
        ALGOS.forEach((algo, k) => {
            const a = sample('legacy', algo, over, trials, seed + 7 + 101 * k);
            const b = sample('unified', algo, over, trials, seed + 7 + 101 * k + 50);
            const diff = U.mean(b) - U.mean(a), sed = Math.hypot(U.se(a), U.se(b)), z = sed > 0 ? diff / sed : (diff === 0 ? 0 : Infinity);
            const pass = Math.abs(diff) <= Math.max(3 * sed, 1e-6);
            console.log(U.pad(algo, 11), U.rpad(fmt(a), 15), U.rpad(fmt(b), 15), U.rpad(U.f(diff, 3), 9), U.rpad(U.f(sed, 3), 8), U.rpad(U.f(z, 1), 7), pass ? 'PASS' : 'FAIL');
            checks.push(U.check(`T1a ${algo}: mean SINR dB, unified - legacy (K -> inf)`, `${U.f(diff, 3)} dB`, `3*SE = ${U.f(3 * sed, 3)} dB`, pass, `z = ${U.f(z, 1)}`));
        });

        // ---------------------------------------------------------------- T1b
        console.log(`\nT1b  (informational) K = 20 dB, v = 0, tau = 0, ${trials} trials per cell; mean SINR dB ± SE`);
        console.log(U.pad('algorithm', 11), U.rpad('legacy ws', 14), U.rpad('unified ws', 14), U.rpad('legacy sf', 14), U.rpad('unified sf', 14),
            U.rpad('gap_ws', 8), U.rpad('gap_sf', 8), U.rpad('ws - sf', 8), U.rpad('SE', 6), U.rpad('L: sf-ws', 9), U.rpad('U: sf-ws', 9), U.rpad('opt L', 7), U.rpad('opt U', 7), U.rpad('Lws-opt', 8), U.rpad('Uws-opt', 8), U.rpad('Lsf-opt', 8), U.rpad('Usf-opt', 8));
        ALGOS.forEach((algo, k) => {
            const s0 = seed + 9000 + 211 * k;
            const Lw = sample('legacy', algo, { kDb: 20, trainMode: 'withSignal' }, trials, s0), Uw = sample('unified', algo, { kDb: 20, trainMode: 'withSignal' }, trials, s0 + 50);
            const Ls = sample('legacy', algo, { kDb: 20, trainMode: 'signalFree' }, trials, s0), Us = sample('unified', algo, { kDb: 20, trainMode: 'signalFree' }, trials, s0 + 50);
            const gws = U.mean(Uw) - U.mean(Lw), gsf = U.mean(Us) - U.mean(Ls);
            const se = Math.sqrt(U.se(Lw) ** 2 + U.se(Uw) ** 2 + U.se(Ls) ** 2 + U.se(Us) ** 2);
            console.log(U.pad(algo, 11), U.rpad(fmt(Lw), 14), U.rpad(fmt(Uw), 14), U.rpad(fmt(Ls), 14), U.rpad(fmt(Us), 14),
                U.rpad(U.f(gws, 2), 8), U.rpad(U.f(gsf, 2), 8), U.rpad(U.f(gws - gsf, 2), 8), U.rpad(U.f(se, 2), 6),
                U.rpad(U.f(U.mean(Ls) - U.mean(Lw), 2), 9), U.rpad(U.f(U.mean(Us) - U.mean(Uw), 2), 9),
                U.rpad(U.f(U.mean(Lw.opt), 2), 7), U.rpad(U.f(U.mean(Uw.opt), 2), 7), U.rpad(U.f(U.mean(Lw) - U.mean(Lw.opt), 2), 8), U.rpad(U.f(U.mean(Uw) - U.mean(Uw.opt), 2), 8),
                U.rpad(U.f(U.mean(Ls) - U.mean(Ls.opt), 2), 8), U.rpad(U.f(U.mean(Us) - U.mean(Us.opt), 2), 8));
            info.push({ name: `T1b ${algo}`, value: `gap_ws ${U.f(gws, 2)}, gap_sf ${U.f(gsf, 2)}, explained by training ${U.f(gws - gsf, 2)} dB`,
                note: `SINR_opt legacy ${U.f(U.mean(Lw.opt), 2)} / unified ${U.f(U.mean(Uw.opt), 2)} dB; gap to SINR_opt (legacy ws, unified ws, legacy sf, unified sf) = ${U.f(U.mean(Lw) - U.mean(Lw.opt), 2)} / ${U.f(U.mean(Uw) - U.mean(Uw.opt), 2)} / ${U.f(U.mean(Ls) - U.mean(Ls.opt), 2)} / ${U.f(U.mean(Us) - U.mean(Us.opt), 2)} dB; legacy ws ${fmt(Lw)}, unified ws ${fmt(Uw)}, legacy sf ${fmt(Ls)}, unified sf ${fmt(Us)}; SE of a combination ~ ${U.f(se, 2)} dB` });
        });
        console.log('     gap_* = unified - legacy; ws = withSignal (MPDR), sf = signalFree (MVDR); "ws - sf" = part of the gap explained by training data that contain the signal');
        return { id: this.id, title: this.title, checks, info };
    }
};
