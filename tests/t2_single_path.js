'use strict';
/**
 * T2  Single-path invariance. K -> infinity (only the LoS path survives), no angle drift (d_min -> infinity), tau = 0:
 *     changing v must not change R_hat and SINR_inst (only N_ICI/S). Required error < 1e-9.
 *
 *   T2a  literal specification: same seed, same everything, v = 0 / 100 / 300 / 500; compare R_hat and SINR_inst.
 *   T2b  noise-free and jammer-free: R_hat = (1/L) sum |s_n|^2 h_n h_n^H.  Isolates the model from the finite-sample
 *        cross terms (see the note on T2a).
 *   T2c  statistical invariance with noise + jammer: mean SINR_inst and mean R_hat entries for v = 0 vs v = 300 over many trials.
 *   T2d  N_ICI/S of a single path equals 1 - sinc^2(eps) exactly.
 */
const Core = require('../core.js');
const U = require('./_util.js');

const BASE = { model: 'unified', freshRealization: true, kDb: 400, d_min: 1e30, latMs: 0, N: 8, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, calDeg: 0,
    taper: 'NONE', algo: 'SMI', gammaDL: 0.01, mod: 'QAM16' };
const VS = [0, 100, 300, 500];

function oneDraw(sys, v, seed) { Core.setSeed(seed); sys.v = v; sys.snaps = []; sys.computeMath(); return { R: U.cloneR(sys.R_raw), sinr: Math.pow(10, sys.sinrDb / 10), nu: sys.nuICI, eps: sys.eps0 }; }

module.exports = {
    id: 'T2', title: 'single-path invariance (K -> inf, no drift): v must not change R_hat / SINR_inst',
    async run({ trials = 4000, seed = 777 } = {}) {
        const checks = [];
        const sys = Core.createSys(); Core.setSeed(1); sys.rollCal();

        // ---------------------------------------------------------------- T2a: literal
        Object.assign(sys, BASE);
        const d = VS.map(v => oneDraw(sys, v, seed));
        let maxR = 0, maxS = 0;
        for (let k = 1; k < VS.length; k++) { maxR = Math.max(maxR, U.maxAbsDiffR(d[k].R, d[0].R)); maxS = Math.max(maxS, Math.abs(d[k].sinr - d[0].sinr) / d[0].sinr); }
        console.log(`T2a  same seed, v = ${VS.join(' / ')} km/h:  max |dR_hat| = ${U.e(maxR)},  max relative |dSINR_inst| = ${U.e(maxS)}`);
        checks.push(U.check('T2a R_hat unchanged by v (noise + jammer present)', U.e(maxR), '1e-9', maxR < 1e-9,
            'FAILS BY CONSTRUCTION: the sample covariance contains the finite-sample cross terms (1/L) sum s_n n_n^H e^{jw t_n}, (1/L) sum s_n j_n^* e^{j(w1-w2) t_n}; ' +
            'the Doppler phases rotate the random symbols, so for a fixed random draw R_hat changes with v even though its expectation does not (see T2b-T2d)'));
        checks.push(U.check('T2a SINR_inst unchanged by v (noise + jammer present)', U.e(maxS), '1e-9', maxS < 1e-9, 'same cause: the weights are computed from R_hat'));

        // ---------------------------------------------------------------- T2b: noise-free, jammer-free
        Object.assign(sys, BASE, { snr: 400, sir: 400, algo: 'FOURIER' });
        const d2 = VS.map(v => oneDraw(sys, v, seed));
        let maxR2 = 0, maxS2 = 0;
        for (let k = 1; k < VS.length; k++) { maxR2 = Math.max(maxR2, U.maxAbsDiffR(d2[k].R, d2[0].R)); maxS2 = Math.max(maxS2, Math.abs(d2[k].sinr - d2[0].sinr) / d2[0].sinr); }
        console.log(`T2b  noise-free, jammer-free: max |dR_hat| = ${U.e(maxR2)},  max relative |dSINR_inst| (FOURIER) = ${U.e(maxS2)}`);
        checks.push(U.check('T2b R_hat unchanged by v (noise-free, jammer-free)', U.e(maxR2), '1e-9', maxR2 < 1e-9));
        checks.push(U.check('T2b SINR_inst unchanged by v (FOURIER, noise-free, jammer-free)', U.e(maxS2), '1e-9', maxS2 < 1e-9));

        // ---------------------------------------------------------------- T2c: statistical invariance
        Object.assign(sys, BASE);
        const stat = (v, sd) => {
            Core.setSeed(sd); sys.v = v;
            const N = sys.N, sumR = Array.from({ length: N }, () => Array.from({ length: N }, () => ({ r: [], i: [] }))), sinr = [];
            for (let t = 0; t < trials; t++) {
                sys.snaps = []; sys.computeMath(); sinr.push(sys.sinrDb);
                for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) { sumR[i][j].r.push(sys.R_raw[i][j].r); sumR[i][j].i.push(sys.R_raw[i][j].i); }
            }
            return { sinr, sumR };
        };
        const A = stat(0, seed + 1), B = stat(300, seed + 2);
        const zS = (U.mean(B.sinr) - U.mean(A.sinr)) / Math.hypot(U.se(A.sinr), U.se(B.sinr));
        let zMax = 0;
        for (let i = 0; i < sys.N; i++) for (let j = 0; j < sys.N; j++) for (const k of ['r', 'i']) {
            const a = A.sumR[i][j][k], b = B.sumR[i][j][k], s = Math.hypot(U.se(a), U.se(b));
            if (s > 0) zMax = Math.max(zMax, Math.abs((U.mean(b) - U.mean(a)) / s));
        }
        console.log(`T2c  ${trials} trials each, v = 0 vs 300:  mean SINR ${U.f(U.mean(A.sinr), 2)} vs ${U.f(U.mean(B.sinr), 2)} dB (z = ${U.f(zS, 2)}),  max |z| over the ${2 * sys.N * sys.N} R_hat entries = ${U.f(zMax, 2)}`);
        checks.push(U.check('T2c mean SINR_inst, v = 300 vs v = 0 (noise + jammer)', `z = ${U.f(zS, 2)}`, '|z| <= 3', Math.abs(zS) <= 3));
        checks.push(U.check('T2c mean R_hat entries, v = 300 vs v = 0', `max |z| = ${U.f(zMax, 2)}`, `|z| <= 4 (${2 * sys.N * sys.N} comparisons)`, zMax <= 4));

        // ---------------------------------------------------------------- T2d: N_ICI/S of one path
        let maxN = 0;
        Object.assign(sys, BASE);
        for (const v of VS) { const r = oneDraw(sys, v, seed); const x = r.eps; const s = x === 0 ? 1 : Math.sin(Math.PI * x) / (Math.PI * x); maxN = Math.max(maxN, Math.abs(r.nu - (1 - s * s))); }
        console.log(`T2d  single path: max | N_ICI/S - (1 - sinc^2(eps)) | = ${U.e(maxN)}`);
        checks.push(U.check('T2d N_ICI/S = 1 - sinc^2(eps) for a single path', U.e(maxN), '1e-9', maxN < 1e-9));
        return { id: this.id, title: this.title, checks };
    }
};
