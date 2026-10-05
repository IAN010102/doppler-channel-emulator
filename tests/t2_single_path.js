'use strict';
/**
 * T2  Single-path invariance. K -> infinity (only the LoS path survives), no angle drift (d_min -> infinity), tau = 0, withSignal.
 *     The Doppler shift changes only N_ICI/S: the statistics of R_hat and of SINR_inst must not depend on v.
 *
 *   (T2a - same seed, identical R_hat for different v - is removed: impossible by construction, the finite-sample cross terms of the
 *    sample covariance rotate with the Doppler phases.)
 *   T2b  noise-free and jammer-free: R_hat = (1/L) sum |s_n|^2 h_n h_n^H exactly, so it does not depend on v. Required error < 1e-9.
 *   T2c  statistical invariance with noise + jammer, v = 0 vs v = 300, >= 4000 trials each: z-test on the mean SINR_inst (|z| < 3, one comparison) and on
 *        every non-trivial element of R_hat (real and imaginary parts of all N x N entries; the identically-zero imaginary diagonal is skipped).
 *        R_hat criterion: Bonferroni-corrected. With m simultaneous comparisons of (nearly independent) statistics the expected maximum of |z| grows with m
 *        (about sqrt(2 ln m)), so a fixed threshold of 3 is not a controlled test: the family-wise false-alarm rate would be 1 - (1 - 0.0027)^m (~ 53 % for m = 276).
 *        Instead: z_crit = Phi^-1(1 - alpha/(2m)), alpha = 0.01, two-sided, i.e. P(|Z| > z_crit) = alpha/m; pass if max |z| < z_crit. m, z_crit and the
 *        largest |z| are reported. (R_hat is Hermitian, so the entries are not all independent; the correction is then conservative.)
 *   T2d  N_ICI/S of a single path equals 1 - sinc^2(eps) exactly (1e-9).
 */
const Core = require('../core.js');
// z with P(|Z| > z) = p   (bisection on erfc)
function zTwoSided(p) { let lo = 0, hi = 40; for (let i = 0; i < 200; i++) { const mid = (lo + hi) / 2; if (Core.erfc(mid / Math.SQRT2) > p) lo = mid; else hi = mid; } return (lo + hi) / 2; }
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

        // ---------------------------------------------------------------- T2b: noise-free, jammer-free
        Object.assign(sys, BASE, { snr: 400, sir: 400, algo: 'FOURIER' });
        const d2 = VS.map(v => oneDraw(sys, v, seed));
        let maxR2 = 0, maxS2 = 0;
        for (let k = 1; k < VS.length; k++) { maxR2 = Math.max(maxR2, U.maxAbsDiffR(d2[k].R, d2[0].R)); maxS2 = Math.max(maxS2, Math.abs(d2[k].sinr - d2[0].sinr) / d2[0].sinr); }
        console.log(`T2b  noise-free, jammer-free: max |dR_hat| = ${U.e(maxR2)},  max relative |dSINR_inst| (FOURIER) = ${U.e(maxS2)}`);
        checks.push(U.check('T2b R_hat unchanged by v (noise-free, jammer-free)', U.e(maxR2), '1e-9', maxR2 < 1e-9));
        checks.push(U.check('T2b SINR_inst unchanged by v (FOURIER, noise-free, jammer-free)', U.e(maxS2), '1e-9', maxS2 < 1e-9));

        // ---------------------------------------------------------------- T2c: statistical invariance
        // N = 12 so that the R_hat has 2*144 - 12 (identically zero imaginary diagonal) = 276 >= 128 non-trivial values (N = 8 gives only 120)
        Object.assign(sys, BASE, { N: 12 });
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
        let zMax = 0, nEl = 0;
        for (let i = 0; i < sys.N; i++) for (let j = 0; j < sys.N; j++) for (const k of ['r', 'i']) {
            const a = A.sumR[i][j][k], b = B.sumR[i][j][k], s = Math.hypot(U.se(a), U.se(b));
            if (s > 0) { zMax = Math.max(zMax, Math.abs((U.mean(b) - U.mean(a)) / s)); nEl++; }
        }
        console.log(`T2c  ${trials} trials each, v = 0 vs 300:  mean SINR ${U.f(U.mean(A.sinr), 2)} vs ${U.f(U.mean(B.sinr), 2)} dB (z = ${U.f(zS, 2)}),  max |z| over ${nEl} R_hat values (re/im of ${sys.N * sys.N} entries) = ${U.f(zMax, 4)}`);
        checks.push(U.check('T2c mean SINR_inst, v = 300 vs v = 0 (noise + jammer)', `|z| = ${U.f(Math.abs(zS), 2)}`, '|z| < 3', Math.abs(zS) < 3));
const alpha = 0.01, zCrit = zTwoSided(alpha / nEl);
        console.log(`T2c  Bonferroni: m = ${nEl} comparisons, alpha = ${alpha}, z_crit = ${U.f(zCrit, 4)} (fixed-3 reference: family-wise false-alarm rate ${(100 * (1 - Math.pow(1 - 0.0026998, nEl))).toFixed(0)} %),  max |z| = ${U.f(zMax, 4)}`);
        checks.push(U.check('T2c mean R_hat entries, v = 300 vs v = 0 (Bonferroni)', `m = ${nEl}, z_crit = ${U.f(zCrit, 4)}, max |z| = ${U.f(zMax, 4)}`, 'max |z| < z_crit (alpha = 0.01 family-wise)', zMax < zCrit && nEl >= 128));

        // ---------------------------------------------------------------- T2d: N_ICI/S of one path
        let maxN = 0;
        Object.assign(sys, BASE);
        for (const v of VS) { const r = oneDraw(sys, v, seed); const x = r.eps; const s = x === 0 ? 1 : Math.sin(Math.PI * x) / (Math.PI * x); maxN = Math.max(maxN, Math.abs(r.nu - (1 - s * s))); }
        console.log(`T2d  single path: max | N_ICI/S - (1 - sinc^2(eps)) | = ${U.e(maxN)}`);
        checks.push(U.check('T2d N_ICI/S = 1 - sinc^2(eps) for a single path', U.e(maxN), '1e-9', maxN < 1e-9));
        return { id: this.id, title: this.title, checks };
    }
};
