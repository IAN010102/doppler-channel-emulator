'use strict';
/**
 * T12  Small-snapshot region (Commit 15). SMI with R_hat of rank < N: smiSingular = 'pinv' (default) uses the Moore-Penrose pseudo-inverse
 *      (eigenvalues above epsRank * lambda_max), w = R+ a / (a^H R+ a); 'clamp' is the legacy behaviour (pivot clamped to 1e-6).
 *   T12a  R_hat full rank (L >= N): the pseudo-inverse (Core.pinvHermitian) and the direct inverse give the same SMI weights, relative error < 1e-9
 *         (L = 24, 48, 100; 20 realisations each; unified). The default SMI weights at L >= N are the direct-inverse ones (same numbers as before).
 *   T12b  L < N (L = 2, 4, 6) with fixed seeds: the pinv weights satisfy w^H a(theta_hat_1) = 1, error < 1e-9. The constraint holds by the normalisation
 *         w = R+ a / (a^H R+ a) as long as a^H R+ a != 0 (a not orthogonal to the range of R_hat); the test also counts fallbacks.
 *   T12c  (informational) mean SINR vs L = 2 ... 2N (N = 8) for SMI (pinv), DL, BEAMSPACE, MMSE-P; signalFree and withSignal; unified; K = 20 dB, v = 0; 2000
 *         trials per cell with SE. Marks whether SMI shows a jump around L = N (a step of more than 3 dB between neighbouring L, or more than 5 SE).
 *   T12d  DL for L < N does not depend on smiSingular (R_hat + gamma I is non-singular): the weights are identical for 'pinv' and 'clamp' (exactly).
 */
const Core = require('../core.js');
const U = require('./_util.js');
const { Cplx, matMulVec, vecDot } = Core;

const BASE = { N: 8, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, v: 0, latMs: 0, calDeg: 0, taper: 'NONE', kDb: 20, gammaRelDb: 10, trainMode: 'withSignal' };
function mk(model, over, seed) {
    Core.setSeed(seed);
    const s = Core.createSys(); s.rollCal(); Object.assign(s, BASE, over, { model, freshRealization: true });
    s.snaps = []; s.snapKey = ''; s.computeMath(); return s;
}
const cell = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`;

module.exports = {
    id: 'T12', title: 'SMI in the small-snapshot region: pinv vs direct inverse, constraint, SINR vs L (info), DL independent of smiSingular',
    async run({ trials = 2000, seed = 1212 } = {}) {
        const checks = [], info = [];
        // ---- T12a
        for (const L of [24, 48, 100]) {
            let maxRel = 0, maxDef = 0;
            for (let t = 0; t < 20; t++) {
                const s = mk('unified', { algo: 'SMI', L }, seed + t), a = s.steer(s.thTo);
                const pi = Core.pinvHermitian(s.R_raw, Core.CONFIG.epsRank).pinv, num = matMulVec(pi, a), den = vecDot(a, num), wp = num.map(c => Cplx.div(c, den));
                const wi = s.mvdrWeights(s.R_raw, a, true);
                let d = 0, n = 0, d2 = 0; wp.forEach((c, i) => { d += Math.pow(c.r - wi[i].r, 2) + Math.pow(c.i - wi[i].i, 2); n += wi[i].mag2(); d2 += Math.pow(s.weights[i].r - wi[i].r, 2) + Math.pow(s.weights[i].i - wi[i].i, 2); });
                maxRel = Math.max(maxRel, Math.sqrt(d / n)); maxDef = Math.max(maxDef, Math.sqrt(d2 / n));
            }
            console.log(`T12a L = ${L}: max relative |w_pinv - w_inv| / |w_inv| = ${U.e(maxRel)};  default SMI weights vs direct inverse: ${U.e(maxDef)}`);
            checks.push(U.check(`T12a L=${L}: pinv weights equal the direct-inverse weights`, U.e(maxRel), '1e-9', maxRel < 1e-9));
            checks.push(U.check(`T12a L=${L}: default SMI weights (rank N) are the direct-inverse weights`, U.e(maxDef), '1e-12', maxDef < 1e-12));
        }
        // ---- T12b
        for (const L of [2, 4, 6]) {
            let maxErr = 0, fb = 0, n = 50;
            for (let t = 0; t < n; t++) {
                const s = mk('unified', { algo: 'SMI', L }, seed + 100 + t), a = s.steer(s.thTo), g = vecDot(s.weights, a);
                if (s.fallback) fb++; else maxErr = Math.max(maxErr, Math.hypot(g.r - 1, g.i));
            }
            console.log(`T12b L = ${L}: max |w^H a - 1| over ${n} fixed seeds = ${U.e(maxErr)}  (status: RANK-DEFICIENT (pinv); fallbacks: ${fb})`);
            checks.push(U.check(`T12b L=${L}: w^H a(theta_hat_1) = 1 with the pinv weights`, U.e(maxErr), `1e-9 (fallbacks ${fb}/${n})`, maxErr < 1e-9 && fb === 0));
        }
        // ---- T12d
        for (const model of ['legacy', 'unified']) for (const L of [2, 4, 6]) {
            const a = mk(model, { algo: 'DL', L, smiSingular: 'pinv' }, seed + 200), b = mk(model, { algo: 'DL', L, smiSingular: 'clamp' }, seed + 200);
            let d = 0; a.weights.forEach((c, i) => { d = Math.max(d, Math.abs(c.r - b.weights[i].r), Math.abs(c.i - b.weights[i].i)); });
            checks.push(U.check(`T12d ${model} L=${L}: DL weights identical for smiSingular pinv / clamp`, U.e(d), '0 (exact)', d === 0));
        }
        console.log('T12d  DL weights are identical for both smiSingular settings (see the checks).');

        // ---- T12c (informational)
        const Ls = []; for (let L = 2; L <= 16; L++) Ls.push(L);
        for (const tm of ['signalFree', 'withSignal']) {
            console.log(`\nT12c  ${tm}, unified, K = 20 dB, v = 0, ${trials} trials; mean SINR dB ± SE vs L   (SINR_opt ≈ 28.9 dB)`);
            console.log(U.pad('algo', 11), Ls.map(L => U.rpad('L=' + L, 13)).join(''));
            for (const algo of ['SMI', 'DL', 'BEAMSPACE', 'MMSEP']) {
                const cols = Ls.map(L => { const a = []; for (let t = 0; t < trials; t++) { const s = mk('unified', { algo, L, trainMode: tm }, seed + 5000 + t); a.push(s.sinrDb); } return a; });
                console.log(U.pad(algo, 11), cols.map(c => U.rpad(cell(c), 13)).join(''));
                let jump = '';
                for (let i = 0; i + 1 < Ls.length; i++) {
                    const d = cols[i + 1].map((x, k) => x - cols[i][k]), m = U.mean(d), se = U.se(d);
                    if (Math.abs(m) > 3 && Math.abs(m) > 5 * se) jump += ` ${Ls[i]}->${Ls[i + 1]}: ${m > 0 ? '+' : ''}${U.f(m, 2)} dB;`;
                }
                console.log(U.pad('', 11), jump ? 'steps > 3 dB between neighbouring L:' + jump : 'no step > 3 dB between neighbouring L');
                info.push({ name: `T12c ${algo} ${tm}`, value: Ls.map((L, i) => `${L}:${U.f(U.mean(cols[i]), 1)}`).join(' '), note: jump ? 'steps >3 dB:' + jump : 'no step > 3 dB' });
            }
        }
        return { id: this.id, title: this.title, checks, info };
    }
};
