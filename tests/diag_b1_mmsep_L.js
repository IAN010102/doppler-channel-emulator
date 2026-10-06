'use strict';
/**
 * Diagnosis B1 (informational, not part of run_all): MMSE-P around L = N. For L = 2 ... 24 (1000 realisations each, unified, K = 20 dB, v = 0, tau = 0, N = 8):
 *   mean SINR +- SE of MMSE-P, SMI and DL (signalFree and withSignal), median and 90th percentile of kappa(R_hat), median of the smallest non-zero eigenvalue of R_hat
 *   (relative threshold epsRank), median ||w|| of MMSE-P; then, for L = 8, 9, 10: direct inverse vs forced pseudo-inverse of the same R_hat and r_xd (weights and SINR).
 * usage: node tests/diag_b1_mmsep_L.js [trials]
 */
const Core = require('../core.js');
const U = require('./_util.js');
const { Cplx, matMulVec, pinvHermitian } = Core;
const trials = +(process.argv[2] || 1000);
const BASE = { N: 8, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, v: 0, latMs: 0, calDeg: 0, taper: 'NONE', kDb: 20, mod: 'QPSK', gammaRelDb: 10, model: 'unified' };
function mk(over, seed, fresh = true) {
    Core.setSeed(seed);
    const s = Core.createSys(); s.calZ = new Array(16).fill(0); Object.assign(s, BASE, { freshRealization: fresh, trainMode: 'signalFree' }, over);
    s.snaps = []; s.snapKey = ''; if (!fresh) s.newRealization(); s.computeMath(); return s;
}
const q = (a, p) => { const b = a.slice().sort((x, y) => x - y); return b[Math.min(b.length - 1, Math.floor(p * b.length))]; };
const norm = w => Math.sqrt(w.reduce((s, c) => s + c.r * c.r + c.i * c.i, 0));
const cell = a => `${U.f(U.mean(a), 2)}±${U.f(U.se(a), 2)}`;
console.log(`B1: MMSE-P, SMI, DL vs L (N = 8), ${trials} realisations per L, unified, K = 20 dB, v = 0, tau = 0; training modes signalFree / withSignal (MMSE-P ignores it)`);
console.log(U.pad('L', 3), U.rpad('MMSE-P', 13), U.rpad('SMI sF', 13), U.rpad('DL sF', 13), U.rpad('SMI wS', 13), U.rpad('DL wS', 13), U.rpad('kappa med', 11), U.rpad('kappa p90', 11), U.rpad('lmin med', 10), U.rpad('|w| med', 9), U.rpad('rank<N', 7));
for (let L = 2; L <= 24; L++) {
    const A = { mp: [], smiF: [], dlF: [], smiW: [], dlW: [], kap: [], lmin: [], wn: [] }; let deficient = 0;
    for (let t = 0; t < trials; t++) {
        const seed = 9100 + t, p = mk({ algo: 'MMSEP', L }, seed);
        A.mp.push(p.sinrDb); A.kap.push(p.kappaRaw); A.wn.push(norm(p.weights)); if (p.rankR < 8) deficient++;
        const lm = p.eig.filter(e => e > p.epsRank * p.eig[0]); A.lmin.push(lm.length ? lm[lm.length - 1] : 0);
        A.smiF.push(mk({ algo: 'SMI', L, trainMode: 'signalFree' }, seed).sinrDb); A.dlF.push(mk({ algo: 'DL', L, trainMode: 'signalFree' }, seed).sinrDb);
        A.smiW.push(mk({ algo: 'SMI', L, trainMode: 'withSignal' }, seed).sinrDb); A.dlW.push(mk({ algo: 'DL', L, trainMode: 'withSignal' }, seed).sinrDb);
    }
    const kf = A.kap.map(k => Number.isFinite(k) ? k : 1e30);
    console.log(U.pad(L, 3), U.rpad(cell(A.mp), 13), U.rpad(cell(A.smiF), 13), U.rpad(cell(A.dlF), 13), U.rpad(cell(A.smiW), 13), U.rpad(cell(A.dlW), 13), U.rpad(U.e(q(kf, 0.5), 1), 11), U.rpad(U.e(q(kf, 0.9), 1), 11), U.rpad(U.e(q(A.lmin, 0.5), 1), 10), U.rpad(U.f(q(A.wn, 0.5), 2), 9), U.rpad(`${deficient}`, 7));
}
console.log('\nB1b: direct inverse vs forced pseudo-inverse of the same R_hat, r_xd (MMSE-P), 500 realisations per L');
console.log(U.pad('L', 3), U.rpad('max rel |w_inv - w_pinv|', 26), U.rpad('SINR inv', 12), U.rpad('SINR pinv', 12), U.rpad('mean diff', 11));
for (const L of [8, 9, 10]) {
    let worst = 0; const a = [], b = [];
    for (let t = 0; t < 500; t++) {
        const s = mk({ algo: 'MMSEP', L }, 9100 + t, false), N = 8, Klin = 100, paths = s.unifiedPaths(0, Klin), gam = s.gamma();
        const Ts = (1 + Core.CONFIG.cpRatio) / s.scs, tApp = (L - 1) * Ts, Pj = 10, sg = 0.01, thJ = 40 * Math.PI / 180;
        const wI = s.weights.map(c => new Cplx(c.r, c.i)), P = pinvHermitian(s.R_hat, 1e-10).pinv, wP = matMulVec(P, s.rxdHat);
        const dif = wI.map((c, i) => Math.hypot(c.r - wP[i].r, c.i - wP[i].i)); worst = Math.max(worst, Math.sqrt(dif.reduce((x, y) => x + y * y, 0)) / norm(wI));
        s.unifiedMetrics(wI, gam, thJ, tApp, paths, Pj, sg); a.push(10 * Math.log10(s.S / (s.I + s.Nn)));
        s.unifiedMetrics(wP, gam, thJ, tApp, paths, Pj, sg); b.push(10 * Math.log10(s.S / (s.I + s.Nn)));
    }
    const d = a.map((x, i) => x - b[i]);
    console.log(U.pad(L, 3), U.rpad(U.e(worst, 2), 26), U.rpad(cell(a), 12), U.rpad(cell(b), 12), U.rpad(U.e(U.mean(d), 2), 11));
}
