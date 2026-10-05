'use strict';
/**
 * Commit 15, step 1 (report only): what the small-snapshot region L < N actually does today, and the rank / conditioning of R_hat for L = 2, 4, 6, 8 (N = 8).
 * Unified model, K = 20 dB, SNR 20 dB, withSignal, 300 trials per L (median [min, max]). Eigenvalues come from Core.hermitianEigvals (each appears twice).
 *   node tests/diag_c15_paths.js
 */
const Core = require('../core.js');
const U = require('./_util.js');
const trials = 300, N = 8;
const med = a => { const s = [...a].sort((x, y) => x - y); return s[Math.floor(s.length / 2)]; };
console.log(`R_hat of the training window, N = ${N}, ${trials} trials; rank = number of eigenvalues > 1e-10 * lambda_max; kappa_code = what the code reports (infinity if lambda_min <= 1e-12 lambda_max)`);
console.log(U.pad('L', 4), U.pad('rank (min..max)', 16), U.pad('kappa_code', 12), U.pad('lambda_max', 12), U.pad('smallest nonzero eig. [min, median, max]', 44), U.pad('largest "zero" eig.', 20), 'kappa on the nonzero part (median)');
for (const L of [2, 4, 6, 8]) {
    const ranks = [], kinf = [], lmax = [], lnz = [], lz = [], kn = [];
    for (let t = 0; t < trials; t++) {
        Core.setSeed(3000 + t);
        const s = Core.createSys(); s.rollCal(); Object.assign(s, { model: 'unified', algo: 'SMI', L, N, freshRealization: true, kDb: 20, aoaT: 0, aoaJ: 40, v: 0 });
        s.snaps = []; s.computeMath();
        const ev = Core.hermitianEigvals(s.R_raw), e = []; for (let i = ev.length - 1; i >= 0; i -= 2) e.push(Math.max(ev[i], 0));   // descending
        const top = e[0], r = e.filter(x => x > 1e-10 * top).length;
        ranks.push(r); kinf.push(s.kappaRaw === Infinity ? 1 : 0); lmax.push(top);
        lnz.push(e[r - 1]); lz.push(r < N ? e[r] : 0); kn.push(top / e[r - 1]);
    }
    console.log(U.pad(L, 4), U.pad(`${Math.min(...ranks)}..${Math.max(...ranks)}`, 16), U.pad(`inf in ${U.f(100 * U.mean(kinf), 0)}%`, 12), U.pad(U.e(med(lmax), 2), 12),
        U.pad(`[${U.e(Math.min(...lnz), 2)}, ${U.e(med(lnz), 2)}, ${U.e(Math.max(...lnz), 2)}]`, 44), U.pad(U.e(Math.max(...lz), 2), 20), U.e(med(kn), 2));
}

// how often is the absolute pivot threshold (|pivot| < 1e-6, i.e. |pivot|^2 < 1e-12 in invertMatrix) reached although R_hat has full rank (L = N = 8)?
{
    let hit = 0, n = 2000, min = Infinity;
    for (let t = 0; t < n; t++) {
        Core.setSeed(9000 + t);
        const s = Core.createSys(); s.rollCal(); Object.assign(s, { model: 'unified', algo: 'SMI', L: 8, N, freshRealization: true, kDb: 20, aoaT: 0, aoaJ: 40, v: 0 });
        s.snaps = []; s.computeMath();
        const ev = Core.hermitianEigvals(s.R_raw), lmin = Math.max(ev[0], 0); min = Math.min(min, lmin);
        if (lmin < 1e-6) hit++;
    }
    console.log(`\nL = N = 8, ${n} trials: lambda_min < 1e-6 (the scale of the absolute pivot threshold) in ${U.f(100 * hit / n, 1)}% of the trials; smallest lambda_min ${U.e(min, 2)}`);
}
