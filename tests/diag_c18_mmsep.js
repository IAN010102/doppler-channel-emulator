'use strict';
/**
 * Commit 18 diagnosis (report only): why does MMSE-P degrade with speed?  Unified model, K = 20 dB, L = 100, SNR 20, SIR -10, theta2 = 40, 500 trials per cell.
 *   (i)   current MMSE-P: mean SINR and gap to SINR_opt
 *   (ii)  genie version: y_n = x_n conj(s_n) is derotated by the TRUE LoS phase (beta phase excluded: only exp(j phi_LoS(t_n))) before averaging
 *   (iii) LoS phase rotation inside the window (cycles = f_d t_est) and the per-snapshot increment omega = 2 pi f_d T_snap (|omega| < pi required)
 *   (iv)  preview of the estimator proposed in the specification: omega_hat = angle(sum y_n^H y_{n+1}), y'_n = y_n exp(-j omega_hat n), r = mean(y'_n)  (NOT part of the code)
 *   (iv-b) the same with the lag-1 product taken on a spatially filtered scalar z_n (first-stage DL-MVDR filter) - evidence only, not the specified estimator
 *   node tests/diag_c18_mmsep.js [trials]
 */
const Core = require('../core.js');
const U = require('./_util.js');
const { Cplx, invertMatrix, matMulVec } = Core;
const trials = +(process.argv[2] || 500), N = 8;

function weights(snaps, derot, afc, aNom, addTarget) {
    const L = snaps.length, y = [], R = Array.from({ length: N }, () => Array.from({ length: N }, () => new Cplx(0, 0)));
    snaps.forEach((sn, n) => {
        const x = Array.from({ length: N }, (_, i) => addTarget ? new Cplx(sn.rr[i] + sn.tr[i], sn.ri[i] + sn.ti[i]) : new Cplx(sn.rr[i], sn.ri[i]));      // full snapshot (the target is always in the Wiener data; withSignal: rr already contains it)
        for (let m = 0; m < N; m++) for (let k = 0; k < N; k++) R[m][k] = Cplx.add(R[m][k], Cplx.mul(x[m], Cplx.conj(x[k])));
        y.push(x.map(c => Cplx.mul(c, new Cplx(sn.s1r, -sn.s1i))));
    });
    let rot = n => 0;
    if (derot) rot = derot;
    if (afc === 'raw') {
        let sr = 0, si = 0;
        for (let n = 0; n + 1 < L; n++) for (let m = 0; m < N; m++) { const p = Cplx.mul(Cplx.conj(y[n][m]), y[n + 1][m]); sr += p.r; si += p.i; }
        const w = Math.atan2(si, sr); rot = n => w * n;
    } else if (afc === 'filtered') {
        // (iv-b) NOT the specified estimator: the lag-1 product is taken on z_n = w1^H y_n, w1 = (R_hat + 0.1 I)^-1 a_nominal (a first-stage DL-MVDR filter that suppresses the jammer)
        const Rl = R.map((row, i) => row.map((c, j) => new Cplx(c.r / L + (i === j ? 0.1 : 0), c.i / L)));
        const w1 = matMulVec(invertMatrix(Rl, true), aNom), z = y.map(v => v.reduce((a, c, m) => Cplx.add(a, Cplx.mul(Cplx.conj(w1[m]), c)), new Cplx(0, 0)));
        let sr = 0, si = 0;
        for (let n = 0; n + 1 < L; n++) { const p = Cplx.mul(Cplx.conj(z[n]), z[n + 1]); sr += p.r; si += p.i; }
        const w = Math.atan2(si, sr); rot = n => w * n;
    }
    const r = Array.from({ length: N }, () => new Cplx(0, 0));
    y.forEach((v, n) => { const e = new Cplx(Math.cos(-rot(n)), Math.sin(-rot(n))); for (let m = 0; m < N; m++) r[m] = Cplx.add(r[m], Cplx.mul(v[m], e)); });
    const inv = invertMatrix(R.map(row => row.map(c => new Cplx(c.r / L, c.i / L))), true);
    return matMulVec(inv, r.map(c => new Cplx(c.r / L, c.i / L)));
}
const cell = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`;

for (const aoaT of [0, 45]) for (const tm of ['signalFree', 'withSignal']) {
    console.log(`\ntheta1 = ${aoaT}, trainMode = ${tm}   (K = 20 dB, L = 100, ${trials} trials; SINR dB ± SE)`);
    console.log(U.pad('v km/h', 8), U.rpad('cycles', 8), U.rpad('omega rad', 10), U.rpad('(i) MMSE-P', 14), U.rpad('gap to opt', 11), U.rpad('(ii) genie', 14), U.rpad('gap to opt', 11), U.rpad('(iv) spec. est.', 14), U.rpad('gap to opt', 11), U.rpad('(iv-b) filtered', 14), U.rpad('gap to opt', 11));
    for (const v of [0, 30, 100, 300]) {
        const a = [], g = [], p = [], q = [], opt = []; let cyc = 0, om = 0;
        for (let t = 0; t < trials; t++) {
            Core.setSeed(18000 + t);
            const s = Core.createSys(); s.calZ = new Array(16).fill(0);
            Object.assign(s, { model: 'unified', algo: 'MMSEP', trainMode: tm, N, L: 100, aoaT, aoaJ: 40, snr: 20, sir: -10, v, latMs: 0, kDb: 20, calDeg: 0, mod: 'QPSK', freshRealization: false });
            s.newRealization(); s.snaps = []; s.computeMath();
            a.push(s.sinrDb); opt.push(s.sinrOptDb);
            const Klin = 100, thT = aoaT * Math.PI / 180, Ts = (1 + Core.CONFIG.cpRatio) / s.scs, tApp = (s.L - 1) * Ts, paths = s.unifiedPaths(thT, Klin), gam = s.gamma();
            const vms = v / 3.6, lam = s.c / s.fc, thJ = 40 * Math.PI / 180, jam = Math.pow(10, 1), nz = Math.pow(10, -2);
            const phi = n => Core.trackPhase(paths.th0[0], vms, s.d_min, lam, -tApp, n * Ts - tApp);                  // LoS phase of snapshot n
            const metric = w => { s.unifiedMetrics(w, gam, thJ, tApp, paths, jam, nz); return 10 * Math.log10(s.S / (s.I + s.Nn)); };
            g.push(metric(weights(s.snaps, phi, false, null, tm === 'signalFree'))); p.push(metric(weights(s.snaps, null, 'raw', null, tm === 'signalFree'))); q.push(metric(weights(s.snaps, null, 'filtered', s.steer(thT), tm === 'signalFree')));
            if (t === 0) { cyc = (phi(s.L - 1) - phi(0)) / (2 * Math.PI); om = 2 * Math.PI * paths.fd[0] * Ts; }
        }
        const mo = U.mean(opt);
        console.log(U.pad(v, 8), U.rpad(U.f(cyc, 2), 8), U.rpad(U.f(om, 3), 10), U.rpad(cell(a), 14), U.rpad(U.f(U.mean(a) - mo, 2), 11), U.rpad(cell(g), 14), U.rpad(U.f(U.mean(g) - mo, 2), 11), U.rpad(cell(p), 14), U.rpad(U.f(U.mean(p) - mo, 2), 11), U.rpad(cell(q), 14), U.rpad(U.f(U.mean(q) - mo, 2), 11));
    }
}
