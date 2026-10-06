'use strict';
/**
 * Diagnosis B4 (report only): frequency-offset compensation for MMSE-P. Unified model, L = 100, SNR 20, SIR -10, theta2 = 40, QPSK symbols, 500 trials per cell,
 * K = 20 dB and K -> infinity (400 dB), theta1 = 0 and 45 deg, training mode signalFree and withSignal (MMSE-P always keeps the target in R_hat, so the training mode does not matter).
 *   (i)   MMSE-P as it is (no compensation): mean SINR and gap to SINR_opt
 *   (ii)  genie: y_n = x_n conj(s_n) is derotated by the TRUE integrated LoS phase of snapshot n before averaging
 *   (iii) LoS rotation inside the window (cycles) and the mean phase increment per snapshot omega = (phi(L-1) - phi(0))/(L-1)  (|omega| < pi required)
 *   (iv)  lag-1 estimator of the specification: omega_hat = angle(sum_n y_n^H y_{n+1}), y'_n = y_n exp(-j omega_hat n), r = mean(y'_n): RMS error of omega_hat, SINR
 *   (v)   periodogram: omega_hat = argmax_omega || sum_n y_n exp(-j omega n) ||^2 over a grid on (-pi, pi) with spacing pi/(8L), refined by a parabola through the peak bin and its two
 *         neighbours: RMS error of omega_hat, SINR with the derotation by omega_hat
 * R_hat is not derotated (the outer product x x^H is invariant to a common phase).
 * usage: node tests/diag_b4_afc.js [trials]
 */
const Core = require('../core.js');
const U = require('./_util.js');
const { Cplx, invertMatrix, matMulVec } = Core;
const trials = +(process.argv[2] || 500), N = 8, L = 100;

// y_n = x_n conj(s_n); returns { y: [ [Cplx x N] x L ], R: sample covariance (N x N Cplx) }
function build(snaps, addTarget) {   // addTarget: signalFree training window (rr holds interference + noise only, the target part tr is added back); withSignal: rr already contains the target
    const y = [], R = Array.from({ length: N }, () => Array.from({ length: N }, () => new Cplx(0, 0)));
    snaps.forEach(sn => {
        const x = Array.from({ length: N }, (_, i) => addTarget ? new Cplx(sn.rr[i] + sn.tr[i], sn.ri[i] + sn.ti[i]) : new Cplx(sn.rr[i], sn.ri[i]));
        for (let m = 0; m < N; m++) for (let k = 0; k < N; k++) R[m][k] = Cplx.add(R[m][k], Cplx.mul(x[m], Cplx.conj(x[k])));
        y.push(x.map(c => Cplx.mul(c, new Cplx(sn.s1r, -sn.s1i))));
    });
    return { y, R: R.map(row => row.map(c => new Cplx(c.r / snaps.length, c.i / snaps.length))) };
}
function solve(y, R, rot) {                              // w = R^-1 mean_n( y_n exp(-j rot(n)) )
    const r = Array.from({ length: N }, () => new Cplx(0, 0));
    y.forEach((v, n) => { const e = new Cplx(Math.cos(-rot(n)), Math.sin(-rot(n))); for (let m = 0; m < N; m++) r[m] = Cplx.add(r[m], Cplx.mul(v[m], e)); });
    return matMulVec(invertMatrix(R, true), r.map(c => new Cplx(c.r / y.length, c.i / y.length)));
}
function lag1(y) {
    let sr = 0, si = 0;
    for (let n = 0; n + 1 < y.length; n++) for (let m = 0; m < N; m++) { const p = Cplx.mul(Cplx.conj(y[n][m]), y[n + 1][m]); sr += p.r; si += p.i; }
    return Math.atan2(si, sr);
}
function periodogram(y) {
    const Ln = y.length, step = Math.PI / (8 * Ln), K = Math.round(2 * Math.PI / step), P = new Float64Array(K), yr = y.map(v => v.map(c => c.r)), yi = y.map(v => v.map(c => c.i));
    const cr = new Float64Array(Ln), ci = new Float64Array(Ln);
    for (let k = 0; k < K; k++) {
        const w = -Math.PI + k * step; for (let n = 0; n < Ln; n++) { cr[n] = Math.cos(w * n); ci[n] = -Math.sin(w * n); }       // exp(-j w n)
        let p = 0;
        for (let m = 0; m < N; m++) { let ar = 0, ai = 0; for (let n = 0; n < Ln; n++) { ar += yr[n][m] * cr[n] - yi[n][m] * ci[n]; ai += yr[n][m] * ci[n] + yi[n][m] * cr[n]; } p += ar * ar + ai * ai; }
        P[k] = p;
    }
    let kb = 0; for (let k = 1; k < K; k++) if (P[k] > P[kb]) kb = k;
    const a = P[(kb - 1 + K) % K], b = P[kb], c = P[(kb + 1) % K], den = a - 2 * b + c, d = den !== 0 ? 0.5 * (a - c) / den : 0;
    return -Math.PI + (kb + d) * step;
}
const cell = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`, rms = a => Math.sqrt(a.reduce((s, x) => s + x * x, 0) / a.length);
const store = {};
for (const kDb of [20, 400]) for (const aoaT of [0, 45]) for (const tm of ['signalFree', 'withSignal']) {
    console.log(`\nK = ${kDb === 400 ? 'inf (400 dB)' : kDb + ' dB'}, theta1 = ${aoaT}, trainMode = ${tm}   (L = ${L}, ${trials} trials; SINR dB ± SE; omega in rad/snapshot)`);
    console.log(U.pad('v km/h', 7), U.rpad('cycles', 7), U.rpad('omega', 8), U.rpad('(i) MMSE-P', 13), U.rpad('gap opt', 8), U.rpad('(ii) genie', 13), U.rpad('gap opt', 8), U.rpad('(iv) lag-1', 13), U.rpad('gap opt', 8), U.rpad('RMS err', 9), U.rpad('(v) period.', 13), U.rpad('gap opt', 8), U.rpad('RMS err', 9));
    for (const v of [0, 30, 100, 300]) {
        const A = [], G = [], P1 = [], P2 = [], O = [], e1 = [], e2 = []; let cyc = 0, om = 0;
        for (let t = 0; t < trials; t++) {
            Core.setSeed(18000 + t);
            const s = Core.createSys(); s.calZ = new Array(16).fill(0);
            Object.assign(s, { model: 'unified', algo: 'MMSEP', trainMode: tm, N, L, aoaT, aoaJ: 40, snr: 20, sir: -10, v, latMs: 0, kDb, calDeg: 0, mod: 'QPSK', freshRealization: false });
            s.newRealization(); s.snaps = []; s.computeMath();
            A.push(s.sinrDb); O.push(s.sinrOptDb);
            const Klin = Math.pow(10, kDb / 10), thT = aoaT * Math.PI / 180, Ts = (1 + Core.CONFIG.cpRatio) / s.scs, tApp = (L - 1) * Ts, paths = s.unifiedPaths(thT, Klin), gam = s.gamma();
            const vms = v / 3.6, lam = s.c / s.fc, thJ = 40 * Math.PI / 180, jam = 10, nz = 0.01;
            const phi = n => Core.trackPhase(paths.th0[0], vms, s.d_min, lam, -tApp, n * Ts - tApp);
            const metric = w => { s.unifiedMetrics(w, gam, thJ, tApp, paths, jam, nz); return 10 * Math.log10(s.S / (s.I + s.Nn)); };
            const { y, R } = build(s.snaps, tm === 'signalFree'), omTrue = (phi(L - 1) - phi(0)) / (L - 1);
            G.push(metric(solve(y, R, phi)));
            const w1 = lag1(y), w2 = periodogram(y);
            P1.push(metric(solve(y, R, n => w1 * n))); P2.push(metric(solve(y, R, n => w2 * n)));
            e1.push(w1 - omTrue); e2.push(w2 - omTrue);
            if (t === 0) { cyc = (phi(L - 1) - phi(0)) / (2 * Math.PI); om = omTrue; }
        }
        if (kDb === 20 && tm === 'signalFree') store[`${aoaT}/${v}`] = { A, G, P2 };
        const mo = U.mean(O);
        console.log(U.pad(v, 7), U.rpad(U.f(cyc, 2), 7), U.rpad(U.f(om, 4), 8), U.rpad(cell(A), 13), U.rpad(U.f(U.mean(A) - mo, 2), 8), U.rpad(cell(G), 13), U.rpad(U.f(U.mean(G) - mo, 2), 8), U.rpad(cell(P1), 13), U.rpad(U.f(U.mean(P1) - mo, 2), 8), U.rpad(U.e(rms(e1), 2), 9), U.rpad(cell(P2), 13), U.rpad(U.f(U.mean(P2) - mo, 2), 8), U.rpad(U.e(rms(e2), 2), 9));
    }
}

// the acceptance criteria of the specification of Commit 18 (T15a, T15b), evaluated on the periodogram variant (informational; nothing is implemented)
console.log('\nCriteria of T15 evaluated with the periodogram estimator (K = 20 dB, L = 100, signalFree = withSignal for MMSE-P)');
for (const aoaT of [0, 45]) {
    const c0 = store[`${aoaT}/0`], c3 = store[`${aoaT}/300`], d0 = c0.P2.map((x, i) => x - c0.A[i]), up = c3.P2.map((x, i) => x - c3.A[i]), gap = c3.G.map((x, i) => x - c3.P2[i]);
    console.log(`theta1 = ${aoaT}: T15a  v = 0, on - off = ${U.f(U.mean(d0), 2)} dB (SE ${U.f(U.se(d0), 2)}, |z| = ${U.f(Math.abs(U.mean(d0)) / U.se(d0), 1)}) -> ${Math.abs(U.mean(d0)) <= 3 * U.se(d0) ? 'would PASS' : 'would FAIL'};  T15b  v = 300: on - off = ${U.f(U.mean(up), 2)} dB (>= 10: ${U.mean(up) >= 10 ? 'yes' : 'no'}), genie - on = ${U.f(U.mean(gap), 2)} dB (< 3: ${U.mean(gap) < 3 ? 'yes' : 'no'})`);
}
