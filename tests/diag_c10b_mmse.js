'use strict';
/**
 * Commit 10(b) (investigation, informational): where does MMSE's cross-correlation vector r_xd come from, and what would a data-driven r_xd give?
 *
 *   Current code (core.js, computeMath, MMSE branch):  r_xd = P_s * a(theta_hat_1)  -- the NOMINAL steering vector at the estimated angle
 *   (`const rxd = nomT.map(c => new Cplx(c.r * sigPow, c.i * sigPow))`), w = R_hat^-1 r_xd.  This is a "model-based" MMSE (known direction, perfectly
 *   correlated pilot with the nominal vector); it never looks at the true h, and it does not use the data to estimate r_xd.
 *
 *   Variants compared here (K = 20 dB, v = 0, tau = 0, SNR 20, SIR -10, theta1 = 0, theta2 = 40, L = 100), mean SINR dB ± SE over the trials:
 *     M0   current MMSE (model-based r_xd = P_s a(theta_1)), R_hat as in the code (the target is always in R_hat for MMSE)
 *     M1   data-driven:  r_xd = (1/L) sum x_n conj(s_n),  x_n = FULL received snapshot (with target), s_n = the known target symbol; R_hat from the same x_n
 *     M2   literal reading with the training data of the chosen trainMode:  x_n = snapshot of the mode (signalFree: target removed)
 *          -> in signalFree r_xd = (1/L) sum x_n conj(s_n) contains no target, so w is meaningless (reported to show exactly that)
 *   Because MMSE keeps the target in its data under both trainModes (PARAMS.md section 6), M0 and M1 are identical for the two modes; only M2 differs.
 *   Reference only: nothing here replaces the existing MMSE or is enabled by default.
 *   node tests/diag_c10b_mmse.js [trials]
 */
const Core = require('../core.js');
const U = require('./_util.js');
const { Cplx, invertMatrix, matMulVec, quadForm } = Core;
const trials = +(process.argv[2] || 2000);

const BASE = { N: 8, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, v: 0, latMs: 0, calDeg: 0, taper: 'NONE', kDb: 20, algo: 'MMSE' };

function dataWeights(snaps, useTarget) {
    const N = snaps[0].rr.length, L = snaps.length;
    const R = Array.from({ length: N }, () => Array.from({ length: N }, () => new Cplx(0, 0))), r = Array.from({ length: N }, () => new Cplx(0, 0));
    for (const sn of snaps) {
        const x = Array.from({ length: N }, (_, i) => new Cplx(sn.rr[i] + (useTarget ? sn.tr[i] : 0), sn.ri[i] + (useTarget ? sn.ti[i] : 0)));
        // this snapshot's trainMode may already contain the target (rr); useTarget adds it back only when it was removed (signalFree)
        for (let m = 0; m < N; m++) {
            r[m] = Cplx.add(r[m], Cplx.mul(x[m], new Cplx(sn.s1r, -sn.s1i)));
            for (let n = 0; n < N; n++) R[m][n] = Cplx.add(R[m][n], Cplx.mul(x[m], Cplx.conj(x[n])));
        }
    }
    const Rn = R.map(row => row.map(c => new Cplx(c.r / L, c.i / L))), rn = r.map(c => new Cplx(c.r / L, c.i / L));
    const inv = invertMatrix(Rn, true); return matMulVec(inv, rn);
}

function evaluate(sys, w, uni, paths, tApp) {
    const Klin = Math.pow(10, sys.kDb / 10), cL2 = Klin / (Klin + 1), cD2 = 1 / (Klin + 1), jam = Math.pow(10, -sys.sir / 10), nz = Math.pow(10, -sys.snr / 10);
    const gam = sys.gamma(), thT = sys.aoaT * Math.PI / 180, thJ = sys.aoaJ * Math.PI / 180;
    if (uni) { sys.unifiedMetrics(w, gam, thJ, tApp, paths, jam, nz); return 10 * Math.log10(sys.S / (sys.I + sys.Nn)); }
    const S = quadForm(w, sys.trueCov(thT, 1, gam, cL2, cD2)), I = quadForm(w, sys.trueCov(thJ, jam, gam, cL2, cD2)), Nn = w.reduce((a, c) => a + c.mag2(), 0) * nz;
    return 10 * Math.log10(S / (I + Nn));
}

const cell = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`;
console.log(`MMSE variants, K = 20 dB, v = 0, tau = 0, ${trials} trials (mean SINR dB ± SE)`);
console.log(U.pad('model', 9), U.pad('trainMode', 11), U.rpad('M0 model-based', 16), U.rpad('M1 data-driven', 16), U.rpad('M2 literal', 16), U.rpad('M1 - M0', 9));
for (const model of ['legacy', 'unified']) for (const train of ['withSignal', 'signalFree']) {
    const m0 = [], m1 = [], m2 = [];
    const uni = model === 'unified';
    for (let t = 0; t < trials; t++) {
        Core.setSeed(7000 + t);
        const s = Core.createSys(); s.calZ = new Array(16).fill(0); Object.assign(s, BASE, { model, trainMode: train, freshRealization: false });
        if (uni) s.newRealization();
        s.snaps = []; s.snapKey = ''; s.computeMath();
        const Klin = Math.pow(10, s.kDb / 10), Tsnap = (1 + Core.CONFIG.cpRatio) / s.scs, tApp = (s.L - 1) * Tsnap;
        const paths = uni ? s.unifiedPaths(s.aoaT * Math.PI / 180, Klin) : null;
        m0.push(s.sinrDb);                                                       // computeMath's own number (model-based MMSE)
        // snaps hold rr (mode data) and tr (target part): M1 needs the full snapshot, M2 the mode's snapshot
        const full = dataWeights(s.snaps, train === 'signalFree');               // signalFree: add the target back
        m1.push(evaluate(s, full, uni, paths, tApp));
        const lit = dataWeights(s.snaps, false);
        m2.push(evaluate(s, lit, uni, paths, tApp));
    }
    console.log(U.pad(model, 9), U.pad(train, 11), U.rpad(cell(m0), 16), U.rpad(cell(m1), 16), U.rpad(cell(m2), 16), U.rpad(U.f(U.mean(m1) - U.mean(m0), 2), 9));
}
