'use strict';
/**
 * Commit 10(a) ablation (investigation, informational): why is BEAMSPACE (signalFree, K = 20 dB, v = 0, tau = 0) 2.9 dB higher in the unified model
 * than in the legacy model?  Also run for SMI and DL for comparison.
 *
 *   Definitions of "SINR" for the SAME weights w:
 *     B  = realised channel:  |w^H h|^2 / (P_j |w^H a_J|^2 + sigma^2 |w|^2)               (unified native; h = this trial's channel vector)
 *     A  = expected covariance:  w^H R_t w / (w^H R_j w + sigma^2 |w|^2)                  (legacy native; R_t from Sys.trueCov, K, sigma_theta)
 *     A' = A with a SINGLE-PATH jammer  (R_j = P_j |w^H Gamma a_J|^2, the unified jammer model)
 *     A_los = A' with only the LoS part of the target counted as signal; the target's diffuse part is counted as interference (item 4 of the task)
 *   Rows:
 *     L0   legacy weights (Rician jammer in the snapshots), definition A                  (= what the UI shows, legacy)
 *     L1   legacy with the jammer's diffuse part switched off in the snapshots, A'        (jammer = single path, as in the unified model)
 *     L2   L1 with A_los
 *     U0   unified weights, definition B                                                  (= what the UI shows, unified)
 *     U1   unified weights, definition A'  (expected-covariance SINR of the unified weights)
 *     U2   unified weights, A_los
 *     U3   unified weights, definition B, but averaged in the LINEAR domain and shown in dB (removes the Jensen effect of averaging dB)
 *   Beam selection (BEAMSPACE only): the beam set chosen by the algorithm in both models, per trial; fraction of identical sets.
 *   Same seeds in both models; everything paired as far as the models allow.
 *   node tests/diag_c10a_beamspace.js [trials]
 */
const Core = require('../core.js');
const U = require('./_util.js');
const { quadForm } = Core;
const trials = +(process.argv[2] || 2000);

const BASE = { N: 8, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, v: 0, latMs: 0, calDeg: 0, taper: 'NONE', gammaDL: 0.01, gammaRelDb: 0, freshRealization: true };

function metrics(sys, w, kind) {
    const Klin = Math.pow(10, sys.kDb / 10), cL2 = Klin / (Klin + 1), cD2 = 1 / (Klin + 1), jam = Math.pow(10, -sys.sir / 10), nz = Math.pow(10, -sys.snr / 10);
    const gam = sys.gamma(), thT = sys.aoaT * Math.PI / 180, thJ = sys.aoaJ * Math.PI / 180;
    const Nn = w.reduce((a, c) => a + c.mag2(), 0) * nz;
    const St = quadForm(w, sys.trueCov(thT, 1, gam, cL2, cD2));
    const Sl = quadForm(w, sys.trueCov(thT, 1, gam, cL2, 0)), Sd = quadForm(w, sys.trueCov(thT, 1, gam, 0, cD2));
    const Ij1 = quadForm(w, sys.trueCov(thJ, jam, gam, 1, 0));                       // single-path jammer
    const IjR = quadForm(w, sys.trueCov(thJ, jam, gam, cL2, cD2));                    // Rician jammer (legacy)
    const IjP = quadForm(w, sys.trueCov(thJ, jam, gam, cL2, 0));                      // legacy jammer with the diffuse part removed (power cL2)
    return { A: St / (IjR + Nn), Ap: St / (Ij1 + Nn), Alos: Sl / (Ij1 + Nn + Sd), AplainL: St / (IjP + Nn), AplainLlos: Sl / (IjP + Nn + Sd) };
}

function run(algo, kDb, train) {
    const over = Object.assign({}, BASE, { algo, kDb, trainMode: train });
    const out = { L0: [], L1: [], L2: [], U0: [], U1: [], U2: [], U0lin: [], binsL: [], binsU: [] };
    // legacy
    for (let t = 0; t < trials; t++) {
        for (const variant of ['L0', 'L1']) {
            Core.setSeed(5000 + t);
            const s = Core.createSys(); s.calZ = new Array(16).fill(0); Object.assign(s, over, { model: 'legacy' });
            if (variant === 'L1') { const orig = s.diffuse.bind(s); let n = 0; s.diffuse = (rad, gr, gi) => { orig(rad, gr, gi); if (n++ % 2 === 1) { gr.fill(0); gi.fill(0); } }; }   // 2nd call per snapshot = jammer
            s.snaps = []; s.snapKey = ''; s.computeMath();
            const m = metrics(s, s.weights);
            if (variant === 'L0') { out.L0.push(10 * Math.log10(m.A)); out.binsL.push(s.bsBins.join(',')); }
            else { out.L1.push(10 * Math.log10(m.AplainL)); out.L2.push(10 * Math.log10(m.AplainLlos)); }
        }
    }
    // unified
    for (let t = 0; t < trials; t++) {
        Core.setSeed(5000 + t);
        const s = Core.createSys(); s.calZ = new Array(16).fill(0); Object.assign(s, over, { model: 'unified' });
        s.snaps = []; s.computeMath();
        const m = metrics(s, s.weights);
        out.U0.push(s.sinrDb); out.U0lin.push(Math.pow(10, s.sinrDb / 10)); out.U1.push(10 * Math.log10(m.Ap)); out.U2.push(10 * Math.log10(m.Alos)); out.binsU.push(s.bsBins.join(','));
    }
    return out;
}
const cell = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`;

for (const [label, kDb] of [['K = 20 dB', 20], ['K -> inf (400 dB)', 400]]) for (const train of ['signalFree', 'withSignal']) {
    console.log(`\n=== ${label}, ${train}, v = 0, tau = 0, ${trials} trials (mean SINR dB ± SE) ===`);
    console.log(U.pad('algo', 10), ['L0', 'L1', 'L2', 'U0', 'U1', 'U2', 'U0 (lin. mean)', 'U0 - L0', 'U1 - L1', 'U2 - L2', 'bins equal'].map((x, i) => U.rpad(x, i < 6 ? 14 : (i === 6 ? 15 : 10))).join(''));
    for (const algo of ['BEAMSPACE', 'SMI', 'DL']) {
        const r = run(algo, kDb, train);
        const linDb = 10 * Math.log10(U.mean(r.U0lin));
        const same = algo === 'BEAMSPACE' ? r.binsL.filter((b, i) => b === r.binsU[i]).length / trials : NaN;
        console.log(U.pad(algo, 10), [r.L0, r.L1, r.L2, r.U0, r.U1, r.U2].map(a => U.rpad(cell(a), 14)).join(''), U.rpad(U.f(linDb, 2), 15),
            U.rpad(U.f(U.mean(r.U0) - U.mean(r.L0), 2), 10), U.rpad(U.f(U.mean(r.U1) - U.mean(r.L1), 2), 10), U.rpad(U.f(U.mean(r.U2) - U.mean(r.L2), 2), 10), U.rpad(Number.isFinite(same) ? U.f(same, 3) : '-', 10));
        if (algo === 'BEAMSPACE') {
            const cnt = (arr) => { const m = {}; arr.forEach(b => m[b] = (m[b] || 0) + 1); return Object.entries(m).sort((a, b) => b[1] - a[1]).slice(0, 4).map(([k, v]) => `{${k}}:${(100 * v / arr.length).toFixed(0)}%`).join('  '); };
            console.log('   beam sets  legacy: ' + cnt(r.binsL) + '    unified: ' + cnt(r.binsU));
        }
    }
}
