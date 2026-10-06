'use strict';
/**
 * Diagnostic for the T1 gap (informational, not part of run_all): why the unified model gives a much lower SINR than the legacy
 * model for the adaptive weights at v = 0, K = 20 dB.
 *
 * In the unified model the diffuse part of the channel is a fixed vector h = sum_i beta_i a(theta_i) during a trial. The adaptive
 * weights are trained on data that contain this very h, but their distortionless constraint is placed on the nominal a(theta_1).
 * When h differs from a(theta_1) (a few degrees in the N-dim space) and the array SNR is high, the MPDR/Capon design cancels the
 * signal (look-direction mismatch / self-nulling). With the oracle constraint on h the same R_hat gives a high SINR.
 * The legacy model re-draws the diffuse part every snapshot, so no fixed h can be nulled.
 *
 *   node tests/diag_t1_mismatch.js
 */
const Core = require('../core.js');
const { Cplx, invertMatrix, matMulVec, vecDot } = Core;
const U = require('./_util.js');

Core.setSeed(99);
const S = Core.createSys(); S.rollCal();
Object.assign(S, { model: 'unified', freshRealization: false, algo: 'SMI', N: 8, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, v: 0, latMs: 0, kDb: 20, calDeg: 0 });
const N = S.N, trials = 600;
const rep = [], nominal = [], oracle = [], angle = [];
for (let t = 0; t < trials; t++) {
    S.newRealization(); S.snaps = []; S.computeMath();                    // one realisation, used for training AND evaluation
    const paths = S.unifiedPaths(0, Math.pow(10, S.kDb / 10));
    const h = Array.from({ length: N }, () => new Cplx(0, 0));
    for (let i = 0; i < paths.P; i++) {
        const k = -2 * Math.PI * S.d_lambda * Math.sin(paths.th0[i]);
        for (let e = 0; e < N; e++) h[e] = Cplx.add(h[e], Cplx.mul(new Cplx(paths.br[i], paths.bi[i]), new Cplx(Math.cos(k * e), Math.sin(k * e))));
    }
    const aN = S.steer(0), aj = S.steer(40 * Math.PI / 180), inv = invertMatrix(S.R_raw, true);
    const sinr = w => 10 * Math.log10(vecDot(w, h).mag2() / (10 * vecDot(w, aj).mag2() + w.reduce((s, c) => s + c.mag2(), 0) * 0.01));
    const design = v => { const n = matMulVec(inv, v), d = vecDot(v, n); return n.map(c => Cplx.div(c, d)); };
    rep.push(S.sinrDb); nominal.push(sinr(design(aN))); oracle.push(sinr(design(h)));
    const hn = Math.sqrt(h.reduce((s, c) => s + c.mag2(), 0));
    angle.push(Math.acos(Math.min(1, vecDot(aN, h).mag() / (hn * Math.sqrt(N)))) * 180 / Math.PI);
}
console.log(`v = 0, K = 20 dB, SMI-MVDR, L = 100, ${trials} trials (same realisation for training and evaluation)`);
console.log(`sinrDb reported by computeMath:                 ${U.f(U.mean(rep), 2)} dB`);
console.log(`re-computed, constraint on a(theta_1):          ${U.f(U.mean(nominal), 2)} dB   (must equal the line above)`);
console.log(`re-computed, ORACLE constraint on h(t_app):     ${U.f(U.mean(oracle), 2)} dB`);
console.log(`mean angle between a(theta_1) and h:            ${U.f(U.mean(angle), 1)} deg`);
const ok = Math.abs(U.mean(rep) - U.mean(nominal)) < 1e-9 && U.mean(oracle) > U.mean(nominal) + 10;
console.log(ok ? 'consistent: the gap is caused by look-direction mismatch under static fading' : 'UNEXPECTED');
process.exit(ok ? 0 : 1);
