'use strict';
/**
 * T17  Pointing error delta_theta (Commit 20): the weights use a(theta_hat_1 + delta_theta) (FOURIER, SMI, DL, BEAMSPACE with its beam choice, MMSE-M). MMSE-P does not use a
 *      nominal steering vector (r_xd = (1/L) sum x_n conj(s_n) comes from the data), so it is not affected. SINR / EVM are always evaluated on the true channel.
 *   T17a  delta_theta = 0 (set explicitly) reproduces the golden snapshot bit for bit (same cases as T14).
 *   T17b  Convergence check (v = 0, tau = 0, K = 20 dB, SMI, QPSK symbols, 100 realisations per cell, delta_theta in {0, 1, 3, 5} deg): the population solution built DIRECTLY from the path
 *         model, R_exp = P_s h h^H (withSignal only) + P_j g g^H + sigma^2 I  (h = sum_i beta_i Gamma a(theta_i) of the realisation, g = Gamma a(theta_2); v = 0: no Doppler phase),
 *         w = R_exp^-1 a_assumed / (a_assumed^H R_exp^-1 a_assumed), is compared with the sample SMI of L = 20000 snapshots of the same realisation: the mean over realisations of the
 *         paired SINR difference must be within 3 SE. This verifies that the sample pipeline converges to the population solution (it does not test the physics of the model).
 *   T17c  signalFree, SMI (L = 100): mean SINR must decrease with |delta_theta| (|delta_theta| = 0, 1, 2, 3, 5, 7, 10); a step violates this only if the paired increase exceeds one SE.
 *   T17d  (informational) SINR vs delta_theta (-5 ... 5 deg, 0.5 deg steps), SMI / DL / BEAMSPACE / MMSE-M / MMSE-P, signalFree and withSignal, L = 100 and L = 12, 1000 trials, SE.
 *         Also K = 40 dB (nearly pure LoS) for SMI withSignal vs signalFree to show where the self-nulling driven by the pointing error becomes significant.
 *   T17e  (informational) withSignal, delta_theta = 3 deg, L = 100: DL gamma_rel sweep -10 ... +30 dB, best value marked.
 */
const Core = require('../core.js');
const U = require('./_util.js');
const { Cplx, invertMatrix, matMulVec, vecDot } = Core;
const { compute } = require('./golden_cases.js');
const fs = require('fs'), path = require('path');

const BASE = { N: 8, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, v: 0, latMs: 0, calDeg: 0, taper: 'NONE', kDb: 20, mod: 'QPSK', gammaRelDb: 10 };
function mk(over, seed, keepFresh = true) {
    Core.setSeed(seed);
    const s = Core.createSys(); s.calZ = new Array(16).fill(0); Object.assign(s, BASE, { model: 'unified', freshRealization: keepFresh, trainMode: 'signalFree' }, over);
    return s;
}
function run(over, seed) { const s = mk(over, seed); s.snaps = []; s.snapKey = ''; s.computeMath(); return s.sinrDb; }
const cell = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`;

module.exports = {
    id: 'T17', title: 'pointing error: delta_theta = 0 unchanged; sample SMI converges to the population solution; monotone loss; sweeps (info)',
    async run({ trials = 1000, seed = 1717 } = {}) {
        const checks = [], info = [];
        // ---- T17a
        const gold = JSON.parse(fs.readFileSync(path.join(__dirname, 'golden', 'default_snapshot.json'), 'utf8')), now = compute({ pointErrDeg: 0 });
        let worst = 0;
        for (const k of Object.keys(gold.cases)) gold.cases[k].forEach((g, i) => { for (const f of ['sinrDb', 'evm', 'sinrOptDb', 'ser']) worst = Math.max(worst, Math.abs(now.cases[k][i][f] - g[f])); });
        console.log(`T17a  delta_theta = 0 explicit: max |difference| to the golden snapshot over ${Object.keys(gold.cases).length} cases = ${U.e(worst)}`);
        checks.push(U.check('T17a delta_theta = 0: identical to the golden snapshot', U.e(worst), '0 (bit-identical)', worst === 0));

        // ---- T17b
        console.log('\nT17b  population solution from the path model vs sample SMI (L = 20000), K = 20 dB, v = 0, 100 realisations per cell; paired SINR difference (sample - population) [dB]');
        console.log(U.pad('trainMode', 11), U.pad('delta deg', 10), U.rpad('population', 14), U.rpad('sample L=20000', 16), U.rpad('mean diff', 11), U.rpad('SE', 8), 'result');
        for (const tm of ['signalFree', 'withSignal']) for (const dth of [0, 1, 3, 5]) {
            const pop = [], smp = [];
            for (let t = 0; t < 100; t++) {
                const s = mk({ algo: 'SMI', L: 20000, trainMode: tm, pointErrDeg: dth, smiSingular: 'pinv' }, seed + t, false); s.newRealization(); s.snaps = []; s.snapKey = ''; s.computeMath();
                smp.push(s.sinrDb);
                const N = 8, Klin = 100, paths = s.unifiedPaths(0, Klin), gam = s.gamma();
                const h = Array.from({ length: N }, () => new Cplx(0, 0));
                for (let i = 0; i < paths.P; i++) { const k = -2 * Math.PI * 0.5 * Math.sin(paths.th0[i]); for (let e = 0; e < N; e++) h[e] = Cplx.add(h[e], Cplx.mul(new Cplx(paths.br[i], paths.bi[i]), new Cplx(Math.cos(k * e), Math.sin(k * e)))); }
                const hG = h.map((c, e) => Cplx.mul(gam[e], c)), kj = -2 * Math.PI * 0.5 * Math.sin(40 * Math.PI / 180), g = Array.from({ length: N }, (_, e) => Cplx.mul(gam[e], new Cplx(Math.cos(kj * e), Math.sin(kj * e))));
                const Pj = 10, sg = 0.01, R = Array.from({ length: N }, (_, m) => Array.from({ length: N }, (_, n) => {
                    let x = Cplx.mul(g[m], Cplx.conj(g[n])); x = new Cplx(x.r * Pj, x.i * Pj);
                    if (tm === 'withSignal') x = Cplx.add(x, Cplx.mul(hG[m], Cplx.conj(hG[n])));
                    return m === n ? new Cplx(x.r + sg, x.i) : x; }));
                const a = s.steer(dth * Math.PI / 180), inv = invertMatrix(R, true), num = matMulVec(inv, a), den = vecDot(a, num), w = num.map(c => Cplx.div(c, den));
                const Ts = (1 + Core.CONFIG.cpRatio) / s.scs; s.unifiedMetrics(w, gam, 40 * Math.PI / 180, (s.L - 1) * Ts, paths, Pj, sg);
                pop.push(10 * Math.log10(s.S / (s.I + s.Nn)));
            }
            const d = smp.map((x, i) => x - pop[i]), m = U.mean(d), se = U.se(d), pass = Math.abs(m) <= 3 * se;
            console.log(U.pad(tm, 11), U.pad(dth, 10), U.rpad(cell(pop), 14), U.rpad(cell(smp), 16), U.rpad(U.f(m, 4), 11), U.rpad(U.f(se, 4), 8), pass ? 'PASS' : 'FAIL');
            checks.push(U.check(`T17b ${tm}, delta = ${dth} deg: sample SMI (L = 20000) vs population solution`, `${U.f(m, 4)} dB (SE ${U.f(se, 4)})`, '|mean diff| <= 3 SE', pass));
        }

        // ---- T17c
        const ds = [0, 1, 2, 3, 5, 7, 10], cols = ds.map(d => { const a = []; for (let t = 0; t < trials; t++) a.push(run({ algo: 'SMI', L: 100, pointErrDeg: d }, seed + 3000 + t)); return a; });
        let ok = true, worstTxt = '';
        console.log('\nT17c  signalFree SMI, L = 100: mean SINR vs |delta_theta|: ' + ds.map((d, i) => `${d}°: ${cell(cols[i])}`).join(' | '));
        for (let i = 0; i + 1 < ds.length; i++) { const d = cols[i + 1].map((x, k) => x - cols[i][k]), m = U.mean(d), se = U.se(d); if (m > se) { ok = false; worstTxt += ` ${ds[i]}->${ds[i + 1]}: ${U.f(m, 3)}/${U.f(se, 3)}`; } }
        checks.push(U.check('T17c signalFree SMI: mean SINR decreases with |delta_theta|', `${U.f(U.mean(cols[0]), 2)} -> ${U.f(U.mean(cols[ds.length - 1]), 2)} dB`, 'each step: increase <= SE', ok, worstTxt ? 'violations:' + worstTxt : ''));

        // ---- T17d (informational)
        const grid = []; for (let d = -5; d <= 5.0001; d += 0.5) grid.push(Math.round(d * 10) / 10);
        for (const L of [100, 12]) for (const tm of ['signalFree', 'withSignal']) {
            console.log(`\nT17d  ${tm}, L = ${L}, ${trials} trials; mean SINR dB ± SE vs delta_theta`);
            console.log(U.pad('algo', 10), grid.map(d => U.rpad(d, 12)).join(''));
            const table = {};
            for (const algo of ['SMI', 'DL', 'BEAMSPACE', 'MMSE', 'MMSEP']) {
                const row = grid.map(d => { const a = []; for (let t = 0; t < trials; t++) a.push(run({ algo, L, trainMode: tm, pointErrDeg: d }, seed + 6000 + t)); return a; });
                table[algo] = row.map(U.mean);
                console.log(U.pad(algo === 'MMSE' ? 'MMSE-M' : algo === 'MMSEP' ? 'MMSE-P' : algo, 10), row.map(c => U.rpad(`${U.f(U.mean(c), 1)}±${U.f(U.se(c), 1)}`, 12)).join(''));
            }
            info.push({ name: `T17d ${tm} L=${L}`, value: ['SMI', 'DL', 'BEAMSPACE', 'MMSE', 'MMSEP'].map(a => `${a}: ${U.f(table[a][grid.indexOf(0)], 1)} (0°) / ${U.f(table[a][grid.indexOf(3)], 1)} (3°) / ${U.f(table[a][grid.indexOf(5)], 1)} (5°)`).join('; ') });
        }
        console.log('\nT17d  K = 40 dB (nearly pure LoS), SMI, L = 100: signalFree vs withSignal, and the gap');
        console.log(U.pad('delta deg', 10), U.rpad('signalFree', 14), U.rpad('withSignal', 14), U.rpad('gap', 8));
        let onset = null;
        for (const d of [0, 0.5, 1, 1.5, 2, 3, 4, 5]) {
            const a = [], b = []; for (let t = 0; t < trials; t++) { a.push(run({ algo: 'SMI', L: 100, kDb: 40, pointErrDeg: d, trainMode: 'signalFree' }, seed + 8000 + t)); b.push(run({ algo: 'SMI', L: 100, kDb: 40, pointErrDeg: d, trainMode: 'withSignal' }, seed + 8000 + t)); }
            const gap = U.mean(a) - U.mean(b); if (onset === null && gap > 3) onset = d;
            console.log(U.pad(d, 10), U.rpad(cell(a), 14), U.rpad(cell(b), 14), U.rpad(U.f(gap, 2), 8));
        }
        info.push({ name: 'T17d self-nulling onset (K = 40 dB, SMI, L = 100)', value: onset === null ? 'gap never exceeded 3 dB in 0 ... 5 deg' : `smallest delta_theta with (signalFree - withSignal) > 3 dB: ${onset} deg` });

        // ---- T17e (informational)
        const gs = []; for (let g = -10; g <= 30; g += 2) gs.push(g);
        const means = gs.map(g => { const a = []; for (let t = 0; t < trials; t++) a.push(run({ algo: 'DL', L: 100, trainMode: 'withSignal', pointErrDeg: 3, gammaRelDb: g }, seed + 9000 + t)); return U.mean(a); });
        let best = 0; means.forEach((m, i) => { if (m > means[best]) best = i; });
        console.log('\nT17e  withSignal, delta_theta = 3 deg, L = 100, DL: gamma_rel [dB] -> SINR [dB]: ' + gs.map((g, i) => `${g}: ${U.f(means[i], 2)}`).join(' | ') + `   best gamma_rel = ${gs[best]} dB`);
        info.push({ name: 'T17e DL gamma_rel sweep (withSignal, delta = 3 deg)', value: `best gamma_rel ${gs[best]} dB (${U.f(means[best], 2)} dB); +10 dB: ${U.f(means[gs.indexOf(10)], 2)} dB` });
        return { id: this.id, title: this.title, checks, info };
    }
};
