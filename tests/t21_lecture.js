'use strict';
/**
 * T21  Lecture baseline (Commit 22b): covSource 'theory', raw EVM, angleSource 'music'.
 * Scenario E0: N = 8, d = 0.5 lambda, theta1 = -20, theta2 = 30, SNR 30 dB, Es = 1, interferer power = Es (SIR 0 dB), L = 1000 snapshots, K -> infinity (kDb = 400), v = 0, no phase mismatch.
 *   T21a  E0: MMSE-M (theory) and SMI withSignal (theory, the lecture's MVDR B) give weights with the same direction (residual after removing the common scalar < 1e-9) and
 *         w_MMSE = Es (a^H R^-1 a) w_MVDR_B (relative error < 1e-9). [The factor in the specification, sigma_s^2 / (a^H R^-1 a), is the one between ... see the output: it is not the ratio
 *         of these two weight vectors; the checked relation follows from w_MMSE = Es R^-1 a and w_B = R^-1 a / (a^H R^-1 a).]
 *   T21b  E0: the SINR of MMSE-M (theory) is within 0.1 dB of SINR_opt (50 realisations; the largest difference is checked).
 *   T21c  one user, no interferer, no noise: the raw EVM of MVDR is 0 (< 1e-9); for MMSE the bias |g' - 1| equals |1 - Es a^H R^-1 a| (computed with an explicit inverse of the theory covariance,
 *         1e-9) and the closed form sigma^2/(N + sigma^2) (K -> infinity: R = a a^H + sigma^2 I, so a^H R^-1 a = N/(N + sigma^2)).
 *   T21d  (informational) angleSource 'music': SINR vs L (4 ... 100) and vs SNR for SMI / DL / BEAMSPACE / MMSE-M against the nominal angle; distribution of the MUSIC angle errors; relation between
 *         the error and the SINR loss.
 *   T21e  T21d needs the DOA code of T18; skipped (and reported) if it is not available.
 */
const Core = require('../core.js');
const U = require('./_util.js');
const { Cplx, invertMatrix, matMulVec, vecDot } = Core;

const E0 = { N: 8, aoaT: -20, aoaJ: 30, snr: 30, sir: 0, L: 1000, kDb: 400, v: 0, latMs: 0, calDeg: 0, mod: 'QPSK', pointErrDeg: 0, model: 'unified', freshRealization: true, taper: 'NONE', jamWave: 'gaussian' };
function run(over, seed) { Core.setSeed(seed); const s = Core.createSys(); s.calZ = new Array(16).fill(0); Object.assign(s, E0, over); s.snaps = []; s.computeMath(); return s; }
const nrm = w => Math.sqrt(w.reduce((s, c) => s + c.mag2(), 0));
const cell = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`;
const pct = (a, q) => { const s = [...a].sort((x, y) => x - y); return s[Math.min(s.length - 1, Math.floor(q * s.length))]; };

module.exports = {
    id: 'T21', title: 'lecture baseline: MMSE (theory R) = MVDR B direction, SINR at SINR_opt, raw EVM bias, MUSIC angle source (info)',
    async run({ trials = 300, seed = 2121 } = {}) {
        const checks = [], info = [];
        // ---------------------------------------------------------------- T21a
        const m = run({ algo: 'MMSE', covSource: 'theory', trainMode: 'withSignal' }, seed), b = run({ algo: 'SMI', covSource: 'theory', trainMode: 'withSignal' }, seed);
        const u = m.weights.map(c => new Cplx(c.r / nrm(m.weights), c.i / nrm(m.weights))), v = b.weights.map(c => new Cplx(c.r / nrm(b.weights), c.i / nrm(b.weights)));
        const proj = vecDot(v, u), resid = nrm(u.map((c, i) => Cplx.sub(c, Cplx.mul(v[i], proj))));
        const a = m.steer(m.thTo), inv = invertMatrix(b.R_raw, false), ar = vecDot(a, matMulVec(inv, a)).r, k = ar * 1;       // a^H R^-1 a (Es = 1)
        let d = 0, nn = 0; m.weights.forEach((c, i) => { const t = Cplx.mul(new Cplx(k, 0), b.weights[i]); d += Math.pow(c.r - t.r, 2) + Math.pow(c.i - t.i, 2); nn += c.mag2(); });
        const relErr = Math.sqrt(d / nn), specFactor = 1 / ar;
        console.log(`T21a  E0: direction residual ${U.e(resid)}, |w_MMSE - Es (a^H R^-1 a) w_B| / |w_MMSE| = ${U.e(relErr)}; a^H R^-1 a = ${U.f(ar, 4)} (Es / (a^H R^-1 a) = ${U.f(specFactor, 6)} is the factor in the other direction: w_B = w_MMSE / (Es a^H R^-1 a))`);
        checks.push(U.check('T21a E0: MMSE-M (theory) and SMI withSignal (theory) have the same weight direction', U.e(resid), '< 1e-9', resid < 1e-9));
        checks.push(U.check('T21a E0: w_MMSE = Es (a^H R^-1 a) w_MVDR_B', U.e(relErr), '< 1e-9', relErr < 1e-9));
        // ---------------------------------------------------------------- T21b
        let worst = 0; const gaps = [];
        for (let t = 0; t < 50; t++) { const s = run({ algo: 'MMSE', covSource: 'theory', trainMode: 'withSignal' }, seed + 10 + t); gaps.push(s.sinrDb - s.sinrOptDb); worst = Math.max(worst, Math.abs(s.sinrDb - s.sinrOptDb)); }
        console.log(`T21b  E0, 50 realisations: SINR(MMSE-M theory) - SINR_opt: mean ${U.e(U.mean(gaps))}, largest |difference| ${U.e(worst)} dB (SINR_opt = ${U.f(m.sinrOptDb, 2)} dB)`);
        checks.push(U.check('T21b E0: MMSE-M (theory) SINR vs SINR_opt', `${U.e(worst)} dB`, '< 0.1 dB', worst < 0.1));
        // ---------------------------------------------------------------- T21c
        const one = { aoaT: 10, sir: 400, kDb: 400, L: 100 };
        const mv = run(Object.assign({ algo: 'SMI', covSource: 'theory', trainMode: 'signalFree', snr: 400 }, one), seed + 1);
        console.log(`T21c  one user, no interferer, no noise: MVDR raw EVM = ${U.e(mv.evmRaw)}, |g' - 1| = ${U.e(Math.hypot(mv.gRawR - 1, mv.gRawI))}`);
        checks.push(U.check('T21c MVDR raw EVM = 0 (one user, no interferer, no noise)', U.e(mv.evmRaw), '< 1e-9', mv.evmRaw < 1e-9));
        const snr = 60, s2 = Math.pow(10, -snr / 10), mm = run(Object.assign({ algo: 'MMSE', covSource: 'theory', trainMode: 'withSignal', snr }, one), seed + 2);
        const aa = mm.steer(mm.thTo), biasWeb = Math.hypot(mm.gRawR - 1, mm.gRawI), biasR = Math.abs(1 - vecDot(aa, matMulVec(invertMatrix(mm.R_raw, false), aa)).r), biasClosed = s2 / (mm.N + s2);
        console.log(`T21c  MMSE (SNR 60 dB, one user): |g' - 1| = ${U.e(biasWeb)};  |1 - Es a^H R^-1 a| = ${U.e(biasR)};  sigma^2/(N + sigma^2) = ${U.e(biasClosed)};  raw EVM = ${U.f(100 * mm.evmRaw, 4)} %`);
        checks.push(U.check('T21c MMSE bias |g\' - 1| = |1 - Es a^H R^-1 a| (explicit inverse of R)', `${U.e(Math.abs(biasWeb - biasR))}`, '< 1e-9', Math.abs(biasWeb - biasR) < 1e-9));
        checks.push(U.check('T21c MMSE bias = sigma^2/(N + sigma^2) (closed form, K -> infinity)', `${U.e(Math.abs(biasWeb - biasClosed))}`, '< 1e-9', Math.abs(biasWeb - biasClosed) < 1e-9));

        // ---------------------------------------------------------------- T21d / T21e
        if (typeof Core.doaSpectra !== 'function') { console.log('T21e  the DOA code (T18) is not available: T21d skipped'); info.push({ name: 'T21d', value: 'skipped (no DOA code)' }); return { id: this.id, title: this.title, checks, info }; }
        const D = { aoaT: 0, aoaJ: 40, snr: 20, sir: -10, kDb: 20, L: 100, mod: 'QAM16' };
        const sinr = (algo, src, tm, over, sd) => { const s = run(Object.assign({}, D, { algo, trainMode: tm, angleSource: src }, over), sd); return s; };
        console.log(`\nT21d  angleSource = music (SMI, DL, BEAMSPACE, MMSE-M), ${trials} trials, K = 20 dB, theta1 = 0, theta2 = 40; mean SINR dB (nominal angle / MUSIC angle)`);
        for (const tm of ['signalFree', 'withSignal']) {
            console.log(`  ${tm}: SINR vs L`); console.log('  ' + U.pad('algo', 11) + [4, 8, 12, 24, 48, 100].map(L => U.rpad('L=' + L, 18)).join(''));
            for (const algo of ['SMI', 'DL', 'BEAMSPACE', 'MMSE']) {
                const cols = [4, 8, 12, 24, 48, 100].map(L => { const a = [], c = []; for (let t = 0; t < trials; t++) { a.push(sinr(algo, 'true', tm, { L }, seed + 500 + t).sinrDb); c.push(sinr(algo, 'music', tm, { L }, seed + 500 + t).sinrDb); } return [U.mean(a), U.mean(c)]; });
                console.log('  ' + U.pad(algo === 'MMSE' ? 'MMSE-M' : algo, 11) + cols.map(x => U.rpad(`${U.f(x[0], 1)} / ${U.f(x[1], 1)}`, 18)).join(''));
                info.push({ name: `T21d SINR vs L, ${algo}, ${tm}`, value: [4, 8, 12, 24, 48, 100].map((L, i) => `L=${L}: ${U.f(cols[i][0], 1)} / ${U.f(cols[i][1], 1)}`).join('  ') });
            }
            console.log(`  ${tm}: SINR vs SNR (L = 100)`); console.log('  ' + U.pad('algo', 11) + [0, 10, 20, 30, 40].map(x => U.rpad('SNR ' + x, 18)).join(''));
            for (const algo of ['SMI', 'DL', 'BEAMSPACE', 'MMSE']) {
                const cols = [0, 10, 20, 30, 40].map(sn => { const a = [], c = []; for (let t = 0; t < trials; t++) { a.push(sinr(algo, 'true', tm, { snr: sn }, seed + 900 + t).sinrDb); c.push(sinr(algo, 'music', tm, { snr: sn }, seed + 900 + t).sinrDb); } return [U.mean(a), U.mean(c)]; });
                console.log('  ' + U.pad(algo === 'MMSE' ? 'MMSE-M' : algo, 11) + cols.map(x => U.rpad(`${U.f(x[0], 1)} / ${U.f(x[1], 1)}`, 18)).join(''));
            }
        }
        // distribution of the MUSIC errors and the relation to the SINR loss (SMI, signalFree)
        for (const L of [100, 12]) {
            const eT = [], eJ = [], loss = [];
            for (let t = 0; t < trials * 2; t++) { const s = sinr('SMI', 'music', 'signalFree', { L }, seed + 1300 + t), n = sinr('SMI', 'true', 'signalFree', { L }, seed + 1300 + t); if (Number.isFinite(s.musicErrT)) { eT.push(Math.abs(s.musicErrT)); eJ.push(Math.abs(s.musicErrJ)); loss.push(n.sinrDb - s.sinrDb); } }
            const bins = [[0, 0.01], [0.01, 0.3], [0.3, 0.6], [0.6, 90]].map(([lo, hi]) => { const ls = loss.filter((_, i) => eT[i] >= lo && eT[i] < hi); return `${lo}-${hi}°: n=${ls.length}, mean loss ${ls.length ? U.f(U.mean(ls), 2) : '-'} dB`; });
            console.log(`  L = ${L}: |error theta1| median ${U.f(pct(eT, 0.5), 2)}°, 90 % ${U.f(pct(eT, 0.9), 2)}°, 99 % ${U.f(pct(eT, 0.99), 2)}°; |error theta2| median ${U.f(pct(eJ, 0.5), 2)}°, 99 % ${U.f(pct(eJ, 0.99), 2)}°; SINR loss (nominal - MUSIC) by error bin: ${bins.join(' | ')}`);
            info.push({ name: `T21d MUSIC errors, L=${L}`, value: `|err theta1| median ${U.f(pct(eT, 0.5), 2)}, 90% ${U.f(pct(eT, 0.9), 2)}, 99% ${U.f(pct(eT, 0.99), 2)} deg; ${bins.join(' | ')}` });
        }
        return { id: this.id, title: this.title, checks, info };
    }
};
