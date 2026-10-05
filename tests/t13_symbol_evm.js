'use strict';
/**
 * T13  Symbol-level EVM (Commit 16), unified model. Core.symbolLevel: Ns random symbols of the selected modulation through the output gain g = w^H h, the residual
 *      interferer (unit-power QPSK x sqrt(I), random phase), the noise w^H n ~ CN(0, Nn) and a Gaussian ICI term of power |g|^2 N_ICI/S; the receiver divides by g
 *      (perfect channel estimate, g known). EVM_meas = sqrt(sum|s_hat - s|^2 / sum|s|^2).
 *   T13a  four settings (signalFree / withSignal  x  v = 0 / 300 km/h), SMI, 16-QAM, K = 20 dB, theta1 = 0, theta2 = 40, L = 100, 300 channel realisations each,
 *         Ns = 4000 per realisation. Criterion: |mean over realisations of (EVM_meas/EVM_analytic - 1)| <= 3 SE of that mean. The relative difference and SE are reported.
 *         NOTE: the interferer, noise and ICI terms are generated from the same S, I, N, nu the analytic EVM uses, so this verifies the arithmetic of the chain
 *         (scaling, normalisation, sample estimators), not the physics of the ICI model (that is T4/T5).
 *   T13b  no interferer, no noise, v = 0, single path: EVM_meas < 1e-9.
 *   T13c  Monte-Carlo SER (decisions on s_hat) vs the closed-form SER of the UI (Gaussian error, per-axis sigma = EVM/sqrt(2)), 16-QAM and 64-QAM,
 *         check (Commit 19): Gaussian interferer, points with >= 200 expected errors: |relative difference| < 10 %; the other rows are informational.
 *         EVM = 5 ... 15 %, 4e5 symbols per point; once with a Gaussian-only error (noise) and once with an interferer-dominated error (QPSK interferer, same EVM).
 *         A deviation is flagged only where the closed form predicts >= 100 errors (otherwise 'n.s.').
 */
const Core = require('../core.js');
const U = require('./_util.js');

const BASE = { model: 'unified', freshRealization: true, algo: 'SMI', N: 8, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, latMs: 0, calDeg: 0, taper: 'NONE', kDb: 20, mod: 'QAM16', gammaRelDb: 10, d_min: 30 };

module.exports = {
    id: 'T13', title: 'symbol-level EVM vs the analytic sqrt(1/SINR + N_ICI/S); no-impairment limit; SER closed form (info)',
    async run({ realisations = 300, seed = 1313 } = {}) {
        const checks = [], info = [];
        console.log(`T13a  SMI, 16-QAM, K = 20 dB, ${realisations} realisations x 4000 symbols per setting`);
        console.log(U.pad('setting', 24), U.rpad('EVM analytic %', 15), U.rpad('EVM measured %', 15), U.rpad('mean rel. diff', 15), U.rpad('SE', 10), U.rpad('|diff|/SE', 10), 'result');
        for (const tm of ['signalFree', 'withSignal']) for (const v of [0, 300]) {
            const rel = [], an = [], me = [];
            for (let t = 0; t < realisations; t++) {
                Core.setSeed(seed + t);
                const s = Core.createSys(); s.calZ = new Array(16).fill(0); Object.assign(s, BASE, { trainMode: tm, v });
                s.snaps = []; s.computeMath();
                const r = s.symbolEvm(4000, 0);
                rel.push(r.evm / s.evm - 1); an.push(s.evm); me.push(r.evm);
            }
            const m = U.mean(rel), se = U.se(rel), z = Math.abs(m) / se, pass = Math.abs(m) <= 3 * se;
            console.log(U.pad(`${tm}, v = ${v}`, 24), U.rpad(U.f(100 * U.mean(an), 3), 15), U.rpad(U.f(100 * U.mean(me), 3), 15), U.rpad((100 * m).toFixed(4) + ' %', 15), U.rpad((100 * se).toFixed(4) + ' %', 10), U.rpad(z.toFixed(2), 10), pass ? 'PASS' : 'FAIL');
            checks.push(U.check(`T13a ${tm}, v = ${v}: mean EVM_meas / EVM_analytic - 1`, `${(100 * m).toFixed(4)} % (SE ${(100 * se).toFixed(4)} %)`, '|mean| <= 3 SE', pass, `analytic ${U.f(100 * U.mean(an), 3)} %, measured ${U.f(100 * U.mean(me), 3)} %`));
        }
        // ---- T13b
        let worst = 0;
        for (const mod of ['QPSK', 'QAM16', 'QAM64']) {
            Core.setSeed(seed + 99);
            const s = Core.createSys(); s.calZ = new Array(16).fill(0); Object.assign(s, BASE, { kDb: 400, snr: 400, sir: 400, v: 0, mod, algo: 'FOURIER' });
            s.snaps = []; s.computeMath(); worst = Math.max(worst, s.symbolEvm(4000, 0).evm);
        }
        console.log(`T13b  no interferer, no noise, v = 0, single path: max EVM_meas over QPSK / 16-QAM / 64-QAM = ${U.e(worst)}`);
        checks.push(U.check('T13b ideal chain: EVM_meas', U.e(worst), '< 1e-9', worst < 1e-9));

        // ---- T13c
        console.log('\nT13c  Monte-Carlo SER vs closed form (Gaussian interferer / Gaussian noise: 2e6 symbols per point; QPSK interferer: 4e5, informational)');
        console.log(U.pad('mod', 8), U.pad('error type', 22), U.pad('EVM %', 7), U.rpad('SER MC', 11), U.rpad('SER closed', 11), U.rpad('rel. diff', 10), 'flag');
        for (const mod of ['QAM16', 'QAM64']) {
            const mi = Core.modInfo(mod);
            for (const kind of ['Gaussian interferer', 'Gaussian (noise)', 'QPSK interferer']) for (const evm of [0.05, 0.075, 0.10, 0.125, 0.15]) {
                Core.setSeed(seed + 7);
                const Ns = kind.startsWith('QPSK') ? 400000 : 2000000, e2 = evm * evm;
                const o = kind.endsWith('(noise)') ? { S: 1, I: 0, Nn: e2, nu: 0 } : { S: 1, I: e2, Nn: 0, nu: 0, jam: kind.startsWith('QPSK') ? 'qpsk' : 'gaussian' };
                const r = Core.symbolLevel(Object.assign({ mod, Ns, keep: 0 }, o));
                const sig = r.evm / Math.SQRT2, pAxis = Math.min(1, 2 * (1 - 1 / mi.ax) * Core.qfunc(mi.a / sig)), closed = 1 - Math.pow(1 - pAxis, 2);
                const rd = closed > 0 ? (r.ser - closed) / closed : NaN, expected = closed * Ns;
                const gate = kind === 'Gaussian interferer' && expected >= 200;              // PASS/FAIL only where >= 200 errors are expected and the interferer is Gaussian
                const flag = gate ? (Math.abs(rd) < 0.10 ? 'check PASS' : 'check FAIL') : (expected < 100 ? 'n.s. (expected < 100 errors)' : (Math.abs(rd) > 0.2 ? 'DEVIATION > 20 %' : ''));
                console.log(U.pad(mod, 8), U.pad(kind, 22), U.pad((100 * r.evm).toFixed(2), 7), U.rpad(U.e(r.ser, 2), 11), U.rpad(U.e(closed, 2), 11), U.rpad((100 * rd).toFixed(1) + ' %', 10), flag);
                if (gate) checks.push(U.check(`T13c ${mod} Gaussian interferer, EVM ${(100 * evm).toFixed(1)} %: MC SER vs closed form`, `${(100 * rd).toFixed(1)} % (MC ${U.e(r.ser, 2)}, closed ${U.e(closed, 2)})`, '|rel. diff| < 10 % (>= 200 expected errors)', Math.abs(rd) < 0.10));
                else info.push({ name: `T13c ${mod} ${kind} EVM ${(100 * evm).toFixed(1)}%`, value: `MC ${U.e(r.ser, 2)} vs closed ${U.e(closed, 2)}`, note: `rel. diff ${(100 * rd).toFixed(1)} %${flag ? ', ' + flag : ''}` });
            }
        }
        return { id: this.id, title: this.title, checks, info };
    }
};
