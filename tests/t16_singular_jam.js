'use strict';
/**
 * T16  Singular handling and jammer waveform (Commit 19).
 *   T16a  The pivot test of invertMatrix is now relative (|pivot| < 1e-12 max|diag|); before it was absolute (|pivot| < 1e-6, kept as Core.invertMatrixAbs).
 *         Unified, SMI, L = N = 8 (and L = 9, 12, 24, 100), 2000 realisations each: the fraction of trials in which the OLD threshold acted is reported; for every trial in
 *         which it did not act, the new inverse must be bit-for-bit identical to the old one (every element). The trials where the old rule acted are listed separately
 *         (with the fraction in which the new rule still acts).
 *   T16b  MMSE-P (and MMSE-M) for L < N (L = 2, 4, 6) with the pseudo-inverse: the residual |R_hat w - r_xd| / |r_xd| of MMSE-P is < 1e-9 (r_xd = (1/L) sum x_n conj(s_n)
 *         lies in the range of R_hat, so R_hat w = r_xd holds exactly on the non-zero eigenspace).
 *   T16c  (informational) SINR vs L = 2 ... 16 for MMSE-P, MMSE-M and SMI (unified, K = 20 dB, v = 0, withSignal training for SMI = signalFree), 1000 trials: is the +3.8 dB jump of
 *         MMSE-P between L = 8 and 9 gone?
 *   (T13c with the Gaussian jammer is in t13_symbol_evm.js.)
 */
const Core = require('../core.js');
const U = require('./_util.js');
const { Cplx, matMulVec, vecDot } = Core;

const BASE = { N: 8, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, v: 0, latMs: 0, calDeg: 0, taper: 'NONE', kDb: 20, mod: 'QPSK', gammaRelDb: 10 };
function mk(over, seed) {
    Core.setSeed(seed);
    const s = Core.createSys(); s.calZ = new Array(16).fill(0); Object.assign(s, BASE, { model: 'unified', freshRealization: true, trainMode: 'withSignal' }, over);
    s.snaps = []; s.snapKey = ''; s.computeMath(); return s;
}
const cell = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`;

module.exports = {
    id: 'T16', title: 'relative pivot threshold vs the old absolute one; MMSE with the pseudo-inverse for L < N; SINR vs L (info)',
    async run({ trials = 2000, seed = 1616 } = {}) {
        const checks = [], info = [];
        console.log('T16a  SMI, unified, K = 20 dB, SNR 20 dB: R_hat inverted with the old (absolute) and the new (relative) pivot threshold');
        console.log(U.pad('L', 5), U.rpad('old threshold acted', 20), U.rpad('...and new acts too', 20), U.rpad('not acted: identical', 22), U.rpad('not acted: differing', 22));
        for (const L of [8, 9, 12, 24, 100]) {
            let clampedOld = 0, alsoNew = 0, same = 0, diff = 0;
            for (let t = 0; t < trials; t++) {
                const s = mk({ algo: 'SMI', L, smiSingular: 'clamp' }, seed + t), R = s.R_raw;
                const inf = { clamped: false }, a = Core.invertMatrixAbs(R, true, inf), b = Core.invertMatrix(R, true);
                if (inf.clamped) {
                    clampedOld++;
                    let nc = false; { const probe = Core.invertMatrix(R, false); nc = probe === null; } if (nc) alsoNew++;
                } else {
                    let eq = true; for (let i = 0; i < 8 && eq; i++) for (let j = 0; j < 8; j++) if (a[i][j].r !== b[i][j].r || a[i][j].i !== b[i][j].i) { eq = false; break; }
                    if (eq) same++; else diff++;
                }
            }
            console.log(U.pad(L, 5), U.rpad(`${clampedOld}/${trials} (${U.f(100 * clampedOld / trials, 2)} %)`, 20), U.rpad(`${alsoNew}`, 20), U.rpad(`${same}`, 22), U.rpad(`${diff}`, 22));
            checks.push(U.check(`T16a L=${L}: trials where the old threshold did not act give a bit-identical inverse`, `${same} identical, ${diff} differing`, '0 differing', diff === 0, `old threshold acted in ${clampedOld}/${trials} trials (${U.f(100 * clampedOld / trials, 2)} %), new rule also in ${alsoNew}`));
            if (L === 8) info.push({ name: 'T16a L = N = 8: old absolute threshold acted', value: `${clampedOld}/${trials} = ${U.f(100 * clampedOld / trials, 2)} %`, note: `new relative rule acts in ${alsoNew} of them` });
        }

        for (const algo of ['MMSEP', 'MMSE']) for (const L of [2, 4, 6]) {
            let worst = 0;
            for (let t = 0; t < 50; t++) {
                const s = mk({ algo, L }, seed + 500 + t);
                let rxd;
                if (algo === 'MMSEP') rxd = s.rxdHat; else { const a = s.steer(s.thTo); rxd = a.map(c => new Cplx(c.r, c.i)); }
                const Rw = matMulVec(s.R_raw, s.weights); let d = 0, n = 0;
                Rw.forEach((c, i) => { d += Math.pow(c.r - rxd[i].r, 2) + Math.pow(c.i - rxd[i].i, 2); n += rxd[i].mag2(); });
                worst = Math.max(worst, Math.sqrt(d / n));
            }
            if (algo === 'MMSEP') checks.push(U.check(`T16b MMSE-P L=${L}: |R_hat w - r_xd| / |r_xd|`, U.e(worst), '< 1e-9', worst < 1e-9));
            else info.push({ name: `T16b MMSE-M L=${L}: |R_hat w - r_xd| / |r_xd|`, value: U.e(worst), note: 'r_xd = P_s a is generally NOT in the range of R_hat for L < N, so a residual is expected (the pinv solution is the projection)' });
            console.log(`T16b ${algo === 'MMSEP' ? 'MMSE-P' : 'MMSE-M'} L = ${L}: max relative residual |R_hat w - r_xd| / |r_xd| over 50 fixed seeds = ${U.e(worst)}`);
        }

        console.log(`\nT16c  unified, K = 20 dB, v = 0, ${Math.min(trials, 1000)} trials; mean SINR dB ± SE vs L`);
        const Ls = []; for (let L = 2; L <= 16; L++) Ls.push(L);
        console.log(U.pad('algo', 14), Ls.map(L => U.rpad('L=' + L, 13)).join(''));
        for (const [algo, label, tm] of [['MMSEP', 'MMSE-P', 'withSignal'], ['MMSE', 'MMSE-M', 'signalFree'], ['SMI', 'SMI signalFree', 'signalFree']]) {
            const cols = Ls.map(L => { const a = []; for (let t = 0; t < Math.min(trials, 1000); t++) a.push(mk({ algo, L, trainMode: tm }, seed + 9000 + t).sinrDb); return a; });
            console.log(U.pad(label, 14), cols.map(c => U.rpad(cell(c), 13)).join(''));
            const i8 = Ls.indexOf(8), d = cols[i8 + 1].map((x, k) => x - cols[i8][k]);
            info.push({ name: `T16c ${label}: step L=8 -> 9`, value: `${U.f(U.mean(d), 2)} ± ${U.f(U.se(d), 2)} dB`, note: 'paired difference' });
            console.log(U.pad('', 14), `step L = 8 -> 9: ${U.f(U.mean(d), 2)} ± ${U.f(U.se(d), 2)} dB`);
        }
        return { id: this.id, title: this.title, checks, info };
    }
};
