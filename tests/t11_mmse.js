'use strict';
/**
 * T11  The two Wiener variants.  MMSE-M (key 'MMSE'): r_xd = P_s a(theta_hat_1), model-based.  MMSE-P (key 'MMSEP'): r_xd = (1/L) sum x_n conj(s_n), pilot-trained.
 *   T11a  MMSE-M is unchanged: K = 20 dB, v = 0, tau = 0, L = 100, withSignal, 2000 trials (same seeds as the investigation of Commit 10):
 *         legacy 12.03 +- 0.04 dB, unified -5.63 +- 0.08 dB (reference values of the previous round). Criterion: |mean - reference| <= 3 sqrt(SE^2 + SE_ref^2).
 *   T11b  MMSE-P, unified, same scenario: mean (SINR - SINR_opt) > -1 dB, for both training-data settings (MMSE-P ignores trainMode).
 *   T11c  MMSE-P, unified, L = 4, 8, 12, 24, 48, 100: mean SINR must rise monotonically with L; a step violates this only if the paired decrease exceeds one SE
 *         of that paired difference. The whole curve is printed (informational), with the gap to SINR_opt.
 */
const Core = require('../core.js');
const U = require('./_util.js');

const BASE = { N: 8, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, v: 0, latMs: 0, calDeg: 0, taper: 'NONE', kDb: 20, trainMode: 'withSignal' };
function run(model, algo, over, trials, seed0) {
    const out = []; out.opt = [];
    for (let t = 0; t < trials; t++) {
        Core.setSeed(seed0 + t);
        const s = Core.createSys(); s.calZ = new Array(16).fill(0);
        Object.assign(s, BASE, over, { model, algo, freshRealization: true });
        s.snaps = []; s.snapKey = ''; s.computeMath(); out.push(s.sinrDb); out.opt.push(s.sinrOptDb);
    }
    return out;
}
const cell = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`;

module.exports = {
    id: 'T11', title: 'MMSE-M unchanged; MMSE-P within 1 dB of SINR_opt (K = 20 dB, v = 0); MMSE-P rises with L',
    async run({ trials = 2000, seed = 7000 } = {}) {
        const checks = [], info = [];
        for (const [model, ref, refSe] of [['legacy', 12.03, 0.04], ['unified', -5.63, 0.08]]) {
            const a = run(model, 'MMSE', {}, trials, seed), m = U.mean(a), se = U.se(a), tol = 3 * Math.hypot(se, refSe);
            console.log(`T11a ${model} MMSE-M: ${cell(a)} dB (reference ${ref} ± ${refSe}), |diff| = ${U.f(Math.abs(m - ref), 3)}, tolerance ${U.f(tol, 3)}`);
            checks.push(U.check(`T11a MMSE-M ${model}: mean SINR equals the previous value`, `${U.f(m, 2)} dB`, `${ref} ± 3 sqrt(SE² + ${refSe}²) = ${U.f(tol, 2)}`, Math.abs(m - ref) <= tol));
        }
        for (const tm of ['withSignal', 'signalFree']) {
            const a = run('unified', 'MMSEP', { trainMode: tm }, trials, seed), gap = U.mean(a) - U.mean(a.opt);
            console.log(`T11b unified MMSE-P (${tm}): ${cell(a)} dB, SINR_opt ${U.f(U.mean(a.opt), 2)} dB, gap ${U.f(gap, 2)} dB`);
            checks.push(U.check(`T11b MMSE-P unified (${tm}): SINR - SINR_opt`, `${U.f(gap, 2)} dB`, '> -1 dB', gap > -1));
        }
        const Ls = [4, 8, 12, 24, 48, 100], cols = Ls.map(L => run('unified', 'MMSEP', { L }, trials, seed));
        console.log('\nT11c  MMSE-P, unified, K = 20 dB, v = 0, withSignal: mean SINR vs L');
        console.log(U.pad('L', 6), U.rpad('SINR dB', 14), U.rpad('SINR_opt', 9), U.rpad('gap', 8), U.rpad('step vs previous (mean d / SE d)', 34));
        let ok = true, worst = '';
        Ls.forEach((L, i) => {
            let step = '';
            if (i > 0) { const d = cols[i].map((x, k) => x - cols[i - 1][k]), m = U.mean(d), s = U.se(d); step = `${U.f(m, 3)} / ${U.f(s, 3)}`; if (m < -s) { ok = false; worst += ` ${Ls[i - 1]}->${L}: ${U.f(m, 3)}/${U.f(s, 3)}`; } }
            console.log(U.pad(L, 6), U.rpad(cell(cols[i]), 14), U.rpad(U.f(U.mean(cols[i].opt), 2), 9), U.rpad(U.f(U.mean(cols[i]) - U.mean(cols[i].opt), 2), 8), U.rpad(step, 34));
        });
        checks.push(U.check('T11c MMSE-P unified: mean SINR rises with L', `${U.f(U.mean(cols[0]), 2)} -> ${U.f(U.mean(cols[cols.length - 1]), 2)} dB`, 'each step: decrease <= SE(d)', ok, worst ? 'violations:' + worst : ''));
        info.push({ name: 'T11c curve (MMSE-P, unified)', value: Ls.map((L, i) => `L=${L}: ${U.f(U.mean(cols[i]), 2)}`).join(' / ') + ' dB' });
        return { id: this.id, title: this.title, checks, info };
    }
};
