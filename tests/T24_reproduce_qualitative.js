'use strict';
/**
 * T24  Gopala & Slock 重現的「質性」檢查（只檢查有統計依據的質性結論；數值不要求吻合文獻，也不為了吻合而調參數）。
 *      設定同 tests/reproduce_gopala_slock.js：N = 1024、Nr = 2、SNR 0–35 dB（步長 5）、4 路徑、各 200 次通道實現；
 *      到達角兩種做法：固定一組（F1）與每次 U(−90°, 90°) 隨機抽取。所有比較都是「成對」（同一個通道實現）的差，單位 dB。
 *   T24a  最佳波束成形（含 CP 法對頻）的 SINR ≥ MRC 的 SINR，且 ≥ 單天線的 SINR：每個 SNR 點，成對差的平均 ≥ −3 SE（≤ 3 SE 的違反容忍）。
 *   T24b  SNR 高的區間（25、30、35 dB），CP 法對頻的最佳 BF 比不對頻者高：成對差的平均 > 0 且 > 3 SE；差距回報。
 *   T24c  （資訊性，不判定）Estimated BF（R_yy 與 Ĥ，Ns = N，CP 法估計對頻）與 Optimal BF 的 SINR 差，以及用真實 H 取代 Ĥ 的診斷版。
 * 判準在看結果之前訂定，不放寬。
 */
const R = require('./reproduce_gopala_slock.js');
const U = require('./_util.js');

module.exports = {
    id: 'T24', title: 'Gopala & Slock reproduction, qualitative checks: optimal >= MRC and single antenna; CP demodulation helps at high SNR; estimated BF gap (info)',
    async run({ draws = R.CFG.draws } = {}) {
        const checks = [], info = [], modes = [['fixed F1', { mode: 'fixed', set: 'F1' }], ['random angles', { mode: 'random' }]];
        for (const [label, m] of modes) {
            const res = R.computeCurves(Object.assign({ draws }, m)), C = res.curves, snrs = res.snrs, diff = (a, b, i) => C[a][i].map((x, k) => x - C[b][i][k]);
            console.log(`\nT24  ${label}, ${draws} realisations: paired differences in dB (mean ± SE)`);
            console.log(U.pad('SNR', 5), U.rpad('opt - MRC', 16), U.rpad('opt - single', 16), U.rpad('CP - no demod', 16), U.rpad('est - opt', 16), U.rpad('est(true H) - opt', 19));
            let okA1 = true, okA2 = true, worstA = '';
            snrs.forEach((s, i) => {
                const dm = diff('optCp', 'mrc', i), ds = diff('optCp', 'single', i), dc = diff('optCp', 'optNone', i), de = diff('est', 'optCp', i), dt = diff('estTrueH', 'optCp', i);
                const cell = d => `${U.f(U.mean(d), 2)} ± ${U.f(U.se(d), 2)}`;
                console.log(U.pad(s, 5), U.rpad(cell(dm), 16), U.rpad(cell(ds), 16), U.rpad(cell(dc), 16), U.rpad(cell(de), 16), U.rpad(cell(dt), 19));
                if (U.mean(dm) < -3 * U.se(dm)) { okA1 = false; worstA += ` MRC@${s}`; } if (U.mean(ds) < -3 * U.se(ds)) { okA2 = false; worstA += ` single@${s}`; }
                info.push({ name: `T24c ${label}, SNR ${s} dB`, value: `Estimated BF - Optimal BF ${cell(de)} dB; with the true H instead of H_hat ${cell(dt)} dB` });
            });
            checks.push(U.check(`T24a ${label}: optimal BF (CP demod) >= MRC at every SNR`, okA1 ? 'holds' : 'violated at' + worstA, 'mean paired difference >= -3 SE', okA1));
            checks.push(U.check(`T24a ${label}: optimal BF (CP demod) >= single antenna at every SNR`, okA2 ? 'holds' : 'violated at' + worstA, 'mean paired difference >= -3 SE', okA2));
            for (const s of [25, 30, 35]) {
                const i = snrs.indexOf(s), d = diff('optCp', 'optNone', i), m2 = U.mean(d), se = U.se(d);
                checks.push(U.check(`T24b ${label}, SNR ${s} dB: optimal BF with CP demodulation above the one without`, `${U.f(m2, 3)} dB (SE ${U.f(se, 3)}, ${U.f(m2 / se, 1)} SE)`, 'mean paired difference > 0 and > 3 SE', m2 > 0 && m2 > 3 * se));
            }
        }
        return { id: this.id, title: this.title, checks, info };
    }
};
