'use strict';
/**
 * T25  導頻殘差協方差估計器（ofdm_ici.js: pilotResidualWeights）。判準在看結果之前訂定，不放寬；失敗時回報原因。
 *   設定：N = 1024、文獻參數（4 路徑）、隨機到達角、CP 法由含雜訊樣本估計 ε̂、Nr ∈ {2, 4, 8}、SNR ∈ {10, 20, 30, 35} dB、每格 200 次實現；
 *   損失 = Optimal（同一個 ε̂）的 SINR − 該方法的 SINR（dB，成對平均）。
 *   T25a  Ĥ = 真實 H（genie）、|P| 很大（導頻間距 1，即全部 1024 個子載波都是導頻）、無加載：Pilot-residual 與 Optimal 的 SINR 差 < 1 dB
 *         （每個 Nr、每個 SNR 的平均損失）。
 *   T25b  SNR ≥ 30 dB（30 與 35 dB）：Pilot-residual（Ĥ 由導頻估計、導頻間距 12、m = 1、γ_rel = +10 dB、不修正）比 Estimated（含訊號的樣本 R_yy）
 *         至少高 3 dB（每個 Nr、每個 SNR：平均 SINR 差 ≥ 3 dB）。
 *   T25c  （資訊性，單獨列出）其他方法與兩個對照的差：Ĥ 改用真實 H、殘差乘以 n/(n−1)；m 的影響；Ledoit–Wolf 收縮。
 */
const O = require('../ofdm_ici.js');
const U = require('./_util.js');
const R = require('./reproduce_gopala_slock.js');
const P = require('./pilot_residual_curves.js');

const N = R.CFG.N, Ncp = R.CFG.Ncp, dB = x => 10 * Math.log10(x);
module.exports = {
    id: 'T25', title: 'pilot-residual covariance estimator: close to optimal with the true H and dense pilots; at least 3 dB above the estimator with the signal in R_yy at high SNR; controls (info)',
    async run({ draws = 200, seed = 2525 } = {}) {
        const checks = [], info = [], NRS = [2, 4, 8], SNRS = [10, 20, 30, 35];
        // ---------------------------------------------------------------- T25a
        console.log(`T25a  H_hat = true H, all ${N} subcarriers are pilots, no loading: loss of the pilot-residual weights against Optimal [dB], ${draws} realisations`);
        console.log(U.pad('Nr', 4) + SNRS.map(s => U.rpad(`SNR ${s} dB`, 18)).join(''));
        for (const Nr of NRS) {
            const cells = [];
            for (const snr of SNRS) {
                const a = [];
                for (let t = 0; t < draws; t++) {
                    const rp = O.makeRng(seed + 7919 * t), ra = O.makeRng(seed + 500000 + 7919 * t), rs = O.makeRng(seed + 1000000 + 7919 * t), ch = O.paperChannel(R.pickAngles('random', ra), rp, { Nr, N }), sn2 = O.snrToNoise(ch, snr, 1);
                    const sym = O.makeSymbol(ch, { rng: rs, sigmaX2: 1, sigmaN2: sn2, Ncp, pilotSpacing: 1 }), eh = O.epsCpEstimate(sym), Y = O.demodulate(sym, eh), s = O.system(ch, eh, 1, sn2);
                    a.push(dB(O.optimal(s.H, s.R, Nr, 1).sinr) - dB(O.sinrOfWeights(O.pilotResidualWeights([{ sym, Y }], { Hforce: s.H }).G, s.H, s.R, Nr, 1)));
                }
                cells.push(U.rpad(`${U.f(U.mean(a), 3)} ± ${U.f(U.se(a), 3)}`, 18));
                checks.push(U.check(`T25a Nr = ${Nr}, SNR ${snr} dB: pilot-residual (true H, ${N} pilots) against Optimal`, `${U.f(U.mean(a), 3)} dB (SE ${U.f(U.se(a), 3)})`, 'mean loss < 1 dB', U.mean(a) < 1));
            }
            console.log(U.pad(Nr, 4) + cells.join(''));
        }
        // ---------------------------------------------------------------- T25b / T25c
        console.log(`\nT25b  pilot spacing 12, m = 1, gamma_rel = +10 dB, H_hat from the pilots: SINR gain of the pilot-residual weights over the estimator with the signal in R_yy [dB], ${draws} realisations`);
        console.log(U.pad('Nr', 4) + [30, 35].map(s => U.rpad(`SNR ${s} dB`, 20)).join(''));
        for (const Nr of NRS) {
            const r = P.computeNr(Nr, draws, { seed, snrs: SNRS, ms: [1, 8] }), cells = [];
            const [E, PR10, LW, PRg, PRc, PR0, PRn] = [P.METHODS[0], P.METHODS[3], P.METHODS[4], P.METHODS[5], P.METHODS[6], P.METHODS[2], P.METHODS[1]];
            for (const snr of [30, 35]) {
                const i = SNRS.indexOf(snr), d = r.acc[i][E][0].map((x, k) => x - r.acc[i][PR10][0][k]);
                cells.push(U.rpad(`${U.f(U.mean(d), 2)} ± ${U.f(U.se(d), 2)}`, 20));
                checks.push(U.check(`T25b Nr = ${Nr}, SNR ${snr} dB: pilot-residual above Estimated (with signal)`, `${U.f(U.mean(d), 2)} dB (SE ${U.f(U.se(d), 2)})`, 'mean SINR gain >= 3 dB', U.mean(d) >= 3));
            }
            console.log(U.pad(Nr, 4) + cells.join(''));
            // T25c: controls and variants (information)
            SNRS.forEach((snr, i) => {
                const diff = (a, b, k = 0) => r.acc[i][a][k].map((x, j) => x - r.acc[i][b][k][j]), cell = d => `${U.f(U.mean(d), 3)} ± ${U.f(U.se(d), 3)}`;
                info.push({ name: `T25c Nr = ${Nr}, SNR ${snr} dB, m = 1: losses against Optimal [dB]`, value: `Estimated ${cell(r.acc[i][E][0])}; PR none ${cell(r.acc[i][PRn][0])}; PR 0 dB ${cell(r.acc[i][PR0][0])}; PR +10 dB ${cell(r.acc[i][PR10][0])}; PR Ledoit-Wolf ${cell(r.acc[i][LW][0])}; genie signal-free ${cell(r.acc[i][P.METHODS[7]][0])}; AD-MF ${cell(r.acc[i][P.METHODS[8]][0])}` });
                info.push({ name: `T25c Nr = ${Nr}, SNR ${snr} dB: controls at gamma_rel +10 dB (loss difference [dB], paired)`, value: `uncorrected - true H (i): ${cell(diff(PR10, PRg))}; uncorrected - corrected n/(n-1) (ii): ${cell(diff(PR10, PRc))}` });
                info.push({ name: `T25c Nr = ${Nr}, SNR ${snr} dB: loss of PR +10 dB and PR Ledoit-Wolf against m`, value: `PR +10 dB m = 1: ${cell(r.acc[i][PR10][0])}, m = 8: ${cell(r.acc[i][PR10][1])}; Ledoit-Wolf m = 1: ${cell(r.acc[i][LW][0])}, m = 8: ${cell(r.acc[i][LW][1])}; Estimated m = 1: ${cell(r.acc[i][E][0])}, m = 8: ${cell(r.acc[i][E][1])}` });
            });
        }
        return { id: this.id, title: this.title, checks, info };
    }
};
