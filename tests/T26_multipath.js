'use strict';
/**
 * T26  頻率選擇性（多路徑）模型的驗證（ofdm_ici_multipath.js）。判準在看結果之前訂定，不調整。
 *   T26a  確定性：無雜訊時，時域模擬（含非整數延遲 15.36 與 38.4 個取樣點）的 FFT 輸出 Y_l 等於解析式 Σ_m h(l,m) X_m（M1），
 *         在 3 個子載波、3 根天線、有延遲與延遲為 0 兩種通道下，最大誤差 / 最大 |Y| < 1e-9。
 *   T26b  延遲全為 0 時，mpSystem 的 H_l 與 B_l 等於 ofdm_ici.js 既有的 vecAt(iciVectors) 與 iciCov（相對誤差 < 1e-9）。
 *   T26c  最大延遲（2.5 μs = 38.4 個取樣點）小於 CP（72 個取樣點）。
 *   T26d  統計：固定通道、隨機資料與雜訊的 200 個符元，子載波 l 的殘差 Y_l − H_l X_l 的平均功率（對天線加總）
 *         等於 trace(B_l) + Nr σn²（解析），z = (模擬平均 − 解析) / SE，|z| < 4（3 個子載波各一個檢查；SNR 20 dB、Nr = 2）。
 *   T26e  （資訊性）延遲帶來的通道頻率選擇性：相鄰導頻（間距 12 個子載波）的 |H_{k+1} − H_k|² / |H_k|² 的平均，平坦與多路徑各一個值。
 */
const O = require('../ofdm_ici.js');
const M = require('../ofdm_ici_multipath.js');
const U = require('./_util.js');
const R = require('./reproduce_gopala_slock.js');

const N = R.CFG.N, Ncp = R.CFG.Ncp;
module.exports = {
    id: 'T26', title: 'multipath (frequency-selective) ICI model: time-domain simulation equals the analytic h(l,m); reduces to the flat model at zero delay; ICI + noise power matches',
    async run({ draws = 200, seed = 2626 } = {}) {
        const checks = [], info = [];
        // ---------------------------------------------------------------- T26a / T26b / T26c
        let worstRel = 0, worstH = 0, worstB = 0;
        for (const flat of [false, true]) {
            const Nr = 3, ch = M.multipathChannel(R.pickAngles('random', O.makeRng(seed)), O.makeRng(seed + 1), { Nr, N, flat }), eps = 0.03, ls = [7, 300, 901];
            const sym = M.makeSymbolMP(ch, { rng: O.makeRng(seed + 2), sigmaX2: 1, sigmaN2: 0, Ncp }), Y = O.demodulate(sym, eps), sys = M.mpSystem(ch, eps, ls, true), q = [0, 0];
            let worst = 0, scale = 0;
            for (const s of sys) for (let p = 0; p < Nr; p++) {
                let re = 0, im = 0;
                for (let m = 0; m < N; m++) {
                    let hr = 0, hi = 0;
                    for (const pa of ch.paths) { O.Qc(((m - s.l) % N + N) % N + pa.eps - eps, N, q); const ph = -2 * Math.PI * m * pa.d / N, c = Math.cos(ph), sn = Math.sin(ph), wr = q[0] * c - q[1] * sn, wi = q[0] * sn + q[1] * c, cr = pa.Are * wr - pa.Aim * wi, ci = pa.Are * wi + pa.Aim * wr; hr += cr * pa.a.re[p] - ci * pa.a.im[p]; hi += cr * pa.a.im[p] + ci * pa.a.re[p]; }
                    re += hr * sym.Xre[m] - hi * sym.Xim[m]; im += hr * sym.Xim[m] + hi * sym.Xre[m];
                }
                worst = Math.max(worst, Math.hypot(re - Y[p].re[s.l], im - Y[p].im[s.l])); scale = Math.max(scale, Math.hypot(Y[p].re[s.l], Y[p].im[s.l]));
            }
            worstRel = Math.max(worstRel, worst / scale);
            if (flat) {      // T26b：延遲為 0 → 與既有的平坦模型相同
                const V = O.iciVectors(ch, eps);
                for (const s of sys) {
                    const H0 = O.vecAt(V, 0), B0 = O.iciCov(V, null, s.l);
                    worstH = Math.max(worstH, Math.hypot(...Array.from(H0.re, (x, i) => x - s.H.re[i]), ...Array.from(H0.im, (x, i) => x - s.H.im[i])) / Math.hypot(...H0.re, ...H0.im));
                    let dn = 0, bn = 0; for (let i = 0; i < B0.length; i++) { dn += (B0[i] - s.B[i]) ** 2; bn += B0[i] ** 2; } worstB = Math.max(worstB, Math.sqrt(dn / bn));
                }
            }
        }
        checks.push(U.check('T26a time-domain simulation (fractional delays) equals sum_m h(l,m) X_m (noise-free)', `largest relative error ${U.e(worstRel)}`, '< 1e-9', worstRel < 1e-9));
        checks.push(U.check('T26b zero delay: H_l equals vecAt(iciVectors) of the flat model', `relative error ${U.e(worstH)}`, '< 1e-9', worstH < 1e-9));
        checks.push(U.check('T26b zero delay: B_l equals iciCov of the flat model', `relative error ${U.e(worstB)}`, '< 1e-9', worstB < 1e-9));
        const chm = M.multipathChannel(R.pickAngles('random', O.makeRng(seed)), O.makeRng(seed + 1), { Nr: 2, N });
        checks.push(U.check('T26c largest delay shorter than the CP', `${U.f(chm.maxDelaySamples, 2)} samples (CP ${Ncp})`, '< CP', chm.maxDelaySamples < Ncp));
        // ---------------------------------------------------------------- T26d
        {
            const Nr = 2, ch = chm, sn2 = O.snrToNoise(ch, 20, 1), eps = 0.02, ls = [100, 500, 850], sys = M.mpSystem(ch, eps, ls, true), u = ls.map(() => []);
            for (let t = 0; t < draws; t++) {
                const sym = M.makeSymbolMP(ch, { rng: O.makeRng(seed + 100 + 7919 * t), sigmaX2: 1, sigmaN2: sn2, Ncp }), Y = O.demodulate(sym, eps);
                sys.forEach((s, k) => { let pw = 0; for (let p = 0; p < Nr; p++) { const l = s.l, rr = Y[p].re[l] - (sym.Xre[l] * s.H.re[p] - sym.Xim[l] * s.H.im[p]), ri = Y[p].im[l] - (sym.Xre[l] * s.H.im[p] + sym.Xim[l] * s.H.re[p]); pw += rr * rr + ri * ri; } u[k].push(pw); });
            }
            sys.forEach((s, k) => {
                let tr = 0; for (let p = 0; p < Nr; p++) tr += s.B[2 * (p * Nr + p)];
                const ref = tr + Nr * sn2, z = (U.mean(u[k]) - ref) / U.se(u[k]);
                checks.push(U.check(`T26d subcarrier ${s.l}: ICI + noise power (${draws} symbols) against trace(B_l) + Nr sigma_n^2`, `simulated ${U.e(U.mean(u[k]))}, analytic ${U.e(ref)}, z = ${U.f(z, 2)}`, '|z| < 4', Math.abs(z) < 4));
            });
        }
        // ---------------------------------------------------------------- T26e（資訊性）
        for (const flat of [true, false]) {
            const a = [];
            for (let t = 0; t < 50; t++) {
                const ch = M.multipathChannel(R.pickAngles('random', O.makeRng(seed + 7919 * t)), O.makeRng(seed + 500000 + t), { Nr: 2, N, flat }), Hs = M.mpSystem(ch, 0, [300, 312], false).map(x => x.H);
                let dn = 0, hn = 0; for (let p = 0; p < 2; p++) { dn += (Hs[1].re[p] - Hs[0].re[p]) ** 2 + (Hs[1].im[p] - Hs[0].im[p]) ** 2; hn += Hs[0].re[p] ** 2 + Hs[0].im[p] ** 2; } a.push(dn / hn);
            }
            info.push({ name: `T26e ${flat ? 'flat (delay 0)' : 'multipath'}: |H(l+12) - H(l)|^2 / |H(l)|^2 between neighbouring pilots (mean over 50 channels)`, value: `${U.f(U.mean(a), 3)} (SE ${U.f(U.se(a), 3)})` });
        }
        return { id: this.id, title: this.title, checks, info };
    }
};
