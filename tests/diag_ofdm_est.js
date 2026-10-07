'use strict';
/**
 * 診斷（資訊性）：為什麼 Estimated BF（R_yy 含訊號、Ĥ 由導頻估計）在高 SNR 比 Optimal BF 低很多？  node tests/diag_ofdm_est.js [draws]
 * 隨機到達角、Nr = 2、N = 1024、CP 法估計對頻（與 reproduce_gopala_slock.js 相同的通道與種子）。比較（成對、dB，平均 ± SE）：
 *   est            R_yy（含訊號）、Ĥ（導頻）                     —— 文獻的實務接收機
 *   est, true H    R_yy（含訊號）、真實 H
 *   est, R_u       訊號已扣除的樣本協方差 R_u = (1/N) Σ (Y_l − X_l Ĥ)(…)ᴴ、Ĥ（導頻）
 *   est, R_u, true H  R_u（以真實 H 扣除）、真實 H
 * 若 est 與 est,true H 的差主要來自「R_yy 含訊號」（MPDR 自我抵消 + 樣本數有限），則 R_u 版本應接近 Optimal。
 */
const O = require('../ofdm_ici.js'), U = require('./_util.js'), R = require('./reproduce_gopala_slock.js');
const draws = +(process.argv[2] || 200), snrs = [10, 20, 30, 35], base = R.CFG.seed, N = 1024, Nr = 2, Ncp = R.CFG.Ncp, dB = x => 10 * Math.log10(x);
const acc = {}; for (const k of ['est', 'estTrue', 'estRu', 'estRuTrue']) acc[k] = snrs.map(() => []);
for (let t = 0; t < draws; t++) {
    const rp = O.makeRng(base + 7919 * t), ra = O.makeRng(base + 500000 + 7919 * t), rs = O.makeRng(base + 1000000 + 7919 * t), ch = O.paperChannel(R.pickAngles('random', ra), rp, { Nr, N });
    snrs.forEach((snr, i) => {
        const sn2 = O.snrToNoise(ch, snr, 1), sym = O.makeSymbol(ch, { rng: rs, sigmaX2: 1, sigmaN2: sn2, Ncp, pilotSpacing: 12 }), eh = O.epsCpEstimate(sym), Y = O.demodulate(sym, eh), pw = O.practicalWeights(sym, Y, { Ns: N }), s2 = O.system(ch, eh, 1, sn2);
        const ev = G => dB(O.sinrOfWeights(G, s2.H, s2.R, Nr, 1)), opt = dB(O.optimal(s2.H, s2.R, Nr, 1).sinr);
        acc.est[i].push(ev(pw.G) - opt); acc.estTrue[i].push(ev(O.csolve(pw.Ryy, s2.H, Nr)) - opt);
        const Ru = (Hx) => { const M = new Float64Array(2 * Nr * Nr); for (let l = 0; l < N; l++) { const e = { re: new Float64Array(Nr), im: new Float64Array(Nr) }; for (let p = 0; p < Nr; p++) { e.re[p] = Y[p].re[l] - (sym.Xre[l] * Hx.re[p] - sym.Xim[l] * Hx.im[p]); e.im[p] = Y[p].im[l] - (sym.Xre[l] * Hx.im[p] + sym.Xim[l] * Hx.re[p]); } for (let a = 0; a < Nr; a++) for (let b = 0; b < Nr; b++) { M[2 * (a * Nr + b)] += (e.re[a] * e.re[b] + e.im[a] * e.im[b]) / N; M[2 * (a * Nr + b) + 1] += (e.im[a] * e.re[b] - e.re[a] * e.im[b]) / N; } } return M; };
        acc.estRu[i].push(ev(O.csolve(Ru(pw.Hhat), pw.Hhat, Nr)) - opt); acc.estRuTrue[i].push(ev(O.csolve(Ru(s2.H), s2.H, Nr)) - opt);
    });
}
console.log(`Estimated BF minus Optimal BF [dB], random angles, ${draws} realisations, N = 1024, Nr = 2 (paired; mean ± SE)`);
console.log(U.pad('SNR', 5), U.rpad('est (R_yy, H_hat)', 19), U.rpad('est, true H', 15), U.rpad('est, R_u, H_hat', 17), U.rpad('est, R_u, true H', 18));
snrs.forEach((s, i) => console.log(U.pad(s, 5), ...['est', 'estTrue', 'estRu', 'estRuTrue'].map((k, j) => U.rpad(`${U.f(U.mean(acc[k][i]), 2)} ± ${U.f(U.se(acc[k][i]), 2)}`, [19, 15, 17, 18][j]))));
