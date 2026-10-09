'use strict';
/**
 * C2 多路徑穩健度（資訊性；T26 引用這裡的函式做驗證）。   node tests/multipath_pilots.js [--draws 500]
 * 通道：3 個分接，延遲 [0, 1, 2.5] μs（= 0、15.36、38.4 個取樣點，最大延遲 < CP 72 個取樣點 = 4.69 μs），分接功率 [0, −3, −6] dB（正規化總功率 1），
 *   每個分接帶文獻的 4 條都卜勒路徑（共 12 條子路徑，各自隨機相位；到達角每次 U(−90°, 90°)，4 個角度被 3 個分接共用）。假設 A20–A24。
 * 對照「平坦」：同樣 12 條子路徑、同樣的相位與角度，但延遲全為 0（成對、同一組資料符元與雜訊）。
 * 接收機（兩種通道用同一個流程）：原始導頻估計 Ĥ_k；子載波 l 的 Ĥ_l 由相鄰兩個導頻的原始估計線性內插；
 *   全域導頻殘差 d_k = Y_k − Ĥ_k^{LOO} X_k，Ĥ_k^{LOO} = 兩個相鄰導頻原始估計的平均（leave-one-out；排除第一與最後一個導頻）；G_l = (R_res + γI)⁻¹ Ĥ_l；
 *   ε̂ 由 CP 法估計（所有權重在同一個 ε̂ 下與真實 H_l、R_l 比較）；損失 = Optimal_l − 該方法，對 8 個評估子載波（兩個導頻之間）平均。
 * 內插誤差佔殘差的比例（導頻 j = 1…n−2）：e_j = H_j − Ĥ_j^{LOO}（H_j 為解析的真實值，genie 參考），Σ|e|² / Σ|d|²；
 *   其中 curvature = H_j − (H_{j−1} + H_{j+1})/2（只由通道的頻率選擇性造成，確定性），其餘 = e − curvature（相鄰導頻的雜訊加 ICI 經平均後的貢獻）。
 * 輸出 data/multipath_pilots.csv。
 */
const fs = require('fs'), path = require('path');
const O = require('../ofdm_ici.js');
const M = require('../ofdm_ici_multipath.js');
const U = require('./_util.js');
const R = require('./reproduce_gopala_slock.js');

const N = R.CFG.N, Ncp = R.CFG.Ncp, NRS = [2, 4, 8], SNRS = [20, 30], dB = x => 10 * Math.log10(x);
const METHODS = ['PR none', 'PR gamma_rel 0 dB', 'PR gamma_rel +10 dB', 'PR Ledoit-Wolf', 'PR 0 dB, constant H_hat (pilot mean)', 'PR 0 dB, genie H (true H_k and H_l)'];
const EVAL_LS = [0, 1, 2, 3, 4, 5, 6, 7].map(j => 12 * Math.round((j + 0.5) * N / 8 / 12) + 6);
const SHARES = ['share total', 'share curvature', 'share noise+ICI of neighbours'];

function addScaled(S, g, Nr) { const R2 = Float64Array.from(S); for (let p = 0; p < Nr; p++) R2[2 * (p * Nr + p)] += g; return R2; }
/** 一個通道實現、一個 SNR、一個通道型態：回傳 { loss: {method: 數值}, share: [三個比例] } */
function evalOne(ch, snr, rs) {
    const Nr = ch.Nr, sn2 = O.snrToNoise(ch, snr, 1), sym = M.makeSymbolMP(ch, { rng: rs, sigmaX2: 1, sigmaN2: sn2, Ncp, pilotSpacing: 12 }), eh = O.epsCpEstimate(sym), Y = O.demodulate(sym, eh);
    const sys = M.mpSystem(ch, eh, EVAL_LS, true), Hp = M.mpSystem(ch, eh, sym.pilots, false).map(x => x.H), raw = M.rawPilotH(sym, Y);
    const lo = M.residuals(sym, Y, raw), lg = M.residuals(sym, Y, raw, Hp);
    const meanH = { re: new Float64Array(Nr), im: new Float64Array(Nr) }; for (const h of raw) for (let p = 0; p < Nr; p++) { meanH.re[p] += h.re[p] / raw.length; meanH.im[p] += h.im[p] / raw.length; }
    // 常數 Ĥ 的殘差：d_k = Y_k − mean(Ĥ) X_k
    const S0 = new Float64Array(2 * Nr * Nr), dv = { re: new Float64Array(Nr), im: new Float64Array(Nr) };
    for (let j = 1; j < sym.pilots.length - 1; j++) { const l = sym.pilots[j]; for (let p = 0; p < Nr; p++) { dv.re[p] = Y[p].re[l] - (sym.Xre[l] * meanH.re[p] - sym.Xim[l] * meanH.im[p]); dv.im[p] = Y[p].im[l] - (sym.Xre[l] * meanH.im[p] + sym.Xim[l] * meanH.re[p]); } O.addOuter(S0, dv, 1 / lo.n, Nr); }
    const lw = O.ledoitWolfShrink(lo.vecs, lo.S, Nr).S, loss = {}; for (const m of METHODS) loss[m] = 0;
    for (const s of sys) {
        const R = O.rMatrix(s.B, Nr, 1, sn2), optDb = dB(O.optimal(s.H, R, Nr, 1).sinr), Hl = M.interp(raw, sym.pilots, s.l), ev = (Rm, Hh) => optDb - dB(O.sinrOfWeights(O.csolve(Rm, Hh, Nr), s.H, R, Nr, 1));
        loss[METHODS[0]] += ev(lo.S, Hl); loss[METHODS[1]] += ev(addScaled(lo.S, sn2, Nr), Hl); loss[METHODS[2]] += ev(addScaled(lo.S, 10 * sn2, Nr), Hl); loss[METHODS[3]] += ev(lw, Hl);
        loss[METHODS[4]] += ev(addScaled(S0, sn2, Nr), meanH); loss[METHODS[5]] += ev(addScaled(lg.S, sn2, Nr), s.H);
    }
    for (const m of METHODS) loss[m] /= sys.length;
    let sd = 0, se = 0, sc = 0, sr = 0;
    for (let j = 1; j < sym.pilots.length - 1; j++) {
        const l = sym.pilots[j], hl = M.lerpMid(raw[j - 1], raw[j + 1]), cv = M.lerpMid(Hp[j - 1], Hp[j + 1]);
        for (let p = 0; p < Nr; p++) {
            const er = Hp[j].re[p] - hl.re[p], ei = Hp[j].im[p] - hl.im[p], cr = Hp[j].re[p] - cv.re[p], ci = Hp[j].im[p] - cv.im[p];
            se += er * er + ei * ei; sc += cr * cr + ci * ci; sr += (er - cr) * (er - cr) + (ei - ci) * (ei - ci);
        }
    }
    sd = M.traceOf(lo.S, Nr) * lo.n;
    return { loss, share: [se / sd, sc / sd, sr / sd] };
}
/** Nr 一個值：回傳 acc[channel][snrIdx] = { method: [draws], share: [[...],[...],[...]] }，channel ∈ {multipath, flat} */
function computeNr(Nr, draws, o = {}) {
    const snrs = o.snrs || SNRS, base = o.seed || R.CFG.seed, mk = () => snrs.map(() => { const a = { share: SHARES.map(() => []) }; for (const m of METHODS) a[m] = []; return a; }), acc = { multipath: mk(), flat: mk() };
    for (let t = 0; t < draws; t++) {
        const ang = R.pickAngles('random', O.makeRng(base + 500000 + 7919 * t));
        for (const type of ['multipath', 'flat']) {
            const ch = M.multipathChannel(ang, O.makeRng(base + 7919 * t), { Nr, N, flat: type === 'flat' });
            snrs.forEach((snr, i) => { const r = evalOne(ch, snr, O.makeRng(base + 1000000 + 7919 * t)); for (const m of METHODS) acc[type][i][m].push(r.loss[m]); r.share.forEach((v, k) => acc[type][i].share[k].push(v)); });
        }
    }
    return { acc, snrs };
}

if (require.main === module) {
    const argv = process.argv.slice(2), draws = argv.includes('--draws') ? +argv[argv.indexOf('--draws') + 1] : 500, rows = [], t0 = Date.now(), cell = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`;
    for (const Nr of NRS) {
        const r = computeNr(Nr, draws);
        r.snrs.forEach((snr, i) => {
            console.log(`\nNr = ${Nr}, SNR ${snr} dB, ${draws} realisations (${((Date.now() - t0) / 1000).toFixed(0)} s): loss against Optimal [dB] (mean ± SE over realisations; paired flat vs multipath)`);
            console.log(U.pad('method', 42) + U.rpad('flat (delay 0)', 18) + U.rpad('multipath', 18) + U.rpad('multipath - flat', 20));
            for (const m of METHODS) {
                const a = r.acc.flat[i][m], b = r.acc.multipath[i][m], d = b.map((x, k) => x - a[k]);
                console.log(U.pad(m, 42) + U.rpad(cell(a), 18) + U.rpad(cell(b), 18) + U.rpad(cell(d), 20));
                rows.push([Nr, snr, 'flat', m, U.mean(a).toFixed(4), U.se(a).toFixed(4), draws], [Nr, snr, 'multipath', m, U.mean(b).toFixed(4), U.se(b).toFixed(4), draws]);
            }
            console.log(U.pad('interpolation-error share of the residual', 42) + SHARES.map((s, k) => `${s}: flat ${cell(r.acc.flat[i].share[k])}, multipath ${cell(r.acc.multipath[i].share[k])}`).join('; '));
            SHARES.forEach((s, k) => rows.push([Nr, snr, 'flat', s, U.mean(r.acc.flat[i].share[k]).toFixed(5), U.se(r.acc.flat[i].share[k]).toFixed(5), draws], [Nr, snr, 'multipath', s, U.mean(r.acc.multipath[i].share[k]).toFixed(5), U.se(r.acc.multipath[i].share[k]).toFixed(5), draws]));
        });
    }
    const hdr = ['Pilot-residual weights in a flat versus a frequency-selective (3-tap) channel; loss against Optimal per subcarrier [dB], mean over 8 evaluation subcarriers and the realisations',
        `schema_version 1; generated by tests/multipath_pilots.js; seed ${R.CFG.seed}; ${draws} realisations per row; paired (same phases, angles, data symbols and noise); epsilon from the CP method`,
        'multipath: taps at 0 / 1 / 2.5 us (0 / 15.36 / 38.4 samples, CP 72 samples), tap power 0 / -3 / -6 dB, each tap carries the 4 Doppler paths of the literature (12 sub-paths, independent phases); flat: same 12 sub-paths with all delays 0',
        'receiver: H_hat_l by linear interpolation over frequency of the raw pilot estimates; global residual covariance from leave-one-out residuals (H_hat_k = mean of the two neighbouring raw estimates); gamma = gamma_rel sigma_n^2; constant H_hat = mean over all pilots',
        'rows "share ...": interpolation-error power / residual power at the pilots (value column = mean share, not dB): total = |H - H_hat_LOO|^2, curvature = |H - mean(H of the two neighbours)|^2 (true channel, no noise), rest = noise + ICI of the neighbours'];
    fs.writeFileSync(path.join(__dirname, '..', 'data', 'multipath_pilots.csv'), hdr.map(l => '# ' + l).join('\n') + '\nNr,snr_db,channel,method,mean,se,n\n' + rows.map(x => x.join(',')).join('\n') + '\n');
}
module.exports = { METHODS, SHARES, NRS, SNRS, EVAL_LS, evalOne, computeNr };
