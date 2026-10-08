'use strict';
/**
 * 導頻殘差協方差估計器的比較（資訊性；T25 引用這裡的函式做 PASS/FAIL）。   node tests/pilot_residual_curves.js [--draws 500] [--nr 2,4,8]
 * 設定：Fig3 類（N = 1024、導頻間距 12、文獻參數、隨機到達角、CP 法由含雜訊樣本估計 ε̂；所有權重在同一個 ε̂ 下與真實 H、R 比較），Nr ∈ {2, 4, 8}，
 *   SNR 10、20、30、35 dB，每格 ≥ 500 次實現，m 個獨立 OFDM 符元（m ∈ {1, 2, 4, 8}，同一個通道，ε̂ 取第一個符元的 CP 估計）。
 * 方法（損失 = Optimal 的 SINR − 該方法的 SINR，dB，正值較差；成對）：
 *   Estimated              樣本 R_yy（含訊號，1024·m 個樣本）＋導頻 Ĥ（前一輪的現行做法）
 *   PR none / 0 dB / +10 dB 導頻殘差協方差 R_res，G = (R_res + γ I)⁻¹ Ĥ，γ = γ_rel σn²：無加載、γ_rel = 0 dB（γ = σn²）、+10 dB
 *   PR LW                  導頻殘差協方差的 Ledoit–Wolf 線性收縮（ofdm_ici.js 的 ledoitWolfShrink，公式與出處在那裡）
 *   PR +10 dB, true H      對照 (i)：殘差與權重都用真實 H（genie）
 *   PR +10 dB, corrected   對照 (ii)：殘差乘以 n/(n−1)（n = 導頻總數）
 *   genie signal-free      所有子載波、扣掉真實 H 訊號的樣本協方差（只含 ICI 加雜訊）＋ Ĥ
 *   AD-MF                  角度域匹配濾波（genie 基準，與 m 無關）
 * 輸出 data/pilot_residual_nr{Nr}.csv（mean 與 SE）。
 */
const fs = require('fs'), path = require('path');
const O = require('../ofdm_ici.js');
const U = require('./_util.js');
const R = require('./reproduce_gopala_slock.js');

const N = R.CFG.N, Ncp = R.CFG.Ncp, SNRS = [10, 20, 30, 35], MS = [1, 2, 4, 8], dB = x => 10 * Math.log10(x);
const METHODS = ['Estimated (sample R_yy with signal)', 'PR no loading', 'PR gamma_rel 0 dB', 'PR gamma_rel +10 dB', 'PR Ledoit-Wolf', 'PR +10 dB with true H (genie)', 'PR +10 dB corrected n/(n-1)', 'Genie signal-free sample covariance', 'AD-MF (genie)'];

/** 一個通道實現、一個 SNR、m 個符元：回傳 { opt (dB), loss: { method: [m 的各值] } } */
function evalDraw(ch, snr, rs, o = {}) {
    const Nr = ch.Nr, sn2 = O.snrToNoise(ch, snr, 1), ms = o.ms || MS, mmax = Math.max(...ms), pilotSpacing = o.pilotSpacing === undefined ? 12 : o.pilotSpacing;
    const items = []; let eh = 0, sys = null, optDb = 0;
    const loss = {}; for (const m of METHODS) loss[m] = [];
    const Rsum = new Float64Array(2 * Nr * Nr), Rusum = new Float64Array(2 * Nr * Nr), Hsum = { re: new Float64Array(Nr), im: new Float64Array(Nr) }; let nPil = 0, cnt = 0;
    for (let k = 1; k <= mmax; k++) {
        const sym = O.makeSymbol(ch, { rng: rs, sigmaX2: 1, sigmaN2: sn2, Ncp, pilotSpacing });
        if (k === 1) { eh = O.epsCpEstimate(sym); sys = O.system(ch, eh, 1, sn2); optDb = dB(O.optimal(sys.H, sys.R, Nr, 1).sinr); }
        const Y = O.demodulate(sym, eh); items.push({ sym, Y });
        const y = { re: new Float64Array(Nr), im: new Float64Array(Nr) }, e = { re: new Float64Array(Nr), im: new Float64Array(Nr) };
        for (let l = 0; l < N; l++) {
            for (let p = 0; p < Nr; p++) { y.re[p] = Y[p].re[l]; y.im[p] = Y[p].im[l]; e.re[p] = y.re[p] - (sym.Xre[l] * sys.H.re[p] - sym.Xim[l] * sys.H.im[p]); e.im[p] = y.im[p] - (sym.Xre[l] * sys.H.im[p] + sym.Xim[l] * sys.H.re[p]); }
            O.addOuter(Rsum, y, 1, Nr); O.addOuter(Rusum, e, 1, Nr);
        }
        cnt++;
        for (const l of sym.pilots) { const d = sym.Xre[l] * sym.Xre[l] + sym.Xim[l] * sym.Xim[l]; for (let p = 0; p < Nr; p++) { Hsum.re[p] += (Y[p].re[l] * sym.Xre[l] + Y[p].im[l] * sym.Xim[l]) / d; Hsum.im[p] += (Y[p].im[l] * sym.Xre[l] - Y[p].re[l] * sym.Xim[l]) / d; } nPil++; }
        const mi = ms.indexOf(k); if (mi < 0) continue;
        const ev = G => optDb - dB(O.sinrOfWeights(G, sys.H, sys.R, Nr, 1)), part = items.slice(0, k), Hm = { re: Hsum.re.map(x => x / nPil), im: Hsum.im.map(x => x / nPil) };
        loss[METHODS[0]].push(ev(O.csolve(Rsum.map(x => x / (N * cnt)), Hm, Nr)));
        loss[METHODS[1]].push(ev(O.pilotResidualWeights(part, {}).G)); loss[METHODS[2]].push(ev(O.pilotResidualWeights(part, { gammaRel: 1 }).G)); loss[METHODS[3]].push(ev(O.pilotResidualWeights(part, { gammaRel: 10 }).G));
        loss[METHODS[4]].push(ev(O.pilotResidualWeights(part, { shrink: 'lw' }).G));
        loss[METHODS[5]].push(ev(O.pilotResidualWeights(part, { gammaRel: 10, Hforce: sys.H }).G)); loss[METHODS[6]].push(ev(O.pilotResidualWeights(part, { gammaRel: 10, correct: true }).G));
        loss[METHODS[7]].push(ev(O.csolve(Rusum.map(x => x / (N * cnt)), Hm, Nr)));
    }
    const ad = dB(O.angleDomainMf(ch, 1, sn2).sinr); ms.forEach(() => loss[METHODS[8]].push(optDb - ad));
    return { opt: optDb, loss };
}
/** Nr 一個值：回傳 acc[snr][method][mIdx] = [draws 的損失]，optSinr[snr] = [] */
function computeNr(Nr, draws, o = {}) {
    const snrs = o.snrs || SNRS, ms = o.ms || MS, base = o.seed || R.CFG.seed, acc = snrs.map(() => { const a = {}; for (const m of METHODS) a[m] = ms.map(() => []); return a; }), opt = snrs.map(() => []);
    for (let t = 0; t < draws; t++) {
        const rp = O.makeRng(base + 7919 * t), ra = O.makeRng(base + 500000 + 7919 * t), rs = O.makeRng(base + 1000000 + 7919 * t), ch = O.paperChannel(R.pickAngles('random', ra), rp, { Nr, N });
        snrs.forEach((snr, i) => { const r = evalDraw(ch, snr, rs, o); opt[i].push(r.opt); for (const m of METHODS) r.loss[m].forEach((v, k) => acc[i][m][k].push(v)); });
    }
    return { acc, opt, snrs, ms };
}

if (require.main === module) {
    const argv = process.argv.slice(2), draws = argv.includes('--draws') ? +argv[argv.indexOf('--draws') + 1] : 500, nrs = argv.includes('--nr') ? argv[argv.indexOf('--nr') + 1].split(',').map(Number) : [2, 4, 8], t0 = Date.now();
    for (const Nr of nrs) {
        const r = computeNr(Nr, draws), rows = [], cell = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`;
        console.log(`\nNr = ${Nr}, ${draws} realisations (${((Date.now() - t0) / 1000).toFixed(0)} s): loss against Optimal [dB], m = 1`);
        console.log(U.pad('method', 40) + r.snrs.map(s => U.rpad(`SNR ${s} dB`, 16)).join(''));
        console.log(U.pad('mean SINR_opt [dB]', 40) + r.snrs.map((s, i) => U.rpad(U.f(U.mean(r.opt[i]), 2), 16)).join(''));
        for (const m of METHODS) console.log(U.pad(m, 40) + r.snrs.map((s, i) => U.rpad(cell(r.acc[i][m][0]), 16)).join(''));
        console.log(`loss against m (SNR 30 dB):`); console.log(U.pad('method', 40) + r.ms.map(m => U.rpad(`m = ${m}`, 16)).join(''));
        const i30 = r.snrs.indexOf(30); for (const m of METHODS) console.log(U.pad(m, 40) + r.ms.map((mm, k) => U.rpad(cell(r.acc[i30][m][k]), 16)).join(''));
        r.snrs.forEach((s, i) => { for (const m of METHODS) r.ms.forEach((mm, k) => { const a = r.acc[i][m][k]; rows.push([Nr, s, mm, m, U.mean(a).toFixed(4), U.se(a).toFixed(4), (U.mean(r.opt[i]) - U.mean(a)).toFixed(4), draws]); }); });
        const hdr = [`Pilot-residual covariance estimator versus the other weights; Nr = ${Nr}`, `schema_version 1; generated by tests/pilot_residual_curves.js; seed ${R.CFG.seed}; ${draws} realisations per row; paired; loss = SINR of the optimal weights (at the CP-estimated epsilon) minus the SINR of the method [dB]`,
            'N 1024, delta_f 15 kHz, pilot spacing 12 (86 pilots per symbol), Ncp 72, 4 paths (Doppler 1080 / -1080 / 758 / 220 Hz, relative power 0 / 0 / -11 / -0.7 dB), random phases and random angles U(-90, 90) deg (A3, A4), SNR as A1',
            'm independent OFDM symbols of the same channel; epsilon_hat from the first symbol (CP method, noisy samples); H_hat averaged over the pilots of all m symbols; PR = pilot residual R_res = (1/n) sum d d^H, d = Y_k - H_hat X_k; gamma = gamma_rel sigma_n^2',
            'PR Ledoit-Wolf: linear shrinkage of R_res toward mu I (Ledoit & Wolf 2004, J. Multivariate Anal. 88(2):365-411; complex-data version of the same formulas); mean_sinr_db = mean SINR of the method (= mean SINR_opt - mean loss)'];
        fs.writeFileSync(path.join(__dirname, '..', 'data', `pilot_residual_nr${Nr}.csv`), hdr.map(l => '# ' + l).join('\n') + '\nNr,snr_db,m,method,mean_loss_db,se_db,mean_sinr_db,n\n' + rows.map(x => x.join(',')).join('\n') + '\n');
    }
}
module.exports = { METHODS, SNRS, MS, evalDraw, computeNr };
