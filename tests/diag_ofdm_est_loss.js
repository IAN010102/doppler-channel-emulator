'use strict';
/**
 * 診斷（資訊性，不進 run_all；不改任何既有行為）：Estimated BF 的損失來源。   node tests/diag_ofdm_est_loss.js [draws]
 * 設定同 Fig3 類：Nr = 2、N = 1024、文獻參數（4 路徑、都卜勒 [1080, −1080, 758, 220] Hz、相對強度 [0, 0, −11, −0.7] dB、導頻間距 12）、
 * 到達角每次 U(−90°, 90°) 隨機、CP 法由含雜訊樣本估計 ε̂（所有權重在同一個 ε̂ 下，與「真實」H、R 比較），SNR = 10、20、30、35 dB，≥ 200 次實現。
 * 六種權重（都用同一個通道、同一個符元；成對）：
 *   (1) Optimal                      G = R⁻¹ H
 *   (2) Estimated（現行）            樣本 R_yy（含訊號）＋導頻 Ĥ
 *   (3) 樣本 R_yy（含訊號）＋真實 H
 *   (4) 母體 R_yy = σx² H Hᴴ + R（解析式，不抽樣）＋ Ĥ
 *   (5) 母體 R_yy ＋真實 H            理論上與 (1) 完全相同（矩陣反轉引理）：差 < 1e-9 dB，否則是實作錯誤
 *   (6) 不含訊號的樣本協方差 R_u = (1/N) Σ (Y_l − X_l Ĥ)(Y_l − X_l Ĥ)ᴴ（扣除訊號用導頻 Ĥ，只含 ICI 加雜訊）＋ Ĥ
 * 樣本量掃描：m 個獨立 OFDM 符元（同一個通道、各自的資料與雜訊；ε̂ 取第一個符元的 CP 估計，m 個符元都用它對頻；Ĥ 對 m 個符元的所有導頻平均）
 *   的 R_yy 平均，樣本數 1024·m，m = 1、4、16、64：輸出 (2) 與 (3) 對 Optimal 的 SINR 損失（dB）隨 m 的變化，附標準誤。
 */
const fs = require('fs'), path = require('path');
const O = require('../ofdm_ici.js'), U = require('./_util.js'), R = require('./reproduce_gopala_slock.js');
const draws = +(process.argv[2] || 200), snrs = [10, 20, 30, 35], ms = [1, 4, 16, 64], base = R.CFG.seed, N = 1024, Nr = 2, Ncp = R.CFG.Ncp, dB = x => 10 * Math.log10(x);
const out = [];
const log = s => { out.push(s); console.log(s); };

function outerAdd(M, v, w) { for (let i = 0; i < Nr; i++) for (let j = 0; j < Nr; j++) { M[2 * (i * Nr + j)] += w * (v.re[i] * v.re[j] + v.im[i] * v.im[j]); M[2 * (i * Nr + j) + 1] += w * (v.im[i] * v.re[j] - v.re[i] * v.im[j]); } }
const names = ['(1) Optimal', '(2) sample R_yy + H_hat (current)', '(3) sample R_yy + true H', '(4) population R_yy + H_hat', '(5) population R_yy + true H', '(6) signal-free sample R_u + H_hat'];
const loss = names.map(() => snrs.map(() => [])), optSinr = snrs.map(() => []), chk5 = []; const sweep = {}; for (const k of ['2', '3']) sweep[k] = snrs.map(() => ms.map(() => []));
const sinrIn = snrs.map(() => []);                      // 損失與輸出 SINR 的關係：每個實現的 Optimal SINR（dB）與 (3) 的損失

for (let t = 0; t < draws; t++) {
    const rp = O.makeRng(base + 7919 * t), ra = O.makeRng(base + 500000 + 7919 * t), rs = O.makeRng(base + 1000000 + 7919 * t), rs2 = O.makeRng(base + 3000000 + 7919 * t);
    const ch = O.paperChannel(R.pickAngles('random', ra), rp, { Nr, N });
    snrs.forEach((snr, si) => {
        const sn2 = O.snrToNoise(ch, snr, 1), sym = O.makeSymbol(ch, { rng: rs, sigmaX2: 1, sigmaN2: sn2, Ncp, pilotSpacing: 12 }), eh = O.epsCpEstimate(sym), Y = O.demodulate(sym, eh), s = O.system(ch, eh, 1, sn2), H = s.H;
        const ev = G => dB(O.sinrOfWeights(G, H, s.R, Nr, 1)), opt = O.optimal(H, s.R, Nr, 1), optDb = dB(opt.sinr);
        const pw = O.practicalWeights(sym, Y, { Ns: N });
        const Rpop = Float64Array.from(s.R); outerAdd(Rpop, H, 1);                                                     // R + σx² H Hᴴ
        const Ru = new Float64Array(2 * Nr * Nr), e = { re: new Float64Array(Nr), im: new Float64Array(Nr) };
        for (let l = 0; l < N; l++) { for (let p = 0; p < Nr; p++) { e.re[p] = Y[p].re[l] - (sym.Xre[l] * pw.Hhat.re[p] - sym.Xim[l] * pw.Hhat.im[p]); e.im[p] = Y[p].im[l] - (sym.Xre[l] * pw.Hhat.im[p] + sym.Xim[l] * pw.Hhat.re[p]); } outerAdd(Ru, e, 1 / N); }
        const w = [opt.G, pw.G, O.csolve(pw.Ryy, H, Nr), O.csolve(Rpop, pw.Hhat, Nr), O.csolve(Rpop, H, Nr), O.csolve(Ru, pw.Hhat, Nr)];
        w.forEach((G, k) => loss[k][si].push(optDb - ev(G))); optSinr[si].push(optDb); chk5.push(Math.abs(optDb - ev(w[4])));
        sinrIn[si].push([optDb, optDb - ev(w[2])]);
        // ---- 樣本量掃描：同一個通道、m 個符元（第一個就是上面的符元）
        const Rsum = new Float64Array(2 * Nr * Nr), Hsum = { re: new Float64Array(Nr), im: new Float64Array(Nr) }; let nPilot = 0, cnt = 0;
        const addSymbol = (sy, YY) => {
            const y = { re: new Float64Array(Nr), im: new Float64Array(Nr) };
            for (let l = 0; l < N; l++) { for (let p = 0; p < Nr; p++) { y.re[p] = YY[p].re[l]; y.im[p] = YY[p].im[l]; } outerAdd(Rsum, y, 1); }
            for (const l of sy.pilots) { const d = sy.Xre[l] * sy.Xre[l] + sy.Xim[l] * sy.Xim[l]; for (let p = 0; p < Nr; p++) { Hsum.re[p] += (YY[p].re[l] * sy.Xre[l] + YY[p].im[l] * sy.Xim[l]) / d; Hsum.im[p] += (YY[p].im[l] * sy.Xre[l] - YY[p].re[l] * sy.Xim[l]) / d; } nPilot++; }
            cnt++;
        };
        addSymbol(sym, Y);
        for (let k = 1; k <= ms[ms.length - 1]; k++) {
            if (k > 1) { const sy = O.makeSymbol(ch, { rng: rs2, sigmaX2: 1, sigmaN2: sn2, Ncp, pilotSpacing: 12 }); addSymbol(sy, O.demodulate(sy, eh)); }
            const mi = ms.indexOf(k);
            if (mi >= 0) {
                const Rm = Float64Array.from(Rsum).map(x => x / (N * cnt)), Hm = { re: Hsum.re.map(x => x / nPilot), im: Hsum.im.map(x => x / nPilot) };
                sweep['2'][si][mi].push(optDb - ev(O.csolve(Rm, Hm, Nr))); sweep['3'][si][mi].push(optDb - ev(O.csolve(Rm, H, Nr)));
            }
        }
    });
}
const cell = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`;
log(`Estimated BF diagnosis: loss against Optimal BF [dB] (positive = worse), N = 1024, Nr = 2, random angles, ${draws} realisations, eps_hat from the CP method, mean ± SE (paired)`);
log(U.pad('weights', 38) + snrs.map((s, i) => U.rpad(`SNR ${s} dB`, 16)).join(''));
log(U.pad('mean SINR_opt [dB]', 38) + snrs.map((s, i) => U.rpad(U.f(U.mean(optSinr[i]), 2), 16)).join(''));
names.forEach((n, k) => log(U.pad(n, 38) + snrs.map((s, i) => U.rpad(cell(loss[k][i]), 16)).join('')));
const mx = Math.max(...chk5); log(`\ncheck (5) against (1): largest |difference| over all ${chk5.length} cases = ${U.e(mx)} dB  ->  ${mx < 1e-9 ? 'PASS (< 1e-9 dB)' : 'FAIL: possible implementation problem'}`);
log(`\nSample-size sweep: m independent OFDM symbols (1024 m samples for R_yy), loss against Optimal [dB], mean ± SE`);
log(U.pad('weights', 28) + U.pad('SNR', 8) + ms.map(m => U.rpad(`m = ${m}`, 16)).join(''));
for (const k of ['2', '3']) snrs.forEach((s, si) => log(U.pad(k === '2' ? '(2) sample R_yy + H_hat' : '(3) sample R_yy + true H', 28) + U.pad(s + ' dB', 8) + ms.map((m, mi) => U.rpad(cell(sweep[k][si][mi]), 16)).join('')));
log('\nloss of (3) against the ratio SINR_opt / (N m)  (the self-nulling loss of a sample matrix with the signal grows with this ratio):');
log(U.pad('SNR', 8) + U.rpad('mean SINR_opt [dB]', 20) + ms.map(m => U.rpad(`m = ${m}: ratio -> loss`, 24)).join(''));
snrs.forEach((s, si) => log(U.pad(s + ' dB', 8) + U.rpad(U.f(U.mean(optSinr[si]), 2), 20) + ms.map((m, mi) => U.rpad(`${U.f(U.mean(sinrIn[si].map(x => Math.pow(10, x[0] / 10))) / (N * m), 3)} -> ${U.f(U.mean(sweep['3'][si][mi]), 2)}`, 24)).join('')));
fs.writeFileSync(path.join(__dirname, '..', 'docs', 'diagnostics', 'ofdm_est_loss_detail.txt'), out.join('\n') + '\n');
