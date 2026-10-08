'use strict';
/**
 * 選做（只回報，不下結論）：用保護頻帶子載波估計「ICI 加雜訊」協方差。   node tests/guard_band_study.js [draws]
 * N = 1024，每端 104 個保護子載波（占用 816 個），Nr = 4，文獻參數，隨機角度，CP 法 ε̂；SNR 20、30、35 dB；預設 200 次實現。
 * (a) 保護頻帶內的 ICI 功率（tr B_l / Nr，解析式）相對於頻帶內部（l = 512）的比例：保護頻帶只收到來自占用子載波的干擾，且離占用區越遠越小。
 * (b) 以保護頻帶子載波的樣本協方差 R_g = (1/|G|) Σ_{l∈G} Y_l Y_lᴴ（沒有所需訊號，只含 ICI 加雜訊）做權重 G = R_g⁻¹ Ĥ，在頻帶內部的子載波（l = 512，
 *     用解析的 R_l 與 H）算 SINR，對照該子載波的 Optimal，以及導頻殘差法（導頻間距 12、γ_rel = +10 dB）。損失 = Optimal − 該方法（dB，成對）。
 */
const fs = require('fs'), path = require('path');
const O = require('../ofdm_ici.js'), U = require('./_util.js'), R = require('./reproduce_gopala_slock.js');
const draws = +(process.argv[2] || 200), snrs = [20, 30, 35], N = 1024, Nr = 4, Ncp = R.CFG.Ncp, base = R.CFG.seed, dB = x => 10 * Math.log10(x), G = 104, occ = O.maskGuard(N, G, G), l0 = 512;
const out = [], log = s => { out.push(s); console.log(s); };
const guardIdx = []; for (let l = 0; l < N; l++) if (!occ[l]) guardIdx.push(l);
const ratio = { near: [], mid: [], deep: [] }, loss = { guard: snrs.map(() => []), pr: snrs.map(() => []) };
for (let t = 0; t < draws; t++) {
    const rp = O.makeRng(base + 7919 * t), ra = O.makeRng(base + 500000 + 7919 * t), rs = O.makeRng(base + 1000000 + 7919 * t), ch = O.paperChannel(R.pickAngles('random', ra), rp, { Nr, N });
    snrs.forEach((snr, si) => {
        const sn2 = O.snrToNoise(ch, snr, 1), sym = O.makeSymbol(ch, { rng: rs, sigmaX2: 1, sigmaN2: sn2, Ncp, pilotSpacing: 12, occ }), eh = O.epsCpEstimate(sym), Y = O.demodulate(sym, eh), V = O.iciVectors(ch, eh), H = O.vecAt(V, 0);
        const tr = B => { let s = 0; for (let p = 0; p < Nr; p++) s += B[2 * (p * Nr + p)]; return s / Nr; };
        const Bint = O.iciCov(V, occ, l0), Rint = O.rMatrix(Bint, Nr, 1, sn2), opt = dB(O.optimal(H, Rint, Nr, 1).sinr);
        if (si === 0) { const t0 = tr(Bint); ratio.near.push(tr(O.iciCov(V, occ, G - 1)) / t0); ratio.mid.push(tr(O.iciCov(V, occ, 52)) / t0); ratio.deep.push(tr(O.iciCov(V, occ, 2)) / t0); }   // 保護頻帶最靠近占用區 / 中間 / 最外側
        // R_g：保護子載波的樣本協方差
        const Rg = new Float64Array(2 * Nr * Nr), y = { re: new Float64Array(Nr), im: new Float64Array(Nr) };
        for (const l of guardIdx) { for (let p = 0; p < Nr; p++) { y.re[p] = Y[p].re[l]; y.im[p] = Y[p].im[l]; } O.addOuter(Rg, y, 1 / guardIdx.length, Nr); }
        const pr = O.pilotResidualWeights([{ sym, Y }], { gammaRel: 10 }), Gg = O.csolve(Rg, pr.Hhat, Nr);
        loss.guard[si].push(opt - dB(O.sinrOfWeights(Gg, H, Rint, Nr, 1))); loss.pr[si].push(opt - dB(O.sinrOfWeights(pr.G, H, Rint, Nr, 1)));
    });
}
const cell = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`;
log(`Guard-band study (information only): N = 1024, ${G} guard subcarriers at each end (${N - 2 * G} occupied), Nr = 4, ${draws} realisations`);
log(`(a) ICI power tr(B_l)/Nr at a guard subcarrier relative to the band interior (l = ${l0}): nearest to the band (l = ${G - 1}): ${U.f(U.mean(ratio.near), 3)} ± ${U.f(U.se(ratio.near), 3)}; middle of the guard (l = 52): ${U.f(U.mean(ratio.mid), 3)} ± ${U.f(U.se(ratio.mid), 3)}; outermost (l = 2): ${U.f(U.mean(ratio.deep), 3)} ± ${U.f(U.se(ratio.deep), 3)}`);
log('(b) loss against the optimal weights at the interior subcarrier [dB] (mean ± SE, paired):');
log(U.pad('method', 44) + snrs.map(s => U.rpad(`SNR ${s} dB`, 16)).join(''));
log(U.pad('weights from the guard-band covariance R_g', 44) + snrs.map((s, i) => U.rpad(cell(loss.guard[i]), 16)).join(''));
log(U.pad('pilot residual (spacing 12; +10 dB loading)', 44) + snrs.map((s, i) => U.rpad(cell(loss.pr[i]), 16)).join(''));
fs.writeFileSync(path.join(__dirname, '..', 'docs', 'diagnostics', 'guard_band_study.txt'), out.join('\n') + '\n');
