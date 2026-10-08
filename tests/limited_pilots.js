'use strict';
/**
 * C1 導頻樣本數有限（資訊性；T27 引用這裡的函式做 PASS/FAIL）。   node tests/limited_pilots.js [--draws 500]
 * 導頻殘差協方差用導頻子載波的子集（P_eff ∈ {2, 4, 8, 16, 32, 85}；每個符元有 86 個導頻，間距 12：從中均勻取子集，P_eff = 85 時丟掉其中一個），m = 1，
 * Nr ∈ {2, 4, 8}，SNR 20 與 30 dB，每格 500 次實現，隨機到達角，文獻參數，CP 法估計 ε̂（所有權重在同一個 ε̂ 下與真實 H、R 比較）。
 * 方法（損失 = Optimal − 該方法，dB，成對）：
 *   Estimated（含訊號）  R_yy 取同一個子集的 P_eff 個子載波（含訊號），Ĥ 由子集導頻估計
 *   PR none / 0 dB / +10 dB  導頻殘差 R_res（子集），G = (R_res + γ I)⁻¹ Ĥ，γ = γ_rel σn²
 *   PR Ledoit–Wolf       R_res 的 Ledoit–Wolf 收縮
 *   genie signal-free    殘差與權重都用真實 H、無加載（只含 ICI 加雜訊的樣本協方差，同樣只有 P_eff 個樣本）
 * P_eff < Nr 時協方差的秩 ≤ P_eff < Nr：沒有加載或收縮的方法用偽逆（Moore–Penrose，只用高於 epsRank·λmax 的特徵值，epsRank = 1e-10，與 core.js 相同）。
 * 回報：加載或收縮開始優於無加載的 P_eff（每個 Nr、SNR：最佳的加載／收縮版本相對於無加載的成對損失差 > 2 SE 的最大 P_eff；判法事先訂定）。
 * 輸出 data/limited_pilots.csv。
 */
const fs = require('fs'), path = require('path');
const O = require('../ofdm_ici.js');
const U = require('./_util.js');
const R = require('./reproduce_gopala_slock.js');

const N = R.CFG.N, Ncp = R.CFG.Ncp, PEFF = [2, 4, 8, 16, 32, 85], NRS = [2, 4, 8], SNRS = [20, 30], dB = x => 10 * Math.log10(x);
const METHODS = ['Estimated (signal in R_yy)', 'PR none', 'PR gamma_rel 0 dB', 'PR gamma_rel +10 dB', 'PR Ledoit-Wolf', 'genie signal-free (true H)'];

function subsetOf(pilots, P) { const out = []; for (let j = 0; j < P; j++) out.push(pilots[Math.floor(j * pilots.length / P)]); return out; }
function evalDraw(ch, snr, rs, o = {}) {
    const Nr = ch.Nr, sn2 = O.snrToNoise(ch, snr, 1), peffs = o.peffs || PEFF, sym = O.makeSymbol(ch, { rng: rs, sigmaX2: 1, sigmaN2: sn2, Ncp, pilotSpacing: 12 });
    const eh = O.epsCpEstimate(sym), Y = O.demodulate(sym, eh), sys = O.system(ch, eh, 1, sn2), optDb = dB(O.optimal(sys.H, sys.R, Nr, 1).sinr), loss = {}; for (const m of METHODS) loss[m] = [];
    for (const P of peffs) {
        const sub = subsetOf(sym.pilots, P), sym2 = Object.assign({}, sym, { pilots: sub }), items = [{ sym: sym2, Y }], ev = G => optDb - dB(O.sinrOfWeights(G, sys.H, sys.R, Nr, 1));
        const pr0 = O.pilotResidualWeights(items, {});
        const Ryy = new Float64Array(2 * Nr * Nr), y = { re: new Float64Array(Nr), im: new Float64Array(Nr) };
        for (const l of sub) { for (let p = 0; p < Nr; p++) { y.re[p] = Y[p].re[l]; y.im[p] = Y[p].im[l]; } O.addOuter(Ryy, y, 1 / P, Nr); }
        let Ge = P < Nr ? null : O.csolve(Ryy, pr0.Hhat, Nr); if (!Ge) Ge = O.pinvSolveHermitian(Ryy, pr0.Hhat, Nr, 1e-10);
        loss[METHODS[0]].push(ev(Ge)); loss[METHODS[1]].push(ev(pr0.G)); loss[METHODS[2]].push(ev(O.pilotResidualWeights(items, { gammaRel: 1 }).G)); loss[METHODS[3]].push(ev(O.pilotResidualWeights(items, { gammaRel: 10 }).G));
        loss[METHODS[4]].push(ev(O.pilotResidualWeights(items, { shrink: 'lw' }).G)); loss[METHODS[5]].push(ev(O.pilotResidualWeights(items, { Hforce: sys.H }).G));
    }
    return { opt: optDb, loss };
}
function computeNr(Nr, draws, o = {}) {
    const snrs = o.snrs || SNRS, peffs = o.peffs || PEFF, base = o.seed || R.CFG.seed, acc = snrs.map(() => { const a = {}; for (const m of METHODS) a[m] = peffs.map(() => []); return a; });
    for (let t = 0; t < draws; t++) {
        const rp = O.makeRng(base + 7919 * t), ra = O.makeRng(base + 500000 + 7919 * t), rs = O.makeRng(base + 1000000 + 7919 * t), ch = O.paperChannel(R.pickAngles('random', ra), rp, { Nr, N });
        snrs.forEach((snr, i) => { const r = evalDraw(ch, snr, rs, { peffs }); for (const m of METHODS) r.loss[m].forEach((v, k) => acc[i][m][k].push(v)); });
    }
    return { acc, snrs, peffs };
}
/** 加載或收縮開始優於無加載的 P_eff：每個 P_eff 取三個加載／收縮版本中平均損失最小的一個，與無加載的成對差 > 2 SE 稱為「較好」；回傳較好的最大 P_eff（無則 null）與逐點表 */
function criticalPeff(r, si) {
    const none = r.acc[si][METHODS[1]], rows = [];
    r.peffs.forEach((P, k) => {
        const cand = [2, 3, 4].map(j => ({ name: METHODS[j], a: r.acc[si][METHODS[j]][k] })).sort((x, y) => U.mean(x.a) - U.mean(y.a))[0], d = none[k].map((x, i) => x - cand.a[i]);
        rows.push({ P, best: cand.name, diff: U.mean(d), se: U.se(d), better: U.mean(d) > 2 * U.se(d) });
    });
    const better = rows.filter(x => x.better).map(x => x.P);
    return { rows, critical: better.length ? Math.max(...better) : null, allBetter: rows.every(x => x.better) };
}

if (require.main === module) {
    const argv = process.argv.slice(2), draws = argv.includes('--draws') ? +argv[argv.indexOf('--draws') + 1] : 500, rowsCsv = [], t0 = Date.now();
    for (const Nr of NRS) {
        const r = computeNr(Nr, draws), cell = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`;
        r.snrs.forEach((snr, si) => {
            console.log(`\nNr = ${Nr}, SNR ${snr} dB, ${draws} realisations (${((Date.now() - t0) / 1000).toFixed(0)} s): loss against Optimal [dB]`);
            console.log(U.pad('method', 30) + r.peffs.map(P => U.rpad(`P_eff = ${P}`, 15)).join(''));
            for (const m of METHODS) { console.log(U.pad(m, 30) + r.peffs.map((P, k) => U.rpad(cell(r.acc[si][m][k]), 15)).join('')); r.peffs.forEach((P, k) => rowsCsv.push([Nr, snr, P, m, U.mean(r.acc[si][m][k]).toFixed(4), U.se(r.acc[si][m][k]).toFixed(4), draws])); }
            const c = criticalPeff(r, si);
            console.log('loaded / shrunk (best of PR 0 dB, +10 dB, Ledoit-Wolf) against PR none, paired difference of the loss (positive = loading better): ' + c.rows.map(x => `P_eff ${x.P}: ${U.f(x.diff, 2)} ± ${U.f(x.se, 2)} (${x.best.replace('PR ', '')})${x.better ? ' *' : ''}`).join('; '));
            console.log(`  -> loading / shrinkage is better than no loading (> 2 SE) up to P_eff = ${c.critical === null ? 'none' : c.critical}${c.allBetter ? ' (at every P_eff of the grid)' : ''}`);
        });
    }
    const hdr = ['Limited number of pilot samples for the pilot-residual covariance: loss against Optimal [dB]', `schema_version 1; generated by tests/limited_pilots.js; seed ${R.CFG.seed}; ${draws} realisations per row; paired; m = 1; random angles; epsilon from the CP method; H_hat from the same subset of pilots`,
        'P_eff = number of pilot subcarriers used (uniform subset of the 86 pilots of a symbol; 85 drops one); P_eff < Nr: rank-deficient covariance, methods without loading or shrinkage use the pseudo-inverse (relative threshold 1e-10)',
        'Estimated = sample R_yy over the same P_eff subcarriers (signal included); PR = pilot residual; gamma = gamma_rel sigma_n^2; genie signal-free = true H in residual and weights'];
    fs.writeFileSync(path.join(__dirname, '..', 'data', 'limited_pilots.csv'), hdr.map(l => '# ' + l).join('\n') + '\nNr,snr_db,p_eff,method,mean_loss_db,se_db,n\n' + rowsCsv.map(x => x.join(',')).join('\n') + '\n');
}
module.exports = { METHODS, PEFF, NRS, SNRS, evalDraw, computeNr, criticalPeff };
