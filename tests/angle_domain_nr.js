'use strict';
/**
 * 角度域匹配濾波（AD-MF）基準與陣列大小 Nr 軸（資訊性，不判定）。   node tests/angle_domain_nr.js [--draws 200]
 * AD-MF（ofdm_ici.js 的 angleDomainMf）是 genie 基準：對每條路徑以真實到達角做空間匹配濾波、再以該路徑真實的都卜勒對頻、最後 MRC 合併；
 * 不做 DoA 估計、不做都卜勒估計、不考慮 ICI 與路徑間干擾（IPI）。它代表「大陣列」下角度域處理的上限型基準，不是可實作的接收機。
 * 比較對象：最佳自適應合併 Optimal（全陣列 G = R⁻¹ H，對頻用 CP 法的解析期望值）。
 * 設定：文獻設定（N = 1024、Δf = 15 kHz、4 條路徑、都卜勒 [1080, −1080, 758, 220] Hz、相對強度 [0, 0, −11, −0.7] dB），Nr ∈ {2, 4, 8, 16, 32}，
 *   到達角：固定一組 F1 = 0°、25°、−35°、55°（只抽隨機相位）與每次 U(−90°, 90°) 隨機各 200 次；SNR 20 與 30 dB。
 * EVM（假設 A13）= 100/√SINR %（把干擾當雜訊、不計增益偏差）。IPI = 其他路徑洩漏到分支 p 的功率相對於該路徑自己的訊號（對路徑與實現平均），以 dB 表示。
 * 輸出 data/angle_domain_nr.csv；並找出 AD-MF 與 Optimal 相交的 Nr（平均 SINR 相等處，在 Nr 的對數軸上線性內插；另以 Nr = 16…32 步長 2 細查）。
 * 這是在本模型下的結果（ICI 與 IPI 的模型、平坦通道、隨機角度的假設，見 docs/ofdm_ici_assumptions.md）。
 */
const fs = require('fs'), path = require('path');
const O = require('../ofdm_ici.js');
const U = require('./_util.js');
const R = require('./reproduce_gopala_slock.js');

const NRS = [2, 4, 8, 16, 32], SNRS = [20, 30], dB = x => 10 * Math.log10(x), Ncp = R.CFG.Ncp;
function point(mode, Nr, draws, seed = R.CFG.seed) {
    const res = SNRS.map(() => ({ opt: [], ad: [], ipi: [] }));
    for (let t = 0; t < draws; t++) {
        const rp = O.makeRng(seed + 7919 * t), ra = O.makeRng(seed + 500000 + 7919 * t), ch = O.paperChannel(R.pickAngles(mode, ra, 'F1'), rp, { Nr, N: R.CFG.N });
        const e = O.epsCpAnalytic(ch, Ncp), sys = O.system(ch, e, 1, 0);
        SNRS.forEach((snr, i) => {
            const sn2 = O.snrToNoise(ch, snr, 1), ad = O.angleDomainMf(ch, 1, sn2);
            res[i].opt.push(dB(O.optimal(sys.H, O.rMatrix(sys.B, Nr, 1, sn2), Nr, 1).sinr)); res[i].ad.push(dB(ad.sinr)); res[i].ipi.push(ad.ipiRatio);
        });
    }
    return res;
}
const cross = (xs, a, b) => {               // a, b: 平均 SINR（dB）對 xs（Nr）；回傳 a − b 由負變正的位置（對數軸線性內插），沒有則 null
    for (let i = 0; i + 1 < xs.length; i++) { const d0 = a[i] - b[i], d1 = a[i + 1] - b[i + 1]; if (d0 < 0 && d1 >= 0) { const f = -d0 / (d1 - d0); return Math.pow(2, Math.log2(xs[i]) + f * (Math.log2(xs[i + 1]) - Math.log2(xs[i]))); } }
    return null;
};

if (require.main === module) {
    const argv = process.argv.slice(2), draws = argv.includes('--draws') ? +argv[argv.indexOf('--draws') + 1] : 200, rows = [], summary = [];
    // 單一路徑的自我檢查：AD-MF 等於 Nr × SNR
    { let worst = 0; for (const Nr of [2, 8, 32]) { const ch = O.makeChannel({ N: 1024, Nr, paths: [{ eps: 0.07, thetaDeg: 20, powerDb: 0, phase: 0.4 }] }); worst = Math.max(worst, Math.abs(O.angleDomainMf(ch, 1, O.snrToNoise(ch, 20, 1)).sinr / (Nr * 100) - 1)); } console.log(`self-check: single path, AD-MF SINR = Nr x SNR, largest relative difference ${U.e(worst)}  ->  ${worst < 1e-9 ? 'PASS' : 'FAIL'}`); }
    for (const mode of ['fixed', 'random']) {
        const per = NRS.map(Nr => point(mode, Nr, draws));
        console.log(`\nangles ${mode}, ${draws} realisations: mean SINR [dB] (± SE), EVM [%], IPI [dB]`);
        console.log(U.pad('SNR', 5), U.pad('Nr', 4), U.rpad('Optimal (CP demod)', 22), U.rpad('AD-MF (genie)', 20), U.rpad('EVM opt %', 11), U.rpad('EVM AD-MF %', 12), U.rpad('IPI dB', 8));
        SNRS.forEach((snr, i) => {
            NRS.forEach((Nr, k) => {
                const r = per[k][i], o = U.mean(r.opt), a = U.mean(r.ad), ipi = dB(U.mean(r.ipi));
                console.log(U.pad(snr, 5), U.pad(Nr, 4), U.rpad(`${U.f(o, 2)} ± ${U.f(U.se(r.opt), 2)}`, 22), U.rpad(`${U.f(a, 2)} ± ${U.f(U.se(r.ad), 2)}`, 20), U.rpad(U.f(100 / Math.sqrt(Math.pow(10, o / 10)), 3), 11), U.rpad(U.f(100 / Math.sqrt(Math.pow(10, a / 10)), 3), 12), U.rpad(U.f(ipi, 1), 8));
                rows.push([mode, snr, Nr, 'Optimal (CP demod)', o.toFixed(4), U.se(r.opt).toFixed(4), (100 / Math.sqrt(Math.pow(10, o / 10))).toFixed(4), '', draws]);
                rows.push([mode, snr, Nr, 'AD-MF (genie)', a.toFixed(4), U.se(r.ad).toFixed(4), (100 / Math.sqrt(Math.pow(10, a / 10))).toFixed(4), ipi.toFixed(2), draws]);
            });
            const x = cross(NRS, per.map(r => U.mean(r[i].ad)), per.map(r => U.mean(r[i].opt)));
            // 細查 Nr = 16 … 32（步長 2）
            const fine = []; for (let Nr = 16; Nr <= 32; Nr += 2) { const r = point(mode, Nr, draws)[i]; fine.push([Nr, U.mean(r.ad) - U.mean(r.opt)]); }
            const first = fine.find(f => f[1] >= 0);
            summary.push(`angles ${mode}, SNR ${snr} dB: AD-MF reaches the optimal adaptive combining at Nr ~ ${x === null ? 'no crossing within Nr <= 32' : U.f(x, 1) + ' (log-linear interpolation between the grid points)'}; finer search (step 2): ${first ? 'first Nr with AD-MF >= Optimal = ' + first[0] : 'none up to 32'}   [differences AD-MF - Optimal in dB: ${fine.map(f => `${f[0]}: ${U.f(f[1], 2)}`).join(', ')}]`);
        });
    }
    console.log('\ncrossing of AD-MF with the optimal adaptive combining (in this model):'); summary.forEach(s => console.log('  ' + s));
    const hdr = ['AD-MF (angle-domain matched filter) baseline versus array size Nr: SINR and EVM; genie baseline (true angles and Dopplers, MRC across the path branches; no DoA / parameter estimation; ICI and inter-path interference are not mitigated)',
        `schema_version 1; generated by tests/angle_domain_nr.js; seed ${R.CFG.seed}; ${draws} realisations per row; N 1024, delta_f 15 kHz, 4 paths (Doppler 1080 / -1080 / 758 / 220 Hz, relative power 0 / 0 / -11 / -0.7 dB), random phases (A3)`,
        'angle_mode fixed = F1 (0, 25, -35, 55 deg), random = U(-90, 90) deg per realisation (A4); SNR = sigma_x^2 sum|A_i|^2 / sigma_n^2 per antenna (A1)',
        'Optimal (CP demod) = G = R^-1 H with the CP-method epsilon (expected value); EVM = 100 / sqrt(SINR) percent (assumption A13); ipi_db = leakage of the other paths into a path branch relative to its own signal (mean over paths and realisations)',
        'mean and standard error (SE) of the per-realisation SINR in dB'];
    fs.writeFileSync(path.join(__dirname, '..', 'data', 'angle_domain_nr.csv'), hdr.map(l => '# ' + l).join('\n') + '\nangle_mode,snr_db,nr,curve,mean_sinr_db,se_db,evm_pct,ipi_db,n\n' + rows.map(r => r.join(',')).join('\n') + '\n');
    console.log('written data/angle_domain_nr.csv');
}
module.exports = { NRS, point, cross };
