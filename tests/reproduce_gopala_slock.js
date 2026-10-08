'use strict';
/**
 * 重現 Gopala & Slock (EURECOM 2016) 的接收端 ICI 感知波束成形結果（資訊性為主）；輸出 data/reproduce_*.csv。
 *   node tests/reproduce_gopala_slock.js [--draws 200] [--fast]
 * 設定：N = 1024、Δf = 15 kHz、Nr = 2、SNR 0–35 dB（步長 5）、4 條路徑（LoS 加 3 條）、都卜勒 [1080, −1080, 758, 220] Hz、相對強度 [0, 0, −11, −0.7] dB、
 * 導頻間距 12、平坦頻率通道。到達角文獻未寫明：(i) 固定一組角度（200 次只抽隨機相位）與 (ii) 每次從 U(−90°, 90°) 隨機抽取（200 次）各跑一遍。
 * 全部假設見 docs/ofdm_ici_assumptions.md。不調任何參數去貼近文獻。
 *   Fig1 類：最佳波束成形（含 CP 法對頻）、MRC（不考慮 ICI、不對頻）、單天線（不對頻）。
 *   Fig3 類：Optimal BF、Estimated BF（R_yy 與 Ĥ，CP 法對頻）、Approx BF（線性近似，CP 法對頻）、Approx BF 不對頻、以 LoS 都卜勒對頻（Optimal 與 Approx）；
 *            「多重對頻」（可選做）未實作。
 *   Fig2 類：對頻頻率估計方法比較（CP 法解析期望值、CP 法由含雜訊樣本估計、對 SINR 的窮舉、對訊號功率的窮舉；文獻式 24 的簡化近似式：
 *            提示詞沒有給出式子，未實作），10 組通道實現。
 * 每點附標準誤（SE）。CSV 檔頭為 # 開頭的註解（參數、種子、次數、假設、schema 版本）。
 */
const fs = require('fs'), path = require('path');
const O = require('../ofdm_ici.js');
const U = require('./_util.js');

const CFG = { N: 1024, df: 15e3, Nr: 2, Ncp: 72, snrs: [0, 5, 10, 15, 20, 25, 30, 35], draws: 200, seed: 20160101, fixedAngles: [0, 25, -35, 55], pilotSpacing: 12, schema: 1 };
const ANGLE_SETS = { F1: [0, 25, -35, 55], F2: [10, -20, 40, -60], F3: [-5, 5, 15, 45] };

function pickAngles(mode, rng, set, range = 90) {
    if (mode === 'fixed') return ANGLE_SETS[set || 'F1'];
    return [0, 1, 2, 3].map(() => -range + 2 * range * rng.u());
}
const dB = x => 10 * Math.log10(x);

/** 一個通道實現、一組 SNR：回傳各曲線的線性 SINR（陣列，對應 snrs）
 *  o = { snrs, Nr, N, Ncp, losIndex, withEstimated (預設 true), rng (模擬符元用) } */
function evalChannel(ch, o) {
    const { N, Nr } = ch, snrs = o.snrs, Ncp = o.Ncp, losIdx = o.losIndex || 0, sx2 = 1;
    const eCp = O.epsCpAnalytic(ch, Ncp), eLos = O.epsLos(ch, losIdx), E = { none: 0, cp: eCp, los: eLos }, S = {};
    for (const k of Object.keys(E)) S[k] = O.system(ch, E[k], sx2, 0);                      // H 與 B（B 與 σn² 無關）
    const out = { eps: { cp: eCp, los: eLos } }, names = ['optCp', 'optNone', 'optLos', 'mrc', 'single', 'est', 'estTrueH', 'apxCp', 'apxNone', 'apxLos', 'estEps'];
    for (const n of names) out[n] = [];
    snrs.forEach(snr => {
        const sn2 = O.snrToNoise(ch, snr, sx2), R = k => O.rMatrix(S[k].B, Nr, sx2, sn2);
        out.optCp.push(O.optimal(S.cp.H, R('cp'), Nr, sx2).sinr); out.optNone.push(O.optimal(S.none.H, R('none'), Nr, sx2).sinr); out.optLos.push(O.optimal(S.los.H, R('los'), Nr, sx2).sinr);
        out.mrc.push(O.sinrOfWeights(S.none.H, S.none.H, R('none'), Nr, sx2));                // G = H（不對頻、不考慮 ICI）
        out.single.push(sx2 * (S.none.H.re[0] ** 2 + S.none.H.im[0] ** 2) / R('none')[0]);   // 只用天線 0
        for (const [key, k] of [['apxCp', 'cp'], ['apxNone', 'none'], ['apxLos', 'los']]) { const G = O.linearWeights(ch, E[k], sx2, sn2); out[key].push(O.sinrOfWeights(G, S[k].H, R(k), Nr, sx2)); }
        if (o.withEstimated !== false) {
            const sym = O.makeSymbol(ch, { rng: o.rng, sigmaX2: sx2, sigmaN2: sn2, Ncp, pilotSpacing: CFG.pilotSpacing }), eh = O.epsCpEstimate(sym), Y = O.demodulate(sym, eh), pw = O.practicalWeights(sym, Y, { Ns: N });
            const s2 = O.system(ch, eh, sx2, sn2); out.est.push(O.sinrOfWeights(pw.G, s2.H, s2.R, Nr, sx2)); out.estEps.push(eh);
            out.estTrueH.push(O.sinrOfWeights(O.csolve(pw.Ryy, s2.H, Nr), s2.H, s2.R, Nr, sx2));   // 診斷：R_yy 照用，Ĥ 換成真實 H（看導頻估計誤差的影響）
        }
    });
    return out;
}

/** 一個角度設定的 Fig1/Fig3 曲線：draws 次通道實現。回傳 { curves: {name: [snrIdx][draw] (dB)} , epsCp: [draw], ... } */
function computeCurves(o) {
    const draws = o.draws || CFG.draws, snrs = o.snrs || CFG.snrs, Nr = o.Nr || CFG.Nr, N = o.N || CFG.N, Ncp = o.Ncp !== undefined ? o.Ncp : CFG.Ncp, base = o.seed || CFG.seed;
    const names = ['optCp', 'optNone', 'optLos', 'mrc', 'single', 'est', 'estTrueH', 'apxCp', 'apxNone', 'apxLos'];
    const curves = {}; names.forEach(n => { curves[n] = snrs.map(() => []); });
    const epsCp = [], epsLos = [], estEps = snrs.map(() => []);
    for (let t = 0; t < draws; t++) {
        // 共同亂數（common random numbers）：相位 rngPhase 只依實現編號 t；角度 rngAngle 與模擬符元的 rngSym 各自獨立，
        // 所以同一個 t 在不同情境（角度、LoS 對應、CP 長度、都卜勒縮放、SNR 清單）下是同一組相位，比較時成對
        const rngPhase = O.makeRng(base + 7919 * t), rngAngle = O.makeRng(base + 500000 + 7919 * t), rngSym = O.makeRng(base + 1000000 + 7919 * t);
        const ang = pickAngles(o.mode, rngAngle, o.set, o.range), ch = O.paperChannel(ang, rngPhase, { Nr, N, dopplerScale: o.dopplerScale });
        const r = evalChannel(ch, { snrs, Nr, N, Ncp, losIndex: o.losIndex, rng: rngSym, withEstimated: o.withEstimated });
        names.forEach(n => { if (r[n].length) r[n].forEach((v, k) => curves[n][k].push(dB(v))); });
        epsCp.push(r.eps.cp); epsLos.push(r.eps.los); if (r.estEps.length) r.estEps.forEach((v, k) => estEps[k].push(v));
    }
    return { curves, epsCp, epsLos, estEps, snrs, draws };
}

/** Fig2 類：10 組通道（隨機角度）的對頻估計比較，SNR = snrDb */
function computeFig2(o = {}) {
    const snrDb = o.snrDb || 20, n = o.n || 10, rng = O.makeRng((o.seed || CFG.seed) + 33), rows = [];
    for (let t = 0; t < n; t++) {
        const ang = pickAngles('random', rng), ch = O.paperChannel(ang, rng), sn2 = O.snrToNoise(ch, snrDb, 1), sx2 = 1;
        const sym = O.makeSymbol(ch, { rng, sigmaX2: sx2, sigmaN2: sn2, Ncp: CFG.Ncp, pilotSpacing: CFG.pilotSpacing });
        const g = O.epsGenie(ch, sx2, sn2, { objective: 'sinr' }), pw = O.epsGenie(ch, sx2, sn2, { objective: 'power' });
        const ests = { 'CP (expected; noise free)': O.epsCpAnalytic(ch, CFG.Ncp), 'CP (from noisy samples)': O.epsCpEstimate(sym), 'exhaustive search on SINR (genie)': g.eps, 'exhaustive search on signal power': pw.eps, 'LoS Doppler': O.epsLos(ch, 0), 'none': 0 };
        for (const [m, e] of Object.entries(ests)) rows.push({ draw: t, method: m, eps: e, sinrDb: dB(O.sinrOptimalAt(ch, e, sx2, sn2)), angles: ang.map(a => a.toFixed(1)).join('/') });
    }
    return rows;
}

/** 假設的敏感度：改變一個未明訂的設定，其餘不變；SNR 20 與 35 dB；回傳 [{scenario, curve, snr, mean, se}] */
function computeSensitivity(draws) {
    const sc = [['A4 angles: random U(-90, 90) (main)', { mode: 'random' }], ['A4 angles: random U(-60, 60)', { mode: 'random', range: 60 }], ['A4 angles: fixed F1 = 0, 25, -35, 55', { mode: 'fixed', set: 'F1' }],
        ['A4 angles: fixed F2 = 10, -20, 40, -60', { mode: 'fixed', set: 'F2' }], ['A4 angles: fixed F3 = -5, 5, 15, 45', { mode: 'fixed', set: 'F3' }],
        ['A2 LoS = path 2 (-1080 Hz) instead of path 1', { mode: 'random', losIndex: 1 }], ['A5 CP length 14 % (Ncp = 144) instead of 7 %', { mode: 'random', Ncp: 144 }],
        ['A6 Doppler x0.5 (all paths)', { mode: 'random', dopplerScale: 0.5 }], ['A6 Doppler x2 (all paths)', { mode: 'random', dopplerScale: 2 }]], rows = [];
    for (const [name, o] of sc) {
        const r = computeCurves(Object.assign({ draws, snrs: [20, 35] }, o));
        for (const c of ['optCp', 'optNone', 'optLos', 'mrc', 'single', 'est', 'apxCp']) r.snrs.forEach((s2, i) => { const t = stat(r.curves[c][i]); rows.push({ scenario: name, curve: c, snr: s2, mean: t.m, se: t.se }); });
        rows.push({ scenario: name, curve: 'mean |eps_CP|', snr: '', mean: U.mean(r.epsCp.map(Math.abs)), se: U.se(r.epsCp.map(Math.abs)) });
    }
    return rows;
}

const stat = a => ({ m: U.mean(a), se: U.se(a) });
function writeCsv(file, header, cols, rows) { fs.mkdirSync(path.dirname(file), { recursive: true }); fs.writeFileSync(file, header.map(l => '# ' + l).join('\n') + '\n' + cols.join(',') + '\n' + rows.map(r => r.join(',')).join('\n') + '\n'); }
const hdr = (title, extra = []) => [title, `schema_version ${CFG.schema}; generated by tests/reproduce_gopala_slock.js; seed ${CFG.seed}`, `N ${CFG.N}, delta_f ${CFG.df} Hz, Nr ${CFG.Nr}, Ncp ${CFG.Ncp} (7 %), pilot spacing ${CFG.pilotSpacing}, sigma_x^2 1, QPSK`,
    'paths (Doppler Hz / relative power dB): 1080 / 0 (LoS, assumption A2), -1080 / 0, 758 / -11, 220 / -0.7; random phases (A3); SNR = sigma_x^2 sum|A_i|^2 / sigma_n^2 per antenna (A1)', 'angles: not given in the paper: fixed set F1 = 0, 25, -35, 55 deg, or drawn from U(-90, 90) deg for every realisation (A4)', ...extra,
    'every point: mean of the per-realisation SINR in dB and its standard error (SE = sample standard deviation / sqrt(n)); see docs/ofdm_ici_assumptions.md'];

if (require.main === module) {
    const argv = process.argv.slice(2), draws = argv.includes('--draws') ? +argv[argv.indexOf('--draws') + 1] : CFG.draws, root = path.join(__dirname, '..', 'data');
    const modes = [{ mode: 'fixed', set: 'F1' }, { mode: 'random' }];
    const rows1 = [], rows3 = [], t0 = Date.now(); const res = {};
    for (const m of modes) {
        res[m.mode] = computeCurves(Object.assign({ draws }, m)); console.log(`angles ${m.mode}: ${draws} draws done (${((Date.now() - t0) / 1000).toFixed(0)} s)`);
        const R = res[m.mode];
        const f1 = { optCp: 'Optimal BF (CP demod)', mrc: 'MRC (no ICI model; no demod)', single: 'Single antenna (no demod)' }, f3 = { optCp: 'Optimal BF (CP demod)', est: 'Estimated BF (R_yy and H_hat; CP demod)', estTrueH: 'Estimated BF with the true H instead of H_hat (diagnostic)', apxCp: 'Approx BF (linear; CP demod)', apxNone: 'Approx BF (linear; no demod)', optNone: 'Optimal BF (no demod)', optLos: 'Optimal BF (LoS Doppler demod)', apxLos: 'Approx BF (linear; LoS Doppler demod)' };
        for (const [k, name] of Object.entries(f1)) R.snrs.forEach((s, i) => { const t = stat(R.curves[k][i]); rows1.push([m.mode, s, name, t.m.toFixed(4), t.se.toFixed(4), R.curves[k][i].length]); });
        for (const [k, name] of Object.entries(f3)) R.snrs.forEach((s, i) => { const t = stat(R.curves[k][i]); rows3.push([m.mode, s, name, t.m.toFixed(4), t.se.toFixed(4), R.curves[k][i].length]); });
    }
    writeCsv(path.join(root, 'reproduce_fig1.csv'), hdr('Fig. 1 type: SINR versus SNR, Nr = 2'), ['angle_mode', 'snr_db', 'curve', 'mean_sinr_db', 'se_db', 'n'], rows1);
    writeCsv(path.join(root, 'reproduce_fig3.csv'), hdr('Fig. 3 type: SINR versus SNR, Nr = 2, Optimal / Estimated / Approx beamforming with and without demodulation'), ['angle_mode', 'snr_db', 'curve', 'mean_sinr_db', 'se_db', 'n'], rows3);
    const f2 = computeFig2();
    writeCsv(path.join(root, 'reproduce_fig2.csv'), hdr('Fig. 2 type: demodulation frequency (eps) chosen by different methods, 10 channel realisations at SNR 20 dB; the simplified approximation of eq. (24) is not implemented (formula not given)'), ['draw', 'method', 'eps', 'sinr_opt_db_with_this_eps', 'angles_deg'], f2.map(r => [r.draw, r.method, r.eps.toFixed(5), r.sinrDb.toFixed(3), r.angles]));
    const sens = computeSensitivity(draws);
    writeCsv(path.join(root, 'reproduce_sensitivity.csv'), hdr('Sensitivity of the results to the assumptions A2, A4, A5, A6 (SNR 20 and 35 dB); one setting changed at a time'), ['scenario', 'curve', 'snr_db', 'mean', 'se', 'n'], sens.map(r => [r.scenario.replace(/,/g, ';'), r.curve, r.snr, r.mean.toFixed(4), r.se.toFixed(4), draws]));
    console.log('written data/reproduce_fig1.csv, reproduce_fig2.csv, reproduce_fig3.csv, reproduce_sensitivity.csv');
}
module.exports = { CFG, ANGLE_SETS, computeCurves, computeFig2, computeSensitivity, evalChannel, pickAngles, stat };
