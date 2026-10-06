'use strict';
/**
 * node tests/experiment_measure.js [--write]   Measures the numbers quoted in the experiment texts (fixed seeds, unified model) and prints them;
 * with --write the result is stored in data/experiment_numbers.json. `measure()` is also used by T19d (the stored numbers must be reproduced).
 */
const fs = require('fs'), path = require('path');
const Core = require('../core.js');
const X = require('../experiments.js');

function setup(over, seed) {
    Core.setSeed(seed);
    const s = Core.createSys(); s.calZ = new Array(16).fill(0);
    Object.assign(s, X.BASE, over, { freshRealization: true });
    return s;
}
const mean = a => a.reduce((x, y) => x + y, 0) / a.length;
const r2 = x => Math.round(x * 100) / 100;
function sinr(over, trials, seed0) { let a = 0, rk = 0; for (let t = 0; t < trials; t++) { const s = setup(over, seed0 + t); s.snaps = []; s.computeMath(); a += s.sinrDb; rk = s.rankR; } return { sinr: r2(a / trials), rank: rk }; }
function evmRms(over, trials, seed0) { let e2 = 0, ici = 0; for (let t = 0; t < trials; t++) { const s = setup(over, seed0 + t); s.snaps = []; s.computeMath(); e2 += s.evm * s.evm; ici = Math.sqrt(s.nuIciFloor); } return { evm: Math.sqrt(e2 / trials), ici }; }

function measureE6(trials) {
    // E6: speed -> pointing deviation (pointingMode 'mobility'); SINR of four series per speed, delta_theta_eff, its ratio to the 3 dB beamwidth, ICI floor, and the speed at which
    // the ICI floor N_ICI/S exceeds the aging loss of SMI signalFree (1/SINR(v) - 1/SINR(0), linear) -- both are contributions to the EVM^2 = 1/SINR + N_ICI/S
    const out = {};
    out.E6 = { byV: {}, cross_ici_over_aging_kmh: null };
    const e6 = X.byId('E6'), p6 = X.params('E6'), res6 = X.runCurve(Core, 'E6', '', Math.min(trials, 100));
    const s0 = res6.series[0].ys[0];
    res6.xs.forEach((v, i) => {
        const lin = db => Math.pow(10, db / 10), aging = 1 / lin(res6.series[0].ys[i]) - 1 / lin(s0), ici = lin(res6.aux.iciDb[i]);
        if (out.E6.cross_ici_over_aging_kmh === null && v > 0 && ici > aging) out.E6.cross_ici_over_aging_kmh = v;
        if ([0, 50, 100, 200, 300].includes(v)) out.E6.byV[v] = { dEffT: r2(res6.aux.dEffT[i]), dEffJ: r2(res6.aux.dEffJ[i]), ratio_pct: r2(100 * Math.abs(res6.aux.dEffT[i]) / res6.aux.bw3[i]), bw3: r2(res6.aux.bw3[i]), iciDb: r2(res6.aux.iciDb[i]), agingDb: r2(10 * Math.log10(Math.max(aging, 1e-30))),
            smiSF: r2(res6.series[0].ys[i]), smiWS: r2(res6.series[1].ys[i]), dlSF: r2(res6.series[2].ys[i]), mmseP: r2(res6.series[3].ys[i]) };
    });
    out.E6.peak_smiWS = { v: res6.xs[res6.series[1].ys.indexOf(Math.max(...res6.series[1].ys))], sinr: r2(Math.max(...res6.series[1].ys)) };
    // the same with tau = 0: no angle deviation, only the Doppler inside the training window
    const t0 = sinr(Object.assign({}, X.params('E6'), { algo: 'SMI', trainMode: 'signalFree', latMs: 0, v: 300 }), trials, 5000), w0 = sinr(Object.assign({}, X.params('E6'), { algo: 'SMI', trainMode: 'withSignal', latMs: 0, v: 300 }), trials, 5000), w1 = sinr(Object.assign({}, X.params('E6'), { algo: 'SMI', trainMode: 'withSignal', latMs: 0, v: 100 }), trials, 5000);
    out.E6.tau0 = { smiSF_300: t0.sinr, smiWS_300: w0.sinr, smiWS_100: w1.sinr };
    return out.E6;
}

function measure(trials = 300) {
    const out = {};
    // E0: the lecture baseline (at most 200 realisations of L = 1000 snapshots each)
    out.E0 = {};
    for (const r of X.byId('E0').rows) { const s = { sinr: 0, opt: 0, raw: 0, bias: 0 }, n = Math.min(trials, 200); for (let t = 0; t < n; t++) { const sy = setup(Object.assign({}, X.byId('E0').variants[''], r.over), 700 + t); sy.snaps = []; sy.computeMath(); s.sinr += sy.sinrDb; s.opt += sy.sinrOptDb; s.raw += sy.evmRaw; s.bias += Math.hypot(sy.gRawR - 1, sy.gRawI); } out.E0[r.label] = { sinr: r2(s.sinr / n), opt: r2(s.opt / n), evmRaw_pct: r2(100 * s.raw / n), bias: Math.round(1e6 * s.bias / n) / 1e6 }; }
    // E1: EVM (rms, as in the sweep) of SMI vs speed and the ICI floor, for both (fc, df) pairs; the speed at which the SMI EVM crosses the 16-QAM SER = 1e-3 threshold (13.14 %)
    const thr = 0.1314;
    for (const key of ['28', '5']) {
        const row = {}; let cross = null;
        for (let v = 0; v <= 300; v += 10) {
            const r = evmRms(Object.assign({}, X.byId('E1').variants[key], { v, algo: 'SMI' }), trials, 1000 + v);
            if (v === 0 || v === 100 || v === 200 || v === 300) row[v] = { evm_pct: r2(100 * r.evm), ici_pct: r2(100 * r.ici) };
            if (cross === null && r.evm > thr) cross = v;
        }
        out['E1_' + key] = { byV: row, crossing_speed_kmh: cross };
    }
    // E2: L = 4, N = 8
    out.E2 = {};
    for (const a of ['SMI', 'DL', 'BEAMSPACE']) out.E2[a] = sinr(Object.assign({}, X.byId('E2').variants[''], { algo: a }), 1000, 2000);
    // E3: withSignal, delta = 3 deg vs delta = 0 (and signalFree SMI as reference)
    out.E3 = { ref_signalFree_SMI_d0: sinr({ algo: 'SMI', trainMode: 'signalFree', pointErrDeg: 0 }, 1000, 3000).sinr };
    for (const a of ['SMI', 'DL', 'MMSEP']) out.E3[a] = { d0: sinr({ algo: a, trainMode: 'withSignal', pointErrDeg: 0 }, 1000, 3000).sinr, d3: sinr({ algo: a, trainMode: 'withSignal', pointErrDeg: 3 }, 1000, 3000).sinr };
    // E4: SMI vs delta_theta, both training modes
    out.E4 = {};
    for (const tm of ['signalFree', 'withSignal']) { out.E4[tm] = {}; for (const d of [0, 1, 3, 5]) out.E4[tm][d] = sinr({ algo: 'SMI', trainMode: tm, pointErrDeg: d }, 1000, 4000).sinr; }
    // E5: aging
    out.E5 = {};
    for (const dm of [30, 5]) { out.E5['d' + dm] = {}; for (const a of ['SMI', 'DL']) { out.E5['d' + dm][a] = {}; for (const tau of [0, 5, 10]) out.E5['d' + dm][a][tau] = sinr(Object.assign({}, X.byId('E5').variants[''], { d_min: dm, latMs: tau, algo: a }), 1000, 5000).sinr; } }
    out.E6 = measureE6(trials);
    return out;
}

if (require.main === module) {
    const nTr = +(process.argv.find(a => /^\d+$/.test(a)) || 300), onlyE6 = process.argv.includes('--only-e6');
    const res = onlyE6 ? Object.assign(JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'data', 'experiment_numbers.json'), 'utf8')), { E6: measureE6(nTr) }) : measure(nTr);
    console.log(JSON.stringify(res, null, 1));
    if (process.argv.includes('--write')) { const f = path.join(__dirname, '..', 'data', 'experiment_numbers.json'); fs.mkdirSync(path.dirname(f), { recursive: true }); fs.writeFileSync(f, JSON.stringify(res, null, 1) + '\n'); console.log('written', f);
        const g = path.join(__dirname, '..', 'experiment_numbers.js');            // the same numbers for the page (no fetch needed)
        fs.writeFileSync(g, '/* generated by tests/experiment_measure.js --write (fixed seeds); T19d re-measures and compares */\n(function (root, f) { if (typeof module === "object" && module.exports) module.exports = f(); else root.EXP_NUMBERS = f(); }(typeof self !== "undefined" ? self : this, function () { return ' + JSON.stringify(res) + '; }));\n'); console.log('written', g);
    }
}
module.exports = { measure, measureE6 };
