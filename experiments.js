/* experiments.js -- the five one-click experiment presets (DOM-free, shared by the page and the Node tests).
 * Each preset is a set of Sys fields plus what is run afterwards (a sweep, a table or a curve). The texts live in i18n.js; the numbers quoted in the texts are measured by
 * tests/experiment_measure.js (fixed seeds) and stored in experiment_numbers.js / data/experiment_numbers.json (T19d measures them again). */
(function (root, factory) {
    if (typeof module === 'object' && module.exports) module.exports = factory(); else root.Experiments = factory();
}(typeof self !== 'undefined' ? self : this, function () {
    'use strict';
    // the page defaults of the simulator: every experiment starts from them
    const BASE = { N: 8, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, v: 0, kDb: 20, latMs: 0, calDeg: 0, taper: 'NONE', mod: 'QAM16', pointErrDeg: 0, d_min: 30,
        fc: 5e9, scs: 15e3, model: 'unified', trainMode: 'signalFree', jamWave: 'gaussian', smiSingular: 'pinv', gammaRelDb: 10, algo: 'SMI', covSource: 'sample', angleSource: 'true' };
    const EXPERIMENTS = [
        // E0: the lecture baseline (N = 8, d = 0.5 lambda, theta1 = -20, theta2 = 30, SNR 30 dB, Es = 1, interferer power = Es, L = 1000 snapshots, K = 400 dB = infinity: no diffuse paths)
        { id: 'E0', variants: { '': { N: 8, aoaT: -20, aoaJ: 30, snr: 30, sir: 0, L: 1000, kDb: 400, v: 0, pointErrDeg: 0, mod: 'QPSK', trainMode: 'withSignal', covSource: 'theory' } }, defaultVariant: '', action: 'table', algo: 'MMSE',
          rows: [{ label: 'MMSE-M (theory R)', over: { algo: 'MMSE', covSource: 'theory' } }, { label: 'MMSE-M (sample R)', over: { algo: 'MMSE', covSource: 'sample' } },
                 { label: 'SMI withSignal (theory R) = MVDR B', over: { algo: 'SMI', trainMode: 'withSignal', covSource: 'theory' } }, { label: 'SMI withSignal (sample R)', over: { algo: 'SMI', trainMode: 'withSignal', covSource: 'sample' } },
                 { label: 'SMI signalFree (theory R) = MVDR A', over: { algo: 'SMI', trainMode: 'signalFree', covSource: 'theory' } }, { label: 'SMI signalFree (sample R)', over: { algo: 'SMI', trainMode: 'signalFree', covSource: 'sample' } }] },
        { id: 'E1', variants: { '28': { fc: 28e9, scs: 120e3, v: 300 }, '5': { fc: 5e9, scs: 15e3, v: 300 } }, defaultVariant: '28', action: 'bench', algo: 'SMI' },
        { id: 'E2', variants: { '': { L: 4 } }, defaultVariant: '', action: 'table', algos: ['SMI', 'DL', 'BEAMSPACE'], algo: 'SMI' },
        { id: 'E3', variants: { '': { trainMode: 'withSignal', pointErrDeg: 3 } }, defaultVariant: '', action: 'table', algos: ['SMI', 'DL', 'MMSEP'], algo: 'SMI' },
        { id: 'E4', variants: { '': { trainMode: 'signalFree', pointErrDeg: 3 } }, defaultVariant: '', action: 'curve', algos: ['SMI'], algo: 'SMI',
          curve: { x: 'pointErrDeg', from: -5, to: 5, step: 0.5, unit: '°', series: [{ label: 'signalFree', over: { trainMode: 'signalFree' } }, { label: 'withSignal', over: { trainMode: 'withSignal' } }] } },
        { id: 'E5', variants: { '': { v: 300, latMs: 10, d_min: 30, aoaT: 0, aoaJ: -30 } }, defaultVariant: '', action: 'curve', algos: ['SMI', 'DL'], algo: 'SMI',
          curve: { x: 'latMs', from: 0, to: 10, step: 1, unit: 'ms', series: [{ label: 'SMI', over: { algo: 'SMI' } }, { label: 'DL', over: { algo: 'DL' } }] } }
    ];
    const byId = id => EXPERIMENTS.find(e => e.id === id);
    // the parameters an experiment sets (BASE + variant + the experiment's algorithm)
    function params(id, variant) { const e = byId(id); return Object.assign({}, BASE, e.variants[variant === undefined || variant === null ? e.defaultVariant : variant], { algo: e.algo }); }
    function applyTo(sys, id, variant) { Object.assign(sys, params(id, variant)); sys.snaps = []; sys.snapKey = ''; sys.real = null; sys.dirty = true; return sys; }
    // computations behind "run" (Core is passed in so that the page and the tests use the same code); trials: realisations per point
    function sinrOf(Core, p, over, seed, trials) {
        let sum = 0, rank = 0, opt = 0, raw = 0;
        for (let t = 0; t < trials; t++) {
            Core.setSeed(seed + t);
            const s = Core.createSys(); s.calZ = new Array(16).fill(0); Object.assign(s, p, over, { freshRealization: true }); s.snaps = []; s.computeMath(); sum += s.sinrDb; rank = s.rankR; opt += s.sinrOptDb; raw += s.evmRaw;
        }
        return { sinr: sum / trials, rank, opt: opt / trials, evmRaw: raw / trials };
    }
    function runTable(Core, id, variant, trials) {
        const e = byId(id), p = params(id, variant);
        if (e.rows) return e.rows.map(r => Object.assign({ algo: r.over.algo, label: r.label }, sinrOf(Core, p, r.over, 777, trials)));
        return e.algos.map(a => Object.assign({ algo: a }, sinrOf(Core, p, { algo: a }, 777, trials)));
    }
    function runCurve(Core, id, variant, trials) {
        const e = byId(id), p = params(id, variant), c = e.curve, xs = []; for (let x = c.from; x <= c.to + 1e-9; x += c.step) xs.push(Math.round(x * 1000) / 1000);
        return { x: c.x, unit: c.unit, xs, series: c.series.map(sr => ({ label: sr.label, ys: xs.map(x => sinrOf(Core, p, Object.assign({}, sr.over, { [c.x]: x }), 888, trials).sinr) })) };
    }
    return { BASE, EXPERIMENTS, byId, params, applyTo, runTable, runCurve, sinrOf };
}));
