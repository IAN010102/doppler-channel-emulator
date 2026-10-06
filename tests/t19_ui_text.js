'use strict';
/**
 * T19  Texts, tooltips and the five experiment presets (Commit 22).
 *   T19a  i18n tables: every Chinese phrase found in index.html has an English entry (missing ones are listed); no empty string in PHRASES / TIPS / UI / EXP; the tooltip tables have both languages
 *         and every tooltip is <= 40 characters; every tooltip id exists in the page; UI and EXP have the same keys in both languages (functions are called with the measured numbers).
 *   T19b  Each preset, applied to a fresh system (Experiments.applyTo, the code the page uses), sets exactly the expected parameter values (fixed table below; everything else = page defaults).
 *   T19c  After applying, the corresponding computation runs without an exception and returns finite numbers (sweep point / table / curve), and the page code path (ExpUI.run in the stubbed page)
 *         does not throw.
 *   T19d  The numbers quoted in the experiment texts (experiment_numbers.js) are reproduced by a new measurement with the same seeds (0.01 dB / 0.01 %), and the phenomena the texts state hold
 *         (stated per experiment in the output). The texts are written only after these hold.
 */
const fs = require('fs'), vm = require('vm'), path = require('path');
const root = path.join(__dirname, '..');
const U = require('./_util.js');
const Core = require('../core.js'), I18N = require('../i18n.js'), X = require('../experiments.js'), NUM = require('../experiment_numbers.js');
const { measure } = require('./experiment_measure.js');
const html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');

const CJKR = /[一-鿿　-〿＀-￯]+/g;
function flat(o, p = '') { return Object.entries(o).flatMap(([k, v]) => (v && typeof v === 'object' && !Array.isArray(v)) ? flat(v, p + k + '.') : [[p + k, v]]); }
function leaves(o, p = '') { return Object.entries(o).flatMap(([k, v]) => (v && typeof v === 'object') ? leaves(v, p + k + '.') : [[p + k, v]]); }

module.exports = {
    id: 'T19', title: 'i18n coverage and tooltips; experiment presets (parameters, runs, measured numbers and stated phenomena)',
    async run({ trials = 300 } = {}) {
        const checks = [], info = [];
        // ---------------------------------------------------------------- T19a
        const have = new Set(I18N.PHRASES.map(p => p[0])), runs = [...new Set(html.match(CJKR))], missing = runs.filter(r => !have.has(r));
        console.log(`T19a  ${runs.length} Chinese phrases in index.html, ${I18N.PHRASES.length} entries; missing: ${missing.length ? missing.join(' | ') : 'none'}`);
        checks.push(U.check('T19a every Chinese phrase of the page has an English entry', `${missing.length} missing`, '0 missing', missing.length === 0, missing.slice(0, 20).join(' | ')));
        const emptyPh = I18N.PHRASES.filter(p => !p[0] || p[1] === ''), dup = I18N.PHRASES.length - have.size;
        checks.push(U.check('T19a no empty phrase / duplicate key in PHRASES', `${emptyPh.length} empty, ${dup} duplicate`, '0 / 0', emptyPh.length === 0 && dup === 0));
        const tipBad = [], tipId = [];
        const ctlKeys = [...html.matchAll(/key: '([A-Za-z_]+)',\s+zh:/g)].map(m => m[1]);
        for (const [id, t] of Object.entries(I18N.TIPS)) {
            if (!t.zh || !t.en) tipBad.push(id + ': empty'); else if (t.zh.length > 40 || t.en.length > 40) tipBad.push(`${id}: ${t.zh.length}/${t.en.length} chars`);
            const present = id.startsWith('ctl-') ? ctlKeys.includes(id.slice(4)) : (html.includes(`id="${id}"`) || html.includes(`'${id}'`));
            if (!present) tipId.push(id);
        }
        console.log(`T19a  tooltips: ${Object.keys(I18N.TIPS).length}; too long / empty: ${tipBad.join('; ') || 'none'}; ids not found in the page: ${tipId.join(', ') || 'none'}`);
        checks.push(U.check('T19a tooltips: both languages, <= 40 characters', `${tipBad.length} problems`, '0', tipBad.length === 0, tipBad.join('; ')));
        checks.push(U.check('T19a every tooltip id exists in the page', `${tipId.length} missing`, '0', tipId.length === 0, tipId.join(', ')));
        // UI / EXP: same keys in both languages, no empty strings
        const uiZ = flat(I18N.UI.zh).map(x => x[0]), uiE = flat(I18N.UI.en).map(x => x[0]);
        const uiMiss = uiZ.filter(k => !uiE.includes(k)).concat(uiE.filter(k => !uiZ.includes(k)));
        const expProblems = [];
        for (const e of X.EXPERIMENTS) {
            const t = I18N.EXP[e.id]; if (!t) { expProblems.push(e.id + ': no texts'); continue; }
            for (const f of ['title', 'look', 'why', 'real']) for (const L of ['zh', 'en']) {
                let s = t[f] && t[f][L]; if (typeof s === 'function') s = s(NUM, e.defaultVariant);
                if (typeof s !== 'string' || !s.trim()) expProblems.push(`${e.id}.${f}.${L}`);
                if (typeof s === 'string' && /undefined|NaN|\[object/.test(s)) expProblems.push(`${e.id}.${f}.${L}: bad number`);
            }
            for (const [k, v] of Object.entries(e.variants)) { const lab = (t.variantLabel || {})[k]; if (Object.keys(e.variants).length > 1 && (!lab || !lab.zh || !lab.en)) expProblems.push(`${e.id} variant ${k} label`); }
        }
        const uiEmpty = flat(I18N.UI.zh).concat(flat(I18N.UI.en)).filter(([k, v]) => typeof v === 'string' && !v.trim() && !/apply$/.test(k)).map(x => x[0]);
        console.log(`T19a  UI keys differing between languages: ${uiMiss.join(', ') || 'none'}; empty UI strings: ${uiEmpty.join(', ') || 'none'}; experiment text problems: ${expProblems.join(', ') || 'none'}`);
        checks.push(U.check('T19a UI strings: same keys in both languages, none empty', `${uiMiss.length} missing, ${uiEmpty.length} empty`, '0 / 0', uiMiss.length === 0 && uiEmpty.length === 0));
        checks.push(U.check('T19a experiment texts: title/look/why/real in both languages, numbers filled in', `${expProblems.length} problems`, '0', expProblems.length === 0, expProblems.join(', ')));

        // ---------------------------------------------------------------- T19b
        const D = X.BASE, EXPECT = {
            'E0/': { N: 8, aoaT: -20, aoaJ: 30, snr: 30, sir: 0, L: 1000, kDb: 400, v: 0, pointErrDeg: 0, mod: 'QPSK', trainMode: 'withSignal', covSource: 'theory', algo: 'MMSE' },
            'E1/28': { fc: 28e9, scs: 120e3, v: 300 }, 'E1/5': { fc: 5e9, scs: 15e3, v: 300 }, 'E2/': { L: 4 }, 'E3/': { trainMode: 'withSignal', pointErrDeg: 3 },
            'E4/': { trainMode: 'signalFree', pointErrDeg: 3 }, 'E5/': { v: 300, latMs: 10, d_min: 30, aoaT: 0, aoaJ: -30 },
            'E6/': { v: 300, latMs: 10, d_min: 10, aoaT: 20, aoaJ: -30, pointingMode: 'mobility' }
        };
        for (const [key, exp] of Object.entries(EXPECT)) {
            const [id, variant] = key.split('/'), s = Core.createSys(); X.applyTo(s, id, variant);
            const want = Object.assign({}, { N: 8, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, v: 0, kDb: 20, latMs: 0, calDeg: 0, taper: 'NONE', mod: 'QAM16', pointErrDeg: 0, d_min: 30, fc: 5e9, scs: 15e3,
                model: 'unified', trainMode: 'signalFree', jamWave: 'gaussian', smiSingular: 'pinv', gammaRelDb: 10, algo: 'SMI', covSource: 'sample', angleSource: 'true', pointingMode: 'manual' }, exp);
            const bad = Object.keys(want).filter(k => s[k] !== want[k]).map(k => `${k}: ${s[k]} (expected ${want[k]})`);
            checks.push(U.check(`T19b ${id}${variant ? ' (' + variant + ')' : ''}: parameters after applying`, bad.length ? bad.join('; ') : `${Object.keys(want).length} values as expected`, 'all equal to the table', bad.length === 0));
        }

        // ---------------------------------------------------------------- T19c
        let cErr = [];
        const finite = a => a.every(Number.isFinite);
        for (const e of X.EXPERIMENTS) {
            try {
                const v = e.defaultVariant;
                if (e.action === 'table') { const r = X.runTable(Core, e.id, v, 4); if (!finite(r.map(x => x.sinr))) cErr.push(e.id + ': non-finite'); }
                else if (e.action === 'curve') { const r = X.runCurve(Core, e.id, v, 2); if (!r.series.every(s => finite(s.ys) && s.ys.length === r.xs.length)) cErr.push(e.id + ': non-finite'); }
                else { const r = X.sinrOf(Core, X.params(e.id, v), { v: 150 }, 5, 3); if (!Number.isFinite(r.sinr)) cErr.push(e.id + ': non-finite'); }
            } catch (err) { cErr.push(`${e.id}: ${err.message}`); }
        }
        // page code path
        try {
            const ctx = makePage();
            for (const e of X.EXPERIMENTS) ctx.ev(`ExpUI.select('${e.id}'); ExpUI.run();`);
            await new Promise(r => setTimeout(r, 600));
            ctx.ev("document.getElementById('btn-lang') && 0");
        } catch (err) { cErr.push('page: ' + err.message); }
        console.log(`T19c  computations of the 5 presets: ${cErr.length ? cErr.join('; ') : 'ran without exception, all values finite'}`);
        checks.push(U.check('T19c presets run without exceptions (library and page code path)', cErr.length ? cErr.join('; ') : 'ok', 'no exception, finite values', cErr.length === 0));

        // ---------------------------------------------------------------- T19d
        const m = measure(trials), A = leaves(m), B = Object.fromEntries(leaves(NUM));
        let worst = 0, worstKey = '';
        for (const [k, v] of A) { if (v === null || B[k] === null) { if (v !== B[k]) { worst = Infinity; worstKey = k; } continue; } const d = Math.abs(v - B[k]); if (d > worst) { worst = d; worstKey = k; } }
        console.log(`T19d  re-measured ${A.length} numbers: largest difference to experiment_numbers.js = ${U.e(worst)} (${worstKey})`);
        checks.push(U.check('T19d numbers quoted in the texts are reproduced by a new measurement', U.e(worst), '<= 0.01', worst <= 0.01, worstKey));
        const E0 = m.E0, th = E0['MMSE-M (theory R)'], tB = E0['SMI withSignal (theory R) = MVDR B'], tA = E0['SMI signalFree (theory R) = MVDR A'], sW = E0['SMI withSignal (sample R)'], sM = E0['MMSE-M (sample R)'], sF = E0['SMI signalFree (sample R)'];
        const claims = [
            ['E0: MMSE-M, MVDR B and MVDR A with the theoretical covariance are within 0.1 dB of SINR_opt', [th, tB, tA].every(r => Math.abs(r.sinr - r.opt) < 0.1)],
            ['E0: with the sample covariance and the target in the training data (MMSE-M, SMI) the SINR is more than 10 dB below the theory', th.sinr - sM.sinr > 10 && th.sinr - sW.sinr > 10],
            ['E0: SMI trained without the target (sample covariance) stays within 0.5 dB of the theory', th.sinr - sF.sinr < 0.5],
            ['E0: the output scaling bias is 0 for MVDR and > 0 but < 1e-3 for MMSE with the theoretical covariance', tB.bias < 1e-9 && tA.bias < 1e-9 && th.bias > 0 && th.bias < 1e-3],
            ['E0: the raw EVM with the theoretical covariance is below 2 %', th.evmRaw_pct < 2 && tB.evmRaw_pct < 2],
            ['E1: 28 GHz/120 kHz stays below the 16-QAM threshold (13.14 %) up to 300 km/h', m.E1_28.crossing_speed_kmh === null],
            ['E1: 5 GHz/15 kHz crosses 13.14 % at some speed <= 300 km/h', m.E1_5.crossing_speed_kmh !== null],
            ['E1: EVM rises with speed and stays within 1.5 points of the ICI limit at 300 km/h (both pairs)', [m.E1_28, m.E1_5].every(r => r.byV[0].evm_pct < r.byV[100].evm_pct && r.byV[100].evm_pct < r.byV[200].evm_pct && r.byV[200].evm_pct < r.byV[300].evm_pct && r.byV[300].evm_pct - r.byV[300].ici_pct < 1.5)],
            ['E2: rank of R_hat = 4 < N = 8', m.E2.SMI.rank === 4],
            ['E2: DL and BEAMSPACE are above SMI at L = 4', m.E2.DL.sinr > m.E2.SMI.sinr && m.E2.BEAMSPACE.sinr > m.E2.SMI.sinr],
            ['E3: SMI and DL drop by more than 10 dB from 0 to 3 degrees, MMSE-P changes by < 0.05 dB', m.E3.SMI.d0 - m.E3.SMI.d3 > 10 && m.E3.DL.d0 - m.E3.DL.d3 > 10 && Math.abs(m.E3.MMSEP.d0 - m.E3.MMSEP.d3) < 0.05],
            ['E3: SMI without the signal in the training data is far above SMI with it (0 degrees)', m.E3.ref_signalFree_SMI_d0 - m.E3.SMI.d0 > 20],
            ['E4: signalFree falls slowly (0 > 3 > 5 degrees, total < 5 dB), withSignal falls steeply (0 > 1 > 3 degrees, > 10 dB)', m.E4.signalFree[0] > m.E4.signalFree[3] && m.E4.signalFree[3] > m.E4.signalFree[5] && m.E4.signalFree[0] - m.E4.signalFree[5] < 5 && m.E4.withSignal[0] > m.E4.withSignal[1] && m.E4.withSignal[1] > m.E4.withSignal[3] && m.E4.withSignal[0] - m.E4.withSignal[3] > 10],
            ['E5: SINR falls with tau for SMI and DL (0 > 5 > 10 ms) at d_min = 30 m and 5 m', ['d30', 'd5'].every(d => ['SMI', 'DL'].every(a => m.E5[d][a][0] > m.E5[d][a][5] && m.E5[d][a][5] > m.E5[d][a][10]))],
            ['E6: delta_theta_eff grows with the speed and stays below 10 % of the 3 dB beamwidth', m.E6.byV[0].dEffT === 0 && m.E6.byV[100].dEffT < m.E6.byV[200].dEffT && m.E6.byV[200].dEffT < m.E6.byV[300].dEffT && m.E6.byV[300].ratio_pct < 10],
            ['E6: SMI (signalFree) and DL (signalFree) fall from 0 to 300 km/h, MMSE-P drops by more than 20 dB', m.E6.byV[0].smiSF > m.E6.byV[300].smiSF && m.E6.byV[0].dlSF > m.E6.byV[300].dlSF && m.E6.byV[0].mmseP - m.E6.byV[300].mmseP > 20],
            ['E6: SMI withSignal rises first (peak strictly between 0 and 300 km/h, > 5 dB above its 0 km/h value) and is lower at 300 km/h than at the peak', m.E6.peak_smiWS.v > 0 && m.E6.peak_smiWS.v < 300 && m.E6.peak_smiWS.sinr - m.E6.byV[0].smiWS > 5 && m.E6.byV[300].smiWS < m.E6.peak_smiWS.sinr],
            ['E6: the ICI floor is above the SMI (signalFree) aging loss from the first non-zero speed of the sweep (25 km/h)', m.E6.cross_ici_over_aging_kmh === 25],
            ['E6: with tau = 0, SMI signalFree at 300 km/h is within 0.5 dB of its 0 km/h value and SMI withSignal is > 5 dB above its 0 km/h value', Math.abs(m.E6.tau0.smiSF_300 - m.E6.byV[0].smiSF) < 0.5 && m.E6.tau0.smiWS_300 - m.E6.byV[0].smiWS > 5]
        ];
        for (const [txt, ok] of claims) { console.log(`T19d  ${ok ? 'holds ' : 'FAILS '} ${txt}`); checks.push(U.check('T19d stated phenomenon: ' + txt, ok ? 'holds' : 'does not hold', 'holds', ok)); }
        return { id: this.id, title: this.title, checks, info };
    }
};

// the page in a stubbed DOM (same stubs as ui_smoke.js)
function makePage() {
    function stubEl(id) {
        const L = {}, el = { id, innerHTML: '', innerText: '', textContent: '', className: '', style: {}, dataset: {}, children: [], disabled: false, title: '', value: '', clientWidth: 600, clientHeight: 300, parentElement: null, scrollTop: 0, scrollHeight: 0,
            classList: { toggle() {}, add() {}, remove() {}, contains() { return false; } },
            addEventListener(t, f) { (L[t] = L[t] || []).push(f); }, appendChild(c) { c.parentElement = el; el.children.push(c); return c; }, insertBefore(c) { c.parentElement = el; el.children.unshift(c); return c; },
            removeChild() {}, get firstChild() { return el.children[0]; }, closest() { return null; }, querySelector() { return Object.assign(stubEl(), { parentElement: el }); },
            setAttribute(k, v) { el[k] = v; }, getAttribute(k) { return el[k] === undefined ? null : el[k]; }, querySelectorAll() { return []; },
            getContext() { return new Proxy({}, { get: (t, p) => p in t ? t[p] : () => {}, set: (t, p, v) => { t[p] = v; return true; } }); } };
        return el;
    }
    const els = {}, inline = html.slice(html.lastIndexOf('<script>') + 8, html.lastIndexOf('</script>'));
    const ctx = vm.createContext({ document: { getElementById: id => els[id] || (els[id] = stubEl(id)), createElement: () => stubEl(), body: stubEl(), activeElement: null, querySelectorAll: () => [], documentElement: {}, createTreeWalker: () => ({ nextNode: () => null }) },
        window: { devicePixelRatio: 1 }, performance, console, setTimeout, clearTimeout, URL: {}, Blob: function () {}, requestAnimationFrame: () => {} });
    for (const f of ['core.js', 'i18n.js', 'experiment_numbers.js', 'experiments.js']) vm.runInContext(fs.readFileSync(path.join(root, f), 'utf8'), ctx);
    vm.runInContext(inline, ctx);
    return { ev: code => vm.runInContext(code, ctx) };
}
