'use strict';
/**
 * T10  Seeds and CSV (Commit 12), run against the real page script (vm context with a stubbed DOM, as in ui_smoke.js).
 *   T10a  same seed and parameters -> two independent page loads produce byte-identical CSV, except the timestamp row (one row, its own line).
 *   T10b  a different seed -> different CSV (sweep table and the seed row).
 *   T10c  any single sweep point re-run alone (Bench.point with its point index / algorithm index) reproduces the number in the sweep.
 *   T10d  schema: csv_schema_version = 2 and every required key is present (model, train_mode, gamma_rel_dB, d_min, M, K, sigma_theta, sigma_phi, P_s, tau,
 *         theta1, theta2, fc_GHz, delta_f_kHz, N, L, seed, point_seed, max_angle_drift_deg, SINR_opt_dB, sweep_info with train=, algorithm, algorithm_legacy).
 *   Run for the legacy and the unified model.
 */
const fs = require('fs'), vm = require('vm'), path = require('path');
const U = require('./_util.js');
const root = path.join(__dirname, '..');
const html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');
const inline = html.slice(html.lastIndexOf('<script>') + 8, html.lastIndexOf('</script>'));
const coreSrc = fs.readFileSync(path.join(root, 'core.js'), 'utf8');

function stubEl(id) {
    const L = {}, el = { id, innerHTML: '', innerText: '', textContent: '', className: '', style: {}, dataset: {}, children: [], disabled: false, title: '', value: '',
        clientWidth: 600, clientHeight: 300, parentElement: null, scrollTop: 0, scrollHeight: 0,
        classList: { toggle() {}, add() {}, remove() {}, contains() { return false; } },
        addEventListener(t, f) { (L[t] = L[t] || []).push(f); },
        appendChild(c) { c.parentElement = el; el.children.push(c); return c; }, insertBefore(c) { c.parentElement = el; el.children.unshift(c); return c; },
        removeChild() {}, get firstChild() { return el.children[0]; }, closest() { return null },
        querySelector() { return Object.assign(stubEl(), { parentElement: el }); },
        getContext() { return new Proxy({}, { get: (t, p) => p in t ? t[p] : () => {}, set: (t, p, v) => { t[p] = v; return true; } }); } };
    return el;
}
function loadPage() {
    const els = {};
    const ctx = vm.createContext({
        document: { getElementById: id => els[id] || (els[id] = stubEl(id)), createElement: () => stubEl(), body: stubEl(), activeElement: null, querySelectorAll: () => [] },
        window: { devicePixelRatio: 1 }, performance, console, setTimeout, clearTimeout, URL: {}, Blob: function () {}, requestAnimationFrame: () => {}
    });
    vm.runInContext(coreSrc, ctx); vm.runInContext(inline, ctx);
    return code => vm.runInContext(code, ctx);
}
async function makeCsv(model, seed) {
    const ev = loadPage();
    ev(`Sys.model = '${model}'; Sys.isRunning = false; applySeed(${seed}); Sys.resetCache(); setAlgo('SMI'); Sys.computeMath();`);
    ev('Bench.start()');
    while (ev('Bench.running')) await new Promise(r => setTimeout(r, 20));
    return { csv: ev('Exporter.buildCsv()'), ev };
}
const strip = csv => csv.split('\r\n').filter(l => !l.startsWith('meta,timestamp,')).join('\r\n');

module.exports = {
    id: 'T10', title: 'seed reproducibility of the sweep and of the CSV; per-point re-run; CSV schema v2',
    async run() {
        const checks = [], REQ = ['model', 'train_mode', 'gamma_rel_dB', 'd_min', 'M', 'K', 'sigma_theta', 'sigma_phi', 'P_s', 'tau', 'theta1', 'theta2', 'fc_GHz', 'delta_f_kHz', 'N', 'L', 'seed', 'point_seed',
            'max_angle_drift_deg', 'SINR_opt_dB', 'algorithm', 'algorithm_legacy'];
        for (const model of ['legacy', 'unified']) {
            const A = await makeCsv(model, 12345), B = await makeCsv(model, 12345), C = await makeCsv(model, 54321);
            const same = strip(A.csv) === strip(B.csv), tsRows = A.csv.split('\r\n').filter(l => l.startsWith('meta,timestamp,')).length;
            console.log(`T10a [${model}] CSV length ${A.csv.length}, identical apart from the timestamp row: ${same}, timestamp rows: ${tsRows}`);
            checks.push(U.check(`T10a ${model}: same seed -> byte-identical CSV (timestamp row excluded)`, `${same}, ${tsRows} timestamp row`, 'true, 1', same && tsRows === 1));
            const diff = strip(A.csv) !== strip(C.csv);
            checks.push(U.check(`T10b ${model}: different seed -> different CSV`, `${diff}`, 'true', diff));
            // T10c: re-run sweep point index 12 (v = 120 km/h), SMI (algorithm index 2) alone
            const ev = A.ev, ks = ev('Bench.ALGOS.map(a => a[0])'), r = ev('Bench.result'), pi = 12, ai = ks.indexOf('SMI');
            const one = ev(`Bench.point(Bench.shadow('SMI'), ${r.v[pi]}, ${pi}, ${ai})`);
            const dEvm = Math.abs(one.evm - r.curves.SMI.evm[pi]), dSinr = Math.abs(one.sinr - r.curves.SMI.sinr[pi]);
            console.log(`T10c [${model}] point ${pi} (v = ${r.v[pi]} km/h) SMI re-run alone: evm ${one.evm} vs ${r.curves.SMI.evm[pi]}, sinr ${one.sinr} vs ${r.curves.SMI.sinr[pi]}`);
            checks.push(U.check(`T10c ${model}: single sweep point re-run alone gives the same numbers`, `|d evm| = ${U.e(dEvm)}, |d sinr| = ${U.e(dSinr)}`, '0 (bit-identical)', dEvm === 0 && dSinr === 0));
            const miss = REQ.filter(k => !new RegExp(`(^|\r\n)param,${k},`).test(A.csv) && !new RegExp(`(^|\r\n)result,${k},`).test(A.csv));
            const ver = /meta,csv_schema_version,2,/.test(A.csv), sw = /sweep_info,train=(withSignal|signalFree),/.test(A.csv), pcol = /sweep,velocity_kmh,point_seed,/.test(A.csv);
            console.log(`T10d [${model}] schema v2: ${ver}, missing keys: [${miss.join(', ')}], sweep_info has train=: ${sw}, point_seed column: ${pcol}`);
            checks.push(U.check(`T10d ${model}: CSV schema v2 keys present`, `missing: [${miss.join(', ')}]`, 'none missing, version 2, train= and point_seed present', miss.length === 0 && ver && sw && pcol));
        }
        return { id: this.id, title: this.title, checks };
    }
};
