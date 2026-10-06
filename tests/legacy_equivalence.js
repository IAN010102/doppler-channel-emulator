'use strict';
/**
 * Commit-1 check: the refactored core.js must reproduce the pre-refactor index.html (commit 66d29ec) bit for bit
 * when both are driven by the same seeded generator (mulberry32) and the same parameters.
 *
 *   node tests/legacy_equivalence.js [legacy-git-ref]        (default 66d29ec)
 *
 * The legacy page is read straight from git, its inline script is executed inside a vm context with a stubbed DOM, and its
 * Math.random is replaced by Core.mulberry32(seed). The new code uses Core.setSeed(seed). Both consume the generator in the
 * same order, so every output must be identical.
 */
const { execSync } = require('child_process');
const vm = require('vm');
const path = require('path');
const Core = require('../core.js');

const ref = process.argv[2] || '66d29ec';
const html = execSync(`git show ${ref}:index.html`, { cwd: path.join(__dirname, '..'), maxBuffer: 1 << 26 }).toString('utf8');
const script = html.slice(html.indexOf('<script>') + 8, html.lastIndexOf('</script>'));

// ---- minimal DOM stub (enough for the legacy script to load; the physics never touches it)
function stubEl() {
    const el = { innerHTML: '', innerText: '', textContent: '', className: '', style: {}, dataset: {}, children: [], disabled: false, title: '',
        clientWidth: 600, clientHeight: 300, parentElement: null, scrollTop: 0, scrollHeight: 0, value: '',
        classList: { toggle() {}, add() {}, remove() {}, contains() { return false; } },
        addEventListener() {}, appendChild(c) { c.parentElement = el; el.children.push(c); return c; }, insertBefore(c) { el.children.unshift(c); return c; },
        removeChild() {}, get firstChild() { return el.children[0]; }, closest() { return null; },
        querySelector() { const c = stubEl(); c.parentElement = el; return c; }, getContext() { return new Proxy({}, { get: (t, p) => p in t ? t[p] : () => {}, set: (t, p, v) => { t[p] = v; return true; } }); } };
    return el;
}
const els = {};
const ctx = vm.createContext({
    document: { getElementById: id => els[id] || (els[id] = stubEl()), createElement: () => stubEl(), body: stubEl(), activeElement: null, querySelectorAll: () => [] },
    window: { devicePixelRatio: 1 }, performance, console, setTimeout, clearTimeout, requestAnimationFrame: () => {}, URL: {}, Blob: function () {}
});
vm.runInContext(script.replace(/\bconst (Sys)\b/g, 'var $1'), ctx);   // only Sys is needed from the legacy script
const legacy = ctx.Sys;

// ---- scenarios
const base = { N: 8, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, v: 0, calDeg: 0, latMs: 0, kDb: 20, sll: 35, gammaDL: 0.01, taper: 'NONE', mod: 'QPSK', algo: 'SMI', smiSingular: 'clamp', model: 'legacy', trainMode: 'withSignal', jamWave: 'qpsk' };   // QPSK: the old page always sent QPSK target symbols   // clamp = the behaviour of the pre-refactor page (the new default is pinv)
const scenarios = [];
for (const algo of ['FOURIER', 'MMSE', 'SMI', 'DL', 'BEAMSPACE']) scenarios.push({ name: `${algo} default`, p: { algo } });
scenarios.push({ name: 'SMI L=4 (singular, forced inverse)', p: { algo: 'SMI', L: 4 } });
scenarios.push({ name: 'MMSE L=4 (singular)', p: { algo: 'MMSE', L: 4 } });
scenarios.push({ name: 'BEAMSPACE L=2 (< K)', p: { algo: 'BEAMSPACE', L: 2 } });
scenarios.push({ name: 'SMI + Chebyshev (GSC)', p: { algo: 'SMI', taper: 'CHEBYSHEV' } });
scenarios.push({ name: 'DL + Hamming (GSC)', p: { algo: 'DL', taper: 'HAMMING', gammaDL: 0.1 } });
scenarios.push({ name: 'FOURIER + Chebyshev', p: { algo: 'FOURIER', taper: 'CHEBYSHEV' } });
scenarios.push({ name: 'all impairments N=12', p: { algo: 'SMI', N: 12, L: 60, aoaT: 10, aoaJ: -35, calDeg: 10, latMs: 5, v: 300, kDb: 5, snr: 25, sir: -15 } });
scenarios.push({ name: 'DL N=16 L=200 K=-10 dB', p: { algo: 'DL', N: 16, L: 200, kDb: -10, v: 500, latMs: 10 } });
scenarios.push({ name: 'BEAMSPACE N=16 L=5', p: { algo: 'BEAMSPACE', N: 16, L: 5, aoaT: 30, aoaJ: -50 } });

const calZ = Array.from({ length: 16 }, (_, i) => (i === 0 ? 0 : Math.sin(1.7 * i + 0.3) * 1.3));   // fixed calibration draws
const seeds = [1, 20260101, 0xDEADBEEF];

function run(sys, seedFn, p, seed) {
    Object.assign(sys, base, p, { calZ: calZ.slice(), snaps: [], snapKey: '', dirty: true });
    seedFn(seed);
    sys.computeMath();                    // first call: fills the whole window
    const first = snapshot(sys);
    sys.computeMath();                    // second call: sliding-window refresh path (REFRESH)
    return [first, snapshot(sys)];
}
const cplxList = a => a.map(c => [c.r, c.i]);
function snapshot(s) {
    return {
        weights: cplxList(s.weights), pattern: s.pattern.slice(), eig: s.eig.slice(), diagR: s.diagR.slice(),
        R_hat: s.R_hat.map(row => cplxList(row)),
        scalars: { sinrDb: s.sinrDb, evm: s.evm, ser: s.ser, S: s.S, I: s.I, Nn: s.Nn, nuICI: s.nuICI, nIciDb: s.nIciDb, kappa: s.kappa, kappaRaw: s.kappaRaw,
            nullDb: s.nullDb, psll: s.psll, bw3: s.bw3, outGain: s.outGain, fm: s.fm, fd: s.fd, eps0: s.eps0, epsM: s.epsM, isDb: s.isDb, lambdaQ: s.lambdaQ,
            gammaUsed: s.gammaUsed, delta: s.delta, thTo: s.thTo, thJo: s.thJo, dTdeg: s.dTdeg, dJdeg: s.dJdeg },
        status: s.status, collapsed: s.collapsed, bsBins: s.bsBins.slice(), bsAngles: s.bsAngles.slice()
    };
}
function flat(o, out = [], pre = '') {
    if (Array.isArray(o)) o.forEach((v, i) => flat(v, out, pre + '[' + i + ']'));
    else if (o && typeof o === 'object') Object.keys(o).forEach(k => flat(o[k], out, pre + '.' + k));
    else out.push([pre, o]);
    return out;
}

let nVals = 0, nIdentical = 0, maxAbs = 0, worst = '';
const rows = [];
for (const sc of scenarios) {
    let sMax = 0, sId = 0, sN = 0;
    for (const seed of seeds) {
        const a = run(legacy, s => { ctx.__r = Core.mulberry32(s); vm.runInContext('Math.random = __r', ctx); }, sc.p, seed);
        const sysNew = Core.createSys();
        const b = run(sysNew, s => Core.setSeed(s), sc.p, seed);
        for (let k = 0; k < 2; k++) {
            const fa = flat(a[k]), fb = flat(b[k]);
            if (fa.length !== fb.length) throw new Error(`${sc.name}: shape mismatch ${fa.length} vs ${fb.length}`);
            for (let i = 0; i < fa.length; i++) {
                const [ka, va] = fa[i], [kb, vb] = fb[i];
                if (ka !== kb) throw new Error(`${sc.name}: key mismatch ${ka} vs ${kb}`);
                sN++; nVals++;
                if (Object.is(va, vb)) { sId++; nIdentical++; }
                else if (typeof va === 'number' && typeof vb === 'number') {
                    const d = (Number.isFinite(va) && Number.isFinite(vb)) ? Math.abs(va - vb) : Infinity;
                    if (d > sMax) sMax = d;
                    if (d > maxAbs) { maxAbs = d; worst = `${sc.name} ${ka} legacy=${va} new=${vb}`; }
                } else { sMax = Infinity; maxAbs = Infinity; worst = `${sc.name} ${ka} legacy=${va} new=${vb}`; }
            }
        }
    }
    rows.push({ name: sc.name, values: sN, identical: sId, maxAbs: sMax });
}
console.log(`legacy reference: ${ref}   seeds: ${seeds.join(', ')}   (each scenario: full-window call + sliding-window call)`);
console.log('scenario'.padEnd(40), 'values'.padStart(8), 'bit-identical'.padStart(14), 'max |diff|'.padStart(12));
for (const r of rows) console.log(r.name.padEnd(40), String(r.values).padStart(8), String(r.identical).padStart(14), r.maxAbs.toExponential(2).padStart(12));
console.log(`\nTOTAL values compared: ${nVals}, bit-identical: ${nIdentical}, max |diff| = ${maxAbs.toExponential(2)}`);
const pass = maxAbs < 1e-12;
console.log(pass ? 'PASS  (max |diff| < 1e-12)' : `FAIL  worst: ${worst}`);
process.exit(pass ? 0 : 1);
