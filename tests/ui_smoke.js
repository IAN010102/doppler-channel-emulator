'use strict';
/**
 * UI smoke test: loads core.js + the inline script of index.html in a vm context with a stubbed DOM, runs animation frames,
 * switches every algorithm / modulation / pause state and checks that nothing throws and the outputs stay finite.
 *
 *   node tests/ui_smoke.js
 */
const fs = require('fs');
const vm = require('vm');
const path = require('path');

const root = path.join(__dirname, '..');
const html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');
const inline = html.slice(html.lastIndexOf('<script>') + 8, html.lastIndexOf('</script>'));
const coreSrc = fs.readFileSync(path.join(root, 'core.js'), 'utf8');

function stubEl(id) {
    const listeners = {};
    const el = { id, innerHTML: '', innerText: '', textContent: '', className: '', style: {}, dataset: {}, children: [], disabled: false, title: '', value: '',
        clientWidth: 600, clientHeight: 300, parentElement: null, scrollTop: 0, scrollHeight: 0,
        classList: { toggle() {}, add() {}, remove() {}, contains() { return false; } },
        addEventListener(t, f) { (listeners[t] = listeners[t] || []).push(f); },
        appendChild(c) { c.parentElement = el; el.children.push(c); return c; }, insertBefore(c) { c.parentElement = el; el.children.unshift(c); return c; },
        removeChild(c) { el.children.splice(el.children.indexOf(c), 1); }, get firstChild() { return el.children[0]; }, closest() { return null; },
        querySelector() { return Object.assign(stubEl(), { parentElement: el }); },
        getContext() { return new Proxy({}, { get: (t, p) => p in t ? t[p] : () => {}, set: (t, p, v) => { t[p] = v; return true; } }); },
        fire(t, ev) { (listeners[t] || []).forEach(f => f(ev || { target: el })); } };
    return el;
}
const els = {};
let pending = null, now = 1000;
const ctx = vm.createContext({
    document: { getElementById: id => els[id] || (els[id] = stubEl(id)), createElement: () => stubEl(), body: stubEl(), activeElement: null, querySelectorAll: () => [] },
    window: { devicePixelRatio: 1 }, performance, console, setTimeout, clearTimeout, URL: {}, Blob: function () {},
    requestAnimationFrame: f => { pending = f; }
});
vm.runInContext(coreSrc, ctx);                       // UMD: no `module` in the context -> defines the global `Core`
vm.runInContext(inline, ctx);                        // the page's own script (top-level const/function stay reachable below)
const ev = code => vm.runInContext(code, ctx);
const frames = n => { for (let i = 0; i < n; i++) { const f = pending; pending = null; now += 50; if (f) f(now); } };

let fails = 0;
const check = (name, ok, extra = '') => { if (!ok) fails++; console.log((ok ? 'PASS  ' : 'FAIL  ') + name + (extra ? '  ' + extra : '')); };

frames(30);
check('page script loaded, frames ran', ev('Sys.simTime') > 0, 'simTime=' + ev('Sys.simTime').toFixed(2));
const finite = () => ev('[Sys.sinrDb, Sys.evm, Sys.ser, Sys.kappaRaw === Infinity ? 0 : Sys.kappaRaw, Sys.psll].every(Number.isFinite) || Sys.psll === -Infinity');
for (const algo of ['FOURIER', 'MMSE', 'MMSEP', 'SMI', 'DL', 'BEAMSPACE']) {
    for (const L of [100, 4]) {
        ev(`setAlgo('${algo}'); Sys.L = ${L}; Sys.dirty = true;`); frames(4);
        check(`algo ${algo} L=${L}`, finite(), `SINR=${ev('Sys.sinrDb').toFixed(1)} dB  status=${ev('Sys.status')}`);
    }
}
for (const mod of ['QPSK', 'QAM16', 'QAM64']) { ev(`Sys.mod = '${mod}'; Sys.dirty = true;`); frames(3); check(`modulation ${mod}`, finite(), `SER=${ev('Sys.ser').toExponential(1)}`); }
ev('setRunning(false)'); const t0 = ev('Sys.simTime'); frames(10); check('pause freezes the clock', ev('Sys.simTime') === t0);
ev('setRunning(true)'); frames(5); check('resume advances the clock', ev('Sys.simTime') > t0);
ev('Sys.resetCache()'); frames(3); check('reset clears the window and time', ev('Sys.simTime') < 1);
ev('Bench.start()');
check('benchmark sweep starts', ev('Bench.running') === true);
console.log(fails ? `\n${fails} check(s) FAILED` : '\nall UI smoke checks passed');
process.exit(fails ? 1 : 0);
