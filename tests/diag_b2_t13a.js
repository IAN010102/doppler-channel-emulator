'use strict';
/**
 * Diagnosis B2 (informational): the T13a cell signalFree, v = 300 (mean of EVM_meas / EVM_analytic - 1 = +0.1957 % with SE 0.0526 %, |z| = 3.72 at Commit 19 and later; +0.0708 % before Commit 19).
 * Ablation, N realisations per combination (default 400; SMI, K = 20 dB, theta1 = 0, theta2 = 40, L = 100, Ns = 4000 symbols per realisation, seeds 1313 + t):
 *   jammer waveform at symbol level (gaussian | qpsk)  x  transmit modulation (QPSK | 16-QAM)  x  v (0 | 300)  x  ICI term (on | off).
 * "ICI off": the symbol chain is run with nu = 0 and compared with sqrt(1/SINR) (the analytic value without the ICI term).
 * Also: the largest cell (16-QAM, gaussian, v = 300, ICI on) with 2000 realisations, and the dependence on the number of symbols per realisation.
 * usage: node tests/diag_b2_t13a.js [N]
 */
const Core = require('../core.js');
const U = require('./_util.js');
const N = +(process.argv[2] || 400);
const BASE = { model: 'unified', freshRealization: true, algo: 'SMI', N: 8, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, latMs: 0, calDeg: 0, taper: 'NONE', kDb: 20, mod: 'QAM16', gammaRelDb: 10, d_min: 30, trainMode: 'signalFree' };
function cell(over, ici, jam, n, Ns = 4000, seed0 = 1313) {
    const rel = [];
    for (let t = 0; t < n; t++) {
        Core.setSeed(seed0 + t);
        const s = Core.createSys(); s.calZ = new Array(16).fill(0); Object.assign(s, BASE, over); s.snaps = []; s.computeMath();
        const r = Core.symbolLevel({ S: s.S, I: s.I, Nn: s.Nn, nu: ici ? s.nuICI : 0, mod: s.mod, jam, Ns, keep: 0 });
        const an = ici ? s.evm : Math.sqrt((s.I + s.Nn) / s.S);
        rel.push(r.evm / an - 1);
    }
    return { m: U.mean(rel), se: U.se(rel) };
}
// blocks: the T13a cell (16-QAM, gaussian jammer, v = 300, ICI on) in consecutive blocks of 300 realisations (seeds 1313 + 300 k + t); block 0 is the T13a test itself
if (process.argv.includes('--blocks') || process.argv.includes('--only-blocks')) {
    console.log('\nblocks of 300 realisations (seed0 = 1313 + 300 k)');
    const all = [];
    for (let k = 0; k < 10; k++) { const c = cell({ mod: 'QAM16', v: 300 }, true, 'gaussian', 300, 4000, 1313 + 300 * k); all.push(c.m); console.log(`  block ${k}: ${(100 * c.m).toFixed(4)} % ± ${(100 * c.se).toFixed(4)} %  z = ${(c.m / c.se).toFixed(2)}`); }
    console.log(`  mean of the 10 block means ${(100 * U.mean(all)).toFixed(4)} % ± ${(100 * U.se(all)).toFixed(4)} % (from the spread of the blocks)`);
    if (process.argv.includes('--only-blocks')) process.exit(0);
}

console.log(`B2 ablation, ${N} realisations per combination; mean of EVM_meas/EVM_analytic - 1 [%] ± SE [%] (|z|)`);
console.log(U.pad('jammer', 9), U.pad('Tx mod', 7), U.pad('v', 4), U.pad('ICI', 4), U.rpad('rel diff %', 11), U.rpad('SE %', 8), U.rpad('|z|', 6));
for (const jam of ['gaussian', 'qpsk']) for (const mod of ['QPSK', 'QAM16']) for (const v of [0, 300]) for (const ici of [true, false]) {
    const c = cell({ mod, v }, ici, jam, N);
    console.log(U.pad(jam, 9), U.pad(mod, 7), U.pad(v, 4), U.pad(ici ? 'on' : 'off', 4), U.rpad((100 * c.m).toFixed(4), 11), U.rpad((100 * c.se).toFixed(4), 8), U.rpad((Math.abs(c.m) / c.se).toFixed(2), 6));
}
console.log('\nlarge sample: 16-QAM, gaussian, v = 300, ICI on, seeds 1313+t');
for (const n of [300, 2000]) { const c = cell({ mod: 'QAM16', v: 300 }, true, 'gaussian', n); console.log(`  ${n} realisations: ${(100 * c.m).toFixed(4)} % ± ${(100 * c.se).toFixed(4)} % (|z| ${(Math.abs(c.m) / c.se).toFixed(2)})`); }
console.log('\nsymbols per realisation (16-QAM, gaussian, v = 300, ICI on, 400 realisations)');
for (const Ns of [1000, 4000, 16000]) { const c = cell({ mod: 'QAM16', v: 300 }, true, 'gaussian', 400, Ns); console.log(`  Ns = ${Ns}: ${(100 * c.m).toFixed(4)} % ± ${(100 * c.se).toFixed(4)} %`); }
