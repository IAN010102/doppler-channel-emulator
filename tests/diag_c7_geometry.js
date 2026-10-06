'use strict';
/**
 * Commit 7 / Commit 9 before/after (informational): mean SINR of SMI and DL at v = 300 km/h, theta_1 = 45 deg, theta_2 = -30 deg, K = 20 dB, L = 100,
 * d_min = 30 / 50 / 500 m, tau = 0 / 2 ms.
 *   before = a previous core.js passed on the command line (Commit 7: exact angle drift, Doppler phase 2 pi f_d(t_app) t, first order;
 *            git show 043b318:core.js > old_core.js)
 *   after  = current core.js (Commit 9: integrated Doppler phase)
 * Paired trials (same seed for every d_min), 2000 trials. The DL loading is the same absolute value everywhere (0.01 = gamma_rel 0 dB at SNR 20 dB).
 *   node tests/diag_c7_geometry.js <path to the previous core.js> [trials]
 */
const path = require('path');
const U = require('./_util.js');
const NEW = require('../core.js');
const OLDPATH = process.argv[2], trials = +(process.argv[3] || 2000);
const OLD = OLDPATH ? require(path.resolve(OLDPATH)) : null;

function run(Core, which, algo, dmin, tauMs, trainMode, seed) {
    const out = [];
    for (let t = 0; t < trials; t++) {
        Core.setSeed(seed + t);
        const s = Core.createSys(); s.calZ = Array.from({ length: 16 }, () => 0);
        Object.assign(s, { model: 'unified', freshRealization: true, algo, N: 8, L: 100, aoaT: 45, aoaJ: -30, v: 300, snr: 20, sir: -10, kDb: 20, latMs: tauMs, calDeg: 0, gammaDL: 0.01, gammaRelDb: 0, trainMode });
        s.d_min = dmin;
        s.snaps = []; s.computeMath(); out.push(s.sinrDb);
    }
    return out;
}
const cell = (a) => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`;
for (const tau of [0, 2]) {
    console.log(`\ntau = ${tau} ms, v = 300 km/h, theta1 = 45, theta2 = -30, K = 20 dB, L = 100, ${trials} paired trials: mean SINR dB ± SE`);
    console.log(U.pad('version / training', 24), U.pad('algo', 5), U.rpad('d_min=30', 16), U.rpad('d_min=50', 16), U.rpad('d_min=500', 16), U.rpad('50 - 500', 16));
    const rows = [];
    if (OLD) rows.push(['before (Commit 7)', OLD, 'old', 'withSignal'], ['before, signalFree', OLD, 'old', 'signalFree']);
    rows.push(['after, withSignal', NEW, 'new', 'withSignal'], ['after, signalFree', NEW, 'new', 'signalFree']);
    for (const [label, Core, which, tm] of rows) for (const algo of ['SMI', 'DL']) {
        const z = run(Core, which, algo, 30, tau, tm, 1000), a = run(Core, which, algo, 50, tau, tm, 1000), b = run(Core, which, algo, 500, tau, tm, 1000);
        const d = a.map((x, i) => x - b[i]);
        console.log(U.pad(label, 24), U.pad(algo, 5), U.rpad(cell(z), 16), U.rpad(cell(a), 16), U.rpad(cell(b), 16), U.rpad(cell(d), 16));
    }
}
