'use strict';
/**
 * T18  DOA spectra (Commit 21): Capon P = 1/(a^H R^-1 a) (pseudo-inverse if R_hat is rank deficient) and MUSIC (2 assumed sources: target and jammer), 0.5 deg grid.
 *   T18a  no noise (SNR 400 dB), K -> infinity (pure LoS target and jammer), withSignal, L = 200, three well separated pairs (theta1, theta2) = (-20, 40), (10, -50), (-5, 25): the MUSIC peak
 *         nearest to each true angle is within 0.5 deg (the grid step; the maximum possible grid error is 0.25 deg).
 *   T18b  number of peaks = number of sources: withSignal, SNR 40 dB, K -> infinity, L = 200: the MUSIC spectrum has exactly 2 peaks within 30 dB of its maximum (same count reported for Capon).
 *   T18c  (informational) estimation error vs L (4 ... 100) and vs the angle separation (theta1 = 0, theta2 = 5, 10, 20, 40 deg): mean absolute error of the MUSIC and Capon peak nearest to each
 *         true angle, K = 20 dB, SNR 20 dB, withSignal, 500 trials, with SE.
 */
const Core = require('../core.js');
const U = require('./_util.js');

const BASE = { N: 8, L: 200, snr: 400, sir: -10, v: 0, latMs: 0, calDeg: 0, taper: 'NONE', kDb: 400, mod: 'QPSK', algo: 'SMI', model: 'unified', freshRealization: true, trainMode: 'withSignal' };
function mk(over, seed) { Core.setSeed(seed); const s = Core.createSys(); s.calZ = new Array(16).fill(0); Object.assign(s, BASE, over); s.snaps = []; s.computeMath(); return s; }
const cell = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`;

module.exports = {
    id: 'T18', title: 'DOA spectra: MUSIC peak position (noise-free), peak count, error vs L and separation (info)',
    async run({ trials = 500, seed = 1818 } = {}) {
        const checks = [], info = [];
        // ---- T18a
        for (const [t1, t2] of [[-20, 40], [10, -50], [-5, 25]]) {
            const s = mk({ aoaT: t1, aoaJ: t2 }, seed), d = s.computeDoa(), e = d.err;
            const e1 = e.musicT ? Math.abs(e.musicT.deg - t1) : Infinity, e2 = e.musicJ ? Math.abs(e.musicJ.deg - t2) : Infinity;
            console.log(`T18a  (theta1, theta2) = (${t1}, ${t2}): rank ${d.rank}, MUSIC peaks ${d.peaksMusic.map(p => p.deg.toFixed(1)).join(', ')}; errors ${e1.toFixed(2)} / ${e2.toFixed(2)} deg`);
            checks.push(U.check(`T18a (${t1}, ${t2}): MUSIC peak error, target / jammer`, `${e1.toFixed(2)} / ${e2.toFixed(2)} deg`, '< 0.5 deg', e1 < 0.5 && e2 < 0.5));
        }
        // ---- T18b
        for (const [t1, t2] of [[-20, 40], [10, -50]]) {
            const s = mk({ aoaT: t1, aoaJ: t2, snr: 40 }, seed + 1), d = s.computeDoa();
            const nm = d.peaksMusic.filter(p => p.db >= -30).length, nc = d.peaksCapon.filter(p => p.db >= -30).length;
            console.log(`T18b  (${t1}, ${t2}), SNR 40 dB: peaks within 30 dB of the maximum: MUSIC ${nm}, Capon ${nc}`);
            checks.push(U.check(`T18b (${t1}, ${t2}): number of MUSIC peaks = number of sources`, `${nm}`, '2', nm === 2));
            info.push({ name: `T18b (${t1}, ${t2}) Capon peak count`, value: `${nc}`, note: 'Capon with a rank-2 + noise covariance' });
        }
        // ---- T18c (informational)
        console.log(`\nT18c  K = 20 dB, SNR 20 dB, withSignal, ${trials} trials; mean |error| [deg] of the peak nearest to the true angle (MUSIC target / MUSIC jammer / Capon target / Capon jammer)`);
        const f = (p, tru) => p ? Math.abs(p.deg - tru) : 90;
        console.log('vs L (theta1 = 0, theta2 = 40):');
        console.log(U.pad('L', 6), ['MUSIC T', 'MUSIC J', 'Capon T', 'Capon J'].map(x => U.rpad(x, 16)).join(''));
        for (const L of [4, 8, 12, 24, 48, 100]) {
            const a = [[], [], [], []];
            for (let t = 0; t < trials; t++) { const s = mk({ aoaT: 0, aoaJ: 40, L, snr: 20, kDb: 20 }, seed + 100 + t), d = s.computeDoa(), e = d.err; a[0].push(f(e.musicT, 0)); a[1].push(f(e.musicJ, 40)); a[2].push(f(e.caponT, 0)); a[3].push(f(e.caponJ, 40)); }
            console.log(U.pad(L, 6), a.map(x => U.rpad(cell(x), 16)).join(''));
            info.push({ name: `T18c L=${L}`, value: `MUSIC T ${cell(a[0])}, J ${cell(a[1])}; Capon T ${cell(a[2])}, J ${cell(a[3])} deg` });
        }
        console.log('vs separation (theta1 = 0, L = 100):');
        console.log(U.pad('theta2', 8), ['MUSIC T', 'MUSIC J', 'Capon T', 'Capon J'].map(x => U.rpad(x, 16)).join(''));
        for (const t2 of [5, 10, 20, 40]) {
            const a = [[], [], [], []];
            for (let t = 0; t < trials; t++) { const s = mk({ aoaT: 0, aoaJ: t2, L: 100, snr: 20, kDb: 20 }, seed + 200 + t), d = s.computeDoa(), e = d.err; a[0].push(f(e.musicT, 0)); a[1].push(f(e.musicJ, t2)); a[2].push(f(e.caponT, 0)); a[3].push(f(e.caponJ, t2)); }
            console.log(U.pad(t2, 8), a.map(x => U.rpad(cell(x), 16)).join(''));
            info.push({ name: `T18c separation ${t2} deg`, value: `MUSIC T ${cell(a[0])}, J ${cell(a[1])}; Capon T ${cell(a[2])}, J ${cell(a[3])} deg` });
        }
        return { id: this.id, title: this.title, checks, info };
    }
};
