'use strict';
/**
 * T5  Time-domain OFDM verification (independent of core.js: this is a plain OFDM chain).
 *     N_FFT = 64, cyclic prefix, no noise, frequency offset eps*df, unitary FFTs.
 *
 *   T5a  one path with offset eps: measured leakage fraction (power of Y_k not carried by the single-tap term, over the total power
 *        of Y_k, random QPSK on all subcarriers) vs  1 - (sin(pi eps) / (N sin(pi eps/N)))^2.   Required error < 1 %.
 *   T5b  report: finite-N formula vs the sinc^2 used by the UI, eps <= 0.3.
 *   T5c  two paths (different Doppler, different power), >= 2000 random phase pairs: measured leakage fraction vs the
 *        non-coherent approximation  sum_i q_i (1 - |D_i|^2),  q_i = P_i / sum P.   Required error < 3 %.
 *        (The UI formula, with sinc^2 instead of |D|^2, is shown as well.)
 *
 *   The single-tap gain of each realisation is obtained from the same chain with a one-subcarrier probe symbol (no estimation
 *   noise); the residual Y_k - H*X_k is the ICI.
 */
const U = require('./_util.js');
const Core = require('../core.js');

const N = 64;
const CP = Math.round(Core.CONFIG.cpRatio * N);        // 0.07 * 64 = 4.48 -> 4 samples

// ---- tiny radix-2 FFT (unitary)
function fft(re, im, inverse) {
    const n = re.length, outR = new Float64Array(re), outI = new Float64Array(im);
    for (let i = 1, j = 0; i < n; i++) {                                   // bit reversal
        let bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit;
        if (i < j) { [outR[i], outR[j]] = [outR[j], outR[i]]; [outI[i], outI[j]] = [outI[j], outI[i]]; }
    }
    for (let len = 2; len <= n; len <<= 1) {
        const ang = 2 * Math.PI / len * (inverse ? 1 : -1), wr = Math.cos(ang), wi = Math.sin(ang);
        for (let i = 0; i < n; i += len) {
            let cr = 1, ci = 0;
            for (let k = 0; k < len / 2; k++) {
                const a = i + k, b = i + k + len / 2;
                const tr = outR[b] * cr - outI[b] * ci, ti = outR[b] * ci + outI[b] * cr;
                outR[b] = outR[a] - tr; outI[b] = outI[a] - ti; outR[a] += tr; outI[a] += ti;
                const ncr = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = ncr;
            }
        }
    }
    const s = 1 / Math.sqrt(n);
    for (let i = 0; i < n; i++) { outR[i] *= s; outI[i] *= s; }
    return [outR, outI];
}

// one OFDM symbol through the chain.  paths: [{ gr, gi (complex gain), eps }]  ->  Y (frequency domain)
function chain(Xr, Xi, paths) {
    const [xr, xi] = fft(Xr, Xi, true);                                    // IFFT
    const L = N + CP, txr = new Float64Array(L), txi = new Float64Array(L);
    for (let m = 0; m < L; m++) { const src = m < CP ? N - CP + m : m - CP; txr[m] = xr[src]; txi[m] = xi[src]; }
    const rr = new Float64Array(N), ri = new Float64Array(N);
    for (let m = CP; m < L; m++) {                                         // receiver drops the CP; the CFO ramp runs over CP + symbol
        let ar = 0, ai = 0;
        for (const p of paths) {
            const ph = 2 * Math.PI * p.eps * m / N, c = Math.cos(ph), s = Math.sin(ph);
            const gr = p.gr * c - p.gi * s, gi = p.gr * s + p.gi * c;
            ar += gr * txr[m] - gi * txi[m]; ai += gr * txi[m] + gi * txr[m];
        }
        rr[m - CP] = ar; ri[m - CP] = ai;
    }
    return fft(rr, ri, false);
}

const Dsq = eps => { const x = Math.sin(Math.PI * eps) / (N * Math.sin(Math.PI * eps / N)); return x * x; };
const sincSq = eps => { const x = eps === 0 ? 1 : Math.sin(Math.PI * eps) / (Math.PI * eps); return x * x; };

function qpskSymbol(rnd) {
    const Xr = new Float64Array(N), Xi = new Float64Array(N), a = Math.SQRT1_2;
    for (let k = 0; k < N; k++) { Xr[k] = rnd() < 0.5 ? -a : a; Xi[k] = rnd() < 0.5 ? -a : a; }
    return [Xr, Xi];
}
function probeGain(paths) {                      // single-tap gain: one-subcarrier probe at k0, read the same bin
    const k0 = N / 2, Xr = new Float64Array(N), Xi = new Float64Array(N); Xr[k0] = 1;
    const [Yr, Yi] = chain(Xr, Xi, paths); return [Yr[k0], Yi[k0]];
}
// ICI power and total power of Y for one random symbol, given the single-tap gain H
function icOnce(rnd, paths, H) {
    const [Xr, Xi] = qpskSymbol(rnd), [Yr, Yi] = chain(Xr, Xi, paths);
    let ici = 0, tot = 0;
    for (let k = 0; k < N; k++) {
        const er = Yr[k] - (H[0] * Xr[k] - H[1] * Xi[k]), ei = Yi[k] - (H[0] * Xi[k] + H[1] * Xr[k]);
        ici += er * er + ei * ei; tot += Yr[k] * Yr[k] + Yi[k] * Yi[k];
    }
    return [ici, tot];
}

module.exports = {
    id: 'T5', title: 'time-domain OFDM chain: leakage vs 1 - |D|^2 (single path) and the non-coherent two-path approximation',
    async run({ symbols = 4000, realisations = 4000, seed = 20260505 } = {}) {
        const checks = [];
        const rnd = Core.mulberry32(seed);
        console.log(`N_FFT = ${N}, CP = ${CP} samples (cpRatio ${Core.CONFIG.cpRatio} x ${N} = ${(Core.CONFIG.cpRatio * N).toFixed(2)}, rounded), no noise`);

        // ---------------------------------------------------------------- T5a
        console.log(`\nT5a  single path, ${symbols} random OFDM symbols x ${N} subcarriers per eps`);
        console.log(U.pad('eps', 6), U.rpad('measured', 11), U.rpad('1-|D|^2', 11), U.rpad('rel.err', 9), U.rpad('1-sinc^2', 11), U.rpad('sinc vs |D|', 12), 'result');
        for (const eps of [0.02, 0.05, 0.1, 0.15, 0.2, 0.3]) {
            const paths = [{ gr: 1, gi: 0, eps }], H = probeGain(paths);
            let ici = 0, tot = 0;
            for (let s = 0; s < symbols; s++) { const [a, b] = icOnce(rnd, paths, H); ici += a; tot += b; }
            const meas = ici / tot, theory = 1 - Dsq(eps), relErr = Math.abs(meas - theory) / theory;
            const pass = relErr < 0.01;
            console.log(U.pad(eps, 6), U.rpad(U.e(meas), 11), U.rpad(U.e(theory), 11), U.rpad((100 * relErr).toFixed(2) + '%', 9), U.rpad(U.e(1 - sincSq(eps)), 11),
                U.rpad((100 * ((1 - sincSq(eps)) / theory - 1)).toFixed(4) + '%', 12), pass ? 'PASS' : 'FAIL');
            checks.push(U.check(`T5a eps=${eps}: measured leakage vs 1 - |D|^2 (N=${N})`, `${(100 * relErr).toFixed(2)} %`, '< 1 %', pass));
        }

        // ---------------------------------------------------------------- T5b
        console.log('\nT5b  finite-N (N = 64) vs sinc^2 used by the UI: relative difference of the leakage 1 - (.)');
        console.log('     ' + [0.05, 0.1, 0.15, 0.2, 0.25, 0.3].map(e => `eps=${e}: ${(100 * ((1 - sincSq(e)) / (1 - Dsq(e)) - 1)).toFixed(4)}%`).join('   '));

        // ---------------------------------------------------------------- T5c
        const P1 = 1.0, P2 = 0.4, e1 = 0.15, e2 = -0.10, perReal = 4;
        console.log(`\nT5c  two paths: P = ${P1}/${P2}, eps = ${e1}/${e2}, ${realisations} random phase pairs x ${perReal} symbols`);
        let ici = 0, tot = 0; const perR = [];
        for (let r = 0; r < realisations; r++) {
            const f1 = 2 * Math.PI * rnd(), f2 = 2 * Math.PI * rnd();
            const paths = [{ gr: Math.sqrt(P1) * Math.cos(f1), gi: Math.sqrt(P1) * Math.sin(f1), eps: e1 }, { gr: Math.sqrt(P2) * Math.cos(f2), gi: Math.sqrt(P2) * Math.sin(f2), eps: e2 }];
            const H = probeGain(paths); let a = 0, b = 0;
            for (let s = 0; s < perReal; s++) { const [x, y] = icOnce(rnd, paths, H); a += x; b += y; }
            ici += a; tot += b; perR.push(a / b);
        }
        const meas = ici / tot, q1 = P1 / (P1 + P2), q2 = P2 / (P1 + P2);
        const approx = q1 * (1 - Dsq(e1)) + q2 * (1 - Dsq(e2)), approxUI = q1 * (1 - sincSq(e1)) + q2 * (1 - sincSq(e2));
        const relErr = Math.abs(meas - approx) / approx, seRel = U.se(perR) / U.mean(perR);
        const pass = relErr < 0.03;
        console.log(`     measured ${U.e(meas)},  sum q_i (1-|D_i|^2) = ${U.e(approx)}  (rel. err ${(100 * relErr).toFixed(2)} %, sampling SE ~ ${(100 * seRel).toFixed(2)} %),  UI form with sinc^2 = ${U.e(approxUI)}  ${pass ? 'PASS' : 'FAIL'}`);
        checks.push(U.check('T5c two paths: measured leakage vs sum q_i (1-|D_i|^2)', `${(100 * relErr).toFixed(2)} %`, '< 3 %', pass,
            `${realisations} phase pairs; sampling SE ~ ${(100 * seRel).toFixed(2)} %`));
        return { id: this.id, title: this.title, checks };
    }
};
