'use strict';
/* ofdm_ici_multipath.js -- 頻率選擇性（多路徑延遲）的擴充（C2）。不改 ofdm_ici.js 的既有函式；只在 Node 端使用。
 * 模型（假設 A20–A24，見 docs/ofdm_ici_assumptions.md）：路徑 i 有複數振幅 A_i、到達角 θ_i、正規化都卜勒 ε_i、延遲 d_i（取樣點，d_i = τ_i·N·Δf，可為非整數）。
 *   接收訊號  y(t) = Σ_i A_i a(θ_i) exp(j2π ε_i t/N) x(t − d_i)，x(t) = Σ_m X_m exp(j2π m t/N)（CP 內為週期延拓；延遲小於 CP 時不產生符元間干擾，前一個符元的尾巴不模擬）。
 *   對頻 ε 後做 FFT：  Y_l = Σ_m h(l,m) X_m + 雜訊，  h(l,m) = Σ_i A_i a(θ_i) exp(−j2π m d_i/N) Q(m − l + ε_i − ε)    （M1；d_i = 0 時還原成 ofdm_ici.js 的 O2）
 *   因此 H_l = h(l,l)，子載波 l 的 ICI 協方差 B_l = Σ_{m≠l} h(l,m) h(l,m)ᴴ，R_l = σx² B_l + σn² I。 */
const O = require('./ofdm_ici.js');
const TWO_PI = 2 * Math.PI;
const DELAYS_US = [0, 1, 2.5], TAP_DB = [0, -3, -6];

/** 3 個分接（tap）× 文獻的 4 條都卜勒路徑 = 12 條子路徑；每個分接有自己的隨機相位，到達角與都卜勒沿用 4 條路徑。flat = true：同樣 12 條子路徑但延遲全為 0（對照用） */
function multipathChannel(angles, rng, o = {}) {
    const Nr = o.Nr, N = o.N || 1024, df = o.df || 15e3, delays = o.delaysUs || DELAYS_US, tapDb = o.tapDb || TAP_DB, paths = [], dl = [];
    delays.forEach((tau, k) => O.PAPER.dopplerHz.forEach((f, i) => { paths.push({ fHz: f, thetaDeg: angles[i], powerDb: O.PAPER.relPowerDb[i] + tapDb[k] }); dl.push(o.flat ? 0 : tau * 1e-6 * N * df); }));
    const ch = O.makeChannel({ N, Nr, df, rng, paths }); ch.paths.forEach((p, i) => { p.d = dl[i]; }); ch.maxDelaySamples = Math.max(...dl);
    return ch;
}

/** 時域符元：同 ofdm_ici.js 的 makeSymbol，但每條路徑用自己的延遲版本（FFT 相位斜坡，週期延拓） */
function makeSymbolMP(ch, o) {
    const { N, Nr, paths } = ch, rng = o.rng, sx2 = o.sigmaX2 === undefined ? 1 : o.sigmaX2, sn2 = o.sigmaN2, Ncp = o.Ncp, ps = o.pilotSpacing === undefined ? 12 : o.pilotSpacing;
    const Xre = new Float64Array(N), Xim = new Float64Array(N), a = Math.sqrt(sx2 / 2), pilots = [];
    for (let m = 0; m < N; m++) { Xre[m] = rng.u() < 0.5 ? -a : a; Xim[m] = rng.u() < 0.5 ? -a : a; if (ps > 0 && m % ps === 0) pilots.push(m); }
    const cache = new Map();
    const delayed = d => {
        if (!cache.has(d)) { const re = new Float64Array(N), im = new Float64Array(N); for (let m = 0; m < N; m++) { const ph = -TWO_PI * m * d / N, c = Math.cos(ph), s = Math.sin(ph); re[m] = Xre[m] * c - Xim[m] * s; im[m] = Xre[m] * s + Xim[m] * c; } O.fft(re, im, true); cache.set(d, { re, im }); }
        return cache.get(d);
    };
    const T = N + Ncp, rRe = [], rIm = [], sw = Math.sqrt(N * sn2 / 2), xs = paths.map(p => delayed(p.d));
    for (let p = 0; p < Nr; p++) { rRe.push(new Float64Array(T)); rIm.push(new Float64Array(T)); }
    for (let k = 0; k < T; k++) {
        const t = k - Ncp, n = ((t % N) + N) % N, hr = new Float64Array(Nr), hi = new Float64Array(Nr);
        paths.forEach((pa, i) => {
            const ph = TWO_PI * pa.eps * t / N, c = Math.cos(ph), s = Math.sin(ph), ar = pa.Are * c - pa.Aim * s, ai = pa.Are * s + pa.Aim * c, xr = xs[i].re[n], xi = xs[i].im[n], br = ar * xr - ai * xi, bi = ar * xi + ai * xr;
            for (let p = 0; p < Nr; p++) { hr[p] += br * pa.a.re[p] - bi * pa.a.im[p]; hi[p] += br * pa.a.im[p] + bi * pa.a.re[p]; }
        });
        for (let p = 0; p < Nr; p++) { rRe[p][k] = hr[p] + sw * rng.n(); rIm[p][k] = hi[p] + sw * rng.n(); }
    }
    return { N, Ncp, Nr, Xre, Xim, rRe, rIm, pilots, sigmaX2: sx2, sigmaN2: sn2 };
}

/** M1：給定對頻 ε，回傳子載波列表 ls 各自的 H_l；withB = true 時另回傳 B_l（ICI 協方差，不含雜訊） */
function mpSystem(ch, epsDemod, ls, withB) {
    const { N, Nr, paths } = ch, P = paths.length, q = [0, 0], out = [];
    const Hof = l => {
        const re = new Float64Array(Nr), im = new Float64Array(Nr);
        for (const pa of paths) { O.Qc(pa.eps - epsDemod, N, q); const ph = -TWO_PI * l * pa.d / N, c = Math.cos(ph), s = Math.sin(ph), wr = q[0] * c - q[1] * s, wi = q[0] * s + q[1] * c, cr = pa.Are * wr - pa.Aim * wi, ci = pa.Are * wi + pa.Aim * wr; for (let p = 0; p < Nr; p++) { re[p] += cr * pa.a.re[p] - ci * pa.a.im[p]; im[p] += cr * pa.a.im[p] + ci * pa.a.re[p]; } }
        return { re, im };
    };
    if (!withB) return ls.map(l => ({ l, H: Hof(l) }));
    const qt = new Float64Array(2 * N * P), cs = new Float64Array(2 * N * P);   // q[d][i]（d = (m − l) mod N）、exp(−j2π m d_i/N)
    for (let d = 0; d < N; d++) paths.forEach((pa, i) => { O.Qc(d + pa.eps - epsDemod, N, q); qt[2 * (d * P + i)] = q[0]; qt[2 * (d * P + i) + 1] = q[1]; });
    for (let m = 0; m < N; m++) paths.forEach((pa, i) => { const ph = -TWO_PI * m * pa.d / N; cs[2 * (m * P + i)] = Math.cos(ph); cs[2 * (m * P + i) + 1] = Math.sin(ph); });
    const v = { re: new Float64Array(Nr), im: new Float64Array(Nr) };
    for (const l of ls) {
        const B = new Float64Array(2 * Nr * Nr);
        let Hl = null;
        for (let m = 0; m < N; m++) {
            v.re.fill(0); v.im.fill(0); const d = ((m - l) % N + N) % N;
            for (let i = 0; i < P; i++) {
                const pa = paths[i], qr = qt[2 * (d * P + i)], qi = qt[2 * (d * P + i) + 1], cr0 = cs[2 * (m * P + i)], ci0 = cs[2 * (m * P + i) + 1], wr = qr * cr0 - qi * ci0, wi = qr * ci0 + qi * cr0, cr = pa.Are * wr - pa.Aim * wi, ci = pa.Are * wi + pa.Aim * wr;
                for (let p = 0; p < Nr; p++) { v.re[p] += cr * pa.a.re[p] - ci * pa.a.im[p]; v.im[p] += cr * pa.a.im[p] + ci * pa.a.re[p]; }
            }
            if (m === l) { Hl = { re: Float64Array.from(v.re), im: Float64Array.from(v.im) }; continue; }
            O.addOuter(B, v, 1, Nr);
        }
        out.push({ l, H: Hl, B });
    }
    return out;
}

/** 各導頻的原始通道估計 Ĥ_k = Y_k / X_k */
function rawPilotH(sym, Y) {
    const Nr = sym.Nr;
    return sym.pilots.map(l => { const d = sym.Xre[l] * sym.Xre[l] + sym.Xim[l] * sym.Xim[l], re = new Float64Array(Nr), im = new Float64Array(Nr); for (let p = 0; p < Nr; p++) { re[p] = (Y[p].re[l] * sym.Xre[l] + Y[p].im[l] * sym.Xim[l]) / d; im[p] = (Y[p].im[l] * sym.Xre[l] - Y[p].re[l] * sym.Xim[l]) / d; } return { re, im }; });
}
/** 相鄰兩個導頻之間的線性內插（對頻率）；pilots 為升冪，l 落在兩個導頻之間（含端點） */
function interp(raw, pilots, l) {
    let j = 0; while (j + 2 < pilots.length && pilots[j + 1] <= l) j++;
    const w = (l - pilots[j]) / (pilots[j + 1] - pilots[j]), Nr = raw[0].re.length, re = new Float64Array(Nr), im = new Float64Array(Nr);
    for (let p = 0; p < Nr; p++) { re[p] = (1 - w) * raw[j].re[p] + w * raw[j + 1].re[p]; im[p] = (1 - w) * raw[j].im[p] + w * raw[j + 1].im[p]; }
    return { re, im };
}
const lerpMid = (a, b) => ({ re: a.re.map((x, i) => 0.5 * (x + b.re[i])), im: a.im.map((x, i) => 0.5 * (x + b.im[i])) });

/** 全域導頻殘差協方差：d_k = Y_k − Ĥ_k X_k，Ĥ_k = 兩個相鄰導頻原始估計的平均（leave-one-out 內插；排除第一與最後一個導頻）；Hk 給定（每個導頻一個向量）時用它當 Ĥ_k（genie） */
function residuals(sym, Y, raw, Hk) {
    const Nr = sym.Nr, vecs = [], S = new Float64Array(2 * Nr * Nr), n = sym.pilots.length - 2;
    for (let j = 1; j < sym.pilots.length - 1; j++) {
        const l = sym.pilots[j], Hh = Hk ? Hk[j] : lerpMid(raw[j - 1], raw[j + 1]), v = { re: new Float64Array(Nr), im: new Float64Array(Nr) };
        for (let p = 0; p < Nr; p++) { v.re[p] = Y[p].re[l] - (sym.Xre[l] * Hh.re[p] - sym.Xim[l] * Hh.im[p]); v.im[p] = Y[p].im[l] - (sym.Xre[l] * Hh.im[p] + sym.Xim[l] * Hh.re[p]); }
        vecs.push(v); O.addOuter(S, v, 1 / n, Nr);
    }
    return { S, vecs, n };
}
const traceOf = (S, Nr) => { let t = 0; for (let i = 0; i < Nr; i++) t += S[2 * (i * Nr + i)]; return t; };

module.exports = { DELAYS_US, TAP_DB, multipathChannel, makeSymbolMP, mpSystem, rawPilotH, interp, lerpMid, residuals, traceOf };
