/* core.js -- DOM-free physics core of the Doppler / array-processing emulator.
 *
 * Loads in the browser (<script src="core.js">, exposes global `Core`) and in Node (require('./core.js')).
 * No external dependencies. UI, drawing and DOM access live in index.html.
 *
 * Contents: configuration (CONFIG), the single seedable random source (rng), complex / matrix helpers,
 * array manifold + covariance estimation + FOURIER / MMSE / SMI / DL / BEAMSPACE weights, SINR / EVM / SER, ICI formulas.
 */
(function (root, factory) {
    if (typeof module === 'object' && module.exports) module.exports = factory();
    else root.Core = factory();
}(typeof self !== 'undefined' ? self : this, function () {
    'use strict';

    /* ---------------------------------------------------------------- configuration */
    const CONFIG = {
        fc: 5e9,              // carrier frequency (Hz)
        c: 3e8,               // speed of light (m/s)
        scs: 15e3,            // OFDM subcarrier spacing Delta f (Hz)
        N: 8,                 // default number of array elements
        d_lambda: 0.5,        // element spacing in wavelengths
        cpRatio: 0.07,        // cyclic-prefix ratio (T_snap = (1 + cpRatio) / scs) -- used by the unified model
        R_RNG: 30,            // effective range (m): angular rate = v*sin(theta)/R
        SIGMA_ANG_DEG: 10,    // angular spread of the diffuse (NLoS) component (deg)
        M_SCAT: 8,            // scatterers per source per snapshot (legacy model)
        LAMBDA_Q: 30, REL_Q: 10.0,   // quiescent-preserving loading of the tapered (GSC) adaptive path
        REFRESH: 0.2          // fraction of the snapshot window replaced per update (legacy model)
    };

    /* ---------------------------------------------------------------- the single random source */
    // mulberry32: small, fast, 32-bit-state PRNG. All randomness in the core (and the UI scatter plot) goes through rng().
    function mulberry32(a) {
        return function () {
            a |= 0; a = a + 0x6D2B79F5 | 0;
            let t = Math.imul(a ^ a >>> 15, 1 | a);
            t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
            return ((t ^ t >>> 14) >>> 0) / 4294967296;
        };
    }
    let _seed = 0, _gen = null;
    function randomSeed() {
        if (typeof crypto !== 'undefined' && crypto.getRandomValues) { const b = new Uint32Array(1); crypto.getRandomValues(b); return b[0]; }
        return (Math.random() * 4294967296) >>> 0;
    }
    function setSeed(s) { _seed = s >>> 0; _gen = mulberry32(_seed); }
    function getSeed() { return _seed; }
    function rng() { return _gen(); }
    setSeed(randomSeed());      // default: random seed (tests and the UI may override with setSeed)

        /**
         * Complex math library
         */
        class Cplx {
            constructor(r=0, i=0) { this.r = r; this.i = i; }
            static add(a, b) { return new Cplx(a.r + b.r, a.i + b.i); }
            static sub(a, b) { return new Cplx(a.r - b.r, a.i - b.i); }
            static mul(a, b) { return new Cplx(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r); }
            static div(a, b) { const den = b.r*b.r + b.i*b.i; return new Cplx((a.r*b.r + a.i*b.i)/den, (a.i*b.r - a.r*b.i)/den); }
            static conj(a) { return new Cplx(a.r, -a.i); }
            mag2() { return this.r*this.r + this.i*this.i; }
            mag() { return Math.sqrt(this.mag2()); }
        }

        // Gauss-Jordan complex matrix inversion with partial pivoting. Returns null if singular.
        function invertMatrix(M, force = false) {
            const n = M.length;
            const A = M.map(row => row.map(c => new Cplx(c.r, c.i)));
            const I = Array(n).fill(0).map((_, i) => Array(n).fill(0).map((_, j) => new Cplx(i === j ? 1 : 0, 0)));
            for (let i = 0; i < n; i++) {
                let maxRow = i, maxVal = A[i][i].mag2();
                for (let k = i + 1; k < n; k++) if (A[k][i].mag2() > maxVal) { maxVal = A[k][i].mag2(); maxRow = k; }
                if (maxRow !== i) {
                    const tA = A[i]; A[i] = A[maxRow]; A[maxRow] = tA;
                    const tI = I[i]; I[i] = I[maxRow]; I[maxRow] = tI;
                }
                let diag = A[i][i];
                if (diag.mag2() < 1e-12) {                  // singular (e.g. L < N without regularisation)
                    if (!force) return null;
                    diag = A[i][i] = new Cplx(1e-6, 0);     // force: carry on with a tiny pivot -> meaningless, huge weights (collapse)
                }
                for (let j = 0; j < n; j++) { A[i][j] = Cplx.div(A[i][j], diag); I[i][j] = Cplx.div(I[i][j], diag); }
                for (let k = 0; k < n; k++) {
                    if (k === i) continue;
                    const f = A[k][i];
                    for (let j = 0; j < n; j++) {
                        A[k][j] = Cplx.sub(A[k][j], Cplx.mul(f, A[i][j]));
                        I[k][j] = Cplx.sub(I[k][j], Cplx.mul(f, I[i][j]));
                    }
                }
            }
            return I;
        }
        function matMulVec(M, v) { return M.map(row => row.reduce((acc, val, i) => Cplx.add(acc, Cplx.mul(val, v[i])), new Cplx(0,0))); }
        function vecDot(v1, v2) { return v1.reduce((acc, val, i) => Cplx.add(acc, Cplx.mul(Cplx.conj(val), v2[i])), new Cplx(0,0)); } // v1^H v2
        function quadForm(w, R) { // w^H R w
            let s = 0;
            for (let m = 0; m < w.length; m++) for (let n = 0; n < w.length; n++)
                s += Cplx.mul(Cplx.mul(Cplx.conj(w[m]), R[m][n]), w[n]).r;
            return s;
        }

        // Eigenvalues of a real symmetric matrix (cyclic Jacobi). A is modified in place.
        function eigvalsSym(A) {
            const n = A.length;
            let scale = 0;
            for (let i = 0; i < n; i++) scale += A[i][i] * A[i][i];
            scale = scale || 1;
            for (let sweep = 0; sweep < 40; sweep++) {
                let off = 0;
                for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += A[p][q] * A[p][q];
                if (off < 1e-26 * scale) break;
                for (let p = 0; p < n - 1; p++) for (let q = p + 1; q < n; q++) {
                    if (Math.abs(A[p][q]) < 1e-300) continue;
                    const theta = (A[q][q] - A[p][p]) / (2 * A[p][q]);
                    const t = (theta >= 0 ? 1 : -1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
                    const c = 1 / Math.sqrt(t * t + 1), s = t * c;
                    for (let k = 0; k < n; k++) { const akp = A[k][p], akq = A[k][q]; A[k][p] = c * akp - s * akq; A[k][q] = s * akp + c * akq; }
                    for (let k = 0; k < n; k++) { const apk = A[p][k], aqk = A[q][k]; A[p][k] = c * apk - s * aqk; A[q][k] = s * apk + c * aqk; }
                }
            }
            return A.map((row, i) => row[i]);
        }
        // Hermitian H = A + jB  ->  real symmetric [[A,-B],[B,A]] (each eigenvalue appears twice)
        function hermitianEigvals(H) {
            const n = H.length;
            const M = Array.from({ length: 2 * n }, () => new Float64Array(2 * n));
            for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
                const a = H[i][j].r, b = H[i][j].i;
                M[i][j] = a; M[i + n][j + n] = a; M[i][j + n] = -b; M[i + n][j] = b;
            }
            return eigvalsSym(M).sort((x, y) => x - y);
        }
        function sinc(x) { return x === 0 ? 1 : Math.sin(Math.PI * x) / (Math.PI * x); }
        function erfc(x) { // Abramowitz-Stegun 7.1.26
            const z = Math.abs(x), t = 1 / (1 + 0.3275911 * z);
            const y = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
            const r = y * Math.exp(-z * z);
            return x >= 0 ? r : 2 - r;
        }
        function qfunc(x) { return 0.5 * erfc(x / Math.SQRT2); }

        /**
         * Amplitude tapers (real, symmetric, max-normalised)
         */
        function hammingWindow(N) {
            if (N <= 1) return [1];
            return Array.from({ length: N }, (_, n) => 0.54 - 0.46 * Math.cos(2 * Math.PI * n / (N - 1)));
        }
        // Dolph-Chebyshev window (same construction as scipy.signal.windows.chebwin)
        function chebWindow(N, sllDb) {
            if (N <= 1) return [1];
            const order = N - 1, r = Math.pow(10, sllDb / 20);
            const beta = Math.cosh(Math.acosh(r) / order);
            const p = [];
            for (let k = 0; k < N; k++) {
                const x = beta * Math.cos(k * Math.PI / N);
                if (x > 1) p.push(Math.cosh(order * Math.acosh(x)));
                else if (x < -1) p.push((order % 2 ? -1 : 1) * Math.cosh(order * Math.acosh(-x)));
                else p.push(Math.cos(order * Math.acos(x)));
            }
            const even = N % 2 === 0, W = [];
            for (let m = 0; m < N; m++) {
                let s = 0;
                for (let k = 0; k < N; k++) s += p[k] * Math.cos(-2 * Math.PI * k * m / N + (even ? k * Math.PI / N : 0));
                W.push(s);
            }
            const w = [];
            if (even) { const n = N / 2 + 1; for (let i = n - 1; i >= 1; i--) w.push(W[i]); for (let i = 1; i < n; i++) w.push(W[i]); }
            else      { const n = (N + 1) / 2; for (let i = n - 1; i >= 1; i--) w.push(W[i]); for (let i = 0; i < n; i++) w.push(W[i]); }
            const mx = Math.max(...w);
            return w.map(x => x / mx);
        }

        // Square M-QAM with unit average power: levels (2k-(m-1))*a, a = sqrt(3/(2(M-1)))
        const MODS = {
            QPSK:  { name: 'QPSK',   M: 4,  bits: 2 },
            QAM16: { name: '16-QAM', M: 16, bits: 4 },
            QAM64: { name: '64-QAM', M: 64, bits: 6 },
        };
        function modInfo(key) {
            const m = MODS[key] || MODS.QAM16, ax = Math.sqrt(m.M);
            const a = Math.sqrt(3 / (2 * (m.M - 1)));
            return { ...m, ax, a, levels: Array.from({ length: ax }, (_, k) => (2 * k - (ax - 1)) * a) };
        }

        /**
         * System parameters and signal-processing engine
         */
        function createSys() {
        const Sys = {
            // UI parameters
            N: CONFIG.N, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, v: 0,
            calDeg: 0, latMs: 0, kDb: 20, sll: 35,
            algo: 'SMI', gammaDL: 0.01, taper: 'NONE', mod: 'QAM16',
            // lifecycle
            isRunning: true, simTime: 0, forceOnce: false,
            // constants
            fc: CONFIG.fc, c: CONFIG.c, d_lambda: CONFIG.d_lambda, scs: CONFIG.scs,
            R_RNG: CONFIG.R_RNG,                      // effective range (m): angular rate = v*sin(theta)/R
            SIGMA_ANG: CONFIG.SIGMA_ANG_DEG * Math.PI / 180,  // angular spread of the diffuse (NLoS) component
            M_SCAT: CONFIG.M_SCAT,                      // scatterers per source per snapshot
            LAMBDA_Q: CONFIG.LAMBDA_Q, REL_Q: CONFIG.REL_Q,      // quiescent-preserving loading of the adaptive path (see computeMath)
            REFRESH: CONFIG.REFRESH,                   // fraction of the snapshot window replaced per update
            calZ: [], calVer: 0, dirty: true, lastCompute: 0,
            snaps: [], snapKey: '',

            // computed state
            R_hat: [], weights: [], pattern: [], tVec: [], applied: {},
            kappa: 1, status: 'OK', fallback: false, invOk: true, delta: 0, lambdaQ: 0, noisePow: 0,
            outGain: 1, R_raw: [], kappaRaw: 1, gammaUsed: 0, collapsed: false, bsBins: [], bsAngles: [], R_B: [], bsSingular: false,
            eig: [], diagR: [], S: 0, I: 0, Nn: 0, sinrDb: 0, isDb: 0, nuICI: 0, evm: 0, ser: 0,
            fm: 0, fd: 0, eps0: 0, epsM: 0, nIciDb: -Infinity,
            thTo: 0, thJo: 0, dTdeg: 0, dJdeg: 0, nullDb: 0,
            psll: -Infinity, bw3: 0, mainL: 0, mainR: 0, jamInSL: false, jamGainDb: 0,

            rollCal() {
                this.calZ = Array.from({ length: 16 }, () => this.randn());
                this.calZ[0] = 0; // element 0 is the phase reference
                this.calVer++; this.dirty = true;
            },
            randn() { const u = 1 - rng(), v = rng(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); },
            qpsk(power) { const a = Math.sqrt(power / 2); return new Cplx(rng() < 0.5 ? -a : a, rng() < 0.5 ? -a : a); },
            taperWeights() {
                if (this.taper === 'HAMMING') return hammingWindow(this.N);
                if (this.taper === 'CHEBYSHEV') return chebWindow(this.N, this.sll);
                return Array(this.N).fill(1);
            },
            steer(rad) { // nominal a(theta)
                return Array.from({ length: this.N }, (_, n) => {
                    const p = -2 * Math.PI * n * this.d_lambda * Math.sin(rad);
                    return new Cplx(Math.cos(p), Math.sin(p));
                });
            },
            gamma() { // hardware phase mismatch Gamma = diag(e^{j phi_n})
                return Array.from({ length: this.N }, (_, n) => {
                    const ph = this.calDeg * Math.PI / 180 * this.calZ[n];
                    return new Cplx(Math.cos(ph), Math.sin(ph));
                });
            },
            diffuse(rad, gr, gi) { // diffuse (NLoS) component g(l): M scatterers, Gaussian angular spread
                const N = this.N, M = this.M_SCAT, bs = 1 / Math.sqrt(2 * M);
                gr.fill(0); gi.fill(0);
                for (let m = 0; m < M; m++) {
                    const th = rad + this.randn() * this.SIGMA_ANG;
                    const br = this.randn() * bs, bi = this.randn() * bs;
                    const k = -2 * Math.PI * this.d_lambda * Math.sin(th);
                    for (let n = 0; n < N; n++) {
                        const c = Math.cos(k * n), s = Math.sin(k * n);
                        gr[n] += br * c - bi * s; gi[n] += br * s + bi * c;
                    }
                }
            },
            trueCov(rad, pow, gam, cL2, cD2) { // expected covariance of one source
                const N = this.N, R = [];
                for (let m = 0; m < N; m++) {
                    R.push([]);
                    for (let n = 0; n < N; n++) {
                        const k = m - n;
                        const ph = -2 * Math.PI * this.d_lambda * k * Math.sin(rad);
                        const spread = Math.exp(-0.5 * Math.pow(2 * Math.PI * this.d_lambda * k * Math.cos(rad) * this.SIGMA_ANG, 2));
                        const amp = pow * (cL2 + cD2 * spread);
                        R[m].push(Cplx.mul(new Cplx(amp * Math.cos(ph), amp * Math.sin(ph)), Cplx.mul(gam[m], Cplx.conj(gam[n]))));
                    }
                }
                return R;
            },
            blockingMatrix(a) { // orthonormal basis B (as column arrays) of the complement of a: B^H a = 0, B^H B = I
                const N = a.length, cols = [];
                const nrm = Math.sqrt(a.reduce((s, c) => s + c.mag2(), 0));
                const basis = [a.map(c => new Cplx(c.r / nrm, c.i / nrm))];
                for (let k = 0; k < N && cols.length < N - 1; k++) {
                    let v = Array.from({ length: N }, (_, i) => new Cplx(i === k ? 1 : 0, 0));
                    for (const q of basis) {
                        const pr = vecDot(q, v);
                        v = v.map((c, i) => Cplx.sub(c, Cplx.mul(q[i], pr)));
                    }
                    const nv = Math.sqrt(v.reduce((s, c) => s + c.mag2(), 0));
                    if (nv < 1e-3) continue;
                    const u = v.map(c => new Cplx(c.r / nv, c.i / nv));
                    basis.push(u); cols.push(u);
                }
                return cols;
            },

            // w = R^-1 a / (a^H R^-1 a); returns null if the inverse does not exist (unless `force`) or the result is not finite
            mvdrWeights(R, a, force) {
                const inv = invertMatrix(R, force); if (!inv) return null;
                const num = matMulVec(inv, a), den = vecDot(a, num);
                if (!(den.mag2() > 1e-300)) return null;
                const w = num.map(c => Cplx.div(c, den));
                return w.every(c => Number.isFinite(c.r) && Number.isFinite(c.i)) ? w : null;
            },
            // Wiener / MMSE: w = R^-1 r_xd  (not normalised, so it differs from MVDR only by a scalar); null if not finite
            wienerWeights(R, rxd, force) {
                const inv = invertMatrix(R, force); if (!inv) return null;
                const w = matMulVec(inv, rxd);
                return w.every(c => Number.isFinite(c.r) && Number.isFinite(c.i)) ? w : null;
            },
            // Beamspace-MVDR: project onto K = min(3,N) orthonormal DFT beams, invert only the K x K matrix
            //   B = [d_k1 ... d_kK] (N x K), d_k[n] = e^{-j2 pi n k/N}/sqrt(N)  (sin(theta_k) = 2k/N for d = 0.5 lambda)
            //   r_B = B^H r,  R_B = B^H R B,  a_B = B^H a,  w_B = R_B^-1 a_B / (a_B^H R_B^-1 a_B),  w = B w_B
            // beams: the one nearest the target + the (K-1) with the highest estimated power (they cover strong interferers)
            beamspace(aNom, thTo) {
                const N = this.N, K = Math.min(3, N), R = this.R_raw;
                const dft = k => Array.from({ length: N }, (_, n) => {
                    const p = -2 * Math.PI * n * k / N, s = Math.sqrt(N);
                    return new Cplx(Math.cos(p) / s, Math.sin(p) / s);
                });
                const kT = ((Math.round(Math.sin(thTo) * N / 2) % N) + N) % N;
                const pw = [];
                for (let k = 0; k < N; k++) { const d = dft(k); pw.push({ k, p: vecDot(d, matMulVec(R, d)).r }); }
                const others = pw.filter(o => o.k !== kT).sort((a, b) => b.p - a.p).slice(0, K - 1).map(o => o.k);
                const bins = [kT, ...others].sort((a, b) => a - b);
                const B = bins.map(dft);
                const aB = B.map(b => vecDot(b, aNom));
                const RBc = B.map(b => matMulVec(R, b));
                const M = B.map(bi => RBc.map(rb => vecDot(bi, rb)));          // R_B (K x K)
                this.R_B = M; this.bsBins = bins;
                this.bsAngles = bins.map(k => { const kk = k >= N / 2 ? k - N : k; return Math.asin(Math.max(-1, Math.min(1, 2 * kk / N))) * 180 / Math.PI; });
                const ev = hermitianEigvals(M), lmax = ev[ev.length - 1], lmin = Math.max(ev[0], 0);
                this.kappa = (lmin <= 1e-12 * lmax) ? Infinity : lmax / lmin;   // condition number of the small matrix
                this.bsSingular = this.kappa === Infinity;                      // only when L < K
                const inv = invertMatrix(M, true);
                const wBnum = matMulVec(inv, aB), den = vecDot(aB, wBnum);
                if (!(den.mag2() > 1e-300)) return null;
                const wB = wBnum.map(c => Cplx.div(c, den));
                const w = Array.from({ length: N }, (_, n) => {
                    let acc = new Cplx(0, 0);
                    for (let j = 0; j < K; j++) acc = Cplx.add(acc, Cplx.mul(wB[j], B[j][n]));
                    return acc;
                });
                return w.every(c => Number.isFinite(c.r) && Number.isFinite(c.i)) ? w : null;
            },

            computeMath() {
                const N = this.N, L = this.L, D2R = Math.PI / 180;
                const sigPow = 1;
                const noisePow = Math.pow(10, -this.snr / 10);
                const jamPow = Math.pow(10, -this.sir / 10);
                const Klin = Math.pow(10, this.kDb / 10);
                const cL = Math.sqrt(Klin / (Klin + 1)), cD = Math.sqrt(1 / (Klin + 1));
                const thT = this.aoaT * D2R, thJ = this.aoaJ * D2R;

                // Channel aging: weights are estimated at t - tau, applied at t.
                const tau = this.latMs * 1e-3, vms = this.v / 3.6;
                const thTo = thT + vms * Math.sin(thT) * tau / this.R_RNG;
                const thJo = thJ - vms * Math.sin(thJ) * tau / this.R_RNG;
                this.thTo = thTo; this.thJo = thJo;
                this.dTdeg = (thTo - thT) / D2R; this.dJdeg = (thJo - thJ) / D2R;

                const gam = this.gamma();
                const withG = (a) => a.map((c, n) => Cplx.mul(c, gam[n]));

                // ---- snapshot window: r(l) = a1(l) s1(l) + a2(l) s2(l) + n(l)  (sources at the OLD angles)
                // statistics-changing parameters flush the window; otherwise a fraction REFRESH is replaced per update
                const key = [N, this.aoaT, this.aoaJ, this.snr, this.sir, this.calDeg, this.calVer, this.latMs, this.v, this.kDb].join('|');
                if (key !== this.snapKey) { this.snaps = []; this.snapKey = key; }
                const have = this.snaps.length;
                let add = have === 0 ? L : Math.max(0, L - have);
                if (have > 0) add = Math.max(add, Math.ceil(this.REFRESH * L));
                const nomT = this.steer(thTo), nomJ = this.steer(thJo);
                const sd = Math.sqrt(noisePow / 2);
                const gTr = new Float64Array(N), gTi = new Float64Array(N), gJr = new Float64Array(N), gJi = new Float64Array(N);
                for (let l = 0; l < add; l++) {
                    const s1 = this.qpsk(sigPow), s2 = this.qpsk(jamPow);
                    this.diffuse(thTo, gTr, gTi); this.diffuse(thJo, gJr, gJi);
                    const rr = new Float64Array(N), ri = new Float64Array(N);
                    for (let n = 0; n < N; n++) {
                        const t1r = cL * nomT[n].r + cD * gTr[n], t1i = cL * nomT[n].i + cD * gTi[n];
                        const t2r = cL * nomJ[n].r + cD * gJr[n], t2i = cL * nomJ[n].i + cD * gJi[n];
                        const a1r = gam[n].r * t1r - gam[n].i * t1i, a1i = gam[n].r * t1i + gam[n].i * t1r;
                        const a2r = gam[n].r * t2r - gam[n].i * t2i, a2i = gam[n].r * t2i + gam[n].i * t2r;
                        rr[n] = a1r * s1.r - a1i * s1.i + a2r * s2.r - a2i * s2.i + this.randn() * sd;
                        ri[n] = a1r * s1.i + a1i * s1.r + a2r * s2.i + a2i * s2.r + this.randn() * sd;
                    }
                    this.snaps.push({ rr, ri });
                }
                while (this.snaps.length > L) this.snaps.shift();

                const Rr = new Float64Array(N * N), Ri = new Float64Array(N * N);
                for (const { rr, ri } of this.snaps) {
                    for (let m = 0; m < N; m++) for (let n = 0; n < N; n++) { // r r^H
                        Rr[m * N + n] += rr[m] * rr[n] + ri[m] * ri[n];
                        Ri[m * N + n] += ri[m] * rr[n] - rr[m] * ri[n];
                    }
                }
                const Lw = this.snaps.length;
                this.R_hat = Array.from({ length: N }, (_, m) => Array.from({ length: N }, (_, n) =>
                    new Cplx(Rr[m * N + n] / Lw, Ri[m * N + n] / Lw)));
                this.noisePow = noisePow;
                const algo = this.algo;

                // ---- raw sample covariance R_hat (before any regularisation): eigenvalues, diagonal, condition number
                this.R_raw = this.R_hat;
                const evRaw = hermitianEigvals(this.R_raw);
                this.eig = []; for (let i = evRaw.length - 1; i >= 0; i -= 2) this.eig.push(Math.max(evRaw[i], 0));   // each eigenvalue appears twice
                this.diagR = this.R_raw.map((row, m) => row[m].r);
                const lmax = evRaw[evRaw.length - 1], lmin = Math.max(evRaw[0], 0);
                this.kappaRaw = (lmin <= 1e-12 * lmax) ? Infinity : lmax / lmin;

                // ---- DL-MVDR:  R_DL = R_hat + gamma I   (also the matrix shown in panel D)
                this.gammaUsed = algo === 'DL' ? this.gammaDL : 0;
                this.delta = this.gammaUsed;
                if (this.gammaUsed > 0) this.R_hat = this.R_raw.map((row, m) => row.map((c, n) => m === n ? new Cplx(c.r + this.gammaUsed, c.i) : c));
                this.kappa = algo === 'DL' ? (lmax + this.gammaUsed) / (lmin + this.gammaUsed) : this.kappaRaw;   // condition number of the matrix that is inverted

                // ---- quiescent weight  w_q = (t . a) / sum(t)   (so that w_q^H a = 1)
                const t = this.taperWeights(); this.tVec = t;
                const sumT = t.reduce((a, b) => a + b, 0);
                const wq = nomT.map((c, n) => new Cplx(c.r * t[n] / sumT, c.i * t[n] / sumT));

                // ---- beamformer
                this.fallback = false; this.collapsed = false; this.lambdaQ = 0; this.invOk = true;
                this.bsBins = []; this.bsAngles = []; this.R_B = []; this.bsSingular = false;
                const gscTaper = this.taper !== 'NONE' && (algo === 'SMI' || algo === 'DL');
                if (algo === 'FOURIER') {
                    this.weights = wq;
                } else if (algo === 'MMSE') {
                    // MMSE / Wiener:  w = R^-1 r_xd,  r_xd = E[r d*] = P_s a(theta_t)  (known pilot d, perfectly correlated with the target).
                    // Inverted like SMI (no safeguard): an exactly singular R_hat (L < N) collapses the beamformer.
                    const rxd = nomT.map(c => new Cplx(c.r * sigPow, c.i * sigPow));
                    const w = this.wienerWeights(this.R_hat, rxd, true);
                    if (w) this.weights = w; else { this.weights = wq; this.fallback = true; }
                    this.collapsed = this.kappaRaw === Infinity;
                } else if (algo === 'BEAMSPACE') {
                    const w = this.beamspace(nomT, thTo);
                    if (w) this.weights = w; else { this.weights = wq; this.fallback = true; }
                    this.collapsed = this.bsSingular;
                } else if (gscTaper) {
                    // tapered MVDR in generalized-sidelobe-canceller form (quiescent pattern = taper):
                    //   w = w_q - B w_a,  w_a = (B^H R B + lambda I)^-1 B^H R w_q      (R = R_hat, or R_DL for DL-MVDR)
                    const B = this.blockingMatrix(nomT), K = B.length;
                    const RB = B.map(b => matMulVec(this.R_hat, b));
                    const Rwq = matMulVec(this.R_hat, wq);
                    this.lambdaQ = Math.max(this.LAMBDA_Q * noisePow, this.REL_Q * vecDot(wq, Rwq).r);
                    const M = Array.from({ length: K }, (_, i) => Array.from({ length: K }, (_, j) => vecDot(B[i], RB[j])));
                    for (let i = 0; i < K; i++) M[i][i] = Cplx.add(M[i][i], new Cplx(this.lambdaQ, 0));
                    const cvec = B.map(b => vecDot(b, Rwq));
                    const Minv = invertMatrix(M); this.invOk = !!Minv;
                    if (Minv) {
                        const wa = matMulVec(Minv, cvec);
                        this.weights = wq.map((c, n) => {
                            let acc = new Cplx(c.r, c.i);
                            for (let j = 0; j < K; j++) acc = Cplx.sub(acc, Cplx.mul(B[j][n], wa[j]));
                            return acc;
                        });
                    } else { this.weights = wq; this.fallback = true; }
                } else {
                    // plain SMI-MVDR or DL-MVDR:  w = R^-1 a / (a^H R^-1 a).  SMI is inverted *without* any safeguard, so an
                    // exactly singular R_hat (L < N) yields numerically meaningless weights: the beamformer collapses.
                    const w = this.mvdrWeights(this.R_hat, nomT, algo === 'SMI');
                    if (w) this.weights = w; else { this.weights = wq; this.fallback = true; }
                    this.collapsed = (algo === 'SMI' && this.kappaRaw === Infinity);
                }
                if (algo === 'FOURIER') this.status = 'OK';                       // R_hat is not used
                else if (this.collapsed) this.status = 'SINGULAR';
                else if (L < N) this.status = algo === 'BEAMSPACE' ? 'L<N (beamspace)' : ((algo === 'DL' || gscTaper) ? 'L<N (regularized)' : 'L<N');
                else if (this.kappa > 1e6) this.status = 'ILL-COND';
                else this.status = 'OK';

                this.outGain = vecDot(this.weights, nomT).r;     // w^H a: 1 for the distortionless designs, < 1 for the Wiener (MMSE) filter

                // ---- SINR from TRUE covariances at the CURRENT angles (after latency)
                const Rt = this.trueCov(thT, sigPow, gam, cL * cL, cD * cD);
                const Rj = this.trueCov(thJ, jamPow, gam, cL * cL, cD * cD);
                const w = this.weights;
                this.S = quadForm(w, Rt);
                this.I = quadForm(w, Rj);
                this.Nn = w.reduce((a, c) => a + c.mag2(), 0) * noisePow;
                const sinrLin = this.S / (this.I + this.Nn);
                this.sinrDb = 10 * Math.log10(sinrLin);
                this.isDb = 10 * Math.log10(this.I / this.S + 1e-30);

                const gJ2 = vecDot(w, withG(this.steer(thJ))).mag2();
                const gT2 = vecDot(w, withG(this.steer(thT))).mag2();
                this.nullDb = 10 * Math.log10(gJ2 / gT2 + 1e-30);

                // ---- Doppler / ICI
                this.fm = vms * this.fc / this.c;
                this.fd = this.fm * Math.cos(thT);                // fixed heading (+x): fd = fm cos(theta)
                this.eps0 = this.fd / this.scs;
                this.epsM = this.fm / this.scs;
                const los = 1 - Math.pow(sinc(this.eps0), 2);
                let dif = 0; const A = 32;
                for (let k = 0; k < A; k++) dif += 1 - Math.pow(sinc(this.epsM * Math.cos(2 * Math.PI * (k + 0.5) / A)), 2);
                dif /= A;
                this.nuICI = (Klin / (Klin + 1)) * los + (1 / (Klin + 1)) * dif;   // N_ICI / S
                this.nIciDb = 10 * Math.log10(this.nuICI + 1e-30);
                this.evm = Math.sqrt(1 / sinrLin + this.nuICI);

                // ---- symbol error rate of the chosen square M-QAM under Gaussian error (per-axis sigma = EVM/sqrt2)
                const mi = modInfo(this.mod), sig = Math.max(this.evm, 1e-9) / Math.SQRT2;
                const pAxis = Math.min(1, 2 * (1 - 1 / mi.ax) * qfunc(mi.a / sig));
                this.ser = 1 - Math.pow(1 - pAxis, 2);

                // ---- pattern |w^H a(theta)|^2 (nominal array), 1-degree grid, plus lobe metrics
                const g = [];
                for (let deg = -90; deg <= 90; deg++) g.push(vecDot(w, this.steer(deg * D2R)).mag2());
                const gm = Math.max(...g, 1e-30);
                const p = g.map(x => 10 * Math.log10(x / gm + 1e-30));
                this.pattern = p;
                let ip = 0; for (let i = 1; i < p.length; i++) if (p[i] > p[ip]) ip = i;
                let il = ip; while (il > 0 && p[il - 1] <= p[il]) il--;
                let ir = ip; while (ir < p.length - 1 && p[ir + 1] <= p[ir]) ir++;
                this.mainL = il; this.mainR = ir;
                let psll = -Infinity;
                for (let i = 0; i < p.length; i++) if (i < il || i > ir) psll = Math.max(psll, p[i]);
                this.psll = psll;
                let xl = il, xr = ir;
                for (let i = ip; i > il; i--) if (p[i - 1] < -3) { xl = i - 1 + (-3 - p[i - 1]) / (p[i] - p[i - 1]); break; }
                for (let i = ip; i < ir; i++) if (p[i + 1] < -3) { xr = i + (p[i] + 3) / (p[i] - p[i + 1]); break; }
                this.bw3 = xr - xl;
                const jIdx = Math.min(180, Math.max(0, Math.round(this.aoaJ) + 90));
                this.jamInSL = jIdx < il || jIdx > ir;
                this.jamGainDb = p[jIdx];

                // snapshot of the parameters this result was computed with (diagnosis reads these, not the live sliders)
                this.applied = { algo: this.algo, taper: this.taper, gamma: this.gammaUsed, L, N, aoaJ: this.aoaJ, mod: this.mod };
            }
        };
        return Sys;
        }

    return { CONFIG, mulberry32, setSeed, getSeed, randomSeed, rng, Cplx, invertMatrix, matMulVec, vecDot, quadForm, eigvalsSym, hermitianEigvals, sinc, erfc, qfunc, hammingWindow, chebWindow, MODS, modInfo, createSys };
}));
