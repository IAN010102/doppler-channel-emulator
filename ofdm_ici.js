/* ofdm_ici.js -- OFDM 域的 ICI（載波間干擾，inter-carrier interference）與接收端波束成形（receive beamforming）核心模組。
 *
 * 獨立模組：不依賴 DOM、不依賴 core.js，瀏覽器與 Node 都能載入（UMD）。本輪不接到 UI。
 * 參考模型：Gopala & Slock, "Doppler Compensation and Beamforming for High Mobility OFDM Transmissions in Multipath"（EURECOM 2016）。
 * 公式與編號對照見 docs/ofdm_ici_model.md（編號 O1 ... O20 為本專案自訂；文獻式號除式 24 外，提示詞沒有提供，標為「待確認」）。
 * 文獻沒有寫明的設定一律是「假設」，記錄在 docs/ofdm_ici_assumptions.md。
 *
 * 約定（與模擬器其餘部分相同）：
 *   - 均勻線性陣列 ULA（uniform linear array），a_p(θ) = exp(−j2π p (d/λ) sinθ)，p = 0..Nr−1，θ 自陣列法線量起，d/λ = 0.5。
 *   - 複數向量以兩個 Float64Array（re, im）存放；Nr×Nr 的小矩陣以交錯的 Float64Array(2·Nr·Nr)（列優先）存放。
 *   - 子載波間隔 Δf = fs/N；ε_i = f_i/Δf（正規化都卜勒）；時間樣本 t 的通道  h(t) = Σ_i A_i a(θ_i) exp(j2π ε_i t/N)。
 *   - 傳送符元（transmit symbol）X_m，功率 σx²；OFDM 時域樣本 x(t) = Σ_m X_m exp(j2π m t/N)（未縮放的 IDFT）；
 *     接收端 FFT：Y_l = (1/N) Σ_{t=0}^{N−1} y(t) exp(−j2π l t/N)；雜訊在每個子載波的變異數為 σn²（時域樣本變異數 N·σn²）。
 */
(function (root, factory) {
    if (typeof module === 'object' && module.exports) module.exports = factory(); else root.OfdmIci = factory();
}(typeof self !== 'undefined' ? self : this, function () {
    'use strict';
    const PI = Math.PI, TWO_PI = 2 * Math.PI;

    /* ------------------------------------------------------------------ 亂數（mulberry32，本模組自帶一份，不依賴 core.js） */
    function mulberry32(a) {
        return function () {
            a |= 0; a = a + 0x6D2B79F5 | 0;
            let t = Math.imul(a ^ a >>> 15, 1 | a);
            t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
            return ((t ^ t >>> 14) >>> 0) / 4294967296;
        };
    }
    function makeRng(seed) {
        const u = mulberry32(seed >>> 0);
        return { u, n() { return Math.sqrt(-2 * Math.log(1 - u())) * Math.cos(TWO_PI * u()); } };   // n(): 標準常態
    }

    /* ------------------------------------------------------------------ 基本函式 */
    // O1  Q(x) = (1/N) Σ_{k=0}^{N−1} exp(j2π x k / N) = exp(jπ x (N−1)/N) · sin(πx) / (N sin(πx/N))   （封閉式；Q 對 x 的週期為 N）
    function Qc(x, N, out) {
        const y = x - N * Math.round(x / N);                 // 化到 (−N/2, N/2]
        if (Math.abs(y) < 1e-12) { out[0] = 1; out[1] = 0; return out; }
        const s = Math.sin(PI * y) / (N * Math.sin(PI * y / N)), ph = PI * y * (N - 1) / N;
        out[0] = s * Math.cos(ph); out[1] = s * Math.sin(ph); return out;
    }
    function QbruteForce(x, N) { let re = 0, im = 0; for (let k = 0; k < N; k++) { const ph = TWO_PI * x * k / N; re += Math.cos(ph); im += Math.sin(ph); } return [re / N, im / N]; }

    function steer(Nr, thetaRad, dLambda = 0.5) {
        const re = new Float64Array(Nr), im = new Float64Array(Nr);
        for (let p = 0; p < Nr; p++) { const ph = -TWO_PI * p * dLambda * Math.sin(thetaRad); re[p] = Math.cos(ph); im[p] = Math.sin(ph); }
        return { re, im };
    }

    // 未縮放的 radix-2 FFT（就地）；inverse = true 用 exp(+j…)。N 必須是 2 的冪。
    function fft(re, im, inverse) {
        const n = re.length;
        for (let i = 1, j = 0; i < n; i++) {
            let bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit;
            if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
        }
        for (let len = 2; len <= n; len <<= 1) {
            const ang = (inverse ? TWO_PI : -TWO_PI) / len, wr = Math.cos(ang), wi = Math.sin(ang), half = len >> 1;
            for (let i = 0; i < n; i += len) {
                let cr = 1, ci = 0;
                for (let k = 0; k < half; k++) {
                    const a = i + k, b = a + half, tr = re[b] * cr - im[b] * ci, ti = re[b] * ci + im[b] * cr;
                    re[b] = re[a] - tr; im[b] = im[a] - ti; re[a] += tr; im[a] += ti;
                    const nr = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = nr;
                }
            }
        }
    }

    /* ------------------------------------------------------------------ 小型複數線性代數（Nr 很小） */
    // 解 A x = b（A：Nr×Nr 複數矩陣，列優先交錯；b：{re, im}）；高斯消去、部分選主元；奇異時回傳 null
    function csolve(A, b, n) {
        const M = new Float64Array(2 * n * (n + 1));           // 增廣矩陣 [A | b]，每列 n+1 個複數
        const W = n + 1, at = (i, j) => 2 * (i * W + j);
        for (let i = 0; i < n; i++) { for (let j = 0; j < n; j++) { M[at(i, j)] = A[2 * (i * n + j)]; M[at(i, j) + 1] = A[2 * (i * n + j) + 1]; } M[at(i, n)] = b.re[i]; M[at(i, n) + 1] = b.im[i]; }
        let scale = 0; for (let i = 0; i < n; i++) scale = Math.max(scale, Math.hypot(M[at(i, i)], M[at(i, i) + 1]));
        for (let c = 0; c < n; c++) {
            let piv = c, best = -1;
            for (let r = c; r < n; r++) { const m = Math.hypot(M[at(r, c)], M[at(r, c) + 1]); if (m > best) { best = m; piv = r; } }
            if (!(best > 1e-14 * (scale || 1))) return null;
            if (piv !== c) for (let j = 0; j <= n; j++) { for (let q = 0; q < 2; q++) { const t = M[at(c, j) + q]; M[at(c, j) + q] = M[at(piv, j) + q]; M[at(piv, j) + q] = t; } }
            const pr = M[at(c, c)], pi = M[at(c, c) + 1], d = pr * pr + pi * pi;
            for (let r = 0; r < n; r++) {
                if (r === c) continue;
                const fr = M[at(r, c)], fi = M[at(r, c) + 1]; if (fr === 0 && fi === 0) continue;
                const qr = (fr * pr + fi * pi) / d, qi = (fi * pr - fr * pi) / d;          // f / pivot
                for (let j = c; j <= n; j++) { const xr = M[at(c, j)], xi = M[at(c, j) + 1]; M[at(r, j)] -= qr * xr - qi * xi; M[at(r, j) + 1] -= qr * xi + qi * xr; }
            }
        }
        const re = new Float64Array(n), im = new Float64Array(n);
        for (let i = 0; i < n; i++) { const pr = M[at(i, i)], pi = M[at(i, i) + 1], d = pr * pr + pi * pi, br = M[at(i, n)], bi = M[at(i, n) + 1]; re[i] = (br * pr + bi * pi) / d; im[i] = (bi * pr - br * pi) / d; }
        return { re, im };
    }
    // Hermitian 矩陣的偽逆解 A⁺ b（只用高於 epsRank·λmax 的特徵值；以實對稱嵌入 + Jacobi 做特徵分解）
    function pinvSolveHermitian(A, b, n, epsRank = 1e-10) {
        const m = 2 * n, S = Array.from({ length: m }, () => new Float64Array(m)), V = Array.from({ length: m }, (_, i) => { const r = new Float64Array(m); r[i] = 1; return r; });
        for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { const a = A[2 * (i * n + j)], c = A[2 * (i * n + j) + 1]; S[i][j] = a; S[i + n][j + n] = a; S[i][j + n] = -c; S[i + n][j] = c; }
        let scale = 0; for (let i = 0; i < m; i++) scale += S[i][i] * S[i][i]; scale = scale || 1;
        for (let sweep = 0; sweep < 60; sweep++) {
            let off = 0; for (let p = 0; p < m; p++) for (let q = p + 1; q < m; q++) off += S[p][q] * S[p][q];
            if (off < 1e-30 * scale) break;
            for (let p = 0; p < m - 1; p++) for (let q = p + 1; q < m; q++) {
                if (Math.abs(S[p][q]) < 1e-300) continue;
                const th = (S[q][q] - S[p][p]) / (2 * S[p][q]), t = (th >= 0 ? 1 : -1) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
                for (let k = 0; k < m; k++) { const x = S[k][p], y = S[k][q]; S[k][p] = c * x - s * y; S[k][q] = s * x + c * y; }
                for (let k = 0; k < m; k++) { const x = S[p][k], y = S[q][k]; S[p][k] = c * x - s * y; S[q][k] = s * x + c * y; }
                for (let k = 0; k < m; k++) { const x = V[k][p], y = V[k][q]; V[k][p] = c * x - s * y; V[k][q] = s * x + c * y; }
            }
        }
        let lmax = 0; for (let i = 0; i < m; i++) lmax = Math.max(lmax, S[i][i]);
        const bb = new Float64Array(m); for (let i = 0; i < n; i++) { bb[i] = b.re[i]; bb[i + n] = b.im[i]; }
        const x = new Float64Array(m);
        for (let k = 0; k < m; k++) { if (!(S[k][k] > epsRank * lmax)) continue; let dot = 0; for (let i = 0; i < m; i++) dot += V[i][k] * bb[i]; dot /= S[k][k]; for (let i = 0; i < m; i++) x[i] += dot * V[i][k]; }
        return { re: x.slice(0, n), im: x.slice(n) };
    }
    const cdot = (u, v) => { let re = 0, im = 0; for (let i = 0; i < u.re.length; i++) { re += u.re[i] * v.re[i] + u.im[i] * v.im[i]; im += u.re[i] * v.im[i] - u.im[i] * v.re[i]; } return [re, im]; };   // uᴴ v
    function matVec(A, v, n) { const re = new Float64Array(n), im = new Float64Array(n); for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { const ar = A[2 * (i * n + j)], ai = A[2 * (i * n + j) + 1]; re[i] += ar * v.re[j] - ai * v.im[j]; im[i] += ar * v.im[j] + ai * v.re[j]; } return { re, im }; }
    function addOuter(B, v, w, Nr) { for (let i = 0; i < Nr; i++) for (let j = 0; j < Nr; j++) { const k = 2 * (i * Nr + j); B[k] += w * (v.re[i] * v.re[j] + v.im[i] * v.im[j]); B[k + 1] += w * (v.im[i] * v.re[j] - v.re[i] * v.im[j]); } }   // B += w v vᴴ

    /* ------------------------------------------------------------------ 通道 */
    /** 建立通道。o = { N, Nr, df (Hz, 預設 15e3), dLambda (0.5), paths: [{ fHz | eps, thetaDeg, powerDb, phase }], rng, normalize }
     *  路徑功率 powerDb 為相對值；預設正規化使 Σ|A_i|² = 1（SNR 的定義見 snrToNoise）；phase 未給時由 rng 均勻抽取（假設，見 docs/ofdm_ici_assumptions.md）。 */
    function makeChannel(o) {
        const N = o.N, Nr = o.Nr, df = o.df || 15e3, dl = o.dLambda === undefined ? 0.5 : o.dLambda;
        const P = o.paths.map(p => Math.pow(10, (p.powerDb || 0) / 10)), Ptot = P.reduce((a, b) => a + b, 0), norm = o.normalize === false ? 1 : 1 / Ptot;
        const paths = o.paths.map((p, i) => {
            const eps = p.eps !== undefined ? p.eps : p.fHz / df, amp = Math.sqrt(P[i] * norm), ph = p.phase !== undefined ? p.phase : TWO_PI * (o.rng ? o.rng.u() : 0);
            return { eps, fHz: eps * df, thetaDeg: p.thetaDeg, Are: amp * Math.cos(ph), Aim: amp * Math.sin(ph), a: steer(Nr, p.thetaDeg * PI / 180, dl) };
        });
        let tot = 0; for (const p of paths) tot += p.Are * p.Are + p.Aim * p.Aim;
        return { N, Nr, df, dLambda: dl, paths, totalPower: tot };
    }
    /** SNR（每根天線、每個子載波）= σx² Σ|A_i|² / σn²。回傳 σn²。（假設 A1，見假設文件） */
    function snrToNoise(ch, snrDb, sigmaX2 = 1) { return sigmaX2 * ch.totalPower / Math.pow(10, snrDb / 10); }

    // O2/O3  H_ICI(l,m) 只和 d = (m − l) mod N 有關（平坦通道、無延遲擴散）：向量 v_d = Σ_i A_i a(θ_i) Q(d + ε_i − ε)，d = 0..N−1；v_0 = H。
    //        epsDemod = ε：對頻（demodulation）把 ε_i 換成 ε_i − ε（O6）。
    function iciVectors(ch, epsDemod = 0) {
        const { N, Nr, paths } = ch, hr = new Float64Array(N * Nr), hi = new Float64Array(N * Nr), q = [0, 0];
        for (let d = 0; d < N; d++) for (const pa of paths) {
            Qc(d + pa.eps - epsDemod, N, q);
            const cr = pa.Are * q[0] - pa.Aim * q[1], ci = pa.Are * q[1] + pa.Aim * q[0], o = d * Nr;
            for (let p = 0; p < Nr; p++) { hr[o + p] += cr * pa.a.re[p] - ci * pa.a.im[p]; hi[o + p] += cr * pa.a.im[p] + ci * pa.a.re[p]; }
        }
        return { hr, hi, N, Nr };
    }
    const vecAt = (V, d) => ({ re: V.hr.slice(d * V.Nr, (d + 1) * V.Nr), im: V.hi.slice(d * V.Nr, (d + 1) * V.Nr) });
    /** O2  H = Σ_i A_i a(θ_i) Q(ε_i − ε) */
    function Hvec(ch, epsDemod = 0) {
        const { N, Nr, paths } = ch, re = new Float64Array(Nr), im = new Float64Array(Nr), q = [0, 0];
        for (const pa of paths) { Qc(pa.eps - epsDemod, N, q); const cr = pa.Are * q[0] - pa.Aim * q[1], ci = pa.Are * q[1] + pa.Aim * q[0]; for (let p = 0; p < Nr; p++) { re[p] += cr * pa.a.re[p] - ci * pa.a.im[p]; im[p] += cr * pa.a.im[p] + ci * pa.a.re[p]; } }
        return { re, im };
    }
    /** O4  B_l = Σ_{m∈occ, m≠l} H_ICI(l,m) H_ICI(l,m)ᴴ（Nr×Nr）。occ = null：所有子載波皆占用（B 與 l 無關，省略下標）。 */
    function iciCov(V, occ, l = 0) {
        const { N, Nr } = V, B = new Float64Array(2 * Nr * Nr), v = { re: new Float64Array(Nr), im: new Float64Array(Nr) };
        for (let m = 0; m < N; m++) {
            if (m === l || (occ && !occ[m])) continue;
            const d = ((m - l) % N + N) % N;
            for (let p = 0; p < Nr; p++) { v.re[p] = V.hr[d * Nr + p]; v.im[p] = V.hi[d * Nr + p]; }
            addOuter(B, v, 1, Nr);
        }
        return B;
    }
    /** O5  R = σx² B + σn² I */
    function rMatrix(B, Nr, sigmaX2, sigmaN2) { const R = new Float64Array(B.length); for (let i = 0; i < B.length; i++) R[i] = sigmaX2 * B[i]; for (let p = 0; p < Nr; p++) R[2 * (p * Nr + p)] += sigmaN2; return R; }
    /** O7  G_opt = R⁻¹ H，O8  SINR_opt = σx² Hᴴ R⁻¹ H */
    function optimal(H, R, Nr, sigmaX2) { const G = csolve(R, H, Nr), s = cdot(H, G); return { G, sinr: sigmaX2 * s[0] }; }
    /** 任意權重的 SINR = σx² |Gᴴ H|² / (Gᴴ R G)（R：真實的 ICI+雜訊協方差） */
    function sinrOfWeights(G, H, R, Nr, sigmaX2) { const s = cdot(G, H), RG = matVec(R, G, Nr), d = cdot(G, RG); return sigmaX2 * (s[0] * s[0] + s[1] * s[1]) / d[0]; }

    /** 一次算出給定 ε（對頻）下的 H、B、R（所有子載波占用或指定 occ 與子載波 l） */
    function system(ch, epsDemod, sigmaX2, sigmaN2, occ = null, l = 0) {
        const V = iciVectors(ch, epsDemod), H = vecAt(V, 0), B = iciCov(V, occ, l);
        return { H, B, R: rMatrix(B, ch.Nr, sigmaX2, sigmaN2), V };
    }
    /** 各種接收機的 SINR 的基本比較用（真實 H、R；ε = 對頻值） */
    function sinrOptimalAt(ch, epsDemod, sigmaX2, sigmaN2, occ = null, l = 0) { const s = system(ch, epsDemod, sigmaX2, sigmaN2, occ, l); return optimal(s.H, s.R, ch.Nr, sigmaX2).sinr; }
    function sinrMrc(ch, sigmaX2, sigmaN2) { const s = system(ch, 0, sigmaX2, sigmaN2); return sinrOfWeights(s.H, s.H, s.R, ch.Nr, sigmaX2); }                       // G = H（不考慮 ICI、不對頻）
    function sinrSingleAntenna(ch, sigmaX2, sigmaN2, p = 0) { const s = system(ch, 0, sigmaX2, sigmaN2), k = 2 * (p * ch.Nr + p); return sigmaX2 * (s.H.re[p] * s.H.re[p] + s.H.im[p] * s.H.im[p]) / s.R[k]; }

    /* ------------------------------------------------------------------ 對頻的選法 */
    /** (a) 不對頻：0。 (b) LoS：第 losIndex 條路徑的 ε。 */
    const epsNone = () => 0;
    const epsLos = (ch, losIndex = 0) => ch.paths[losIndex].eps;
    /** (c) 循環前綴相關法（CP-based），無雜訊的期望值版本：
     *      C = Σ_{t=−Ncp}^{−1} Σ_p E[y_p(t+N) y_p(t)*] = N σx² Σ_t Σ_p h_p(t+N) h_p(t)*  （CP 讓 x(t+N) = x(t)），ε_CP = angle(C)/(2π)。
     *      這等於選 ε 使對頻後的 Σ y'(t+N) y'(t)* 為實數（O9）。適用範圍：|ε| < 0.5（相位 2π 模糊）；多路徑時得到各路徑的相位加權平均。 */
    function epsCpAnalytic(ch, Ncp) {
        const { N, Nr, paths } = ch; let cr = 0, ci = 0;
        for (let t = -Ncp; t < 0; t++) {
            const h0r = new Float64Array(Nr), h0i = new Float64Array(Nr), h1r = new Float64Array(Nr), h1i = new Float64Array(Nr);
            for (const pa of paths) {
                const p0 = TWO_PI * pa.eps * t / N, p1 = TWO_PI * pa.eps * (t + N) / N, c0 = Math.cos(p0), s0 = Math.sin(p0), c1 = Math.cos(p1), s1 = Math.sin(p1);
                const ar0 = pa.Are * c0 - pa.Aim * s0, ai0 = pa.Are * s0 + pa.Aim * c0, ar1 = pa.Are * c1 - pa.Aim * s1, ai1 = pa.Are * s1 + pa.Aim * c1;
                for (let p = 0; p < Nr; p++) {
                    h0r[p] += ar0 * pa.a.re[p] - ai0 * pa.a.im[p]; h0i[p] += ar0 * pa.a.im[p] + ai0 * pa.a.re[p];
                    h1r[p] += ar1 * pa.a.re[p] - ai1 * pa.a.im[p]; h1i[p] += ar1 * pa.a.im[p] + ai1 * pa.a.re[p];
                }
            }
            for (let p = 0; p < Nr; p++) { cr += h1r[p] * h0r[p] + h1i[p] * h0i[p]; ci += h1i[p] * h0r[p] - h1r[p] * h0i[p]; }     // h(t+N) h(t)*
        }
        return Math.atan2(ci, cr) / TWO_PI;
    }
    /** (c) 有雜訊時：由接收到的時域樣本（含 CP）估計 ε̂ = angle(Σ_{k<Ncp} Σ_p r_p[k+N] r_p[k]*) / (2π)；r_p 長度為 N + Ncp。 */
    function epsCpEstimate(sym) {
        const { N, Ncp, Nr } = sym; let cr = 0, ci = 0;
        for (let p = 0; p < Nr; p++) { const R = sym.rRe[p], I = sym.rIm[p]; for (let k = 0; k < Ncp; k++) { cr += R[k + N] * R[k] + I[k + N] * I[k]; ci += I[k + N] * R[k] - R[k + N] * I[k]; } }
        return Math.atan2(ci, cr) / TWO_PI;
    }
    /** ε 的格點掃描：回傳 [{eps, sinr[], power}]，sinr[] 對應 sigmaN2List。B 與 σn² 無關，所以每個 ε 只算一次 B。 */
    function scanEps(ch, epsList, sigmaX2, sigmaN2List) {
        return epsList.map(e => {
            const V = iciVectors(ch, e), H = vecAt(V, 0), B = iciCov(V, null, 0), power = H.re.reduce((a, x, i) => a + x * x + H.im[i] * H.im[i], 0);
            return { eps: e, power, sinr: sigmaN2List.map(sn => optimal(H, rMatrix(B, ch.Nr, sigmaX2, sn), ch.Nr, sigmaX2).sinr) };
        });
    }
    /** (d) genie：窮舉 ε 取使 SINR_opt 最大者（粗格點 + 黃金分割細化）；objective = 'sinr'（最佳 SINR）或 'power'（訊號功率 ‖H(ε)‖²，文獻圖 2 的比較項）。 */
    function epsGenie(ch, sigmaX2, sigmaN2, o = {}) {
        const objective = o.objective || 'sinr', lo = (o.lo !== undefined ? o.lo : Math.min(...ch.paths.map(p => p.eps)) - 0.05), hi = (o.hi !== undefined ? o.hi : Math.max(...ch.paths.map(p => p.eps)) + 0.05), steps = o.steps || 60;
        const f = e => { const r = scanEps(ch, [e], sigmaX2, [sigmaN2])[0]; return objective === 'power' ? r.power : r.sinr[0]; };
        let bestE = lo, bestV = -Infinity;
        for (let k = 0; k <= steps; k++) { const e = lo + (hi - lo) * k / steps, v = f(e); if (v > bestV) { bestV = v; bestE = e; } }
        let a = Math.max(lo, bestE - (hi - lo) / steps), b = Math.min(hi, bestE + (hi - lo) / steps); const gr = (Math.sqrt(5) - 1) / 2;
        let c = b - gr * (b - a), d = a + gr * (b - a), fc = f(c), fd = f(d);
        for (let it = 0; it < 30; it++) { if (fc > fd) { b = d; d = c; fd = fc; c = b - gr * (b - a); fc = f(c); } else { a = c; c = d; fc = fd; d = a + gr * (b - a); fd = f(d); } }
        const e = (a + b) / 2, v = f(e);
        return v >= bestV ? { eps: e, value: v } : { eps: bestE, value: bestV };
    }

    /* ------------------------------------------------------------------ 線性近似版（通道在符元內線性變化） */
    // O10  h_n = h0 + (n − (N−1)/2) h1：h0 = 符元中心的通道，h1 = 該點的斜率（真實 h1，genie 版：文獻中 h1 的取得方式引用自另一篇，未實作）。
    // O11  Ξ_{k,l} = (1/N) Σ_n (n − (N−1)/2) exp(j2π(k−l)n/N)；Σ_k |Ξ_{k,l}|² = (N²−1)/12（Parseval，與 l 無關，Ξ_{l,l} = 0）。
    function linearModel(ch, epsDemod = 0) {
        const { N, Nr, paths } = ch, c = (N - 1) / 2, h0 = { re: new Float64Array(Nr), im: new Float64Array(Nr) }, h1 = { re: new Float64Array(Nr), im: new Float64Array(Nr) };
        for (const pa of paths) {
            const e = pa.eps - epsDemod, ph = TWO_PI * e * c / N, cr = Math.cos(ph), ci = Math.sin(ph);
            const g0r = pa.Are * cr - pa.Aim * ci, g0i = pa.Are * ci + pa.Aim * cr, k = TWO_PI * e / N;           // g1 = j k g0
            const g1r = -k * g0i, g1i = k * g0r;
            for (let p = 0; p < Nr; p++) {
                h0.re[p] += g0r * pa.a.re[p] - g0i * pa.a.im[p]; h0.im[p] += g0r * pa.a.im[p] + g0i * pa.a.re[p];
                h1.re[p] += g1r * pa.a.re[p] - g1i * pa.a.im[p]; h1.im[p] += g1r * pa.a.im[p] + g1i * pa.a.re[p];
            }
        }
        return { h0, h1, sumXi2: (N * N - 1) / 12 };
    }
    /** O12  G_lin = (σx² S h1 h1ᴴ + σn² I)⁻¹ h0，S = Σ_k |Ξ_{k,l}|² */
    function linearWeights(ch, epsDemod, sigmaX2, sigmaN2) {
        const { Nr } = ch, m = linearModel(ch, epsDemod), R = new Float64Array(2 * Nr * Nr);
        addOuter(R, m.h1, sigmaX2 * m.sumXi2, Nr); for (let p = 0; p < Nr; p++) R[2 * (p * Nr + p)] += sigmaN2;
        return csolve(R, m.h0, Nr);
    }

    /* ------------------------------------------------------------------ 時域模擬與實務接收機 */
    /** 產生一個 OFDM 符元（含 CP）經過時變通道與雜訊後的接收樣本。
     *  o = { rng, sigmaX2, sigmaN2, Ncp, occ (Uint8Array | null), pilotSpacing (12), pilotOffset (0) }
     *  傳送符元：QPSK、功率 σx²（假設）；導頻子載波的符元也是已知的 QPSK。 */
    function makeSymbol(ch, o) {
        const { N, Nr, paths } = ch, rng = o.rng, sx2 = o.sigmaX2 === undefined ? 1 : o.sigmaX2, sn2 = o.sigmaN2, Ncp = o.Ncp, occ = o.occ || null, ps = o.pilotSpacing === undefined ? 12 : o.pilotSpacing, po = o.pilotOffset || 0;
        const Xre = new Float64Array(N), Xim = new Float64Array(N), a = Math.sqrt(sx2 / 2), pilots = [];
        for (let m = 0; m < N; m++) { if (occ && !occ[m]) continue; Xre[m] = rng.u() < 0.5 ? -a : a; Xim[m] = rng.u() < 0.5 ? -a : a; if (ps > 0 && (m - po) % ps === 0 && m >= po) pilots.push(m); }
        const xr = Float64Array.from(Xre), xi = Float64Array.from(Xim); fft(xr, xi, true);                          // x(n) = Σ_m X_m exp(j2π m n/N)
        const T = N + Ncp, rRe = [], rIm = [], sw = Math.sqrt(N * sn2 / 2);
        for (let p = 0; p < Nr; p++) { rRe.push(new Float64Array(T)); rIm.push(new Float64Array(T)); }
        for (let k = 0; k < T; k++) {
            const t = k - Ncp, n = ((t % N) + N) % N;                                                              // CP：x(t) = x(t + N)
            const hr = new Float64Array(Nr), hi = new Float64Array(Nr);
            for (const pa of paths) {
                const ph = TWO_PI * pa.eps * t / N, c = Math.cos(ph), s = Math.sin(ph), ar = pa.Are * c - pa.Aim * s, ai = pa.Are * s + pa.Aim * c;
                for (let p = 0; p < Nr; p++) { hr[p] += ar * pa.a.re[p] - ai * pa.a.im[p]; hi[p] += ar * pa.a.im[p] + ai * pa.a.re[p]; }
            }
            for (let p = 0; p < Nr; p++) { rRe[p][k] = hr[p] * xr[n] - hi[p] * xi[n] + sw * rng.n(); rIm[p][k] = hr[p] * xi[n] + hi[p] * xr[n] + sw * rng.n(); }
        }
        return { N, Ncp, Nr, Xre, Xim, rRe, rIm, pilots, sigmaX2: sx2, sigmaN2: sn2 };
    }
    /** 對頻 ε（把時間樣本 t 乘上 exp(−j2π ε t/N)，t = 0..N−1 為 FFT 窗）後做 FFT：Y_p[l]（O13：Y_l = (1/N) Σ_t y(t) e^{−j2π l t/N}） */
    function demodulate(sym, eps) {
        const { N, Ncp, Nr } = sym, Y = [];
        for (let p = 0; p < Nr; p++) {
            const re = new Float64Array(N), im = new Float64Array(N);
            for (let n = 0; n < N; n++) { const ph = -TWO_PI * eps * n / N, c = Math.cos(ph), s = Math.sin(ph), r = sym.rRe[p][Ncp + n], i = sym.rIm[p][Ncp + n]; re[n] = r * c - i * s; im[n] = r * s + i * c; }
            fft(re, im, false); for (let n = 0; n < N; n++) { re[n] /= N; im[n] /= N; }
            Y.push({ re, im });
        }
        return Y;
    }
    /** 實務接收機：Ĥ 由導頻估計（對各導頻平均，平坦通道），R_yy = (1/Ns) Σ_{l∈S} Y_l Y_lᴴ，G_est = R_yy⁻¹ Ĥ（O14–O16）。
     *  o = { Ns (預設 N), select: 'contiguous' | 'uniform', start (contiguous 的起點，預設 0), reg: 'none' | 'dl', gammaRel (線性倍數，γ = gammaRel·σn²) }
     *  reg = 'none'：Ns ≥ Nr 時直接求逆，否則偽逆；reg = 'dl'：對角加載（diagonal loading）。 */
    function practicalWeights(sym, Y, o = {}) {
        const { N, Nr } = sym, Ns = o.Ns || N, sel = o.select || 'contiguous', start = o.start || 0, reg = o.reg || 'none';
        const Hh = { re: new Float64Array(Nr), im: new Float64Array(Nr) };
        for (const l of sym.pilots) { const d = sym.Xre[l] * sym.Xre[l] + sym.Xim[l] * sym.Xim[l]; for (let p = 0; p < Nr; p++) { Hh.re[p] += (Y[p].re[l] * sym.Xre[l] + Y[p].im[l] * sym.Xim[l]) / d; Hh.im[p] += (Y[p].im[l] * sym.Xre[l] - Y[p].re[l] * sym.Xim[l]) / d; } }
        for (let p = 0; p < Nr; p++) { Hh.re[p] /= sym.pilots.length; Hh.im[p] /= sym.pilots.length; }
        const idx = []; if (sel === 'uniform') for (let k = 0; k < Ns; k++) idx.push(Math.round(k * N / Ns) % N); else for (let k = 0; k < Ns; k++) idx.push((start + k) % N);
        const Ryy = new Float64Array(2 * Nr * Nr), y = { re: new Float64Array(Nr), im: new Float64Array(Nr) };
        for (const l of idx) { for (let p = 0; p < Nr; p++) { y.re[p] = Y[p].re[l]; y.im[p] = Y[p].im[l]; } addOuter(Ryy, y, 1 / idx.length, Nr); }
        if (reg === 'dl') { const gam = (o.gammaRel === undefined ? 1 : o.gammaRel) * sym.sigmaN2; for (let p = 0; p < Nr; p++) Ryy[2 * (p * Nr + p)] += gam; }
        let G = (reg === 'none' && Ns < Nr) ? null : csolve(Ryy, Hh, Nr);
        if (!G) G = pinvSolveHermitian(Ryy, Hh, Nr);
        return { G, Hhat: Hh, Ryy, idx };
    }

    /* ------------------------------------------------------------------ 導頻殘差協方差估計器（pilot-residual covariance estimator） */
    // 只用導頻子載波的殘差估計「ICI 加雜訊」協方差，所以協方差裡沒有所需訊號（避免含訊號訓練的自我抵消）：
    //   d_k = Y_k − Ĥ X_k，k 為導頻子載波（可跨 m 個符元）；  R_res = (1/n) Σ d_k d_kᴴ   （n = 導頻總數）；  G = (R_res + γ I)⁻¹ Ĥ，γ = gammaRel · σn²。
    // 注意：Ĥ 由同一批導頻估計時，殘差偏小（自由度損失）；選項 correct 把 R_res 乘上 n/(n−1)；選項 Hforce 用真實 H（genie）取代 Ĥ。
    // shrink = 'lw'：Ledoit–Wolf 線性收縮（單位矩陣目標），公式見 ledoitWolfShrink。
    /** Ledoit–Wolf (2004), "A well-conditioned estimator for large-dimensional covariance matrices", J. Multivariate Analysis 88(2): 365–411（樣本協方差對 μI 的線性收縮，
     *  已知零平均所以不去平均）：μ = tr(S)/p；d² = ‖S − μI‖²/p；b̄² = (1/n²) Σ_k ‖d_k d_kᴴ − S‖²/p；b² = min(b̄², d²)；a² = d² − b²；S* = (b²/d²) μ I + (a²/d²) S，
     *  其中 ‖A‖² = tr(AAᴴ)（本函式的 /p 對應原文的縮放 Frobenius 範數）。本實作把原文的實數資料直接換成複數資料（外積用共軛轉置），公式不變。 */
    function ledoitWolfShrink(vecs, S, Nr) {
        const n = vecs.length, tr = (() => { let t = 0; for (let i = 0; i < Nr; i++) t += S[2 * (i * Nr + i)]; return t; })(), mu = tr / Nr;
        const fro2 = M => { let t = 0; for (let i = 0; i < M.length; i++) t += M[i] * M[i]; return t; };
        const D = Float64Array.from(S); for (let i = 0; i < Nr; i++) D[2 * (i * Nr + i)] -= mu;
        const d2 = fro2(D) / Nr; let b2bar = 0;
        for (const v of vecs) { const M = new Float64Array(2 * Nr * Nr); addOuter(M, v, 1, Nr); for (let i = 0; i < M.length; i++) M[i] -= S[i]; b2bar += fro2(M) / Nr; }
        b2bar /= n * n; const b2 = Math.min(b2bar, d2), a2 = d2 - b2, out = new Float64Array(S.length);
        for (let i = 0; i < S.length; i++) out[i] = (a2 / d2) * S[i]; for (let i = 0; i < Nr; i++) out[2 * (i * Nr + i)] += (b2 / d2) * mu;
        return { S: out, shrinkage: b2 / d2 };
    }
    /** items = [{ sym, Y }, …]（m 個符元，同一個通道，同一個對頻 ε̂ 下的 Y）；o = { gammaRel (線性倍數；null/undefined = 無加載), shrink: 'lw' | null, correct: bool, Hforce: {re, im} | null } */
    function pilotResidualWeights(items, o = {}) {
        const Nr = items[0].sym.Nr, sn2 = items[0].sym.sigmaN2, H = { re: new Float64Array(Nr), im: new Float64Array(Nr) }; let n = 0;
        for (const { sym, Y } of items) for (const l of sym.pilots) { const d = sym.Xre[l] * sym.Xre[l] + sym.Xim[l] * sym.Xim[l]; for (let p = 0; p < Nr; p++) { H.re[p] += (Y[p].re[l] * sym.Xre[l] + Y[p].im[l] * sym.Xim[l]) / d; H.im[p] += (Y[p].im[l] * sym.Xre[l] - Y[p].re[l] * sym.Xim[l]) / d; } n++; }
        for (let p = 0; p < Nr; p++) { H.re[p] /= n; H.im[p] /= n; }
        const Hu = o.Hforce || H, vecs = [], S = new Float64Array(2 * Nr * Nr);
        for (const { sym, Y } of items) for (const l of sym.pilots) {
            const v = { re: new Float64Array(Nr), im: new Float64Array(Nr) };
            for (let p = 0; p < Nr; p++) { v.re[p] = Y[p].re[l] - (sym.Xre[l] * Hu.re[p] - sym.Xim[l] * Hu.im[p]); v.im[p] = Y[p].im[l] - (sym.Xre[l] * Hu.im[p] + sym.Xim[l] * Hu.re[p]); }
            vecs.push(v); addOuter(S, v, 1 / n, Nr);
        }
        let R = S; if (o.correct) { R = S.map(x => x * n / (n - 1)); }
        let shrinkage = null; if (o.shrink === 'lw') { const r = ledoitWolfShrink(vecs, R, Nr); R = r.S; shrinkage = r.shrinkage; }
        if (o.gammaRel !== undefined && o.gammaRel !== null) { R = Float64Array.from(R); for (let p = 0; p < Nr; p++) R[2 * (p * Nr + p)] += o.gammaRel * sn2; }
        let G = csolve(R, Hu, Nr); if (!G) G = pinvSolveHermitian(R, Hu, Nr);
        return { G, Hhat: H, R, n, shrinkage };
    }

    /* ------------------------------------------------------------------ 其他工具 */
    /** 保護頻帶：回傳占用遮罩（Uint8Array）；nLow、nHigh = 兩端不使用的子載波數 */
    function maskGuard(N, nLow, nHigh) { const m = new Uint8Array(N).fill(1); for (let i = 0; i < nLow; i++) m[i] = 0; for (let i = 0; i < nHigh; i++) m[N - 1 - i] = 0; return m; }
    /** 文獻模擬設定（假設 A2：路徑 1 是 LoS；A3：隨機相位；到達角見 anglesFixed / anglesRandom） */
    const PAPER = { N: 1024, df: 15e3, Nr: 2, dopplerHz: [1080, -1080, 758, 220], relPowerDb: [0, 0, -11, -0.7], pilotSpacing: 12 };
    function paperChannel(angles, rng, o = {}) {
        const Nr = o.Nr || PAPER.Nr, N = o.N || PAPER.N, scale = o.dopplerScale === undefined ? 1 : o.dopplerScale;
        return makeChannel({ N, Nr, df: o.df || PAPER.df, rng, paths: PAPER.dopplerHz.map((f, i) => ({ fHz: f * scale, thetaDeg: angles[i], powerDb: PAPER.relPowerDb[i] })) });
    }
    /* ------------------------------------------------------------------ 角度域匹配濾波接收機（AD-MF，genie 基準；不做 DoA 估計、不做參數估計） */
    // 對每條路徑 p：用真實到達角做空間匹配濾波 u_p(t) = a(θ_p)ᴴ r(t) / Nr，再用該路徑真實的都卜勒 ε_p 對頻（乘 exp(−j2π ε_p t/N)），做 FFT，得到分支輸出 Z_p[l]；
    // 最後把 L 個分支用 MRC 合併（G = h0，h0 是 L 個分支對所需訊號的等效增益，不考慮 ICI 與路徑間干擾）。
    // 分支觀測向量 z[l] = Σ_m X_m h(m−l) + W，h(d)_p = Σ_i A_i c_{p,i} Q(d + ε_i − ε_p)，c_{p,i} = a(θ_p)ᴴ a(θ_i) / Nr；
    // 分支雜訊協方差 Rn[p][q] = σn² c_{p,q} Q(ε_q − ε_p) / Nr（由 u_p 的定義與 FFT 推得）。
    // IPI（inter-path interference，路徑間干擾）：分支 p 在所需子載波上收到其他路徑的洩漏 Σ_{i≠p} A_i c_{p,i} Q(ε_i − ε_p)，相對於自己路徑的訊號 |A_p|²。
    function angleDomainMf(ch, sigmaX2, sigmaN2) {
        const { N, Nr, paths } = ch, L = paths.length, q = [0, 0];
        const cre = [], cim = [];                                              // c[p][i] = a_pᴴ a_i / Nr
        for (let p = 0; p < L; p++) { cre.push(new Float64Array(L)); cim.push(new Float64Array(L)); for (let i = 0; i < L; i++) { let re = 0, im = 0; for (let k = 0; k < Nr; k++) { re += paths[p].a.re[k] * paths[i].a.re[k] + paths[p].a.im[k] * paths[i].a.im[k]; im += paths[p].a.re[k] * paths[i].a.im[k] - paths[p].a.im[k] * paths[i].a.re[k]; } cre[p][i] = re / Nr; cim[p][i] = im / Nr; } }
        const hr = new Float64Array(N * L), hi = new Float64Array(N * L);
        for (let d = 0; d < N; d++) for (let p = 0; p < L; p++) for (let i = 0; i < L; i++) {
            Qc(d + paths[i].eps - paths[p].eps, N, q);
            const ar = paths[i].Are * cre[p][i] - paths[i].Aim * cim[p][i], ai = paths[i].Are * cim[p][i] + paths[i].Aim * cre[p][i];
            hr[d * L + p] += ar * q[0] - ai * q[1]; hi[d * L + p] += ar * q[1] + ai * q[0];
        }
        const B = new Float64Array(2 * L * L), v = { re: new Float64Array(L), im: new Float64Array(L) };
        for (let d = 1; d < N; d++) { for (let p = 0; p < L; p++) { v.re[p] = hr[d * L + p]; v.im[p] = hi[d * L + p]; } addOuter(B, v, 1, L); }
        const Rn = new Float64Array(2 * L * L);
        for (let p = 0; p < L; p++) for (let r = 0; r < L; r++) { Qc(paths[r].eps - paths[p].eps, N, q); Rn[2 * (p * L + r)] = sigmaN2 * (cre[p][r] * q[0] - cim[p][r] * q[1]) / Nr; Rn[2 * (p * L + r) + 1] = sigmaN2 * (cre[p][r] * q[1] + cim[p][r] * q[0]) / Nr; }
        const h0 = { re: hr.slice(0, L), im: hi.slice(0, L) }, Ri = new Float64Array(2 * L * L);
        for (let i = 0; i < Ri.length; i++) Ri[i] = sigmaX2 * B[i] + Rn[i];
        const sinr = sinrOfWeights(h0, h0, Ri, L, sigmaX2);                    // MRC：G = h0
        // IPI：其他路徑在分支 p 的所需子載波上的洩漏功率，相對於自己路徑的訊號 |A_p|²（對 p 平均）
        let ipi = 0;
        for (let p = 0; p < L; p++) { let lr = 0, li = 0; for (let i = 0; i < L; i++) { if (i === p) continue; Qc(paths[i].eps - paths[p].eps, N, q); const ar = paths[i].Are * cre[p][i] - paths[i].Aim * cim[p][i], ai = paths[i].Are * cim[p][i] + paths[i].Aim * cre[p][i]; lr += ar * q[0] - ai * q[1]; li += ar * q[1] + ai * q[0]; } ipi += (lr * lr + li * li) / (paths[p].Are * paths[p].Are + paths[p].Aim * paths[p].Aim); }
        // 分解（資訊性）：只有雜訊時的 SINR（無 ICI、無 IPI 干擾項）
        let sNoise = 0; { const Rn2 = Rn; sNoise = sinrOfWeights(h0, h0, Rn2, L, sigmaX2); }
        return { sinr, ipiRatio: ipi / L, sinrNoiseOnly: sNoise, h0, B, Rn };
    }

    return { PAPER, makeRng, Qc, QbruteForce, steer, fft, csolve, pinvSolveHermitian, makeChannel, snrToNoise, iciVectors, vecAt, Hvec, iciCov, rMatrix, optimal, sinrOfWeights, system,
        sinrOptimalAt, sinrMrc, sinrSingleAntenna, epsNone, epsLos, epsCpAnalytic, epsCpEstimate, scanEps, epsGenie, linearModel, linearWeights, makeSymbol, demodulate, practicalWeights, maskGuard, paperChannel, angleDomainMf, pilotResidualWeights, ledoitWolfShrink, addOuter };
}));
