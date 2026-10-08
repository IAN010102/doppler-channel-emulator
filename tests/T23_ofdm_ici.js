'use strict';
/**
 * T23  ofdm_ici.js: OFDM 域 ICI 與接收端波束成形核心（獨立於 core.js）。
 *   T23z  基本恆等式：Q 的封閉式對暴力加總 < 1e-12；FFT 對 DFT < 1e-10；Σ_d |Q(d+ε)|² = 1（Parseval）< 1e-12；
 *         全部子載波占用時 B_l 與 l 無關（< 1e-12）；保護頻帶時 B_l 隨 l 改變且 trace 不增。
 *   T23a  蒙地卡羅時域驗證（N = 64 與 256 各一組，每組 4000 個符元 ≥ 2000）：隨機 QPSK 符元經時變通道（時域直接逐樣本相乘，與解析公式無關）、
 *         FFT，無雜訊。對選定子載波 l：
 *           ICI 功率：mean_s Σ_p |Y_p[l] − X[l] H_p|² 對解析的 σx² Σ_p B_pp（H 用解析值）；
 *           訊號增益：Ĥ_p = mean_s Y_p[l] conj(X[l]) / σx² 對解析的 H_p，各分量取 |z|。
 *         判定：|差| ≤ 3 SE（SE 由樣本標準差 / √4000 估計）；各組的相對差一併回報。
 *   T23b  單一路徑且對頻選到該路徑的都卜勒：ICI 為零（所有 d ≥ 1 的 |H_ICI| < 1e-9），最佳 SINR 等於無 ICI 的 σx²‖H‖²/σn²（相對差 < 1e-9）。
 *   T23c  單一路徑不對頻：ICI 功率 / 總功率 = 1 − |Q(ε)|²（Parseval，恆等式，< 1e-12）；與既有 T5／UI 的 1 − sinc²(ε) 比較：N = 1024 時
 *         絕對差 < 1e-6（判準事先由展開式 (πε/N)²/3 的量級決定）；N = 64、256、1024 的差異回報（資訊性）。
 *   T23d  （資訊性）文獻設定（4 路徑、Nr = 2 與 4）下，精確 Σ_m H_ICI H_ICIᴴ 的最大特徵值佔總功率比例，以及它與線性近似（秩 1）的相對差，
 *         都卜勒縮放 0.1 ... 3 倍（最大 ε 0.007 ... 0.216）；到達角隨機 200 次。
 *   不通過時不放寬判準，回報原因。
 */
const O = require('../ofdm_ici.js');
const U = require('./_util.js');
const D2R = Math.PI / 180;

function eig2(B, Nr) {            // Hermitian Nr×Nr 的特徵值（Jacobi，沿用模組的偽逆工具的嵌入法）
    const m = 2 * Nr, S = Array.from({ length: m }, () => new Float64Array(m)), V = Array.from({ length: m }, (_, i) => { const r = new Float64Array(m); r[i] = 1; return r; });
    for (let i = 0; i < Nr; i++) for (let j = 0; j < Nr; j++) { const a = B[2 * (i * Nr + j)], c = B[2 * (i * Nr + j) + 1]; S[i][j] = a; S[i + Nr][j + Nr] = a; S[i][j + Nr] = -c; S[i + Nr][j] = c; }
    for (let sw = 0; sw < 60; sw++) {
        let off = 0; for (let p = 0; p < m; p++) for (let q = p + 1; q < m; q++) off += S[p][q] * S[p][q]; if (off < 1e-28) break;
        for (let p = 0; p < m - 1; p++) for (let q = p + 1; q < m; q++) {
            if (Math.abs(S[p][q]) < 1e-300) continue;
            const th = (S[q][q] - S[p][p]) / (2 * S[p][q]), t = (th >= 0 ? 1 : -1) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
            for (let k = 0; k < m; k++) { const x = S[k][p], y = S[k][q]; S[k][p] = c * x - s * y; S[k][q] = s * x + c * y; }
            for (let k = 0; k < m; k++) { const x = S[p][k], y = S[q][k]; S[p][k] = c * x - s * y; S[q][k] = s * x + c * y; }
        }
    }
    return Array.from({ length: m }, (_, i) => S[i][i]).sort((a, b) => b - a).filter((_, i) => i % 2 === 0);
}

module.exports = {
    id: 'T23', title: 'ofdm_ici.js: identities, time-domain Monte Carlo of H and ICI, single path with matching demodulation, T5 formula, rank of the ICI covariance (info)',
    async run({ symbols = 4000, seed = 2323 } = {}) {
        const checks = [], info = [];
        const rng = O.makeRng(seed), q = [0, 0];

        // ---------------------------------------------------------------- T23z
        let wq = 0; for (const N of [8, 64, 257]) for (const x of [-3.7, -0.3, 0, 0.05, 0.5, 2.2, N + 0.1, 5.5, -N + 0.4]) { O.Qc(x, N, q); const b = O.QbruteForce(x, N); wq = Math.max(wq, Math.hypot(q[0] - b[0], q[1] - b[1])); }
        checks.push(U.check('T23z Q(x): closed form vs brute-force sum', U.e(wq), '< 1e-12', wq < 1e-12));
        const nf = 16, fr = Float64Array.from({ length: nf }, () => rng.u()), fi = Float64Array.from({ length: nf }, () => rng.u()), gr = Float64Array.from(fr), gi = Float64Array.from(fi); O.fft(gr, gi, false);
        let wf = 0; for (let k = 0; k < nf; k++) { let a = 0, b = 0; for (let t = 0; t < nf; t++) { const ph = -2 * Math.PI * k * t / nf; a += fr[t] * Math.cos(ph) - fi[t] * Math.sin(ph); b += fr[t] * Math.sin(ph) + fi[t] * Math.cos(ph); } wf = Math.max(wf, Math.hypot(a - gr[k], b - gi[k])); }
        checks.push(U.check('T23z FFT vs direct DFT (N = 16)', U.e(wf), '< 1e-10', wf < 1e-10));
        let wp = 0; for (const N of [64, 1024]) for (const e of [0.013, 0.2, -0.37, 0.5, 1.3]) { let s = 0; for (let d = 0; d < N; d++) { O.Qc(d + e, N, q); s += q[0] * q[0] + q[1] * q[1]; } wp = Math.max(wp, Math.abs(s - 1)); }
        checks.push(U.check('T23z Parseval: sum_d |Q(d + eps)|^2 = 1', U.e(wp), '< 1e-12', wp < 1e-12));
        {
            const ch = O.paperChannel([5, -40, 25, 60], O.makeRng(1), { N: 256 }), V = O.iciVectors(ch, 0);
            const Ba = O.iciCov(V, null, 0), Bb = O.iciCov(V, null, 100); let dB = 0; for (let i = 0; i < Ba.length; i++) dB = Math.max(dB, Math.abs(Ba[i] - Bb[i]));
            checks.push(U.check('T23z all subcarriers occupied: B_l does not depend on l', U.e(dB), '< 1e-12', dB < 1e-12));
            const occ = O.maskGuard(256, 20, 20), Be = O.iciCov(V, occ, 21), Bc = O.iciCov(V, occ, 128), trace = B => { let t = 0; for (let p = 0; p < ch.Nr; p++) t += B[2 * (p * ch.Nr + p)]; return t; };
            checks.push(U.check('T23z guard band: ICI trace at an edge subcarrier differs from the centre and is not larger than with all subcarriers', `edge ${U.f(trace(Be), 6)}, centre ${U.f(trace(Bc), 6)}, all occupied ${U.f(trace(Ba), 6)}`, 'edge != centre; both <= all occupied', Math.abs(trace(Be) - trace(Bc)) > 1e-9 && trace(Be) <= trace(Ba) + 1e-12 && trace(Bc) <= trace(Ba) + 1e-12));
        }

        // ---------------------------------------------------------------- T23a
        console.log(`T23a  time-domain Monte Carlo (noise free), Nr = 2, ${symbols} symbols per N; relative difference of the ICI power and the signal gain with their standard errors`);
        console.log(U.pad('N', 6), U.pad('l', 6), U.rpad('ICI pow analytic', 17), U.rpad('measured', 11), U.rpad('rel diff', 10), U.rpad('SE rel', 9), U.rpad('z', 7), U.rpad('|H| analytic', 13), U.rpad('max |z| of H_hat', 17));
        for (const N of [64, 256]) {
            const ch = O.paperChannel([10, -35, 50, 5], O.makeRng(100 + N), { N }), Nr = ch.Nr, sx2 = 1, V = O.iciVectors(ch, 0), H = O.vecAt(V, 0), B = O.iciCov(V, null, 0);
            let analytic = 0; for (let p = 0; p < Nr; p++) analytic += sx2 * B[2 * (p * Nr + p)];
            const ls = [N / 4, N - 1], acc = ls.map(() => ({ pow: [], gr: Array.from({ length: Nr }, () => []), gi: Array.from({ length: Nr }, () => []) }));
            const r2 = O.makeRng(7000 + N), a = Math.sqrt(sx2 / 2);
            for (let s = 0; s < symbols; s++) {
                const Xr = new Float64Array(N), Xi = new Float64Array(N); for (let m = 0; m < N; m++) { Xr[m] = r2.u() < 0.5 ? -a : a; Xi[m] = r2.u() < 0.5 ? -a : a; }
                const xr = Float64Array.from(Xr), xi = Float64Array.from(Xi); O.fft(xr, xi, true);
                const Yr = [], Yi = [];
                for (let p = 0; p < Nr; p++) {
                    const yr = new Float64Array(N), yi = new Float64Array(N);
                    for (let n = 0; n < N; n++) {                                  // 時域逐樣本：h_p(n) = Σ_i A_i a_p(θ_i) exp(j2π ε_i n / N)，與解析公式無關
                        let hr = 0, hi = 0; for (const pa of ch.paths) { const ph = 2 * Math.PI * pa.eps * n / N, c = Math.cos(ph), sn = Math.sin(ph), ar = pa.Are * c - pa.Aim * sn, ai = pa.Are * sn + pa.Aim * c; hr += ar * pa.a.re[p] - ai * pa.a.im[p]; hi += ar * pa.a.im[p] + ai * pa.a.re[p]; }
                        yr[n] = hr * xr[n] - hi * xi[n]; yi[n] = hr * xi[n] + hi * xr[n];
                    }
                    O.fft(yr, yi, false); for (let n = 0; n < N; n++) { yr[n] /= N; yi[n] /= N; } Yr.push(yr); Yi.push(yi);
                }
                ls.forEach((l, k) => {
                    let pw = 0; for (let p = 0; p < Nr; p++) {
                        const er = Yr[p][l] - (Xr[l] * H.re[p] - Xi[l] * H.im[p]), ei = Yi[p][l] - (Xr[l] * H.im[p] + Xi[l] * H.re[p]); pw += er * er + ei * ei;
                        acc[k].gr[p].push((Yr[p][l] * Xr[l] + Yi[p][l] * Xi[l]) / sx2); acc[k].gi[p].push((Yi[p][l] * Xr[l] - Yr[p][l] * Xi[l]) / sx2);
                    }
                    acc[k].pow.push(pw);
                });
            }
            ls.forEach((l, k) => {
                const m = U.mean(acc[k].pow), se = U.se(acc[k].pow), rel = m / analytic - 1, relSe = se / analytic, z = (m - analytic) / se; let zmax = 0;
                for (let p = 0; p < Nr; p++) { zmax = Math.max(zmax, Math.abs(U.mean(acc[k].gr[p]) - H.re[p]) / U.se(acc[k].gr[p]), Math.abs(U.mean(acc[k].gi[p]) - H.im[p]) / U.se(acc[k].gi[p])); }
                const hn = Math.sqrt(H.re.reduce((s2, x, i) => s2 + x * x + H.im[i] * H.im[i], 0));
                console.log(U.pad(N, 6), U.pad(l, 6), U.rpad(U.e(analytic, 4), 17), U.rpad(U.e(m, 4), 11), U.rpad((100 * rel).toFixed(3) + ' %', 10), U.rpad((100 * relSe).toFixed(3) + ' %', 9), U.rpad(z.toFixed(2), 7), U.rpad(U.f(hn, 4), 13), U.rpad(zmax.toFixed(2), 17));
                checks.push(U.check(`T23a N = ${N}, subcarrier ${l}: ICI power vs sigma_x^2 sum |H_ICI|^2`, `${(100 * rel).toFixed(3)} % (SE ${(100 * relSe).toFixed(3)} %), z = ${z.toFixed(2)}`, '|difference| <= 3 SE', Math.abs(z) <= 3));
                checks.push(U.check(`T23a N = ${N}, subcarrier ${l}: signal gain H_hat vs analytic H (largest |z| of the ${2 * Nr} components)`, `${zmax.toFixed(2)} SE`, '<= 3 SE', zmax <= 3));
            });
        }

        // ---------------------------------------------------------------- T23b
        let wIci = 0, wSinr = 0; const cases = [];
        for (const N of [64, 1024]) for (const Nr of [2, 4]) for (const eps of [0.05, -0.3, 0.72]) {
            const ch = O.makeChannel({ N, Nr, paths: [{ eps, thetaDeg: 20, powerDb: 0, phase: 0.7 }] }), V = O.iciVectors(ch, eps);
            for (let d = 1; d < N; d++) for (let p = 0; p < Nr; p++) wIci = Math.max(wIci, Math.hypot(V.hr[d * Nr + p], V.hi[d * Nr + p]));
            const snr = 12.5, sn2 = O.snrToNoise(ch, snr, 1), s = O.sinrOptimalAt(ch, eps, 1, sn2), want = Nr * Math.pow(10, snr / 10);         // σx²‖H‖²/σn² = Nr · SNR（陣列增益）
            wSinr = Math.max(wSinr, Math.abs(s / want - 1)); cases.push(`N=${N},Nr=${Nr},eps=${eps}`);
        }
        checks.push(U.check('T23b single path, demodulated with its own Doppler: all |H_ICI(d >= 1)|', U.e(wIci), '< 1e-9', wIci < 1e-9, `${cases.length} cases`));
        checks.push(U.check('T23b same: optimal SINR equals the SNR without ICI (array gain Nr x SNR), max relative difference', U.e(wSinr), '< 1e-9', wSinr < 1e-9));

        // ---------------------------------------------------------------- T23c
        console.log('\nT23c  single path, no demodulation: ICI power / total power = 1 - |Q(eps)|^2 against the T5 / UI formula 1 - sinc^2(eps)');
        console.log(U.pad('eps', 7), [64, 256, 1024].map(N => U.rpad(`N=${N}: exact`, 18) + U.rpad('1-sinc^2', 12) + U.rpad('diff', 11)).join(''));
        let wPar = 0, worstTol = 0;
        for (const eps of [0.01, 0.05, 0.1, 0.2, 0.3]) {
            const cells = [], sinc2 = Math.pow(Math.sin(Math.PI * eps) / (Math.PI * eps), 2);
            for (const N of [64, 256, 1024]) {
                const ch = O.makeChannel({ N, Nr: 1, paths: [{ eps, thetaDeg: 0, powerDb: 0, phase: 0 }] }), V = O.iciVectors(ch, 0), B = O.iciCov(V, null, 0), qq = V.hr[0] * V.hr[0] + V.hi[0] * V.hi[0], ratio = B[0] / (B[0] + qq);
                wPar = Math.max(wPar, Math.abs(ratio - (1 - qq))); const diff = ratio - (1 - sinc2); if (N === 1024) worstTol = Math.max(worstTol, Math.abs(diff));
                cells.push(U.rpad(U.e(ratio, 4), 18) + U.rpad(U.e(1 - sinc2, 4), 12) + U.rpad(U.e(diff, 2), 11));
            }
            console.log(U.pad(eps, 7), cells.join(''));
        }
        checks.push(U.check('T23c ICI power / total power = 1 - |Q(eps)|^2 (Parseval identity)', U.e(wPar), '< 1e-12', wPar < 1e-12));
        checks.push(U.check('T23c N = 1024: ICI ratio against 1 - sinc^2(eps) of T5 / the UI (eps 0.01 ... 0.3)', U.e(worstTol), 'absolute difference < 1e-6', worstTol < 1e-6, 'finite-N differences for N = 64, 256, 1024 are printed above'));

        // ---------------------------------------------------------------- T23d (informational)
        console.log('\nT23d  (informational) share of the largest eigenvalue in sum_m H_ICI H_ICI^H (4 paths of the paper setting, random angles, 200 draws) and distance to the rank-1 linear model');
        console.log(U.pad('Nr', 4), U.pad('Doppler scale', 14), U.pad('max eps', 9), U.rpad('lambda_max share %', 20), U.rpad('||B - B_lin||/||B||', 21));
        const r3 = O.makeRng(2323);
        for (const Nr of [2, 4]) for (const sc of [0.1, 0.25, 0.5, 1, 2, 3]) {
            const sh = [], dd = [];
            for (let t = 0; t < 200; t++) {
                const ang = [0, 1, 2, 3].map(() => -90 + 180 * r3.u()), ch = O.paperChannel(ang, r3, { Nr, dopplerScale: sc, N: 1024 }), V = O.iciVectors(ch, 0), B = O.iciCov(V, null, 0), ev = eig2(B, Nr);
                sh.push(100 * ev[0] / ev.reduce((a, b) => a + b, 0));
                const lm = O.linearModel(ch, 0), Bl = new Float64Array(B.length); for (let i = 0; i < Nr; i++) for (let j = 0; j < Nr; j++) { Bl[2 * (i * Nr + j)] = lm.sumXi2 * (lm.h1.re[i] * lm.h1.re[j] + lm.h1.im[i] * lm.h1.im[j]); Bl[2 * (i * Nr + j) + 1] = lm.sumXi2 * (lm.h1.im[i] * lm.h1.re[j] - lm.h1.re[i] * lm.h1.im[j]); }
                let num = 0, den = 0; for (let i = 0; i < B.length; i++) { num += (B[i] - Bl[i]) ** 2; den += B[i] ** 2; } dd.push(Math.sqrt(num / den));
            }
            const maxEps = Math.max(...O.PAPER.dopplerHz.map(f => Math.abs(f * sc / O.PAPER.df)));
            console.log(U.pad(Nr, 4), U.pad(sc, 14), U.pad(maxEps.toFixed(3), 9), U.rpad(`${U.f(U.mean(sh), 2)} ± ${U.f(U.se(sh), 2)}`, 20), U.rpad(`${U.f(U.mean(dd), 4)} ± ${U.f(U.se(dd), 4)}`, 21));
            info.push({ name: `T23d Nr=${Nr}, Doppler x${sc} (max eps ${maxEps.toFixed(3)})`, value: `largest eigenvalue ${U.f(U.mean(sh), 2)} ± ${U.f(U.se(sh), 2)} % of the ICI power; distance to the rank-1 linear model ${U.f(U.mean(dd), 4)} ± ${U.f(U.se(dd), 4)}` });
        }
        return { id: this.id, title: this.title, checks, info };
    }
};
