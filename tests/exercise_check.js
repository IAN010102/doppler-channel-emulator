'use strict';
/**
 * Node-side equivalent of matlab/check_my_work.m (same cases, same rules, same tolerances) with JavaScript ports of the reference solutions,
 * so that the case data and the grading logic can be shown to be self-consistent without MATLAB. Used by T20f.
 *   ex1 a(theta)   vs a_assumed             relative error < 1e-9
 *   ex2 R_hat      vs R_hat_train_web       relative Frobenius error < 1e-9
 *   ex3 MVDR       SINR vs web SMI          abs error < 1e-6 dB and w^H a = 1 (< 1e-9), cases with L >= N only
 *   ex4 SINR       vs web SINR              abs error < 1e-6 dB
 *   ex5 DL         SINR vs web DL           abs error < 1e-6 dB and w^H a = 1 (< 1e-9)
 * Cases: signalFree, v = 0 (10 cases).
 */
const fs = require('fs'), path = require('path');
const Core = require('../core.js');
const { Cplx, invertMatrix, matMulVec, vecDot } = Core;

const cv = z => z.re.map((r, i) => new Cplx(r, z.im[i]));
const cm = z => z.re.map((row, i) => row.map((r, j) => new Cplx(r, z.im[i][j])));
const TOL_REL = 1e-9, TOL_DB = 1e-6;

// ---- the reference solutions (what matlab/solutions/*.m do), in JS
const REF = {
    ex1: (N, d, thDeg) => Array.from({ length: N }, (_, n) => { const p = -2 * Math.PI * n * d * Math.sin(thDeg * Math.PI / 180); return new Cplx(Math.cos(p), Math.sin(p)); }),
    ex2: X => { const N = X.length, L = X[0].length; return Array.from({ length: N }, (_, m) => Array.from({ length: N }, (_, n) => { let a = new Cplx(0, 0); for (let l = 0; l < L; l++) a = Cplx.add(a, Cplx.mul(X[m][l], Cplx.conj(X[n][l]))); return new Cplx(a.r / L, a.i / L); })); },
    ex3: (R, a) => { const x = matMulVec(invertMatrix(R, false), a), den = vecDot(a, x); return x.map(c => Cplx.div(c, den)); },
    ex4: (w, h, aJ, s2) => 10 * Math.log10(vecDot(w, h).mag2() / (vecDot(w, aJ).mag2() + s2 * w.reduce((s, c) => s + c.mag2(), 0))),
    ex5: (R, a, gr, s2) => { const g = gr * s2, Rd = R.map((row, m) => row.map((c, n) => m === n ? new Cplx(c.r + g, c.i) : c)); const x = matMulVec(invertMatrix(Rd, false), a), den = vecDot(a, x); return x.map(c => Cplx.div(c, den)); }
};
// ---- typical mistakes (the ones check_my_work.m recognises); each must be rejected by the grader
const MISTAKES = {
    'ex1: sign of the phase flipped': { ex1: (N, d, th) => REF.ex1(N, d, th).map(c => Cplx.conj(c)) },
    'ex1: forgot 2*pi': { ex1: (N, d, th) => Array.from({ length: N }, (_, n) => { const p = -n * d * Math.sin(th * Math.PI / 180); return new Cplx(Math.cos(p), Math.sin(p)); }) },
    'ex1: degrees not converted to radians': { ex1: (N, d, th) => Array.from({ length: N }, (_, n) => { const p = -2 * Math.PI * n * d * Math.sin(th); return new Cplx(Math.cos(p), Math.sin(p)); }) },
    'ex2: forgot the division by L': { ex2: X => REF.ex2(X).map(r => r.map(c => new Cplx(c.r * X[0].length, c.i * X[0].length))) },
    'ex2: conjugate on the wrong side (R transposed)': { ex2: X => { const R = REF.ex2(X); return R.map((r, m) => r.map((c, n) => R[n][m])); } },
    'ex3: not normalised (w = R^-1 a)': { ex3: (R, a) => matMulVec(invertMatrix(R, false), a) },
    'ex3: R instead of R^-1': { ex3: (R, a) => { const x = matMulVec(R, a), den = vecDot(a, x); return x.map(c => Cplx.div(c, den)); } },
    'ex4: 20*log10 instead of 10*log10': { ex4: (w, h, aJ, s2) => 2 * REF.ex4(w, h, aJ, s2) },
    'ex4: forgot the noise term': { ex4: (w, h, aJ) => 10 * Math.log10(vecDot(w, h).mag2() / vecDot(w, aJ).mag2()) },
    'ex4: forgot ||w||^2 in the noise term': { ex4: (w, h, aJ, s2) => 10 * Math.log10(vecDot(w, h).mag2() / (vecDot(w, aJ).mag2() + s2)) },
    'ex5: gamma without sigma^2': { ex5: (R, a, gr) => { const Rd = R.map((row, m) => row.map((c, n) => m === n ? new Cplx(c.r + gr, c.i) : c)); const x = matMulVec(invertMatrix(Rd, false), a), den = vecDot(a, x); return x.map(c => Cplx.div(c, den)); } },
    'ex5: loading added to the whole matrix': { ex5: (R, a, gr, s2) => { const g = gr * s2, Rd = R.map(row => row.map(c => new Cplx(c.r + g, c.i))); const x = matMulVec(invertMatrix(Rd, true), a), den = vecDot(a, x); return x.map(c => Cplx.div(c, den)); } }
};

function sinrDb(w, h, aJs, s2) { return 10 * Math.log10(vecDot(w, h).mag2() / (vecDot(w, aJs).mag2() + s2 * w.reduce((s, c) => s + c.mag2(), 0))); }
function fro(A, B) { let d = 0, n = 0; A.forEach((r, i) => r.forEach((c, j) => { d += Math.pow(c.r - B[i][j].r, 2) + Math.pow(c.i - B[i][j].i, 2); n += B[i][j].mag2(); })); return Math.sqrt(d / n); }

// funcs: { ex1..ex5 } (missing ones fall back to the reference); returns, per exercise, 'pass' | 'fail' and the largest error
function check(funcs = {}, casesFile = path.join(__dirname, '..', 'data', 'matlab_cases.json')) {
    const J = JSON.parse(fs.readFileSync(casesFile, 'utf8')), f = Object.assign({}, REF, funcs);
    const res = Array.from({ length: 5 }, () => ({ fail: false, err: 0, n: 0, skipped: 0 })); let used = 0;
    for (const c of J.cases) {
        const P = c.params; if (P.trainMode !== 'signalFree' || P.v_kmh !== 0) continue; used++;
        const N = P.N, L = P.L, Xsf = cm(c.X_sf), Rweb = cm(c.R_hat_train_web), a = cv(c.a_assumed), h = cv(c.h), aJs = cv(c.a_J).map(x => new Cplx(x.r * Math.sqrt(P.P_j), x.i * Math.sqrt(P.P_j)));
        const s2 = P.sigma2, sSMI = c.algorithms.SMI.sinr_dB, sDL = c.algorithms.DL.sinr_dB, wWeb = cv(c.algorithms.SMI.w), th = P.theta1_hat_deg + P.delta_theta_deg;
        const rec = (q, e, tol) => { res[q].n++; res[q].err = Math.max(res[q].err, e); if (!(e < tol)) res[q].fail = true; };
        { const o = f.ex1(N, P.d_over_lambda, th); rec(0, o.length === N ? Math.sqrt(o.reduce((s, x, i) => s + Math.pow(x.r - a[i].r, 2) + Math.pow(x.i - a[i].i, 2), 0) / a.reduce((s, x) => s + x.mag2(), 0)) : Infinity, TOL_REL); }
        { const o = f.ex2(Xsf); rec(1, fro(o, Rweb), TOL_REL); }
        const gradeW = (o, ref) => { const e = Math.abs(sinrDb(o, h, aJs, s2) - ref), con = Math.hypot(vecDot(a, o).r - 1, vecDot(a, o).i); return con < 1e-9 ? e : Math.max(e, con); };   // SINR equal AND w^H a = 1
        if (L >= N) { const o = f.ex3(Rweb, a); rec(2, gradeW(o, sSMI), TOL_DB); } else res[2].skipped++;
        { const o = f.ex4(wWeb, h, aJs, s2); rec(3, Math.abs(o - sSMI), TOL_DB); }
        { const o = f.ex5(Rweb, a, Math.pow(10, P.gamma_rel_dB / 10), s2); rec(4, gradeW(o, sDL), TOL_DB); }
    }
    return { used, res, passed: res.filter(r => !r.fail).length };
}
module.exports = { check, REF, MISTAKES };
