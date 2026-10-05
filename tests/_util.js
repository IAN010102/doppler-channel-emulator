'use strict';
// Small helpers shared by the tests (no dependencies).

function mean(a) { let s = 0; for (const x of a) s += x; return s / a.length; }
function variance(a, m = mean(a)) { let s = 0; for (const x of a) s += (x - m) * (x - m); return s / Math.max(1, a.length - 1); }
function se(a) { return Math.sqrt(variance(a) / a.length); }
const pad = (s, n) => String(s).padEnd(n);
const rpad = (s, n) => String(s).padStart(n);
const f = (x, d = 3) => Number.isFinite(x) ? x.toFixed(d) : String(x);
const e = (x, d = 2) => Number.isFinite(x) ? x.toExponential(d) : String(x);
const db10 = x => 10 * Math.log10(x);
const lin10 = x => Math.pow(10, x / 10);

// every test returns { id, title, checks: [{ name, value, tol, pass, note }] }
function check(name, value, tol, pass, note = '') { return { name, value, tol, pass: !!pass, note }; }

function cloneR(R) { return R.map(row => row.map(c => ({ r: c.r, i: c.i }))); }
function maxAbsDiffR(A, B) {
    let m = 0;
    for (let i = 0; i < A.length; i++) for (let j = 0; j < A.length; j++) m = Math.max(m, Math.hypot(A[i][j].r - B[i][j].r, A[i][j].i - B[i][j].i));
    return m;
}

module.exports = { mean, variance, se, pad, rpad, f, e, db10, lin10, check, cloneR, maxAbsDiffR };
