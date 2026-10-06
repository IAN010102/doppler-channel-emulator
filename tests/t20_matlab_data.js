'use strict';
/**
 * T20  The files for the MATLAB scripts (Commits 23-24). The MATLAB scripts themselves cannot run here; these checks cover everything that can be checked in Node.
 *   T20a  tests/make_matlab_cases.js run twice (same seeds): byte-identical JSON.
 *   T20b  The weights stored in data/matlab_cases.json are recomputed from the data stored in the same file (X / X_sf, s, a_assumed, B, parameters) with core.js; they must agree with the stored
 *         values (< 1e-12 relative), and the stored SINR / SINR_opt must agree with the formulas of the JSON header evaluated on the stored h, a_J, w (1e-9 dB). This shows that the export
 *         has no omissions or misplaced entries.
 *   T20c  Every CSV in data/: header columns equal the columns listed for that file in data/README_data.md; no NaN / empty field.
 *   T20d  Every curve point has a mean and a standard error (SE >= 0, finite) and the number of trials is >= the required minimum (>= 500 for the speed curve, >= 1000 otherwise).
 *   T20e  Every exercise file (matlab/exercises/ex1 ... ex5) and every solution file (matlab/solutions/ex*_solution.m) exists; each exercise header has the sections 【生活比喻】, 【輸入輸出】 and
 *         【對應數學式】 (and the exercise has a %%% TODO block, the 「這一關還沒寫」 placeholder message and a 【提示】 section); the solutions have comments; check_my_work.m, selftest.m,
 *         verify_cases.m, mc_independent.m, plot_figures.m and matlab/README.md exist.
 *   T20f  Node-side equivalent of check_my_work.m (tests/exercise_check.js: same cases, rules and tolerances): the reference solutions pass 5 / 5 exercises on all cases; each typical mistake
 *         that the MATLAB grader recognises is rejected (the exercise it belongs to fails, the others still pass).
 */
const fs = require('fs'), path = require('path');
const root = path.join(__dirname, '..');
const U = require('./_util.js');
const Core = require('../core.js');
const { Cplx, invertMatrix, pinvHermitian, hermitianEigvals, matMulVec, vecDot } = Core;
const { build } = require('./make_matlab_cases.js');

const toVec = v => v.re.map((r, i) => new Cplx(r, v.im[i]));
const toMat = M => M.re.map((row, i) => row.map((r, j) => new Cplx(r, M.im[i][j])));

// ---- recomputation of the weights from the case data (the same operations as Sys.computeMath, so that the result can be compared at 1e-12)
function trainingR(X, N, L) {
    const Rr = new Float64Array(N * N), Ri = new Float64Array(N * N);
    for (let l = 0; l < L; l++) for (let m = 0; m < N; m++) for (let n = 0; n < N; n++) {
        const xm = X[m][l], xn = X[n][l];
        Rr[m * N + n] += xm.r * xn.r + xm.i * xn.i; Ri[m * N + n] += xm.i * xn.r - xm.r * xn.i;
    }
    return Array.from({ length: N }, (_, m) => Array.from({ length: N }, (_, n) => new Cplx(Rr[m * N + n] / L, Ri[m * N + n] / L)));
}
function rankOf(R, eps) { const ev = hermitianEigvals(R), lmax = ev[ev.length - 1]; let rk = 0; for (let i = ev.length - 1; i >= 0; i -= 2) if (ev[i] > eps * lmax) rk++; return rk; }
function mvdr(R, a, force) { const inv = invertMatrix(R, force); if (!inv) return null; const num = matMulVec(inv, a), den = vecDot(a, num); return num.map(c => Cplx.div(c, den)); }
function weightsFromCase(c, algo) {
    const P = c.params, N = P.N, L = P.L, eps = P.epsRank, a = toVec(c.a_assumed), sym = toVec(c.s);
    const Xfull = toMat(c.X), Xsf = toMat(c.X_sf), Xtr = P.trainMode === 'signalFree' ? Xsf : Xfull;
    if (algo === 'FOURIER') return a.map(x => new Cplx(x.r * 1 / N, x.i * 1 / N));
    if (algo === 'MMSE-M' || algo === 'MMSE-P') {
        const R = trainingR(Xfull, N, L); let rxd;
        if (algo === 'MMSE-M') rxd = a.map(x => new Cplx(x.r * P.P_s, x.i * P.P_s));
        else { const qr = new Float64Array(N), qi = new Float64Array(N); for (let l = 0; l < L; l++) for (let m = 0; m < N; m++) { const x = Xfull[m][l]; qr[m] += x.r * sym[l].r + x.i * sym[l].i; qi[m] += x.i * sym[l].r - x.r * sym[l].i; } rxd = Array.from({ length: N }, (_, m) => new Cplx(qr[m] / L, qi[m] / L)); }
        if (P.smiSingular === 'pinv' && rankOf(R, eps) < N) return matMulVec(pinvHermitian(R, eps).pinv, rxd);
        return matMulVec(invertMatrix(R, true), rxd);
    }
    const R = trainingR(Xtr, N, L);
    if (algo === 'SMI') {
        if (P.smiSingular === 'pinv' && rankOf(R, eps) < N) { const num = matMulVec(pinvHermitian(R, eps).pinv, a), den = vecDot(a, num); return num.map(x => Cplx.div(x, den)); }
        return mvdr(R, a, true);
    }
    if (algo === 'DL') return mvdr(R.map((row, m) => row.map((x, n) => m === n ? new Cplx(x.r + P.gamma_abs, x.i) : x)), a, false);
    if (algo === 'BEAMSPACE') {
        const Bm = toMat(c.beamspace.B), K = c.beamspace.bins.length, B = Array.from({ length: K }, (_, j) => Array.from({ length: N }, (_, n) => Bm[n][j]));
        const aB = B.map(b => vecDot(b, a)), RBc = B.map(b => matMulVec(R, b)), M = B.map(bi => RBc.map(rb => vecDot(bi, rb)));
        const inv = invertMatrix(M, true), num = matMulVec(inv, aB), den = vecDot(aB, num), wB = num.map(x => Cplx.div(x, den));
        return Array.from({ length: N }, (_, n) => { let acc = new Cplx(0, 0); for (let j = 0; j < K; j++) acc = Cplx.add(acc, Cplx.mul(wB[j], B[j][n])); return acc; });
    }
    throw new Error('unknown algorithm ' + algo);
}
function sinrDb(w, h, aJ, Pj, s2) {
    const wh = vecDot(w, h), wj = vecDot(w, aJ), nw = w.reduce((x, c) => x + c.mag2(), 0);
    return 10 * Math.log10(wh.mag2() / (Pj * wj.mag2() + s2 * nw));
}
function sinrOptDb(h, aJ, Pj, s2) {
    const h2 = h.reduce((x, c) => x + c.mag2(), 0), g2 = aJ.reduce((x, c) => x + c.mag2(), 0), gh = vecDot(aJ, h);
    return 10 * Math.log10((h2 - Pj * gh.mag2() / (s2 + Pj * g2)) / s2);
}

// ---- CSV helpers
function readCsv(file) {
    const lines = fs.readFileSync(file, 'utf8').split(/\r?\n/).filter(l => l.length), comments = lines.filter(l => l.startsWith('#')), rest = lines.filter(l => !l.startsWith('#'));
    return { comments, cols: rest[0].split(','), rows: rest.slice(1).map(l => l.split(',')) };
}
const MIN_TRIALS = { 'curve_evm_vs_velocity.csv': 500, 'curve_sinr_vs_L.csv': 1000, 'curve_sinr_vs_delta_theta.csv': 1000, 'curve_doa_spectrum.csv': 1000 };

module.exports = {
    id: 'T20', title: 'MATLAB data files: reproducible export, weights recomputed from the exported data, CSV columns and standard errors (MATLAB scripts not run here)',
    async run() {
        const checks = [], info = [];
        // ---------------------------------------------------------------- T20a
        const j1 = JSON.stringify(build()), j2 = JSON.stringify(build()), same = j1 === j2;
        console.log(`T20a  make_matlab_cases twice: ${j1.length} characters, identical: ${same}`);
        checks.push(U.check('T20a case export is bit-identical on a second run (same seeds)', `${same}`, 'true', same));
        const onDisk = fs.readFileSync(path.join(root, 'data', 'matlab_cases.json'), 'utf8');
        checks.push(U.check('T20a data/matlab_cases.json is the current output of the generator', `${onDisk === j1}`, 'true', onDisk === j1, 'regenerate with node tests/make_matlab_cases.js'));

        // ---------------------------------------------------------------- T20b
        const J = JSON.parse(onDisk);
        let worstW = 0, worstS = 0, worstO = 0, n = 0;
        for (const c of J.cases) {
            const P = c.params, h = toVec(c.h), aJ = toVec(c.a_J);
            for (const algo of Object.keys(c.algorithms)) {
                const w = weightsFromCase(c, algo), ws = toVec(c.algorithms[algo].w);
                let d = 0, nn = 0; w.forEach((x, i) => { d += Math.pow(x.r - ws[i].r, 2) + Math.pow(x.i - ws[i].i, 2); nn += ws[i].mag2(); });
                worstW = Math.max(worstW, Math.sqrt(d / nn));
                worstS = Math.max(worstS, Math.abs(sinrDb(ws, h, aJ, P.P_j, P.sigma2) - c.algorithms[algo].sinr_dB));
                worstO = Math.max(worstO, Math.abs(sinrOptDb(h, aJ, P.P_j, P.sigma2) - c.algorithms[algo].sinr_opt_dB));
                n++;
            }
        }
        console.log(`T20b  ${J.cases.length} cases x 6 algorithms = ${n}: max relative weight difference ${U.e(worstW)}, max SINR difference ${U.e(worstS)} dB, max SINR_opt difference ${U.e(worstO)} dB`);
        checks.push(U.check('T20b weights recomputed from the exported data equal the stored weights', U.e(worstW), '< 1e-12 (relative)', worstW < 1e-12));
        checks.push(U.check('T20b stored SINR equals the header formula on the stored h, a_J, w', U.e(worstS) + ' dB', '< 1e-9 dB', worstS < 1e-9));
        checks.push(U.check('T20b stored SINR_opt equals h^H R_in^-1 h on the stored h, a_J', U.e(worstO) + ' dB', '< 1e-9 dB', worstO < 1e-9));
        checks.push(U.check('T20b case list: 20 static + 2 dynamic + the lecture baseline E0 = 23 cases', `${J.cases.length}`, '23', J.cases.length === 23 && J.cases[22].id === 'E0_lab_baseline'));

        // ---------------------------------------------------------------- T20c / T20d
        const readme = fs.readFileSync(path.join(root, 'data', 'README_data.md'), 'utf8');
        for (const file of Object.keys(MIN_TRIALS)) {
            const p = path.join(root, 'data', file);
            if (!fs.existsSync(p)) { checks.push(U.check(`T20c ${file}: exists`, 'missing', 'present', false, 'run node tests/make_curves.js')); continue; }
            const csv = readCsv(p), sec = readme.split(/^## /m).find(s => s.startsWith(file)) || '';
            const listed = [...sec.matchAll(/^\| `([A-Za-z_0-9]+)` \|/gm)].map(m => m[1]);
            const same = csv.cols.length === listed.length && csv.cols.every((c, i) => c === listed[i]);
            let bad = 0; csv.rows.forEach(r => { if (r.length !== csv.cols.length || r.some(x => x === '' || /nan|inf/i.test(x))) bad++; });
            const hasSchema = csv.comments.some(l => /csv_schema_version=\d+/.test(l)) && csv.comments.some(l => /seed0=\d+/.test(l)) && csv.comments.some(l => /trials/.test(l));
            console.log(`T20c  ${file}: ${csv.rows.length} rows, columns ${csv.cols.join(' ')}; README lists the same columns: ${same}; bad rows: ${bad}; header comments complete: ${hasSchema}`);
            checks.push(U.check(`T20c ${file}: columns equal the README list; no NaN or empty field; header comments`, `${same ? 'columns equal' : 'columns differ'}, ${bad} bad rows, comments ${hasSchema ? 'ok' : 'incomplete'}`, 'equal, 0 bad rows, complete', same && bad === 0 && hasSchema));
            const seCols = csv.cols.map((c, i) => [c, i]).filter(([c]) => /_se_/.test(c) || /_se$/.test(c) || /se_(dB|pct)/.test(c)), meanCols = csv.cols.filter(c => /mean/.test(c) || /^evm_mean|ici_floor_evm/.test(c));
            const nIdx = csv.cols.indexOf('n_trials');
            let seBad = 0, nMin = Infinity; csv.rows.forEach(r => { seCols.forEach(([, i]) => { const v = parseFloat(r[i]); if (!(Number.isFinite(v) && v >= 0)) seBad++; }); nMin = Math.min(nMin, parseFloat(r[nIdx])); });
            const ok = seCols.length >= 1 && meanCols.length >= 1 && seBad === 0 && nMin >= MIN_TRIALS[file];
            checks.push(U.check(`T20d ${file}: every point has a mean and a standard error >= 0, trials >= ${MIN_TRIALS[file]}`, `${seCols.length} SE column(s), ${seBad} bad SE, min trials ${nMin}`, `SE >= 0, trials >= ${MIN_TRIALS[file]}`, ok));
        }
        // ---------------------------------------------------------------- T20e
        const mdir = path.join(root, 'matlab'), names = ['ex1_steering', 'ex2_covariance', 'ex3_mvdr', 'ex4_sinr', 'ex5_diagload'], miss = [], badHead = [];
        for (const n of names) {
            const fe = path.join(mdir, 'exercises', n + '.m'), fs_ = path.join(mdir, 'solutions', n + '_solution.m');
            if (!fs.existsSync(fe)) { miss.push(fe); continue; } if (!fs.existsSync(fs_)) miss.push(fs_);
            const t = fs.readFileSync(fe, 'utf8'), head = t.split(/\r?\n/).filter(l => l.startsWith('%')).join('\n');
            for (const sec of ['【生活比喻】', '【輸入輸出】', '【對應數學式】', '【這一關要做什麼】', '【提示】']) if (!head.includes(sec)) badHead.push(`${n}: ${sec}`);
            if (!t.includes('%%% TODO')) badHead.push(`${n}: %%% TODO`);
            if (!t.includes('這一關還沒寫')) badHead.push(`${n}: placeholder message`);
            if (!/NaN|nan\(/.test(t)) badHead.push(`${n}: NaN placeholder`);
            if (fs.existsSync(fs_) && (fs.readFileSync(fs_, 'utf8').match(/%/g) || []).length < 4) badHead.push(`${n}_solution: too few comments`);
        }
        for (const f of ['check_my_work.m', 'selftest.m', 'verify_cases.m', 'mc_independent.m', 'plot_figures.m', 'lab_mmse_mvdr_demo.m', 'README.md']) if (!fs.existsSync(path.join(mdir, f))) miss.push('matlab/' + f);
        console.log(`T20e  missing files: ${miss.join(', ') || 'none'}; header problems: ${badHead.join('; ') || 'none'}`);
        checks.push(U.check('T20e exercise, solution and script files exist; exercise headers have the required sections', `${miss.length} missing, ${badHead.length} problems`, '0 / 0', miss.length === 0 && badHead.length === 0, miss.concat(badHead).join('; ')));

        // ---------------------------------------------------------------- T20f
        const EX = require('./exercise_check.js'), ok = EX.check();
        console.log(`T20f  reference solutions on ${ok.used} cases: passed ${ok.passed} / 5; largest errors ${ok.res.map(r => U.e(r.err)).join(' / ')}; skipped (ex3, L < N): ${ok.res[2].skipped}`);
        checks.push(U.check('T20f reference solutions pass 5 / 5 exercises on all cases (Node-side equivalent of check_my_work)', `${ok.passed} / 5 (${ok.used} cases)`, '5 / 5', ok.passed === 5 && ok.used === 10));
        let wrong = [];
        for (const [name, funcs] of Object.entries(EX.MISTAKES)) {
            const r = EX.check(funcs), q = parseInt(Object.keys(funcs)[0].slice(2), 10) - 1;
            const rejected = r.res[q].fail, others = r.res.every((x, i) => i === q || !x.fail);
            console.log(`T20f  mistake "${name}": exercise ${q + 1} ${rejected ? 'rejected' : 'ACCEPTED'}, other exercises ${others ? 'unaffected' : 'AFFECTED'}`);
            if (!rejected || !others) wrong.push(name);
        }
        checks.push(U.check(`T20f the grader rejects each of the ${Object.keys(EX.MISTAKES).length} typical mistakes (and only the exercise it belongs to)`, `${wrong.length} not as expected`, '0', wrong.length === 0, wrong.join('; ')));
        return { id: this.id, title: this.title, checks, info };
    }
};
