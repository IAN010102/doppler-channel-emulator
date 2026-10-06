'use strict';
/**
 * Sensitivity diagnostics for E6 (informational, not part of run_all; no behaviour or criterion is changed).
 *   node tests/diag_s_sensitivity.js s1 [trials]   S1: ICI loss and aging loss on a grid  -> docs/diagnostics/s1_sensitivity.csv / .txt
 *   node tests/diag_s_sensitivity.js s2 [trials]   S2: why DL loses more than SMI at v = 300 (E6 setting)   -> docs/diagnostics/s2_dl_vs_smi.txt
 *   node tests/diag_s_sensitivity.js s3 [trials]   S3: SMI withSignal rise-then-fall: K -> inf / diffuse only / original, target covariance rank  -> docs/diagnostics/s3_withsignal.txt
 * Common setting = E6: N = 8, L = 100, theta1 = 20, theta2 = -30, SNR 20, SIR -10, K = 20 dB (unless stated), unified model, QPSK, signalFree (unless stated).
 * Definitions (PARAMS.md section 21):  SINR_eff = 1/(1/SINR + N_ICI/S);  ICI loss [dB] = 10 log10(SINR_tau0 / SINR_eff_tau0);  aging loss [dB] = SINR(tau = 0) - SINR(tau).
 */
const fs = require('fs'), path = require('path');
const Core = require('../core.js');
const U = require('./_util.js');
const mode = process.argv[2] || 's1', trials = +(process.argv[3] || 500), D2R = Math.PI / 180;
const OUT = path.join(__dirname, '..', 'docs', 'diagnostics');
const BASE = { model: 'unified', N: 8, L: 100, aoaT: 20, aoaJ: -30, snr: 20, sir: -10, kDb: 20, calDeg: 0, taper: 'NONE', mod: 'QPSK', gammaRelDb: 10, d_min: 10, latMs: 10, v: 300, fc: 5e9, scs: 15e3, trainMode: 'signalFree', jamWave: 'gaussian' };
function run(over, seed) {
    Core.setSeed(seed);
    const s = Core.createSys(); s.calZ = new Array(16).fill(0); Object.assign(s, BASE, over, { freshRealization: true }); s.snaps = []; s.snapKey = ''; s.computeMath(); return s;
}
const f2 = x => U.f(x, 2), cellS = a => `${f2(U.mean(a))}±${f2(U.se(a))}`;
const lin = db => Math.pow(10, db / 10), db10 = x => 10 * Math.log10(x);
const write = (name, txt) => { fs.mkdirSync(OUT, { recursive: true }); fs.writeFileSync(path.join(OUT, name), txt); };

if (mode === 's1') {
    const rows = [], lines = [];
    const head = ['fc_GHz', 'scs_kHz', 'd_min_m', 'tau_ms', 'v_kmh', 'algo', 'sinr_tau0_dB', 'sinr_tau0_se', 'sinr_dB', 'sinr_se', 'aging_loss_dB', 'aging_se', 'nici_dB', 'ici_loss_dB', 'ici_se', 'aging_ge_ici'];
    const out = [head.join(',')];
    for (const [fc, scs] of [[5e9, 15e3], [28e9, 120e3]]) for (const dmin of [10, 30, 100, 500]) for (const algo of ['SMI', 'DL']) for (const v of [25, 50, 100, 200, 300]) {
        const base = [], eff = [], nici = [], per = { 1: [], 2: [], 10: [] };
        for (let t = 0; t < trials; t++) {
            const s0 = run({ fc, scs, d_min: dmin, algo, v, latMs: 0 }, 5000 + t);
            base.push(s0.sinrDb); const nu = s0.nuICI; nici.push(db10(nu + 1e-30)); eff.push(db10(1 / (1 / lin(s0.sinrDb) + nu)));
            for (const tau of [1, 2, 10]) per[tau].push(run({ fc, scs, d_min: dmin, algo, v, latMs: tau }, 5000 + t).sinrDb);
        }
        const ici = base.map((b, i) => b - eff[i]);
        for (const tau of [1, 2, 10]) {
            const ag = base.map((b, i) => b - per[tau][i]), ge = U.mean(ag) >= U.mean(ici);
            const r = [fc / 1e9, scs / 1e3, dmin, tau, v, algo, U.mean(base), U.se(base), U.mean(per[tau]), U.se(per[tau]), U.mean(ag), U.se(ag), U.mean(nici), U.mean(ici), U.se(ici), ge ? 1 : 0];
            rows.push({ fc, scs, dmin, tau, v, algo, ag: U.mean(ag), agse: U.se(ag), ici: U.mean(ici), icise: U.se(ici), ge });
            out.push(r.map(x => typeof x === 'number' ? (Number.isInteger(x) ? x : +x.toFixed(5)) : x).join(','));
        }
        process.stderr.write(`s1 ${fc / 1e9} GHz d_min ${dmin} ${algo} v ${v}\n`);
    }
    write('s1_sensitivity.csv', out.join('\n') + '\n');
    lines.push(`S1: ICI loss and aging loss [dB], signalFree, E6 geometry (theta1 20, theta2 -30, K 20 dB, L 100, N 8), ${trials} realisations per cell (paired: the same seeds for tau = 0 and tau > 0); mean ± SE`);
    lines.push('aging = SINR(tau=0) - SINR(tau); ICI = 10log10(SINR_tau0 / SINR_eff_tau0), SINR_eff = 1/(1/SINR + N_ICI/S);  * marks aging >= ICI');
    for (const [fc, scs] of [[5e9, 15e3], [28e9, 120e3]]) for (const algo of ['SMI', 'DL']) {
        lines.push(`\n(fc, df) = (${fc / 1e9} GHz, ${scs / 1e3} kHz), ${algo}${fc === 5e9 && scs === 15e3 ? '   <- the (fc, df) used by E6' : ''}`);
        lines.push(U.pad('d_min m', 8) + U.pad('tau ms', 7) + [25, 50, 100, 200, 300].map(v => U.rpad(`v=${v}: aging / ICI`, 24)).join(''));
        for (const dmin of [10, 30, 100, 500]) for (const tau of [1, 2, 10]) {
            lines.push(U.pad(dmin, 8) + U.pad(tau, 7) + [25, 50, 100, 200, 300].map(v => { const r = rows.find(x => x.fc === fc && x.scs === scs && x.dmin === dmin && x.tau === tau && x.v === v && x.algo === algo); return U.rpad(`${f2(r.ag)} / ${f2(r.ici)}${r.ge ? ' *' : ''}`, 24); }).join(''));
        }
    }
    const ge = rows.filter(r => r.ge);
    lines.push(`\ncells with aging loss >= ICI loss: ${ge.length} of ${rows.length}`);
    if (ge.length) for (const r of ge) lines.push(`  ${r.fc / 1e9} GHz/${r.scs / 1e3} kHz, d_min ${r.dmin} m, tau ${r.tau} ms, ${r.algo}, v = ${r.v} km/h: aging ${f2(r.ag)}±${f2(r.agse)} dB, ICI ${f2(r.ici)}±${f2(r.icise)} dB`);
    else lines.push('  none: in no cell of the grid does the aging loss reach the ICI loss.');
    write('s1_sensitivity.txt', lines.join('\n') + '\n'); console.log(lines.join('\n'));
}

if (mode === 's2') {
    const L = [];
    L.push(`S2: DL versus SMI in the E6 setting (theta1 20, theta2 -30, d_min 10 m, tau 10 ms, K 20 dB, signalFree), ${trials} realisations`);
    L.push('\n(a) SINR [dB] versus gamma_rel (DL), v = 0 (tau = 0 and tau = 10 ms: no angle change at v = 0) and v = 300; loss = SINR(v=0) - SINR(v=300)');
    L.push(U.pad('gamma_rel dB', 13) + U.rpad('SINR v=0', 14) + U.rpad('SINR v=300', 14) + U.rpad('loss', 14));
    let best = null;
    for (let g = -10; g <= 30; g += 2) {
        const a = [], b = [];
        for (let t = 0; t < trials; t++) { a.push(run({ algo: 'DL', v: 0, gammaRelDb: g }, 7000 + t).sinrDb); b.push(run({ algo: 'DL', v: 300, gammaRelDb: g }, 7000 + t).sinrDb); }
        const d = a.map((x, i) => x - b[i]);
        L.push(U.pad(g, 13) + U.rpad(cellS(a), 14) + U.rpad(cellS(b), 14) + U.rpad(cellS(d), 14));
        if (best === null || U.mean(b) > best.b) best = { g, b: U.mean(b), loss: U.mean(d) };
    }
    { const a = [], b = []; for (let t = 0; t < trials; t++) { a.push(run({ algo: 'SMI', v: 0 }, 7000 + t).sinrDb); b.push(run({ algo: 'SMI', v: 300 }, 7000 + t).sinrDb); } const d = a.map((x, i) => x - b[i]); L.push(U.pad('SMI', 13) + U.rpad(cellS(a), 14) + U.rpad(cellS(b), 14) + U.rpad(cellS(d), 14)); }
    L.push(`best gamma_rel at v = 300: ${best.g} dB (${f2(best.b)} dB, loss ${f2(best.loss)} dB); default +10 dB`);
    // (b), (c): beam pattern gain at the jammer angle at the estimation and at the application time
    L.push('\n(b)/(c) at v = 300 and, for reference, v = 0: gain G(theta) = 10log10(|w^H a(theta)|^2 / |w^H a(theta1(t_app))|^2) [dB] at the jammer angle at t_est and at t_app, ||w||, null width, INR per element');
    const Pj = lin(-BASE.sir), sg2 = lin(-BASE.snr), inr = db10(Pj / sg2);
    L.push(`jammer-to-noise ratio per element (INR) = ${f2(inr)} dB (Pj = ${f2(Pj)}, sigma^2 = ${sg2})`);
    L.push(U.pad('algo', 6) + U.pad('v', 5) + U.rpad('d_theta2 deg', 13) + U.rpad('G(t_est) dB', 13) + U.rpad('G(t_app) dB', 13) + U.rpad('G(app)-G(est)', 14) + U.rpad('||w||', 8) + U.rpad('null width deg', 15) + U.rpad('SINR dB', 14));
    for (const algo of ['SMI', 'DL']) for (const v of [0, 300]) {
        const gE = [], gA = [], wn = [], wd = [], sn = [], dd = [];
        for (let t = 0; t < trials; t++) {
            const s = run({ algo, v }, 7000 + t), w = s.weights, th1 = s.aoaT * D2R, thE = s.thJo, thA = s.aoaJ * D2R;
            const pat = th => db10(Core.vecDot(w, s.steer(th)).mag2() + 1e-300), ref = pat(th1);
            gE.push(pat(thE) - ref); gA.push(pat(thA) - ref); wn.push(Math.sqrt(w.reduce((a, c) => a + c.mag2(), 0))); sn.push(s.sinrDb); dd.push((thA - thE) / D2R);
            // null width: contiguous interval around the estimation-time jammer angle where G < -30 dB (0.05 deg grid, limited to +-10 deg)
            let lo = 0, hi = 0; for (let k = 1; k <= 200; k++) { if (pat(thE - k * 0.05 * D2R) - ref < -30) lo = k; else break; } for (let k = 1; k <= 200; k++) { if (pat(thE + k * 0.05 * D2R) - ref < -30) hi = k; else break; }
            wd.push((lo + hi) * 0.05);
        }
        L.push(U.pad(algo, 6) + U.pad(v, 5) + U.rpad(f2(U.mean(dd)), 13) + U.rpad(cellS(gE), 13) + U.rpad(cellS(gA), 13) + U.rpad(f2(U.mean(gA) - U.mean(gE)), 14) + U.rpad(f2(U.mean(wn)), 8) + U.rpad(f2(U.mean(wd)), 15) + U.rpad(cellS(sn), 14));
    }
    write('s2_dl_vs_smi.txt', L.join('\n') + '\n'); console.log(L.join('\n'));
}

if (mode === 's3') {
    const L = [];
    L.push(`S3: SMI withSignal versus speed in the E6 setting (theta1 20, theta2 -30, d_min 10 m, tau 10 ms), ${trials} realisations; mean SINR [dB] ± SE`);
    const variants = [['(i) K -> inf (no diffuse paths)', { kDb: 400 }], ['(ii) diffuse only (no LoS)', { kDb: -400 }], ['(iii) original K = 20 dB', { kDb: 20 }]];
    const vs = []; for (let v = 0; v <= 300; v += 25) vs.push(v);
    L.push(U.pad('variant', 34) + vs.map(v => U.rpad(v, 11)).join(''));
    for (const [name, o] of variants) {
        L.push(U.pad(name, 34) + vs.map(v => { const a = []; for (let t = 0; t < trials; t++) a.push(run(Object.assign({ algo: 'SMI', trainMode: 'withSignal', v }, o), 8000 + t).sinrDb); return U.rpad(cellS(a), 11); }).join(''));
    }
    L.push('\nsame with tau = 0 (no angle drift between estimation and application)');
    L.push(U.pad('variant', 34) + vs.map(v => U.rpad(v, 11)).join(''));
    for (const [name, o] of variants) L.push(U.pad(name, 34) + vs.map(v => { const a = []; for (let t = 0; t < trials; t++) a.push(run(Object.assign({ algo: 'SMI', trainMode: 'withSignal', v, latMs: 0 }, o), 8000 + t).sinrDb); return U.rpad(cellS(a), 11); }).join(''));
    L.push('\neffective rank of the target covariance in the training window: R_t = (1/L) sum_n h_n h_n^H (target part of the snapshots); share of the total power in the 1st, 1st+2nd, 1st+2nd+3rd eigenvalue [%]');
    L.push(U.pad('variant', 34) + vs.map(v => U.rpad(v, 19)).join(''));
    for (const [name, o] of variants) {
        L.push(U.pad(name, 34) + vs.map(v => {
            const sh = [[], [], []];
            for (let t = 0; t < Math.min(trials, 200); t++) {
                const s = run(Object.assign({ algo: 'SMI', trainMode: 'withSignal', v }, o), 8000 + t), N = s.N, Rt = Array.from({ length: 2 * N }, () => new Float64Array(2 * N));
                for (const { tr, ti } of s.snaps) for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) { const a = tr[i] * tr[j] + ti[i] * ti[j], b = ti[i] * tr[j] - tr[i] * ti[j]; Rt[i][j] += a; Rt[i + N][j + N] += a; Rt[i][j + N] += -b; Rt[i + N][j] += b; }
                const ev = Core.eigSymDecomp(Rt).vals.filter((x, k) => k % 2 === 0), tot = ev.reduce((a, b) => a + Math.max(b, 0), 0);
                sh[0].push(100 * ev[0] / tot); sh[1].push(100 * (ev[0] + ev[1]) / tot); sh[2].push(100 * (ev[0] + ev[1] + ev[2]) / tot);
            }
            return U.rpad(sh.map(a => U.f(U.mean(a), 1)).join('/'), 19);
        }).join(''));
    }
    write('s3_withsignal.txt', L.join('\n') + '\n'); console.log(L.join('\n'));
}
