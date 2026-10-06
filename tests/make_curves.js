'use strict';
/**
 * node tests/make_curves.js [--quick]  ->  data/curve_*.csv  (+ data/README_data.md is hand-written and checked by T20c)
 * Report curves for the MATLAB plotting script. Long ("tidy") CSV format, one row per point; every file starts with '#' comment lines (schema version, parameters, seed, trials), then the
 * column header, then the data. Every point has a mean and a standard error of that mean. Trials of a point use the seeds seed0 + t (the same seeds for every point and algorithm: paired).
 * Defaults of the simulator unless stated: unified model, N = 8, L = 100, K = 20 dB, SNR = 20 dB, SIR = -10 dB, theta1 = 0, theta2 = 40, 16-QAM, signalFree training, Gaussian jammer, d_min = 30 m.
 * --quick uses 20 trials per point (development only; the files in data/ are made without it).
 */
const fs = require('fs'), path = require('path');
const Core = require('../core.js');
const quick = process.argv.includes('--quick');
const OUT = path.join(__dirname, '..', 'data');
const SCHEMA = 1;
const ALG = { FOURIER: 'FOURIER', MMSE: 'MMSE-M', MMSEP: 'MMSE-P', SMI: 'SMI', DL: 'DL', BEAMSPACE: 'BEAMSPACE' };

function sysFor(over, seed) {
    Core.setSeed(seed);
    const s = Core.createSys(); s.calZ = new Array(16).fill(0);
    Object.assign(s, { model: 'unified', N: 8, L: 100, aoaT: 0, aoaJ: 40, snr: 20, sir: -10, kDb: 20, v: 0, latMs: 0, calDeg: 0, taper: 'NONE', mod: 'QAM16', d_min: 30, freshRealization: true }, over);
    return s;
}
function stat(a) { const n = a.length, m = a.reduce((x, y) => x + y, 0) / n; let v = 0; for (const x of a) v += (x - m) * (x - m); return { mean: m, se: n > 1 ? Math.sqrt(v / (n - 1) / n) : 0, n }; }
const f = x => Number.isFinite(x) ? Number(x.toPrecision(8)).toString() : 'NaN';
function write(name, header, cols, rows) {
    fs.mkdirSync(OUT, { recursive: true });
    const txt = header.map(l => '# ' + l).join('\n') + '\n' + cols.join(',') + '\n' + rows.map(r => r.join(',')).join('\n') + '\n';
    fs.writeFileSync(path.join(OUT, name), txt); console.log('written', name, rows.length, 'rows');
}

// ------------------------------------------------------------------------------------------------ EVM vs velocity
function evmVsVelocity() {
    const trials = quick ? 20 : 500, rows = [], seed0 = 910000;
    for (const [fc, scs] of [[5e9, 15e3], [28e9, 120e3]]) for (let v = 0; v <= 300; v += 10) for (const k of Object.keys(ALG)) {
        const evm = [], ici = [];
        for (let t = 0; t < trials; t++) { const s = sysFor({ fc, scs, v, algo: k }, seed0 + t); s.snaps = []; s.computeMath(); evm.push(100 * s.evm); ici.push(100 * Math.sqrt(s.nuIciFloor)); }
        const a = stat(evm), b = stat(ici);
        rows.push([fc / 1e9, scs / 1e3, v, ALG[k], f(a.mean), f(a.se), f(b.mean), f(b.se), a.n]);
    }
    write('curve_evm_vs_velocity.csv', [`csv_schema_version=${SCHEMA}`, 'file: EVM and ICI floor versus speed (analytic EVM = sqrt(1/SINR + N_ICI/S), mean over channel realisations of the EVM itself, not rms)',
        `generator: tests/make_curves.js, seed0=${seed0} (trial t uses seed0+t), trials per point=${trials}, unified model, trainMode=signalFree, jamWave=gaussian, smiSingular=pinv`,
        'parameters: N=8, L=100, K=20 dB, SNR=20 dB, SIR=-10 dB, theta1=0, theta2=40 deg, 16-QAM, d_min=30 m, tau=0, delta_theta=0, gamma_rel=+10 dB',
        'ici_floor_evm_pct = sqrt of the ICI floor (no spatial filtering) in percent; it is the same for every algorithm'],
    ['fc_GHz', 'delta_f_kHz', 'velocity_kmh', 'algorithm', 'evm_mean_pct', 'evm_se_pct', 'ici_floor_evm_pct', 'ici_floor_se_pct', 'n_trials'], rows);
}

// ------------------------------------------------------------------------------------------------ SINR vs L
function sinrVsL() {
    const trials = quick ? 20 : 1000, rows = [], seed0 = 920000;
    for (const tm of ['signalFree', 'withSignal']) for (let L = 2; L <= 48; L++) {
        const opt = [];
        for (const k of ['SMI', 'DL', 'BEAMSPACE', 'MMSEP']) {
            const a = [];
            for (let t = 0; t < trials; t++) { const s = sysFor({ L, trainMode: tm, algo: k }, seed0 + t); s.snaps = []; s.computeMath(); a.push(s.sinrDb); if (k === 'SMI') opt.push(s.sinrOptDb); }
            const st = stat(a); rows.push([tm, L, ALG[k], f(st.mean), f(st.se), st.n]);
        }
        const o = stat(opt); rows.push([tm, L, 'SINR_opt', f(o.mean), f(o.se), o.n]);
    }
    write('curve_sinr_vs_L.csv', [`csv_schema_version=${SCHEMA}`, 'file: mean SINR versus the number of snapshots L (v=0, tau=0, delta_theta=0); SINR_opt = genie upper bound (true interference-plus-noise covariance, true channel)',
        `generator: tests/make_curves.js, seed0=${seed0} (trial t uses seed0+t), trials per point=${trials}, unified model, jamWave=gaussian, smiSingular=pinv (L < N: pseudo-inverse)`,
        'parameters: N=8, K=20 dB, SNR=20 dB, SIR=-10 dB, theta1=0, theta2=40 deg, 16-QAM, gamma_rel=+10 dB (DL)', 'MMSE-P ignores the trainMode (its data always contain the target)'],
    ['train_mode', 'L', 'algorithm', 'sinr_mean_dB', 'sinr_se_dB', 'n_trials'], rows);
}

// ------------------------------------------------------------------------------------------------ SINR vs delta_theta
function sinrVsDelta() {
    const trials = quick ? 20 : 1000, rows = [], seed0 = 930000;
    for (const L of [100, 12]) for (const tm of ['signalFree', 'withSignal']) for (let d = -5; d <= 5.0001; d += 0.5) {
        const dd = Math.round(d * 10) / 10, opt = [];
        for (const k of ['SMI', 'DL', 'BEAMSPACE', 'MMSE', 'MMSEP']) {
            const a = [];
            for (let t = 0; t < trials; t++) { const s = sysFor({ L, trainMode: tm, pointErrDeg: dd, algo: k }, seed0 + t); s.snaps = []; s.computeMath(); a.push(s.sinrDb); if (k === 'SMI') opt.push(s.sinrOptDb); }
            const st = stat(a); rows.push([tm, L, dd, ALG[k], f(st.mean), f(st.se), st.n]);
        }
        const o = stat(opt); rows.push([tm, L, dd, 'SINR_opt', f(o.mean), f(o.se), o.n]);
    }
    write('curve_sinr_vs_delta_theta.csv', [`csv_schema_version=${SCHEMA}`, 'file: mean SINR versus the pointing error delta_theta (weights designed for a(theta1 + delta_theta); SINR evaluated on the true channel); MMSE-P does not use the nominal direction',
        `generator: tests/make_curves.js, seed0=${seed0} (trial t uses seed0+t), trials per point=${trials}, unified model, jamWave=gaussian, smiSingular=pinv`,
        'parameters: N=8, K=20 dB, SNR=20 dB, SIR=-10 dB, theta1=0, theta2=40 deg, 16-QAM, v=0, tau=0, gamma_rel=+10 dB (DL)'],
    ['train_mode', 'L', 'delta_theta_deg', 'algorithm', 'sinr_mean_dB', 'sinr_se_dB', 'n_trials'], rows);
}

// ------------------------------------------------------------------------------------------------ DOA spectrum
function doaSpectrum() {
    const trials = quick ? 20 : 1000, seed0 = 940000, T1 = -10, T2 = 40, rows = [];
    const grid = Array.from({ length: 361 }, (_, i) => -90 + 0.5 * i), cap = grid.map(() => []), mus = grid.map(() => []);
    let single = null;
    for (let t = 0; t < trials; t++) {
        const s = sysFor({ aoaT: T1, aoaJ: T2, algo: 'SMI', trainMode: 'withSignal', L: 100 }, seed0 + t); s.snaps = []; s.computeMath(); const d = s.computeDoa();
        d.capon.forEach((v, i) => cap[i].push(v)); d.music.forEach((v, i) => mus[i].push(v)); if (t === 0) single = d;
    }
    grid.forEach((a, i) => { const c = stat(cap[i]), m = stat(mus[i]); rows.push([f(a), f(single.capon[i]), f(single.music[i]), f(c.mean), f(c.se), f(m.mean), f(m.se), trials]); });
    write('curve_doa_spectrum.csv', [`csv_schema_version=${SCHEMA}`, 'file: Capon (MVDR) and MUSIC spatial spectra of the sample covariance of the training window, dB relative to each maximum',
        `generator: tests/make_curves.js, seed0=${seed0} (trial t uses seed0+t), trials=${trials}; the *_single columns are trial 0, the *_mean / *_se columns average the dB values of all trials`,
        `scenario: withSignal training (the spectrum contains target and jammer), theta1=${T1} deg, theta2=${T2} deg, N=8, L=100, K=20 dB, SNR=20 dB, SIR=-10 dB, unified model, MUSIC assumes 2 sources`],
    ['angle_deg', 'capon_single_dB', 'music_single_dB', 'capon_mean_dB', 'capon_se_dB', 'music_mean_dB', 'music_se_dB', 'n_trials'], rows);
}

if (require.main === module) {
    const only = process.argv.filter(a => /^--only=/.test(a)).map(a => a.slice(7));
    const jobs = { evm: evmVsVelocity, L: sinrVsL, delta: sinrVsDelta, doa: doaSpectrum };
    for (const k of Object.keys(jobs)) if (!only.length || only.includes(k)) { const t0 = Date.now(); jobs[k](); console.log(`  (${((Date.now() - t0) / 1000).toFixed(0)} s)`); }
}
module.exports = { evmVsVelocity, sinrVsL, sinrVsDelta, doaSpectrum };
