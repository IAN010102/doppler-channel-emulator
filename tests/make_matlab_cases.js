'use strict';
/**
 * node tests/make_matlab_cases.js  ->  data/matlab_cases.json
 * Fixed comparison cases for the MATLAB scripts (matlab/verify_cases.m, matlab/check_my_work.m, ...). No randomness at read time: every number is stored.
 * Complex values are always split into { re: [...], im: [...] }; matrices are stored column-major as arrays of columns? NO: as arrays of ROWS ([row][col]) in re / im.
 * The same function `build()` is used by T20 (two runs must be bit-identical).
 */
const fs = require('fs'), path = require('path');
const Core = require('../core.js');
const { Cplx } = Core;

const N = 8, SNR = 20, SIR = -10, KDB = 20, THETA1 = 10, THETA2 = 40, GAMMA_REL = 10, CALDEG = 5, MOD = 'QAM16';
const ALGOS = ['FOURIER', 'MMSE', 'MMSEP', 'SMI', 'DL', 'BEAMSPACE'];
const ALGO_NAME = { FOURIER: 'FOURIER', MMSE: 'MMSE-M', MMSEP: 'MMSE-P', SMI: 'SMI', DL: 'DL', BEAMSPACE: 'BEAMSPACE' };

const cvec = v => ({ re: v.map(c => c.r), im: v.map(c => c.i) });
const cmat = M => ({ re: M.map(r => r.map(c => c.r)), im: M.map(r => r.map(c => c.i)) });

function setup(over, seed) {
    Core.setSeed(seed);
    const s = Core.createSys(); s.rollCal();
    Object.assign(s, { N, aoaT: THETA1, aoaJ: THETA2, snr: SNR, sir: SIR, kDb: KDB, gammaRelDb: GAMMA_REL, calDeg: CALDEG, mod: MOD, taper: 'NONE', model: 'unified', freshRealization: true }, over);
    return s;
}

function oneCase(id, p, seed) {
    const T1 = p.th1 !== undefined ? p.th1 : THETA1, T2 = p.th2 !== undefined ? p.th2 : THETA2, SN = p.snr !== undefined ? p.snr : SNR, SI = p.sir !== undefined ? p.sir : SIR,
        KD = p.kDb !== undefined ? p.kDb : KDB, MD = p.mod || MOD, CD = p.calDeg !== undefined ? p.calDeg : CALDEG;
    // snapshots of both training modes share every random draw (same seed): X_full = withSignal data, X_sf = signalFree data
    const run = (algo, trainMode, covSource = 'sample') => { const s = setup({ L: p.L, v: p.v, latMs: p.tauMs, d_min: p.dmin, pointErrDeg: p.delta, trainMode, algo, covSource, aoaT: T1, aoaJ: T2, snr: SN, sir: SI, kDb: KD, mod: MD, calDeg: CD }, seed); s.snaps = []; s.computeMath(); return s; };
    const sW = run('SMI', 'withSignal'), sF = run('SMI', 'signalFree');
    const snap = (s, key) => Array.from({ length: N }, (_, n) => s.snaps.map(sn => new Cplx(sn[key === 'r' ? 'rr' : 'tr'][n], sn[key === 'r' ? 'ri' : 'ti'][n])));
    const toX = s => Array.from({ length: N }, (_, n) => s.snaps.map(sn => new Cplx(sn.rr[n], sn.ri[n])));
    // X = the data that contain the target, exactly as the algorithms see them: the withSignal snapshots, or (signalFree case) the signal-free snapshots plus the stored target part (the Wiener data of Sys.computeMath)
    const Xsf = toX(sF), sym = sW.snaps.map(sn => new Cplx(sn.s1r, sn.s1i));
    const Xfull = p.trainMode === 'withSignal' ? toX(sW) : Array.from({ length: N }, (_, n) => sF.snaps.map(sn => new Cplx(sn.rr[n] + sn.tr[n], sn.ri[n] + sn.ti[n])));
    const cur = run('SMI', p.trainMode);
    // true channel vector and jammer steering vector at t_app (Gamma included; unified model, same expressions as Sys.genieSinr)
    const L = p.L, Ts = (1 + Core.CONFIG.cpRatio) / cur.scs, tApp = (L - 1) * Ts + p.tauMs * 1e-3, paths = (cur.freshRealization = false, cur.unifiedPaths(T1 * Math.PI / 180, Math.pow(10, KD / 10))), gam = cur.gamma();   // freshRealization off: use the realisation of this trial
    const h = Array.from({ length: N }, () => new Cplx(0, 0));
    for (let i = 0; i < paths.P; i++) {
        const k = -2 * Math.PI * 0.5 * Math.sin(paths.th0[i]), ph = Core.trackPhase(paths.th0[i], paths.vms, cur.d_min, paths.lam, -tApp, 0);
        const b = Cplx.mul(new Cplx(paths.br[i], paths.bi[i]), new Cplx(Math.cos(ph), Math.sin(ph)));
        for (let e = 0; e < N; e++) h[e] = Cplx.add(h[e], Cplx.mul(b, new Cplx(Math.cos(k * e), Math.sin(k * e))));
    }
    const hG = h.map((c, e) => Cplx.mul(gam[e], c)), kj = -2 * Math.PI * 0.5 * Math.sin(T2 * Math.PI / 180), aJ = Array.from({ length: N }, (_, e) => Cplx.mul(gam[e], new Cplx(Math.cos(kj * e), Math.sin(kj * e))));
    const aAs = cur.steer(cur.thTo + p.delta * Math.PI / 180);
    const algos = {};
    let bs = null;
    for (const a of ALGOS) {
        const s = run(a, p.trainMode);
        algos[ALGO_NAME[a]] = { w: cvec(s.weights), sinr_dB: s.sinrDb, sinr_opt_dB: s.sinrOptDb, status: s.status, rank_R: s.rankR };
        if (a === 'BEAMSPACE') {
            const bins = s.bsBins.slice(), B = Array.from({ length: N }, (_, n) => bins.map(k => new Cplx(Math.cos(-2 * Math.PI * n * k / N) / Math.sqrt(N), Math.sin(-2 * Math.PI * n * k / N) / Math.sqrt(N))));
            bs = { bins, B: cmat(B) };
        }
    }
    // the lecture baseline case: the theoretical covariances of the lecture (R_r with the target, R_u without it) and the weights of MMSE and of the two MVDR variants built from them
    let lab = null;
    if (p.lab) {
        const mm = run('MMSE', 'withSignal', 'theory'), mb = run('SMI', 'withSignal', 'theory'), ma = run('SMI', 'signalFree', 'theory'), a1 = cur.steer(cur.thTo), a2 = cur.steer(T2 * Math.PI / 180);
        const pack = s => ({ w: cvec(s.weights), sinr_dB: s.sinrDb, sinr_opt_dB: s.sinrOptDb, evm_raw: s.evmRaw, evm_perfect_scaling: s.evm });
        lab = { note: 'Lecture notation: array_vec_est1 = a(theta_user1), a2 = a(theta_user2), Es = 1, sigma = noise power (sigma2), NR = L, rx = X. R_r = Es a1 a1^H + Es a2 a2^H + sigma I (target + interferer + noise, power of the interferer = Es), R_u = Es a2 a2^H + sigma I (without the target). In this case K -> infinity and there is no phase mismatch, so the web theory covariance equals these expressions (the random LoS phase cancels in h h^H).',
            array_vec_est1: cvec(a1), a2: cvec(a2), Es: 1, sigma: Math.pow(10, -SN / 10), NR: L, R_r_theory: cmat(mm.R_raw), R_u_theory: cmat(ma.R_raw),
            mmse_theory: pack(mm), mvdr_B_theory: pack(mb), mvdr_A_theory: pack(ma),
            formulas: { mmse: 'w = R_r^-1 r_rs, r_rs = a1 Es', mvdr_B: 'w = R_r^-1 a1 / (a1^H R_r^-1 a1)', mvdr_A: 'w = R_u^-1 a1 / (a1^H R_u^-1 a1)', relation: 'w_MMSE = Es (a1^H R_r^-1 a1) w_MVDR_B (same direction, a scalar apart)' } };
    }
    return {
        id, seed,
        params: { N, d_over_lambda: 0.5, L, sigma2: Math.pow(10, -SN / 10), SNR_dB: SN, K_dB: KD, P_j: Math.pow(10, -SI / 10), P_s: 1, theta1_hat_deg: cur.thTo * 180 / Math.PI, delta_theta_deg: p.delta, theta2_deg: T2,
            theta1_hat_lab_deg: cur.thTo * 180 / Math.PI, theta2_lab_deg: T2,
            gamma_rel_dB: GAMMA_REL, gamma_abs: Math.pow(10, GAMMA_REL / 10) * Math.pow(10, -SN / 10), model: 'unified', trainMode: p.trainMode, jamWave: cur.jamWave, smiSingular: cur.smiSingular, epsRank: cur.epsRank,
            v_kmh: p.v, tau_ms: p.tauMs, d_min_m: p.dmin, modulation: MD, sigma_phi_deg: CD },
        X: cmat(Xfull), X_sf: cmat(Xsf), s: cvec(sym), h: cvec(hG), a_J: cvec(aJ), a_assumed: cvec(aAs), gamma: cvec(gam),
        R_hat_train_web: cmat(cur.R_raw), beamspace: bs, algorithms: algos, lab
    };
}

function build(baseSeed = 20260610) {
    const cases = []; let idx = 0;
    for (const L of [4, 8, 12, 24, 100]) for (const trainMode of ['signalFree', 'withSignal']) for (const delta of [0, 3]) {
        cases.push(oneCase(`static_L${L}_${trainMode}_d${delta}`, { L, trainMode, delta, v: 0, tauMs: 0, dmin: 30 }, baseSeed + idx++));
    }
    for (const trainMode of ['signalFree', 'withSignal']) cases.push(oneCase(`dynamic_v300_tau2_L100_${trainMode}_d0`, { L: 100, trainMode, delta: 0, v: 300, tauMs: 2, dmin: 30 }, baseSeed + idx++));
    // the lecture baseline (experiment E0): N = 8, theta1 = -20, theta2 = 30, SNR 30 dB, Es = 1, interferer power = Es, L = 1000, K -> infinity, no phase mismatch
    cases.push(oneCase('E0_lab_baseline', { L: 1000, trainMode: 'withSignal', delta: 0, v: 0, tauMs: 0, dmin: 30, th1: -20, th2: 30, snr: 30, sir: 0, kDb: 400, mod: 'QPSK', calDeg: 0, lab: true }, baseSeed + idx++));
    return {
        schema_version: 1,
        description: 'Fixed comparison cases exported by the web simulator (tests/make_matlab_cases.js) for the MATLAB scripts. Nothing here is random at read time.',
        lab_convention: 'The simulator and the lecture use the SAME angle convention: a_n(theta) = exp(-j 2 pi n (d/lambda) sin(theta)), n = 0..N-1, theta from the array broadside (= the heading direction of the vehicle), positive sign as in the lecture. Conversion theta_lab = theta_sim (identity, including the sign). params.theta1_hat_lab_deg / theta2_lab_deg repeat the angles in the lecture convention.',
        complex_format: 'every complex vector/matrix is { re, im } with two real arrays of the same shape; vectors are plain arrays of length N (or L); matrices are arrays of ROWS: M.re[row][col]',
        matrices: 'X, X_sf: N x L (row = antenna, column = snapshot); B: N x K; R_hat_train_web: N x N',
        definitions: {
            R_hat: 'R_hat = (1/L) * X * X^H for the data the algorithm trains on: SMI, DL, BEAMSPACE use X_sf when trainMode = signalFree and X when withSignal; MMSE-M and MMSE-P always use X (the target is part of the Wiener data)',
            snapshots: 'X = h-dependent signal + interferer + noise: x_n = sum_i beta_i Gamma a(theta_i(t_n)) e^{j phi_i(t_n)} s_n + sqrt(P_j) Gamma a(theta_2(t_n)) e^{j phi_2(t_n)} j_n + noise_n; X_sf = the same draws without the target term',
            training_symbols: 's = known target symbols s_n (unit average power), used only by MMSE-P',
            h: 'h = true channel vector at the application time t_app = (L-1) T_snap + tau (Gamma included)',
            a_J: 'a_J = Gamma * a(theta_2) at t_app (jammer, unit modulus entries before Gamma); the interference power is P_j |w^H a_J|^2',
            a_assumed: 'a(theta1_hat + delta_theta), a_n(theta) = exp(-j 2 pi n (d/lambda) sin(theta)), n = 0..N-1 (the nominal steering vector the weights are designed for)',
            sinr: 'SINR = |w^H h|^2 / ( P_j |w^H a_J|^2 + sigma2 * ||w||^2 ),  SINR_dB = 10 log10(SINR); invariant to a scalar factor of w. P_s = 1 is contained in h.',
            sinr_opt: 'SINR_opt = h^H R_in^-1 h with R_in = P_j a_J a_J^H + sigma2 I  (genie upper bound)'
        },
        weight_normalisation: {
            FOURIER: 'w = a_assumed / N (taper NONE): w^H a_assumed = 1',
            'MMSE-M': 'w = R_hat^-1 a_assumed (P_s = 1), NOT normalised (Wiener scale); pseudo-inverse when rank(R_hat) < N',
            'MMSE-P': 'w = R_hat^-1 r_xd, r_xd = (1/L) sum_n x_n conj(s_n), NOT normalised; pseudo-inverse when rank(R_hat) < N',
            SMI: 'w = R_hat^-1 a_assumed / (a_assumed^H R_hat^-1 a_assumed); pseudo-inverse when rank(R_hat) < N (smiSingular = pinv): w = R^+ a / (a^H R^+ a)',
            DL: 'w = (R_hat + gamma I)^-1 a_assumed / (a_assumed^H (R_hat + gamma I)^-1 a_assumed), gamma = gamma_rel * sigma2 (gamma_abs)',
            BEAMSPACE: 'B = [d_k1 d_k2 d_k3], d_k[n] = exp(-j 2 pi n k / N) / sqrt(N) (beamspace.B, beamspace.bins); R_B = B^H R_hat B, a_B = B^H a_assumed; w_B = R_B^-1 a_B / (a_B^H R_B^-1 a_B); w = B w_B'
        },
        pseudo_inverse: 'rank = number of eigenvalues of the Hermitian R_hat above epsRank * lambda_max; R^+ = sum over those eigenpairs of (1/lambda) u u^H',
        cases
    };
}

if (require.main === module) {
    const out = path.join(__dirname, '..', 'data', 'matlab_cases.json');
    fs.mkdirSync(path.dirname(out), { recursive: true });
    const json = build();
    fs.writeFileSync(out, JSON.stringify(json));
    console.log(`written ${out}: ${json.cases.length} cases, ${(fs.statSync(out).size / 1024).toFixed(0)} kB`);
}
module.exports = { build };
