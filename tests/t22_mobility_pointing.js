'use strict';
/**
 * T22  Pointing mode 'mobility' and the effective pointing deviation delta_theta_eff (Commit 22c).
 *      delta_theta_eff = theta1(t_app) - theta1(t_est), t_est = t_app - tau (the time the nominal angle of the weights belongs to). In 'mobility' mode the delta_theta slider is ignored.
 *   T22a  v = 0: delta_theta_eff = 0 (target and jammer) and mobility is bit-identical to manual delta_theta = 0, for all six algorithms and both training modes (with the slider
 *         set to 3 deg, which must be ignored). v = 300 km/h: mobility (slider 3 deg) is bit-identical to manual delta_theta = 0 (the drift is already part of the nominal angle).
 *   T22b  delta_theta_eff against an independent computation from explicit coordinates (receiver at x = v t, ground point at (x_s, +-d_min), theta = atan2(y, x_s - x_r) measured from the
 *         heading), for target and jammer, over a grid of angles, speeds, distances and latencies: < 1e-9 deg.
 *   T22c  (informational) mean SINR of mobility mode at speed v against the manual curve of T17d (v = 0, tau = 0, delta_theta = -delta_theta_eff(v) so that the weights miss the target by
 *         the same angle), SMI / DL, signalFree / withSignal, same seeds; the difference is reported (sources: the jammer drift, the Doppler spread of the diffuse paths,
 *         the aging of the covariance inside the window). The last column repeats the mobility run with tau = 0 (no angle deviation, only the Doppler inside the window).
 */
const Core = require('../core.js');
const U = require('./_util.js');
const D2R = Math.PI / 180;

const BASE = { N: 8, L: 100, aoaT: 20, aoaJ: -30, snr: 20, sir: -10, v: 0, latMs: 0, calDeg: 0, taper: 'NONE', kDb: 20, mod: 'QPSK', gammaRelDb: 10, d_min: 10, model: 'unified' };
function mk(over, seed) {
    Core.setSeed(seed);
    const s = Core.createSys(); s.calZ = new Array(16).fill(0); Object.assign(s, BASE, { freshRealization: true, trainMode: 'signalFree' }, over);
    s.snaps = []; s.snapKey = ''; s.computeMath(); return s;
}
// angle (rad) of a ground point from explicit coordinates: receiver on the x axis at x_r = v t, point at (x0, y) with y = sign(sin th0) d_min and x0 - x_r(t_app) = d_min cot|th0|... (all at t_app)
function angleAt(th0, vms, dmin, dt) {          // angle at t_app + dt
    const s = Math.sin(th0); if (s === 0) return th0;
    const y = (s > 0 ? 1 : -1) * dmin, xrel0 = y / Math.tan(th0);           // x_s - x_r at t_app  (tan(th0) = y / xrel0)
    return Math.atan2(y, xrel0 - vms * dt);
}
const same = (a, b) => a === b || (Number.isNaN(a) && Number.isNaN(b));

module.exports = {
    id: 'T22', title: "pointing mode 'mobility': delta_theta_eff = 0 at v = 0 and identical to manual 0; independent geometry; SINR vs delta_theta_eff (info)",
    async run({ trials = 300, seed = 2222 } = {}) {
        const checks = [], info = [];
        // ---- T22a
        const algos = ['FOURIER', 'MMSE', 'MMSEP', 'SMI', 'DL', 'BEAMSPACE'];
        for (const v of [0, 300]) {
            let worst = 0, n = 0, effMax = 0;
            for (const tm of ['signalFree', 'withSignal']) for (const algo of algos) for (let t = 0; t < 5; t++) {
                const a = mk({ algo, trainMode: tm, v, latMs: v ? 10 : 0, pointingMode: 'mobility', pointErrDeg: 3 }, seed + t), b = mk({ algo, trainMode: tm, v, latMs: v ? 10 : 0, pointingMode: 'manual', pointErrDeg: 0 }, seed + t);
                for (const f of ['sinrDb', 'evm', 'sinrOptDb', 'evmRaw', 'thTo']) { if (!same(a[f], b[f])) worst = Math.max(worst, Math.abs(a[f] - b[f]) || 1); }
                if (!(a.R_hat && b.R_hat)) { /* R_hat not required */ }
                effMax = Math.max(effMax, Math.abs(a.dEffT), Math.abs(a.dEffJ)); n++;
            }
            console.log(`T22a  v = ${v}: ${n} cases mobility(slider 3 deg) vs manual(0): max |difference| = ${U.e(worst)}; max |delta_theta_eff| = ${U.f(effMax, 4)} deg`);
            checks.push(U.check(`T22a v = ${v}: mobility (slider 3 deg) is bit-identical to manual delta_theta = 0`, U.e(worst), '0 (bit-identical)', worst === 0));
            if (v === 0) checks.push(U.check('T22a v = 0: delta_theta_eff = 0 (target and jammer)', U.e(effMax), '0 (exact)', effMax === 0));
            else checks.push(U.check('T22a v = 300, tau = 10 ms, d_min = 10 m: delta_theta_eff is not zero', U.f(effMax, 3) + ' deg', '> 0', effMax > 0));
        }
        // ---- T22b
        let worstB = 0, nB = 0, worstWhere = '';
        for (const aoaT of [-60, -20, 5, 20, 45, 80]) for (const aoaJ of [-30, 40]) for (const v of [30, 100, 300]) for (const dmin of [5, 10, 30, 200]) for (const lat of [1, 10, 50]) {
            const s = mk({ algo: 'SMI', aoaT, aoaJ, v, d_min: dmin, latMs: lat, L: 8 }, seed), vms = v / 3.6, tau = lat * 1e-3;
            const refT = (angleAt(aoaT * D2R, vms, dmin, 0) - angleAt(aoaT * D2R, vms, dmin, -tau)) / D2R, refJ = (angleAt(aoaJ * D2R, vms, dmin, 0) - angleAt(aoaJ * D2R, vms, dmin, -tau)) / D2R;
            const e = Math.max(Math.abs(s.dEffT - refT), Math.abs(s.dEffJ - refJ)); nB++;
            if (e > worstB) { worstB = e; worstWhere = `aoaT ${aoaT}, aoaJ ${aoaJ}, v ${v}, d_min ${dmin}, tau ${lat} ms`; }
        }
        console.log(`T22b  ${nB} settings: max |delta_theta_eff - independent coordinates| = ${U.e(worstB)} deg   (worst: ${worstWhere})`);
        checks.push(U.check('T22b delta_theta_eff equals the explicit-coordinate computation (target and jammer)', U.e(worstB) + ' deg', '< 1e-9 deg', worstB < 1e-9));
        // a target on the track (theta = 0) does not drift
        const s0 = mk({ algo: 'SMI', aoaT: 0, v: 300, latMs: 10 }, seed);
        checks.push(U.check('T22b target on the track (theta1 = 0): no deviation', U.e(Math.abs(s0.dEffT)) + ' deg', '0', s0.dEffT === 0));

        // ---- T22c (informational)
        console.log(`\nT22c  mobility mode vs manual curve (v = 0, tau = 0, delta_theta = -delta_theta_eff), ${trials} trials, tau = 10 ms, d_min = 10 m, theta1 = 20, theta2 = -30; mean SINR dB`);
        console.log(U.pad('algo / training', 20), U.rpad('v km/h', 7), U.rpad('d_eff deg', 10), U.rpad('d_eff/bw3', 10), U.rpad('mobility', 14), U.rpad('manual(-d_eff)', 16), U.rpad('diff (mob-man)', 15), U.rpad('SE', 6), U.rpad('mob, tau=0', 12));
        for (const [algo, tm] of [['SMI', 'signalFree'], ['SMI', 'withSignal'], ['DL', 'signalFree'], ['DL', 'withSignal']]) {
            const row = [];
            for (const v of [100, 200, 300]) {
                const mob = [], man = [], mob0 = []; let dEff = 0, bw = 0;
                for (let t = 0; t < trials; t++) {
                    const a = mk({ algo, trainMode: tm, v, latMs: 10, pointingMode: 'mobility' }, seed + 100 + t); mob.push(a.sinrDb); dEff = a.dEffT; bw = a.bw3;
                    mob0.push(mk({ algo, trainMode: tm, v, latMs: 0, pointingMode: 'mobility' }, seed + 100 + t).sinrDb);
                    man.push(mk({ algo, trainMode: tm, v: 0, latMs: 0, pointingMode: 'manual', pointErrDeg: -dEff }, seed + 100 + t).sinrDb);
                }
                const d = mob.map((x, i) => x - man[i]);
                console.log(U.pad(`${algo} ${tm}`, 20), U.rpad(v, 7), U.rpad(U.f(dEff, 3), 10), U.rpad(U.f(Math.abs(dEff) / bw, 3), 10), U.rpad(U.f(U.mean(mob), 2), 14), U.rpad(U.f(U.mean(man), 2), 16), U.rpad(U.f(U.mean(d), 2), 15), U.rpad(U.f(U.se(d), 2), 6), U.rpad(U.f(U.mean(mob0), 2), 12));
                row.push(`${v}: ${U.f(U.mean(mob), 1)} vs ${U.f(U.mean(man), 1)} (diff ${U.f(U.mean(d), 1)} ± ${U.f(U.se(d), 1)}; mobility with tau = 0: ${U.f(U.mean(mob0), 1)})`);
            }
            info.push({ name: `T22c ${algo} ${tm}`, value: row.join('; ') });
        }
        return { id: this.id, title: this.title, checks, info };
    }
};
