'use strict';
/**
 * T4  H1 check: ICI of the diffuse paths.
 *     For theta_1 in {0, 30, 60, 90} deg the model's own path sampler (Sys.newRealization with M = 4096 paths, Gaussian angular
 *     offsets) estimates  E[1 - sinc^2(eps_m cos(theta_1 + delta))]  as a sample mean; this is compared with the numerical
 *     integral over the same angular distribution (Core.diffuseIciExpectation).
 *     Criterion: |sampled estimate - numerical integral| < 3 x (standard error of the sample mean); the standard error is printed.
 *     (The earlier fixed 1 % tolerance was below the 4096-path sampling error at theta_1 = 90 deg, ~2 %, so it could not be a statistical test.)
 *     The legacy isotropic formula  <1 - sinc^2(eps_m cos(alpha))>_alpha  is printed next to it, with its deviation in dB.
 */
const Core = require('../core.js');
const U = require('./_util.js');

const THETAS = [0, 30, 60, 90];
const VS = [100, 300, 500];

function legacyIsotropic(epsM) {
    let s = 0; const A = 32;
    for (let k = 0; k < A; k++) s += 1 - Math.pow(Core.sinc(epsM * Math.cos(2 * Math.PI * (k + 0.5) / A)), 2);
    return s / A;
}

module.exports = {
    id: 'T4', title: 'H1: diffuse-path ICI, sampled (M = 4096) vs numerical integral; legacy isotropic formula',
    async run({ seed = 20260404, M = 4096 } = {}) {
        const checks = [];
        Core.setSeed(seed);
        const sys = Core.createSys(); sys.rollCal(); sys.M_UNIFIED = M;
        const sigma = sys.SIGMA_ANG;
        console.log(`sigma_theta = ${(sigma * 180 / Math.PI).toFixed(1)} deg (Gaussian), M = ${M} sampled paths, f_c = ${sys.fc / 1e9} GHz, df = ${sys.scs / 1e3} kHz`);
        console.log(U.pad('theta1', 7), U.pad('v km/h', 7), U.rpad('eps_m', 7), U.rpad('integral', 11), U.rpad('sampled', 11), U.rpad('rel.err', 9), U.rpad('rel.SE', 8), U.rpad('err/SE', 6),
            U.rpad('legacy dB', 10), U.rpad('exact dB', 9), U.rpad('legacy-exact', 13), 'result');
        for (const th of THETAS) {
            const re = sys.newRealization();                          // the model's own sampler
            for (const v of VS) {
                const fm = v / 3.6 * sys.fc / sys.c, epsM = fm / sys.scs, thr = th * Math.PI / 180;
                const x = Array.from(re.delta, d => 1 - Math.pow(Core.sinc(epsM * Math.cos(thr + d)), 2));
                const est = U.mean(x), ref = Core.diffuseIciExpectation(thr, epsM, sigma);
                const relErr = Math.abs(est - ref) / ref, relSE = U.se(x) / est;
                const leg = legacyIsotropic(epsM);
                const dB = U.db10, dev = dB(leg) - dB(ref);
                const se = U.se(x), absErr = Math.abs(est - ref), pass = absErr < 3 * se;
                console.log(U.pad(th, 7), U.pad(v, 7), U.rpad(U.f(epsM, 3), 7), U.rpad(U.e(ref), 11), U.rpad(U.e(est), 11), U.rpad((100 * relErr).toFixed(2) + '%', 9), U.rpad((100 * relSE).toFixed(2) + '%', 8), U.rpad((absErr / se).toFixed(2), 6),
                    U.rpad(U.f(dB(leg), 1), 10), U.rpad(U.f(dB(ref), 1), 9), U.rpad((dev >= 0 ? '+' : '') + U.f(dev, 1) + ' dB', 13), pass ? 'PASS' : 'FAIL');
                checks.push(U.check(`T4 theta1=${th} deg, v=${v}: sampled vs integral`, `${(100 * relErr).toFixed(2)} % = ${(absErr / se).toFixed(2)} SE`, `< 3 SE (SE = ${U.e(se)}, ${(100 * relSE).toFixed(2)} %)`, pass,
                    ` legacy isotropic formula is ${dev >= 0 ? '+' : ''}${U.f(dev, 1)} dB off`));
            }
        }
        // convergence information (not a check): the same estimate with 64x more paths
        sys.M_UNIFIED = M * 64; const big = sys.newRealization();
        const th = Math.PI / 2, epsM = 300 / 3.6 * sys.fc / sys.c / sys.scs;
        const xb = Array.from(big.delta, d => 1 - Math.pow(Core.sinc(epsM * Math.cos(th + d)), 2));
        console.log(`info: theta1 = 90 deg, v = 300, M = ${M * 64}: sampled/integral - 1 = ${(100 * (U.mean(xb) / Core.diffuseIciExpectation(th, epsM, sigma) - 1)).toFixed(2)} %`);
        return { id: this.id, title: this.title, checks };
    }
};
