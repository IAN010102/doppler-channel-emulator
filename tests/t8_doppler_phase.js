'use strict';
/**
 * T8  Integrated Doppler phase (Commit 9): Core.trackPhase(th, v, d_min, lambda, ta, tb) = 2 pi [R(ta) - R(tb)] / lambda,  R = d_min/|sin(theta)|,
 *     theta(t) = Core.trackAngle(th, v, d_min, t - t_ref).
 *   T8a  closed form vs numerical integral of 2 pi f_m cos(theta(t)) dt (>= 1e5 points, composite Simpson) for random (theta, v, d_min, ta, tb). Error < 1e-9 rad.
 *   T8b  explicit coordinates: receiver at (v t, 0), source at (x_s, d_min) (or (x_s, -d_min)); the distance difference is computed directly
 *        sqrt((x_s - v ta)^2 + d_min^2) - sqrt((x_s - v tb)^2 + d_min^2) and compared with the phase (and its sign). Error < 1e-9 rad; f_d > 0 while approaching.
 *   T8c  d_min -> very large: phase -> 2 pi f_d t with f_d = f_m cos(theta) (relative error < 1e-6).
 */
const Core = require('../core.js');
const U = require('./_util.js');

module.exports = {
    id: 'T8', title: 'integrated Doppler phase: closed form vs numerical integral, explicit coordinates, d_min -> infinity',
    async run({ seed = 888 } = {}) {
        const checks = [], rnd = Core.mulberry32(seed), D = Math.PI / 180, lam = 3e8 / 5e9;
        const fm = v => v * 5e9 / 3e8;
        // ---- T8a
        let maxA = 0; const cases = 40;
        for (let k = 0; k < cases; k++) {
            const th = (rnd() * 340 - 170) * D, v = rnd() * 140 + 1, dmin = Math.pow(10, 0.5 + 2.5 * rnd()), ta = -rnd() * 0.02, tb = rnd() * 0.02 - 0.01;
            const f = t => 2 * Math.PI * fm(v) * Math.cos(Core.trackAngle(th, v, dmin, t));
            const n = 100000, h = (tb - ta) / n; let s = f(ta) + f(tb);
            for (let i = 1; i < n; i++) s += f(ta + i * h) * (i % 2 ? 4 : 2);
            const num = s * h / 3, cf = Core.trackPhase(th, v, dmin, lam, ta, tb);
            maxA = Math.max(maxA, Math.abs(num - cf));
        }
        console.log(`T8a  ${cases} random (theta, v, d_min, ta, tb), 1e5-point Simpson: max |closed form - numerical integral| = ${U.e(maxA)} rad`);
        checks.push(U.check('T8a closed-form phase vs numerical integral (1e5 points)', U.e(maxA) + ' rad', '1e-9 rad', maxA < 1e-9));

        // ---- T8b
        let maxB = 0, signOk = true;
        for (const deg of [-100, -60, -20, 5, 30, 45, 80, 100, 150]) for (const dmin of [5, 30, 500]) {
            const th = deg * D, v = 300 / 3.6, ds = (Math.sin(th) >= 0 ? 1 : -1) * dmin, xs = ds / Math.tan(th);
            for (const [ta, tb] of [[-0.007, 0], [-0.007, 0.01], [0, 0.0071], [-0.017, -0.01]]) {
                // reference time of th is t_ref = 0 (the angle th holds at t = 0 in these coordinates)
                const R = t => Math.hypot(xs - v * t, ds), direct = 2 * Math.PI * (R(ta) - R(tb)) / lam;
                maxB = Math.max(maxB, Math.abs(Core.trackPhase(th, v, dmin, lam, ta, tb) - direct));
                if (tb > ta && Math.cos(th) > 0.05 && xs - v * ta > 0 && xs - v * tb > 0 && !(Core.trackPhase(th, v, dmin, lam, ta, tb) > 0)) signOk = false;   // approaching: phase grows
            }
        }
        // instantaneous f_d > 0 while approaching: (phi(t+h) - phi(t))/(2 pi h) -> f_m cos(theta) > 0 for theta in (-90, 90)
        let fdOk = true;
        for (const deg of [-80, -45, 10, 45, 80]) { const th = deg * D, v = 300 / 3.6, h = 1e-7; const fd = Core.trackPhase(th, v, 30, lam, 0, h) / (2 * Math.PI * h); if (!(fd > 0) || Math.abs(fd - fm(v) * Math.cos(th)) > 1e-3 * fm(v)) fdOk = false; }
        console.log(`T8b  explicit coordinates: max |phase - 2 pi (R(ta) - R(tb))/lambda| = ${U.e(maxB)} rad; phase grows while approaching: ${signOk}; f_d > 0 and = f_m cos(theta) for |theta| < 90: ${fdOk}`);
        checks.push(U.check('T8b phase equals the explicit-coordinate range difference', U.e(maxB) + ' rad', '1e-9 rad', maxB < 1e-9));
        checks.push(U.check('T8b sign: phase grows (f_d > 0) while approaching', `${signOk && fdOk}`, 'true', signOk && fdOk));

        // ---- T8c   ("large" is d_min >= 1e8 m: the true deviation from 2 pi f_d t is O(v t / d_min), ~1.6e-6 at 1e6 m, so 1e-6 is only meaningful
        //             beyond that; the error per d_min is printed, and its 1/d_min scaling is checked)
        const errAt = dmin => { let m = 0; for (const deg of [-70, -30, 10, 45, 80]) { const th = deg * D, v = 300 / 3.6, t = 0.007, ref = 2 * Math.PI * fm(v) * Math.cos(th) * t; m = Math.max(m, Math.abs(Core.trackPhase(th, v, dmin, lam, 0, t) / ref - 1)); } return m; };
        const row = [1e6, 1e8, 1e9, 1e12, 1e30].map(d => [d, errAt(d)]);
        console.log('T8c  max relative |phase / (2 pi f_d t) - 1| vs d_min: ' + row.map(([d, e]) => `${d.toExponential(0)} m: ${U.e(e)}`).join(',  '));
        const maxC = Math.max(...row.filter(([d]) => d >= 1e8).map(([, e]) => e)), scal = row[0][1] / row[1][1];
        checks.push(U.check('T8c d_min >= 1e8 m: phase -> 2 pi f_d t', U.e(maxC), '< 1e-6', maxC < 1e-6));
        checks.push(U.check('T8c deviation scales as 1/d_min (1e6 vs 1e8 m)', `ratio ${scal.toFixed(1)}`, '100 ± 10', Math.abs(scal - 100) < 10));
        return { id: this.id, title: this.title, checks };
    }
};
