'use strict';
/**
 * T7  Geometry (Commit 7): Core.trackAngle / Core.trackRate against explicit coordinates.
 *     A ground point at (x_s, d_s), a receiver at (v t, 0) moving along +x; theta = atan2(d_s, x_s - v t) is the angle to the heading.
 *     T7a  trackAngle(theta_0, v, d_min, t) equals the coordinate angle for theta_0 in +-{10, 45, 80, 100} deg, several t (including t < 0). Error < 1e-9 rad.
 *     T7b  trackRate equals the numerical time derivative of the coordinate angle. Relative error < 1e-6.
 *     T7c  sign consistency: for theta in (0, 90) the angle grows and f_d = f_m cos(theta) falls while the receiver approaches; for theta < 0 the angle
 *          decreases (more negative) and f_d falls as well; a source on the track (theta = 0) does not drift.
 *     T7d  window diagnostics of the model: with d_min -> infinity the drift read-outs are 0, with d_min = 30 m they are > 0 (unified only).
 */
const Core = require('../core.js');
const U = require('./_util.js');

module.exports = {
    id: 'T7', title: 'straight-track geometry: trackAngle / trackRate vs explicit coordinates, sign consistency, drift read-outs',
    async run() {
        const checks = [], D = Math.PI / 180, v = 300 / 3.6, dmin = 40;
        let maxAng = 0, maxRate = 0;
        for (const deg of [-100, -80, -45, -10, 10, 45, 80, 100]) {
            const th0 = deg * D, ds = Math.sign(Math.sin(th0)) * dmin, xs = ds / Math.tan(th0);
            const ang = t => Math.atan2(ds, xs - v * t);
            for (const t of [-0.5, -0.01, 0, 0.007, 0.3, 1]) maxAng = Math.max(maxAng, Math.abs(Core.trackAngle(th0, v, dmin, t) - ang(t)));
            const h = 1e-6, num = (ang(h) - ang(-h)) / (2 * h);
            maxRate = Math.max(maxRate, Math.abs(num - Core.trackRate(th0, v, dmin)) / Math.abs(num));
        }
        console.log(`T7a  max |trackAngle - coordinate angle| = ${U.e(maxAng)} rad;   T7b  max relative error of trackRate = ${U.e(maxRate)}`);
        checks.push(U.check('T7a trackAngle equals the coordinate geometry', U.e(maxAng) + ' rad', '1e-9 rad', maxAng < 1e-9));
        checks.push(U.check('T7b trackRate equals d(theta)/dt of the coordinate geometry', U.e(maxRate), '1e-6', maxRate < 1e-6));

        const fm = v * 5e9 / 3e8, up = Core.trackAngle(30 * D, v, dmin, 0.01), dn = Core.trackAngle(-30 * D, v, dmin, 0.01);
        const ok = up > 30 * D && fm * Math.cos(up) < fm * Math.cos(30 * D) && dn < -30 * D && fm * Math.cos(dn) < fm * Math.cos(-30 * D) && Core.trackAngle(0, v, dmin, 0.5) === 0;
        console.log(`T7c  theta = +30 -> ${(up / D).toFixed(4)}, theta = -30 -> ${(dn / D).toFixed(4)} after 10 ms;  f_d falls in both cases; theta = 0 stays at 0`);
        checks.push(U.check('T7c drift direction and Doppler sign are consistent', ok ? 'consistent' : 'inconsistent', 'consistent', ok));

        const sys = Core.createSys(); Core.setSeed(3); sys.rollCal();
        const o = { model: 'unified', freshRealization: true, kDb: 20, aoaT: 45, aoaJ: -30, v: 300, L: 100, latMs: 0, algo: 'SMI' };
        Object.assign(sys, o, { d_min: 1e30 }); sys.snaps = []; sys.computeMath(); const a0 = sys.angDrift, p0 = sys.dopPhaseErr;
        Object.assign(sys, o, { d_min: 30 }); sys.snaps = []; sys.computeMath(); const a1 = sys.angDrift, p1 = sys.dopPhaseErr;
        console.log(`T7d  window drift read-outs: d_min = inf -> ${U.e(a0)} deg / ${U.e(p0)} rad;  d_min = 30 m -> ${a1.toFixed(3)} deg / ${p1.toFixed(3)} rad`);
        checks.push(U.check('T7d drift read-outs vanish without drift and are positive with it', `${U.e(a0)} / ${a1.toFixed(3)} deg`, '0 / > 0', a0 < 1e-12 && a1 > 0 && p1 > 0));
        return { id: this.id, title: this.title, checks };
    }
};
