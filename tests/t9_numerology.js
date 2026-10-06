'use strict';
/**
 * T9  Carrier frequency and numerology (Commit 11).
 *   T9a  eps = v * fc * cos(theta) / (c * df) (hand formula) equals the implementation's eps0 (and eps_m = v fc / (c df)), several (v, fc, df, theta), both models. Relative error < 1e-12.
 *   T9b  fc and df both multiplied by 2: eps and the ICI floor are unchanged (< 1e-12).
 *   T9c  df doubled: T_snap halves, so rho = L T_snap B_D halves (B_D does not depend on df); the aging ratio tau B_D does not depend on df and scales with fc.
 *   T9d  fc = 28 GHz, df = 120 kHz, v = 300 km/h: f_m, eps and the ICI floor next to hand values:
 *        f_m = (300/3.6 m/s)(28e9)/(3e8) = 7777.78 Hz;  eps_m = f_m/df = 0.0648148;  K -> inf, theta = 0: floor = 1 - sinc^2(eps) ~ (pi eps)^2/3 - 2 (pi eps)^4/45 = 1.37442e-2 (-18.62 dB).
 */
const Core = require('../core.js');
const U = require('./_util.js');

const sinc2 = x => { const s = x === 0 ? 1 : Math.sin(Math.PI * x) / (Math.PI * x); return s * s; };

function sysWith(over, model) {
    const s = Core.createSys(); Core.setSeed(11); s.rollCal();
    Object.assign(s, { model, freshRealization: false, algo: 'FOURIER', N: 8, L: 100, aoaT: 30, aoaJ: -40, kDb: 20, latMs: 3, d_min: 1e30, calDeg: 0 }, over);
    s.snaps = []; s.snapKey = ''; if (model === 'unified') s.newRealization(); s.computeMath(); return s;
}

module.exports = {
    id: 'T9', title: 'fc and numerology: eps, ICI floor, T_snap, aging ratio scaling; 28 GHz / 120 kHz hand values',
    async run() {
        const checks = [], info = [];
        // ---- T9a
        let maxA = 0;
        for (const model of ['legacy', 'unified']) for (const fc of [3.5e9, 5e9, 28e9]) for (const scs of [15e3, 60e3, 120e3]) for (const v of [30, 300]) for (const th of [0, 30, 75]) {
            const s = sysWith({ fc, scs, v, aoaT: th }, model), eps = v / 3.6 * fc * Math.cos(th * Math.PI / 180) / (3e8 * scs), epsM = v / 3.6 * fc / (3e8 * scs);
            maxA = Math.max(maxA, Math.abs(s.eps0 - eps) / Math.max(Math.abs(eps), 1e-300), Math.abs(s.epsM - epsM) / epsM);
        }
        console.log(`T9a  max relative |eps_impl - v fc cos(theta)/(c df)| over 2 models x 3 fc x 3 df x 2 v x 3 theta = ${U.e(maxA)}`);
        checks.push(U.check('T9a eps equals v fc cos(theta)/(c df)', U.e(maxA), '1e-12', maxA < 1e-12));

        // ---- T9b
        let maxB = 0;
        for (const model of ['legacy', 'unified']) for (const [fc, scs] of [[5e9, 15e3], [3.5e9, 30e3], [28e9, 120e3]]) {
            const a = sysWith({ fc, scs, v: 200 }, model), b = sysWith({ fc: 2 * fc, scs: 2 * scs, v: 200 }, model);
            maxB = Math.max(maxB, Math.abs(a.eps0 - b.eps0), Math.abs(a.epsM - b.epsM), Math.abs(a.nuIciFloor - b.nuIciFloor));
        }
        console.log(`T9b  fc and df doubled: max |change| of eps0, eps_m, ICI floor = ${U.e(maxB)}`);
        checks.push(U.check('T9b fc x2 and df x2 leave eps and the ICI floor unchanged', U.e(maxB), '1e-12', maxB < 1e-12));

        // ---- T9c
        const base = sysWith({ fc: 5e9, scs: 15e3, v: 300 }, 'unified'), dDf = sysWith({ fc: 5e9, scs: 30e3, v: 300 }, 'unified'), dFc = sysWith({ fc: 10e9, scs: 15e3, v: 300 }, 'unified');
        const r1 = dDf.rho / base.rho, r2 = dDf.agingB / base.agingB, r3 = dFc.rho / base.rho, r4 = dFc.agingB / base.agingB;
        console.log(`T9c  df x2: rho ratio ${r1.toFixed(12)} (expect 0.5), tau*B_D ratio ${r2.toFixed(12)} (expect 1);  fc x2: rho ratio ${r3.toFixed(12)} (expect 2), tau*B_D ratio ${r4.toFixed(12)} (expect 2)`);
        const ok = Math.abs(r1 - 0.5) < 1e-12 && Math.abs(r2 - 1) < 1e-12 && Math.abs(r3 - 2) < 1e-12 && Math.abs(r4 - 2) < 1e-12;
        checks.push(U.check('T9c df x2: T_snap halves (rho x0.5), tau*B_D unchanged; fc x2: both x2', `${r1.toFixed(6)} / ${r2.toFixed(6)} / ${r3.toFixed(6)} / ${r4.toFixed(6)}`, 'within 1e-12 of 0.5 / 1 / 2 / 2', ok));
        // absolute check of T_snap: rho = L T_snap B_D, so rho / (L B_D) must equal (1 + cpRatio)/df
        const eT = Math.abs(base.rho / (base.L * base.bD) - (1 + Core.CONFIG.cpRatio) / 15e3) / ((1 + Core.CONFIG.cpRatio) / 15e3);
        checks.push(U.check('T9c T_snap = (1 + CP)/df (from rho / (L B_D))', U.e(eT), '1e-12', eT < 1e-12));

        // ---- T9d
        const s = sysWith({ fc: 28e9, scs: 120e3, v: 300, kDb: 400, aoaT: 0 }, 'unified'), sD = sysWith({ fc: 28e9, scs: 120e3, v: 300, kDb: 20, aoaT: 45 }, 'unified');
        const hand = { fm: 7777.78, eps: 0.0648148, floor: 1.37442e-2 };
        const e1 = Math.abs(s.fm - hand.fm) / hand.fm, e2 = Math.abs(s.epsM - hand.eps) / hand.eps, e3 = Math.abs(s.nuIciFloor - hand.floor) / hand.floor;
        console.log(`T9d  28 GHz, 120 kHz, 300 km/h:  f_m = ${s.fm.toFixed(2)} Hz (hand ${hand.fm}),  eps_m = ${s.epsM.toFixed(7)} (hand ${hand.eps}),  ICI floor (K -> inf, theta = 0) = ${U.e(s.nuIciFloor, 5)} = ${(10 * Math.log10(s.nuIciFloor)).toFixed(2)} dB (hand ${hand.floor} = -18.62 dB)`);
        console.log(`     K = 20 dB, theta = 45 deg: eps0 = ${sD.eps0.toFixed(6)}, ICI floor = ${U.e(sD.nuIciFloor, 4)} (${(10 * Math.log10(sD.nuIciFloor)).toFixed(2)} dB) -> EVM floor ${(100 * Math.sqrt(sD.nuIciFloor)).toFixed(2)} %;  LoS-only at 45 deg: 1 - sinc^2(eps0) = ${U.e(1 - sinc2(sD.eps0), 4)}`);
        checks.push(U.check('T9d f_m at 28 GHz, 300 km/h', `${s.fm.toFixed(2)} Hz`, `hand 7777.78 Hz (1e-6 rel.)`, e1 < 1e-6));
        checks.push(U.check('T9d eps_m at 28 GHz, 120 kHz, 300 km/h', s.epsM.toFixed(7), 'hand 0.0648148 (1e-6 rel.)', e2 < 1e-6));
        checks.push(U.check('T9d ICI floor (K -> inf, theta = 0) vs series hand value', U.e(s.nuIciFloor, 5), 'hand 1.37442e-2 (1e-4 rel.)', e3 < 1e-4));
        return { id: this.id, title: this.title, checks, info };
    }
};
