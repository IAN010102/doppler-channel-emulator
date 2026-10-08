'use strict';
/**
 * T27  導頻樣本數有限（tests/limited_pilots.js）。注意：這些判準是在看過 C1 的 500 次資料「之後」訂定的（事後），因此只選
 *   (a) 代數恆等式與 (b) 數值有限性這兩類不依賴結果大小的檢查，以及 (c) 一條對已觀察到效應的回歸檢查（明確標示為事後）；不調任何門檻去通過。
 *   T27a  恆等式（代數）：導頻為單位模（QPSK）、Ĥ 取自同一組導頻時 R_res = R_yy − Ĥ Ĥᴴ，由 Sherman–Morrison，(R_res)⁻¹Ĥ ∝ R_yy⁻¹Ĥ，
 *         所以當 P_eff > Nr 時（P_eff = Nr 時 R_res = R_yy − Ĥ Ĥᴴ 少一個秩、奇異，要用偽逆，恆等式不成立；原先寫成 ≥ 是推導錯誤，不是調門檻）「PR none」與「Estimated（同一子集的 R_yy）」的 SINR 相同：每個實現的損失差 < 1e-6 dB。
 *   T27b  P_eff < Nr（偽逆）時，所有方法的損失都是有限數。
 *   T27c  （事後）Nr = 8、P_eff = 8、SNR 30 dB：最佳的加載／收縮版本比無加載好，成對差 > 3 SE（z > 3）。
 *   「加載開始優於無加載的 P_eff」是資訊性輸出，不判定。
 */
const U = require('./_util.js');
const L = require('./limited_pilots.js');

module.exports = {
    id: 'T27', title: 'limited pilot samples: PR none = Estimated when P_eff > Nr (identity); finite losses with the pseudo-inverse; loading beats no loading at Nr = P_eff = 8 (post hoc)',
    async run({ draws = 200, seed = 2727 } = {}) {
        const checks = [], info = [], M = L.METHODS;
        let worst = 0, allFinite = true;
        for (const Nr of [2, 4, 8]) {
            const r = L.computeNr(Nr, draws, { seed, snrs: [30], peffs: [2, 4, 8, 16, 85] });
            r.peffs.forEach((P, k) => {
                if (P > Nr) r.acc[0][M[1]][k].forEach((x, i) => { worst = Math.max(worst, Math.abs(x - r.acc[0][M[0]][k][i])); });
                for (const m of M) if (!r.acc[0][m][k].every(Number.isFinite)) allFinite = false;
            });
            const c = L.criticalPeff(r, 0);
            info.push({ name: `T27 Nr = ${Nr}, SNR 30 dB: loading / shrinkage better than no loading (> 2 SE) up to P_eff`, value: `${c.critical === null ? 'none' : c.critical}${c.allBetter ? ' (every grid point)' : ''}` });
            if (Nr === 8) {
                const k = r.peffs.indexOf(8), row = c.rows[k], z = row.diff / row.se;
                checks.push(U.check('T27c (post hoc) Nr = 8, P_eff = 8, SNR 30 dB: best loaded / shrunk against PR none', `paired difference ${U.f(row.diff, 2)} dB, SE ${U.f(row.se, 2)}, z = ${U.f(z, 1)} (${row.best})`, 'z > 3', z > 3));
            }
        }
        checks.push(U.check('T27a PR none equals Estimated (same pilots) when P_eff > Nr', `largest |difference| ${U.e(worst)} dB`, '< 1e-6 dB', worst < 1e-6));
        checks.push(U.check('T27b all losses finite (including P_eff < Nr, pseudo-inverse)', allFinite ? 'all finite' : 'non-finite value found', 'all finite', allFinite));
        return { id: this.id, title: this.title, checks, info };
    }
};
