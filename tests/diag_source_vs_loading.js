'use strict';
/**
 * 診斷（資訊性，不進 run_all）：Estimated 與導頻殘差的差異有多少來自「樣本來源」、多少來自「加載或收縮」。   node tests/diag_source_vs_loading.js [draws]
 * Nr = 2、4、8，SNR 30 dB，m = 1，隨機到達角，CP 法 ε̂，成對。
 *   Estimated（1024 個子載波）= pilot_residual_curves.js 的 Estimated：R_yy 取全部 1024 個子載波（含資料訊號），Ĥ 來自 86 個導頻。
 *   Estimated（86 個導頻）= limited_pilots.js 的 Estimated 在 P_eff = 86：R_yy 取同一批 86 個導頻（與導頻殘差相同的樣本來源），Ĥ 同。
 *   PR none = 導頻殘差無加載（P_eff > Nr 時與 Estimated（86 個導頻）代數上相同）；PR 0 dB / +10 dB / Ledoit-Wolf = 加載或收縮。
 * 樣本來源的效應 = Estimated(1024) − Estimated(86)；加載或收縮的效應 = Estimated(86) − PR（正值 = 加載或收縮較好）。
 */
const A = require('./pilot_residual_curves.js'), B = require('./limited_pilots.js'), U = require('./_util.js');
const draws = +(process.argv[2] || 500), c = a => `${U.f(U.mean(a), 2)} ± ${U.f(U.se(a), 2)}`, d = (x, y) => x.map((v, i) => v - y[i]);
for (const Nr of [2, 4, 8]) {
    const a = A.computeNr(Nr, draws, { snrs: [30], ms: [1] }), b = B.computeNr(Nr, draws, { snrs: [30], peffs: [86] }), M = B.METHODS;
    const E1024 = a.acc[0][A.METHODS[0]][0], E86 = b.acc[0][M[0]][0], N0 = b.acc[0][M[1]][0], P0 = b.acc[0][M[2]][0], P10 = b.acc[0][M[3]][0], LW = b.acc[0][M[4]][0];
    console.log(`Nr ${Nr}: Estimated(1024 subcarriers) ${c(E1024)} | Estimated(86 pilots) ${c(E86)} | PR none ${c(N0)} | PR 0 dB ${c(P0)} | PR +10 dB ${c(P10)} | LW ${c(LW)} || source effect ${c(d(E1024, E86))}; loading effect: 0 dB ${c(d(E86, P0))}, +10 dB ${c(d(E86, P10))}, LW ${c(d(E86, LW))}`);
}
