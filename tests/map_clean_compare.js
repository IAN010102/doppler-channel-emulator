'use strict';
/**
 * 乾淨估計器那片的兩版並列（資訊性，不判定；只讀既有 CSV，不重新模擬，不覆蓋任何舊檔）。   node tests/map_clean_compare.js
 *   舊版：data/map_prototype.csv 的 L_est（導頻殘差，γ_rel = +10 dB，SNR 20 dB）與其 dominant_net。
 *   新版：data/map_v2.csv 的 L_est（導頻殘差，無加載）與 dominant_clean（SNR 20 dB；同一組 704 個格子，同一個種子與實現數）。
 * 兩版的 L_ICI_net 相同（同一條計算路徑），主導因素 = L_ICI_net 與 L_est 中較大者，兩者都 < 0.5 dB 記 '-'；eps_max > 0.5 的格子標記並排除在比例之外。
 * 第三版（Ledoit–Wolf 收縮，data/map_clean_lw.csv，tests/map_clean_lw.js 產生）以新欄位 L_est_lw_db、dominant_lw 附在最後，舊欄位不動。
 * 輸出 data/map_clean_compare.csv 與 docs/diagnostics/map_clean_compare.txt。
 */
const fs = require('fs'), path = require('path');
const U = require('./_util.js');
const rd = f => fs.readFileSync(path.join(__dirname, '..', 'data', f), 'utf8').split('\n').filter(l => l && !l.startsWith('#')).slice(1).map(l => l.split(','));
const dom = (a, b) => Math.max(a, b) < 0.5 ? '-' : (a >= b ? 'I' : 'E');

const old = new Map(rd('map_prototype.csv').map(r => [`${r[0]}|${r[1]}|${r[2]}|${r[3]}`, r]));      // v,fc,df,Nr,eps_max,L_ICI,se,L_est,se,naive,se,admf,dominant,L_ICI_net,dominant_net,aging,n
const lwMap = new Map(rd('map_clean_lw.csv').map(r => [`${r[0]}|${r[1]}|${r[2]}|${r[3]}`, +r[4]]));
const v2 = rd('map_v2.csv').filter(r => +r[4] === 20);                                                // v,fc,df,Nr,snr,eps_max,excl,net,raw,L_est,se,naive,se,dom_clean,dom_naive,n
const rows = [];
for (const r of v2) {
    const o = old.get(`${r[0]}|${r[1]}|${r[2]}|${r[3]}`), net = +r[7], lOld = +o[7], lNew = +r[9];
    const lLw = lwMap.get(`${r[0]}|${r[1]}|${r[2]}|${r[3]}`);
    rows.push({ lLw, dLw: dom(net, lLw), v: r[0], fc: r[1], df: r[2], Nr: r[3], eps: r[5], ex: r[6], net, lOld, lNew, dOld: dom(net, lOld), dNew: dom(net, lNew), netOld: +o[13] });
}
const inc = rows.filter(x => x.ex === '0'), cnt = (a, k) => { const c = { I: 0, E: 0, '-': 0 }; for (const x of a) c[x[k]]++; return c; }, pc = (c, n) => `ICI ${c.I} (${U.f(100 * c.I / n, 1)} %), estimation ${c.E} (${U.f(100 * c.E / n, 1)} %), neither ${c['-']} (${U.f(100 * c['-'] / n, 1)} %)`;
const netDiff = Math.max(...rows.map(x => Math.abs(x.net - x.netOld))), changed = inc.filter(x => x.dOld !== x.dNew).length;
const lines = [`SNR 20 dB, ${rows.length} cells, excluded (eps_max > 0.5): ${rows.length - inc.length}, counted: ${inc.length}`,
    `old  (pilot residual, gamma_rel +10 dB, data/map_prototype.csv): ${pc(cnt(inc, 'dOld'), inc.length)}`,
    `new  (no loading, data/map_v2.csv):                               ${pc(cnt(inc, 'dNew'), inc.length)}`,
    `LW   (Ledoit-Wolf shrinkage, data/map_clean_lw.csv):              ${pc(cnt(inc, 'dLw'), inc.length)}`,
    `cells whose dominant factor differs between old and new: ${changed}; between no loading and LW: ${inc.filter(x => x.dNew !== x.dLw).length}; largest |L_ICI_net(map_v2) - L_ICI_net(map_prototype)| = ${U.e(netDiff)} dB (CSV rounding)`];
for (const Nr of [2, 4, 8, 16]) { const s = inc.filter(x => +x.Nr === Nr), m = k => U.mean(s.map(x => x[k])), mx = k => Math.max(...s.map(x => x[k])); lines.push(`  Nr ${Nr}: old ${pc(cnt(s, 'dOld'), s.length)}; new ${pc(cnt(s, 'dNew'), s.length)}; LW ${pc(cnt(s, 'dLw'), s.length)}; L_est mean old ${U.f(m('lOld'), 2)} / new ${U.f(m('lNew'), 2)} / LW ${U.f(m('lLw'), 2)} dB, max old ${U.f(mx('lOld'), 2)} / new ${U.f(mx('lNew'), 2)} / LW ${U.f(mx('lLw'), 2)} dB`); }
fs.writeFileSync(path.join(__dirname, '..', 'data', 'map_clean_compare.csv'), '# Clean-estimator slice, SNR 20 dB: old pilot residual with gamma_rel +10 dB (data/map_prototype.csv) next to the new one without loading (data/map_v2.csv); both old files are unchanged; the last two columns (Ledoit-Wolf, SNR 20 dB) were added later\n# dominant = larger of L_ICI_net and L_est (I = ICI, E = estimation, - = both < 0.5 dB); excluded = 1 when eps_max > 0.5\nv_kmh,fc_GHz,df_kHz,Nr,eps_max,excluded,L_ICI_net_db,L_est_old_db,L_est_new_db,dominant_old,dominant_new,L_est_lw_db,dominant_lw\n' + rows.map(x => [x.v, x.fc, x.df, x.Nr, x.eps, x.ex, x.net.toFixed(3), x.lOld.toFixed(3), x.lNew.toFixed(3), x.dOld, x.dNew, x.lLw.toFixed(3), x.dLw].join(',')).join('\n') + '\n');
fs.writeFileSync(path.join(__dirname, '..', 'docs', 'diagnostics', 'map_clean_compare.txt'), lines.join('\n') + '\n'); console.log(lines.join('\n'));
