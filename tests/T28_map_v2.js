'use strict';
/**
 * T28  地圖 v2（tests/map_v2.js）的內部一致性。判準在看結果之前訂定，不調整；地圖本身的數字（主導比例、臨界點）是資訊性輸出，不在這裡判定。
 *   T28a  crossing()：合成資料 d = [1, 1, −1, −1] 在 Nr = [2, 4, 8, 16] 上的臨界點 = 2^2.5（誤差 < 1e-9）；d 全為正時回傳 null。
 *   T28b  domOf()：(0.2, 0.3) → '-'；(1, 2) → 'E'；(2, 1) → 'I'；(1, 1) → 'I'（相等歸 ICI）。
 *   T28c  與 map_prototype.csv 交叉驗證：同一格（v = 300 km/h、fc = 7 GHz、Δf = 30 kHz、Nr = 4、SNR 20 dB，200 次實現、同一個種子）重算，
 *         L_est_naive 與 L_ICI（原始定義）和 CSV 的值相差 < 0.0015 dB（CSV 取三位小數）。兩者是同一條計算路徑，這是迴歸檢查，不是獨立驗證。
 *   T28d  data/map_v2.csv：excluded 欄位在每一列都等於 (eps_max > 0.5)；v = 0 的列 L_ICI_net 恰為 0，所有列 L_ICI_net ≥ 0；列數 = 2112。
 */
const fs = require('fs'), path = require('path');
const U = require('./_util.js');
const M = require('./map_v2.js');
const R = require('./reproduce_gopala_slock.js');

const readCsv = f => fs.readFileSync(path.join(__dirname, '..', 'data', f), 'utf8').split('\n').filter(l => l && !l.startsWith('#'));
module.exports = {
    id: 'T28', title: 'map v2: crossing and dominance helpers, cross-check against map_prototype.csv, exclusion flag and baseline consistency of map_v2.csv',
    async run() {
        const checks = [], info = [];
        const x = M.crossing([2, 4, 8, 16], [1, 1, -1, -1], true);
        checks.push(U.check('T28a crossing() of a synthetic sign change', `${U.f(x, 6)}`, `2^2.5 = ${U.f(Math.pow(2, 2.5), 6)} +- 1e-9`, Math.abs(x - Math.pow(2, 2.5)) < 1e-9));
        checks.push(U.check('T28a crossing() without a sign change', String(M.crossing([2, 4, 8, 16], [1, 2, 3, 4], true)), 'null', M.crossing([2, 4, 8, 16], [1, 2, 3, 4], true) === null));
        const dm = [M.domOf(0.2, 0.3), M.domOf(1, 2), M.domOf(2, 1), M.domOf(1, 1)];
        checks.push(U.check('T28b domOf()', dm.join(' '), "- E I I", dm.join(' ') === '- E I I'));
        const pr = readCsv('map_prototype.csv').slice(1).map(l => l.split(',')).find(r => +r[0] === 300 && +r[1] === 7 && +r[2] === 30 && +r[3] === 4), c = M.cell(300, 7, 30, 4, 20, 200, R.CFG.seed);
        checks.push(U.check('T28c cell (300 km/h, 7 GHz, 30 kHz, Nr 4, SNR 20): L_est_naive against map_prototype.csv', `${U.f(c.Lnaive, 4)} vs ${pr[9]}`, '< 0.0015 dB', Math.abs(c.Lnaive - +pr[9]) < 0.0015));
        checks.push(U.check('T28c same cell: L_ICI (raw definition) against map_prototype.csv', `${U.f(c.Lraw, 4)} vs ${pr[5]}`, '< 0.0015 dB', Math.abs(c.Lraw - +pr[5]) < 0.0015));
        const rows = readCsv('map_v2.csv').slice(1).map(l => l.split(',')); // v,fc,df,Nr,snr,eps_max,excluded,L_ICI_net,...
        const badEx = rows.filter(r => (+r[5] > 0.5 ? '1' : '0') !== r[6]).length, badBase = rows.filter(r => +r[0] === 0 && +r[7] !== 0).length, neg = rows.filter(r => +r[7] < 0).length;
        checks.push(U.check('T28d map_v2.csv: number of rows', String(rows.length), '2112', rows.length === 2112));
        checks.push(U.check('T28d excluded flag equals (eps_max > 0.5)', `${badEx} mismatches`, '0', badEx === 0));
        checks.push(U.check('T28d L_ICI_net is 0 at v = 0 and never negative', `${badBase} nonzero at v = 0, ${neg} negative`, '0 and 0', badBase === 0 && neg === 0));
        info.push({ name: 'T28 map v2 cells excluded (eps_max > 0.5)', value: `${rows.filter(r => r[6] === '1').length} of ${rows.length}` });
        return { id: this.id, title: this.title, checks, info };
    }
};
