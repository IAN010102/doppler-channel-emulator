'use strict';
/**
 * node tests/check_report_numbers.js [--suite <path to the saved output of a full `node tests/run_all.js`>] [--out <file>]
 *
 * Read-only check of every number in item 10 of docs/FINAL_REPORT.md (the sentences meant for papers): the number is taken from the report text (regular expression), the
 * source value is read again from its source file, and the two are compared; the allowed difference is the rounding of the last significant digit written in the report
 * (half a unit of the last digit). Sources:
 *   data/experiment_numbers.json, docs/diagnostics/*.{csv,txt}, experiments.js / core.js / the test files (parameters and numbers of realisations, read from their source code),
 *   and, only with --suite, the terminal output of a full run of the test suite (numbers that exist only there: T11b, T17b2, T17c, T17d, T18c).
 * A number without an automatic source is listed as "manual" with the test or script to rerun; nothing is assumed. Numbers of the sentence that no claim covers are listed as "manual (not mapped)".
 * Nothing in the repository is modified (the table is printed, and written to --out when given).
 */
const fs = require('fs'), path = require('path'), os = require('os');
const ROOT = path.join(__dirname, '..');
const rd = f => fs.readFileSync(path.join(ROOT, f), 'utf8').replace(/\r\n/g, '\n');
const argv = process.argv.slice(2), opt = n => { const i = argv.indexOf(n); return i >= 0 ? argv[i + 1] : null; };
const ARCHIVE = path.join(ROOT, 'docs', 'diagnostics', 'suite_output_2026-10-06_198of198.txt');   // the saved output of the complete run (198/198) of 2026-10-06; used for T17 only (T11 and T18 have fresh outputs)
const suitePath = opt('--suite') || (argv.includes('--no-archive') ? null : ARCHIVE);
const suiteText = suitePath ? fs.readFileSync(suitePath, 'utf8').replace(/\r\n/g, '\n') : null;
const freshT11 = rd('docs/diagnostics/fresh_T11_output.txt'), freshT18 = rd('docs/diagnostics/fresh_T18_output.txt');   // outputs of `node` runs of T11 and T18 alone (see report_number_check.txt)
const Core = require('../core.js'), X = require('../experiments.js');
const J = JSON.parse(rd('data/experiment_numbers.json'));

// ------------------------------------------------------------------ the sentences of item 10
const report = rd('docs/FINAL_REPORT.md'), i10 = report.indexOf('## 10. ');
const lines10 = report.slice(i10).split('\n');
const SENT = {};
for (const l of lines10) { const m = /^(\d+)\. /.exec(l); if (m) SENT[+m[1]] = l; else if (l.startsWith('設定縮寫')) SENT[0] = l; }
if (Object.keys(SENT).length !== 12) throw new Error('expected the preamble and 11 sentences in item 10, found ' + Object.keys(SENT).length);

// ------------------------------------------------------------------ source readers
const num = s => parseFloat(String(s).replace(/−/g, '-'));
const dec = s => { const m = /\.(\d+)/.exec(String(s)); return m ? m[1].length : 0; };
function parseBase(file, name = 'BASE') { const m = new RegExp('const ' + name + ' = (\\{[^;]*\\});').exec(rd(file)); return m ? new Function('return ' + m[1])() : null; }
const TB = {}; for (const f of ['t11_mmse', 't17_pointing', 't18_doa', 't21_lecture', 't13_symbol_evm']) { try { TB[f] = parseBase('tests/' + f + '.js'); } catch (e) { TB[f] = null; } }
const XB = X.BASE, XP = (id, k) => X.params(id)[k];
const CFG = k => Core.CONFIG[k];
const s1rows = (() => { const L = rd('docs/diagnostics/s1_sensitivity.csv').trim().split('\n'), h = L[0].split(','); return L.slice(1).map(l => { const c = l.split(','), o = {}; h.forEach((k, i) => { o[k] = isNaN(+c[i]) ? c[i] : +c[i]; }); return o; }); })();
const s1txt = rd('docs/diagnostics/s1_sensitivity.txt'), s3txt = rd('docs/diagnostics/s3_withsignal.txt'), b4txt = rd('docs/diagnostics/b4_afc.txt');
const code = { measure: rd('tests/experiment_measure.js'), diagS: rd('tests/diag_s_sensitivity.js'), t11: rd('tests/t11_mmse.js'), t17: rd('tests/t17_pointing.js'), t18: rd('tests/t18_doa.js'), t21: rd('tests/t21_lecture.js'), t3: rd('tests/t3_aging.js'), core: rd('core.js'), exps: rd('experiments.js') };
const fromCode = (key, re) => { const m = re.exec(code[key]); return m ? num(m[1]) : null; };
const S3 = (() => {   // the first table of s3_withsignal.txt: rows (i) K -> inf, (ii) diffuse only, (iii) original; columns v = 0, 25 ... 300
    const L = s3txt.split('\n'), h = L.findIndex(l => l.startsWith('variant')), vs = L[h].trim().split(/\s+/).slice(1).map(Number), rows = {};
    for (let i = 1; i <= 3; i++) { const l = L[h + i], name = /^\((i+)\)/.exec(l)[1], cells = l.slice(34).trim().split(/\s+/).map(c => parseFloat(c)); rows[name] = cells; }
    return { vs, rows };
})();
function b4block(k, th) { const L = b4txt.split('\n'), i = L.findIndex(l => l.startsWith('K = ' + k) && l.includes('theta1 = ' + th + ',') && l.includes('signalFree')); const out = {}; for (let j = i + 2; j < i + 6; j++) out[+L[j].trim().split(/\s+/)[0]] = L[j].trim(); return out; }
function b4val(k, th, v, col) {   // col: 'mmsep' | 'genie' | 'periodogram' ; the row is "v cycles omega  (i) +- se gap  (ii) +- se gap ..."
    const t = b4block(k, th)[v].replace(/±/g, ' ').trim().split(/\s+/).map(Number);   // v cycles omega i se gap ii se gap ...
    return { mmsep: t[3], genie: t[6] }[col];
}
const suite = (re, g = 1) => { if (suiteText === null) return undefined; const m = re.exec(suiteText); return m ? num(m[g]) : null; };
const fresh11 = (re, g = 1) => { const m = re.exec(freshT11); return m ? num(m[g]) : null; };
// the console tables of T18c: "vs L (theta1 = 0, theta2 = 40):" rows "L  MUSIC T (+-)  MUSIC J  Capon T  Capon J" and "vs separation (theta1 = 0, L = 100):" rows "theta2  ..."
function t18rows(section) { const L = freshT18.split('\n'), i = L.findIndex(l => l.startsWith(section)); const out = []; for (let j = i + 2; j < L.length && /^\d+\s/.test(L[j]); j++) { const t = L[j].replace(/±/g, ' ').trim().split(/\s+/).map(Number); out.push({ key: t[0], musicT: t[1], capT: t[5] }); } return out; }
const SUITEDESC = 'archived terminal output of the complete run of 2026-10-06 (docs/diagnostics/suite_output_2026-10-06_198of198.txt)';

// ------------------------------------------------------------------ claims
const claims = [];
/** claim(sentence, label, regex of the report text (first capture group = the number as written), getter, source text, options) */
function claim(s, label, re, get, src, o = {}) { claims.push({ s, label, re, get, src, o }); }
const SUITE = SUITEDESC, RERUN = (t, min) => `rerun ${t} (about ${min})`;
const needSuite = (t, min) => ({ manualHint: RERUN(t, min) });

// --- preamble (setting abbreviations)
claim(0, 'N', /N = (8)、/, () => CFG('N'), 'core.js CONFIG.N');
claim(0, 'd/lambda', /d = (0\.5)λ/, () => CFG('d_lambda'), 'core.js CONFIG.d_lambda');
claim(0, 'SNR dB', /SNR (20) dB/, () => XB.snr, 'experiments.js BASE.snr');
claim(0, 'SIR dB', /SIR (−10) dB/, () => XB.sir, 'experiments.js BASE.sir');
claim(0, 'theta2 deg', /θ₂ = (40)°/, () => XB.aoaJ, 'experiments.js BASE.aoaJ');
// --- sentence 1
claim(1, 'fc GHz', /fc = (5) GHz/, () => CFG('fc') / 1e9, 'core.js CONFIG.fc');
claim(1, 'df kHz', /Δf = (15) kHz/, () => CFG('scs') / 1e3, 'core.js CONFIG.scs');
claim(1, 'd_min m', /d_min = (30) m/, () => CFG('d_min'), 'core.js CONFIG.d_min');
claim(1, 'L', /L = (100)、/, () => TB.t11_mmse.L, 'tests/t11_mmse.js BASE.L');
claim(1, 'K dB', /K = (20) dB/, () => TB.t11_mmse.kDb, 'tests/t11_mmse.js BASE.kDb');
claim(1, 'v', /v = (0)、/, () => TB.t11_mmse.v, 'tests/t11_mmse.js BASE.v');
claim(1, 'tau ms', /τ = (0)、/, () => TB.t11_mmse.latMs, 'tests/t11_mmse.js BASE.latMs');
claim(1, 'theta1 deg', /θ₁ = (0)°/, () => TB.t11_mmse.aoaT, 'tests/t11_mmse.js BASE.aoaT');
claim(1, 'theta2 deg', /（θ₂ = (40)°、SNR/, () => TB.t11_mmse.aoaJ, 'tests/t11_mmse.js BASE.aoaJ');
claim(1, 'SNR dB', /SNR (20) dB、SIR/, () => TB.t11_mmse.snr, 'tests/t11_mmse.js BASE.snr');
claim(1, 'SIR dB', /SIR (−10) dB、N/, () => TB.t11_mmse.sir, 'tests/t11_mmse.js BASE.sir');
claim(1, 'N', /N = (8)）/, () => TB.t11_mmse.N, 'tests/t11_mmse.js BASE.N');
claim(1, 'MMSE-P SINR dB (T11b)', /都約 (28\.6) dB/, () => fresh11(/T11b unified MMSE-P \(signalFree\): ([\d.]+) ±/), 'fresh run of T11 (docs/diagnostics/fresh_T11_output.txt): "T11b unified MMSE-P (signalFree)"');
claim(1, 'SMI SINR dB (E4, 0 deg)', /SMI：E4 的 0° 點 (28\.59) dB/, () => J.E4.signalFree['0'], 'data/experiment_numbers.json E4.signalFree["0"]');
claim(1, 'gap to SINR_opt dB (T11b)', /SINR_opt 約 (0\.3) dB/, () => { const g = fresh11(/T11b unified MMSE-P \(signalFree\):.*gap (-?[\d.]+) dB/); return g === null ? g : Math.abs(g); }, 'fresh run of T11 (docs/diagnostics/fresh_T11_output.txt): gap of "T11b unified MMSE-P (signalFree)"');
claim(1, 'realisations T11b', /T11b，(2000) 次/, () => fromCode('t11', /trials = (\d+), seed = 7000/), 'tests/t11_mmse.js default trials');
claim(1, 'realisations E4/T17c', /T17c，(1000) 次/, () => fromCode('measure', /E4: SMI vs delta_theta[\s\S]*?pointErrDeg: d\b[^\n]*?, (1000), 4000\)/), 'tests/experiment_measure.js (E4 block, 1000 trials)');
// --- sentence 2
claim(2, 'fc GHz', /fc = (5) GHz/, () => CFG('fc') / 1e9, 'core.js CONFIG.fc');
claim(2, 'df kHz', /Δf = (15) kHz/, () => CFG('scs') / 1e3, 'core.js CONFIG.scs');
claim(2, 'v', /v = (0)、/, () => XP('E3', 'v'), 'experiments.js params(E3).v');
claim(2, 'tau ms', /τ = (0)、/, () => XP('E3', 'latMs'), 'experiments.js params(E3).latMs');
claim(2, 'theta1 deg', /θ₁ = (0)°、/, () => XP('E3', 'aoaT'), 'experiments.js params(E3).aoaT');
claim(2, 'theta2 deg', /θ₂ = (40)°、/, () => XP('E3', 'aoaJ'), 'experiments.js params(E3).aoaJ');
claim(2, 'L', /L = (100)、/, () => XP('E3', 'L'), 'experiments.js params(E3).L');
claim(2, 'K dB', /K = (20) dB；/, () => XP('E3', 'kDb'), 'experiments.js params(E3).kDb');
claim(2, 'delta_theta deg', /δθ = (0)，SMI/, () => null, 'tests/experiment_measure.js E3 block: d0 = pointErrDeg 0', { fileHas: ['measure', /d0: sinr\(\{ algo: a, trainMode: 'withSignal', pointErrDeg: 0 \}/] });
const t17d0 = () => suite(/T17d withSignal L=100   SMI: (-?[\d.]+) \(0°\)/);
claim(2, 'SMI withSignal SINR dB, upper end of the range (largest of the three runs)', /降到約 (−5\.5) 至/, () => { const t = t17d0(); return t === undefined || t === null ? t : Math.max(t, J.E3.SMI.d0, J.E4.withSignal['0']); }, SUITE + ' (T17d) and E3.SMI.d0, E4.withSignal["0"]', needSuite('T17', 'about 21 min'));
claim(2, 'SMI withSignal SINR dB, lower end of the range (smallest of the three runs)', /至 (−5\.7) dB（自我抵消）/, () => { const t = t17d0(); return t === undefined || t === null ? t : Math.min(t, J.E3.SMI.d0, J.E4.withSignal['0']); }, SUITE + ' (T17d) and E3.SMI.d0, E4.withSignal["0"]', needSuite('T17', 'about 21 min'));
claim(2, 'SMI withSignal SINR dB (T17d)', /T17d (−5\.5)、E3/, t17d0, SUITE + ': "T17d withSignal L=100"', needSuite('T17', 'about 21 min'));
claim(2, 'SMI withSignal SINR dB (E3)', /E3 (−5\.56)、E4/, () => J.E3.SMI.d0, 'data/experiment_numbers.json E3.SMI.d0');
claim(2, 'SMI withSignal SINR dB (E4)', /E4 (−5\.72)（T17d/, () => J.E4.withSignal['0'], 'data/experiment_numbers.json E4.withSignal["0"]');
claim(2, 'realisations of the three runs', /三者均為 (1000) 次實現/, () => fromCode('t17', /trials = (\d+), seed = 1717/), 'tests/t17_pointing.js default trials (T17d); E3 and E4: tests/experiment_measure.js (1000)');
const se3 = (() => { const t = rd('docs/diagnostics/withsignal_se_check.txt'), g = k => { const m = new RegExp(k + '.*?: -?[\\d.]+ ± ([\\d.]+)').exec(t); return m ? num(m[1]) : null; }; return [g('E3 setting'), g('E4 setting'), g('T17d setting')]; })();
claim(2, 'SE of the means, smallest dB', /標準誤差約 (0\.10)–0\.11 dB/, () => Math.min(...se3), 'docs/diagnostics/withsignal_se_check.txt');
claim(2, 'SE of the means, largest dB', /標準誤差約 0\.10–(0\.11) dB/, () => Math.max(...se3), 'docs/diagnostics/withsignal_se_check.txt');
claim(2, 'largest difference of the three values dB', /兩兩最大差 (0\.22) dB/, () => { const t = t17d0(); return t === undefined || t === null ? t : Math.max(t, J.E3.SMI.d0, J.E4.withSignal['0']) - Math.min(t, J.E3.SMI.d0, J.E4.withSignal['0']); }, SUITE + ' (T17d) and E3, E4', needSuite('T17', 'about 21 min'));
claim(2, 'standard error of the difference dB (E4 against T17d)', /差的標準誤差（(0\.15) dB）/, () => Math.sqrt(se3[1] * se3[1] + se3[2] * se3[2]), 'docs/diagnostics/withsignal_se_check.txt: sqrt(SE_E4^2 + SE_T17d^2)');
claim(2, 'largest difference / SE of the difference', /的 (1\.5) 倍，彼此/, () => { const t = t17d0(); if (t === undefined || t === null) return t; return (Math.max(t, J.E3.SMI.d0, J.E4.withSignal['0']) - Math.min(t, J.E3.SMI.d0, J.E4.withSignal['0'])) / Math.sqrt(se3[1] * se3[1] + se3[2] * se3[2]); }, 'largest difference (T17d, E3, E4) divided by sqrt(SE_E4^2 + SE_T17d^2)', needSuite('T17', 'about 21 min'));
claim(2, 'pointing error deg', /指向偏差 (3)° 時 SMI/, () => XP('E3', 'pointErrDeg'), 'experiments.js params(E3).pointErrDeg');
claim(2, 'drop SMI withSignal dB (3 deg)', /比 0° 低 (10) dB 以上/, () => J.E3.SMI.d0 - J.E3.SMI.d3, 'E3.SMI.d0 - E3.SMI.d3', { atLeast: true });
claim(2, 'drop DL withSignal dB (3 deg)', /比 0° 低 (10) dB 以上/, () => J.E3.DL.d0 - J.E3.DL.d3, 'E3.DL.d0 - E3.DL.d3', { atLeast: true, second: true });
claim(2, 'MMSE-P change dB', /MMSE-P 不受指向偏差影響/, () => Math.abs(J.E3.MMSEP.d0 - J.E3.MMSEP.d3) < 0.05, 'E3.MMSEP: |d0 - d3| < 0.05 dB', { bool: true });
claim(2, 'signalFree drop dB (3 deg)', /只降約 (0\.6) dB/, () => J.E4.signalFree['0'] - J.E4.signalFree['3'], 'E4.signalFree[0] - E4.signalFree[3]');
claim(2, 'SMI signalFree 0 deg (T17d)', /SMI (28\.6)→28\.0/, () => suite(/T17d signalFree L=100   SMI: (-?[\d.]+) \(0°\)/), SUITE + ': "T17d signalFree L=100"', needSuite('T17', '21 min'));
claim(2, 'SMI signalFree 3 deg (T17d)', /28\.6→(28\.0)，T17d/, () => suite(/T17d signalFree L=100   SMI: -?[\d.]+ \(0°\) \/ (-?[\d.]+) \(3°\)/), SUITE + ': "T17d signalFree L=100"', needSuite('T17', '21 min'));
// --- sentence 3
claim(3, 'N', /N=(8)、/, () => XP('E0', 'N'), 'experiments.js params(E0).N');
claim(3, 'theta1 deg', /θ₁=(−20)°/, () => XP('E0', 'aoaT'), 'experiments.js params(E0).aoaT');
claim(3, 'theta2 deg', /θ₂=(30)°/, () => XP('E0', 'aoaJ'), 'experiments.js params(E0).aoaJ');
claim(3, 'SNR dB', /SNR (30) dB/, () => XP('E0', 'snr'), 'experiments.js params(E0).snr');
claim(3, 'SIR dB', /SIR (0) dB/, () => XP('E0', 'sir'), 'experiments.js params(E0).sir');
claim(3, 'K dB', /（(400) dB）/, () => XP('E0', 'kDb'), 'experiments.js params(E0).kDb');
claim(3, 'v', /v = (0)、/, () => XP('E0', 'v'), 'experiments.js params(E0).v');
claim(3, 'tau ms', /τ = (0)、/, () => XP('E0', 'latMs'), 'experiments.js params(E0).latMs');
claim(3, 'fc GHz', /fc = (5) GHz/, () => XP('E0', 'fc') / 1e9, 'experiments.js params(E0).fc');
claim(3, 'df kHz', /Δf = (15) kHz/, () => XP('E0', 'scs') / 1e3, 'experiments.js params(E0).scs');
claim(3, 'd_min m', /d_min = (30) m/, () => XP('E0', 'd_min'), 'experiments.js params(E0).d_min');
claim(3, 'SINR_opt dB', /SINR_opt（(38\.97) dB）/, () => J.E0['MMSE-M (theory R)'].opt, 'E0["MMSE-M (theory R)"].opt');
claim(3, 'MMSE-M theory SINR dB', /SINR_opt（(38\.97) dB）/, () => J.E0['MMSE-M (theory R)'].sinr, 'E0["MMSE-M (theory R)"].sinr', { second: true });
claim(3, 'MVDR B theory SINR dB', /SINR_opt（(38\.97) dB）/, () => J.E0['SMI withSignal (theory R) = MVDR B'].sinr, 'E0[MVDR B].sinr', { second: true });
claim(3, 'MVDR A theory SINR dB', /SINR_opt（(38\.97) dB）/, () => J.E0['SMI signalFree (theory R) = MVDR A'].sinr, 'E0[MVDR A].sinr', { second: true });
claim(3, 'L', /L=(1000) 的樣本/, () => XP('E0', 'L'), 'experiments.js params(E0).L');
claim(3, 'sample R withSignal SINR dB', /降為 (21\.88) dB/, () => J.E0['SMI withSignal (sample R)'].sinr, 'E0[SMI withSignal (sample R)].sinr');
claim(3, 'MMSE-M sample R SINR dB', /降為 (21\.88) dB/, () => J.E0['MMSE-M (sample R)'].sinr, 'E0[MMSE-M (sample R)].sinr', { second: true });
claim(3, 'SMI signalFree sample R SINR dB', /仍有 (38\.94) dB/, () => J.E0['SMI signalFree (sample R)'].sinr, 'E0[SMI signalFree (sample R)].sinr');
claim(3, 'realisations T21', /預設 (300) 次實現/, () => fromCode('t21', /trials = (\d+), seed = 2121/), 'tests/t21_lecture.js default trials');
claim(3, 'realisations E0', /(200) 次實現的平均/, () => fromCode('measure', /n = Math\.min\(trials, (\d+)\)/), 'tests/experiment_measure.js: E0 uses Math.min(trials, 200)');
claim(3, 'L (again)', /L = (1000)。）/, () => XP('E0', 'L'), 'experiments.js params(E0).L');
// --- sentence 4
claim(4, 'L', /快照數 L=(4) </, () => XP('E2', 'L'), 'experiments.js params(E2).L');
claim(4, 'N', /< N=(8)/, () => XP('E2', 'N'), 'experiments.js params(E2).N');
claim(4, 'rank', /秩只有 (4)，/, () => J.E2.SMI.rank, 'E2.SMI.rank');
claim(4, 'DL SINR dB', /DL（(27\.81) dB）/, () => J.E2.DL.sinr, 'E2.DL.sinr');
claim(4, 'BEAMSPACE SINR dB', /BEAMSPACE（(26\.36) dB）/, () => J.E2.BEAMSPACE.sinr, 'E2.BEAMSPACE.sinr');
claim(4, 'SMI SINR dB', /SMI（(22\.46) dB）/, () => J.E2.SMI.sinr, 'E2.SMI.sinr');
claim(4, 'realisations', /E2，(1000) 次實現/, () => fromCode('measure', /E2: L = 4, N = 8[\s\S]*?\{ algo: a \}\), (1000), 2000\)/), 'tests/experiment_measure.js (E2 block)');
claim(4, 'fc GHz', /fc = (5) GHz/, () => XP('E2', 'fc') / 1e9, 'experiments.js params(E2).fc');
claim(4, 'df kHz', /Δf = (15) kHz/, () => XP('E2', 'scs') / 1e3, 'experiments.js params(E2).scs');
claim(4, 'v', /v = (0)、/, () => XP('E2', 'v'), 'experiments.js params(E2).v');
claim(4, 'tau ms', /τ = (0)、/, () => XP('E2', 'latMs'), 'experiments.js params(E2).latMs');
claim(4, 'theta1 deg', /θ₁ = (0)°、/, () => XP('E2', 'aoaT'), 'experiments.js params(E2).aoaT');
claim(4, 'theta2 deg', /θ₂ = (40)°、/, () => XP('E2', 'aoaJ'), 'experiments.js params(E2).aoaJ');
claim(4, 'K dB', /K = (20) dB、SNR/, () => XP('E2', 'kDb'), 'experiments.js params(E2).kDb');
claim(4, 'SNR dB', /SNR (20) dB、SIR/, () => XP('E2', 'snr'), 'experiments.js params(E2).snr');
claim(4, 'SIR dB', /SIR (−10) dB、signalFree/, () => XP('E2', 'sir'), 'experiments.js params(E2).sir');
claim(4, 'gamma_rel dB', /γ_rel = \+(10) dB/, () => XP('E2', 'gammaRelDb'), 'experiments.js params(E2).gammaRelDb');
// --- sentence 5
claim(5, 'v km/h', /在 v = (300) km\/h/, () => XP('E1', 'v'), 'experiments.js params(E1).v');
claim(5, 'theta1 deg', /θ₁ = (0)°、θ₂/, () => XP('E1', 'aoaT'), 'experiments.js params(E1).aoaT');
claim(5, 'theta2 deg', /θ₂ = (40)°、d_min/, () => XP('E1', 'aoaJ'), 'experiments.js params(E1).aoaJ');
claim(5, 'd_min m', /d_min = (30) m/, () => XP('E1', 'd_min'), 'experiments.js params(E1).d_min');
claim(5, 'tau ms', /τ = (0)、L/, () => XP('E1', 'latMs'), 'experiments.js params(E1).latMs');
claim(5, 'L', /L = (100)、K/, () => XP('E1', 'L'), 'experiments.js params(E1).L');
claim(5, 'K dB', /K = (20) dB（/, () => XP('E1', 'kDb'), 'experiments.js params(E1).kDb');
claim(5, 'SNR dB', /SNR (20) dB、SIR/, () => XP('E1', 'snr'), 'experiments.js params(E1).snr');
claim(5, 'SIR dB', /SIR (−10) dB、signalFree/, () => XP('E1', 'sir'), 'experiments.js params(E1).sir');
claim(5, 'fc GHz (E1 28)', /下，(28) GHz\/120 kHz/, () => XP('E1', 'fc') / 1e9, 'experiments.js params(E1,"28").fc', { variantArgs: ['28'] });
claim(5, 'df kHz (E1 28)', /GHz\/(120) kHz 的 SMI EVM/, () => X.params('E1', '28').scs / 1e3, 'experiments.js params(E1,"28").scs');
claim(5, 'EVM % (28 GHz)', /(12\.31) %、5 GHz/, () => J.E1_28.byV['300'].evm_pct, 'E1_28.byV[300].evm_pct');
claim(5, 'fc GHz (E1 5)', /、(5) GHz\/15 kHz 為/, () => X.params('E1', '5').fc / 1e9, 'experiments.js params(E1,"5").fc');
claim(5, 'df kHz (E1 5)', /GHz\/(15) kHz 為 17/, () => X.params('E1', '5').scs / 1e3, 'experiments.js params(E1,"5").scs');
claim(5, 'EVM % (5 GHz)', /為 (17\.12) %；/, () => J.E1_5.byV['300'].evm_pct, 'E1_5.byV[300].evm_pct');
claim(5, 'crossing speed km/h (5 GHz)', /約 (230) km\/h/, () => J.E1_5.crossing_speed_kmh, 'E1_5.crossing_speed_kmh');
claim(5, 'threshold %', /的 (13\.14) % 門檻/, () => 100 * fromCode('measure', /const thr = ([\d.]+);/), 'tests/experiment_measure.js: const thr (0.1314)');
claim(5, 'no crossing (28 GHz)', /前者在 0–300 km\/h 內沒有越過/, () => J.E1_28.crossing_speed_kmh === null, 'E1_28.crossing_speed_kmh === null', { bool: true });
claim(5, 'range end km/h', /在 0–(300) km\/h 內/, () => 300, 'tests/experiment_measure.js: loop v = 0 ... 300', { codeRe: [/for \(let v = 0; v <= (300); v \+= 10\)/, 'measure'] });
claim(5, 'realisations', /每點 (300) 次實現/, () => fromCode('measure', /function measure\(trials = (\d+)\)/), 'tests/experiment_measure.js: measure(trials = 300), evmRms(..., trials, ...)');
claim(5, 'speed step km/h', /速度步長 (10) km\/h/, () => fromCode('measure', /for \(let v = 0; v <= 300; v \+= (\d+)\)/), 'tests/experiment_measure.js: v += 10');
claim(5, 'crossing resolution km/h', /解析度為 (10) km\/h/, () => fromCode('measure', /for \(let v = 0; v <= 300; v \+= (\d+)\)/), 'tests/experiment_measure.js: v += 10');
// --- sentence 6
claim(6, 'tau ms (from)', /從 (0) 增到/, () => 0, 'E5.d30.SMI keys 0, 5, 10', { codeRe: [/for \(const tau of \[(0), 5, 10\]\)/, 'measure'] });
claim(6, 'tau ms (to)', /增到 (10) ms/, () => XP('E5', 'latMs'), 'experiments.js params(E5).latMs');
claim(6, 'v km/h', /v=(300) km\/h/, () => XP('E5', 'v'), 'experiments.js params(E5).v');
claim(6, 'd_min m', /d_min=(30) m/, () => XP('E5', 'd_min'), 'experiments.js params(E5).d_min (the d30 curve)');
claim(6, 'theta1 deg', /θ₁ = (0)°、θ₂/, () => XP('E5', 'aoaT'), 'experiments.js params(E5).aoaT');
claim(6, 'theta2 deg', /θ₂=(−30)°/, () => XP('E5', 'aoaJ'), 'experiments.js params(E5).aoaJ');
claim(6, 'K dB', /K=(20) dB、L/, () => XP('E5', 'kDb'), 'experiments.js params(E5).kDb');
claim(6, 'L', /L=(100)、fc/, () => XP('E5', 'L'), 'experiments.js params(E5).L');
claim(6, 'fc GHz', /fc = (5) GHz/, () => XP('E5', 'fc') / 1e9, 'experiments.js params(E5).fc');
claim(6, 'df kHz', /Δf = (15) kHz/, () => XP('E5', 'scs') / 1e3, 'experiments.js params(E5).scs');
claim(6, 'SNR dB', /SNR (20) dB、SIR/, () => XP('E5', 'snr'), 'experiments.js params(E5).snr');
claim(6, 'SIR dB', /SIR (−10) dB、signalFree/, () => XP('E5', 'sir'), 'experiments.js params(E5).sir');
claim(6, 'SMI SINR dB at tau 0', /從 (28\.31) 降到/, () => J.E5.d30.SMI['0'], 'E5.d30.SMI[0]');
claim(6, 'SMI SINR dB at tau 10', /降到 (25\.48) dB/, () => J.E5.d30.SMI['10'], 'E5.d30.SMI[10]');
claim(6, 'DL SINR dB at tau 10', /DL 降到 (23\.75) dB/, () => J.E5.d30.DL['10'], 'E5.d30.DL[10]');
claim(6, 'd_min m (second curve)', /離軌道 (5) m 時/, () => 5, 'tests/experiment_measure.js: for (const dm of [30, 5])', { codeRe: [/for \(const dm of \[30, (5)\]\)/, 'measure'] });
claim(6, 'DL SINR dB at tau 0 (5 m)', /DL 從 (26\.19) 降到/, () => J.E5.d5.DL['0'], 'E5.d5.DL[0]');
claim(6, 'DL SINR dB at tau 10 (5 m)', /降到 (17\.17) dB/, () => J.E5.d5.DL['10'], 'E5.d5.DL[10]');
claim(6, 'realisations', /E5，各 (1000) 次實現/, () => fromCode('measure', /for \(const tau of \[0, 5, 10\]\)[^\n]*?, (1000), 5000\)/), 'tests/experiment_measure.js (E5 block)');
claim(6, 'T3a theta1', /K→∞、θ₁ = (0)°、θ₂ = −30°、d_min = 5/, () => null, 'tests/t3_aging.js T3a scenario (comment: "theta_1 = 0", "theta_2 = -30", K -> infinity, d_min in {5, 30})', { fileHas: ['t3', /theta_1 = 0[\s\S]*theta_2 = -30 deg/] });
claim(6, 'T3a theta2', /K→∞、θ₁ = 0°、θ₂ = (−30)°、d_min = 5/, () => null, 'tests/t3_aging.js T3a scenario', { fileHas: ['t3', /theta_2 = -30 deg/] });
claim(6, 'T3a d_min 5', /d_min = (5) 與 30 m/, () => null, 'tests/t3_aging.js T3a: d_min in {5, 30}', { fileHas: ['t3', /d_min in \{5, 30\}/] });
claim(6, 'T3a d_min 30', /d_min = 5 與 (30) m 下/, () => null, 'tests/t3_aging.js T3a: d_min in {5, 30}', { fileHas: ['t3', /d_min in \{5, 30\}/] });
// --- sentence 7
const e6 = X.params('E6');
const grid = k => [...new Set(s1rows.map(r => r[k]))].sort((a, b) => a - b);
const listClaim = (s, label, re, getList, src) => claim(s, label, re, getList, src, { list: true });
claim(7, 'fc GHz', /在 (5) GHz\/15 kHz、d_min/, () => grid('fc_GHz')[0], 'docs/diagnostics/s1_sensitivity.csv column fc_GHz');
claim(7, 'df kHz', /GHz\/(15) kHz、d_min/, () => grid('scs_kHz')[0], 'docs/diagnostics/s1_sensitivity.csv column scs_kHz');
claim(7, 'd_min upper bound m', /d_min ≤ (500) m/, () => Math.max(...grid('d_min_m')), 's1_sensitivity.csv max d_min_m');
claim(7, 'tau upper bound ms', /τ ≤ (10) ms/, () => Math.max(...grid('tau_ms')), 's1_sensitivity.csv max tau_ms');
claim(7, 'v upper bound km/h', /v ≤ (300) km\/h/, () => Math.max(...grid('v_kmh')), 's1_sensitivity.csv max v_kmh');
claim(7, 'fc GHz (second)', /與 (28) GHz\/120 kHz（/, () => grid('fc_GHz')[1], 's1_sensitivity.csv column fc_GHz');
claim(7, 'df kHz (second)', /GHz\/(120) kHz（Δf/, () => grid('scs_kHz')[1], 's1_sensitivity.csv column scs_kHz');
claim(7, 'df kHz (second, again)', /Δf = (120) kHz）/, () => grid('scs_kHz')[1], 's1_sensitivity.csv column scs_kHz');
listClaim(7, 'd_min set m', /d_min ∈ \{([\d, ]+)\} m/, () => grid('d_min_m'), 's1_sensitivity.csv unique d_min_m');
listClaim(7, 'tau set ms', /τ ∈ \{([\d, ]+)\} ms/, () => grid('tau_ms'), 's1_sensitivity.csv unique tau_ms');
listClaim(7, 'v set km/h', /v ∈ \{([\d, ]+)\} km\/h/, () => grid('v_kmh'), 's1_sensitivity.csv unique v_kmh');
claim(7, 'gamma_rel dB', /γ_rel = \+(10) dB/, () => parseBase('tests/diag_s_sensitivity.js').gammaRelDb, 'tests/diag_s_sensitivity.js BASE.gammaRelDb');
claim(7, 'K dB', /K=(20) dB、θ₁/, () => parseBase('tests/diag_s_sensitivity.js').kDb, 'tests/diag_s_sensitivity.js BASE.kDb');
claim(7, 'theta1 deg', /θ₁=(20)°/, () => parseBase('tests/diag_s_sensitivity.js').aoaT, 'tests/diag_s_sensitivity.js BASE.aoaT');
claim(7, 'theta2 deg', /θ₂=(−30)°/, () => parseBase('tests/diag_s_sensitivity.js').aoaJ, 'tests/diag_s_sensitivity.js BASE.aoaJ');
claim(7, 'L', /L=(100)、N/, () => parseBase('tests/diag_s_sensitivity.js').L, 'tests/diag_s_sensitivity.js BASE.L');
claim(7, 'N', /N=(8)、SNR/, () => parseBase('tests/diag_s_sensitivity.js').N, 'tests/diag_s_sensitivity.js BASE.N');
claim(7, 'SNR dB', /SNR (20) dB、SIR/, () => parseBase('tests/diag_s_sensitivity.js').snr, 'tests/diag_s_sensitivity.js BASE.snr');
claim(7, 'SIR dB', /SIR (−10) dB），/, () => parseBase('tests/diag_s_sensitivity.js').sir, 'tests/diag_s_sensitivity.js BASE.sir');
const closest = s1rows.map(r => ({ r, q: r.aging_loss_dB / r.ici_loss_dB })).sort((a, b) => b.q - a.q)[0].r;
claim(7, 'closest cell: fc GHz', /最接近的情況是 (28) GHz\/120 kHz/, () => closest.fc_GHz, 's1_sensitivity.csv: row with the largest aging/ICI ratio');
claim(7, 'closest cell: df kHz', /最接近的情況是 28 GHz\/(120) kHz/, () => closest.scs_kHz, 's1_sensitivity.csv: same row');
claim(7, 'closest cell: d_min m', /kHz、(10) m、10 ms/, () => closest.d_min_m, 's1_sensitivity.csv: same row');
claim(7, 'closest cell: tau ms', /m、(10) ms、300 km\/h：/, () => closest.tau_ms, 's1_sensitivity.csv: same row');
claim(7, 'closest cell: v km/h', /ms、(300) km\/h：老化/, () => closest.v_kmh, 's1_sensitivity.csv: same row');
claim(7, 'aging loss dB', /老化 (5\.53) dB/, () => closest.aging_loss_dB, 's1_sensitivity.csv: aging_loss_dB of that row');
claim(7, 'ICI loss dB', /ICI (9\.91) dB/, () => closest.ici_loss_dB, 's1_sensitivity.csv: ici_loss_dB of that row');
claim(7, 'aging < ICI in every cell', /空間域的老化損失都小於 ICI 損失/, () => s1rows.every(r => r.aging_loss_dB < r.ici_loss_dB), 's1_sensitivity.csv: all rows aging_loss_dB < ici_loss_dB', { bool: true });
claim(7, 'cells', /S1，(240) 格/, () => s1rows.length, 's1_sensitivity.csv: number of rows');
claim(7, 'realisations per cell', /每格 (500) 次實現/, () => { const m = /(\d+) realisations per cell/.exec(s1txt); return m ? +m[1] : null; }, 's1_sensitivity.txt first line');
// --- sentence 8
const s3col = v => S3.vs.indexOf(v);
claim(8, 'K=inf SINR dB at v=0', /v=0 為 (11\.4) dB/, () => S3.rows.i[s3col(0)], 's3_withsignal.txt row (i), v=0');
claim(8, 'K=inf SINR dB at 300', /300 km\/h 為 (6\.7) dB/, () => S3.rows.i[s3col(300)], 's3_withsignal.txt row (i), v=300');
claim(8, 'K=20 SINR dB at v=0', /從 (−5\.5) dB 升到/, () => S3.rows.iii[s3col(0)], 's3_withsignal.txt row (iii), v=0');
claim(8, 'K=20 peak SINR dB', /升到 (7\.2) dB（125/, () => Math.max(...S3.rows.iii), 's3_withsignal.txt row (iii), maximum');
claim(8, 'K=20 peak speed km/h', /dB（(125) km\/h）後/, () => S3.vs[S3.rows.iii.indexOf(Math.max(...S3.rows.iii))], 's3_withsignal.txt row (iii), speed of the maximum');
claim(8, 'K=20 SINR dB at 300', /降到 (4\.6) dB（300/, () => S3.rows.iii[s3col(300)], 's3_withsignal.txt row (iii), v=300');
claim(8, 'speed at the end km/h', /dB（(300) km\/h）。（適用/, () => S3.vs[S3.vs.length - 1], 's3_withsignal.txt last speed column');
claim(8, 'fc GHz', /fc = (5) GHz/, () => e6.fc / 1e9, 'experiments.js params(E6).fc');
claim(8, 'df kHz', /Δf = (15) kHz/, () => e6.scs / 1e3, 'experiments.js params(E6).scs');
claim(8, 'd_min m', /d_min = (10) m/, () => e6.d_min, 'experiments.js params(E6).d_min');
claim(8, 'tau ms', /τ = (10) ms/, () => e6.latMs, 'experiments.js params(E6).latMs');
claim(8, 'theta1 deg', /θ₁ = (20)°/, () => e6.aoaT, 'experiments.js params(E6).aoaT');
claim(8, 'theta2 deg', /θ₂ = (−30)°/, () => e6.aoaJ, 'experiments.js params(E6).aoaJ');
claim(8, 'L', /L = (100)、N/, () => e6.L, 'experiments.js params(E6).L');
claim(8, 'N', /N = (8)、SNR/, () => e6.N, 'experiments.js params(E6).N');
claim(8, 'SNR dB', /SNR (20) dB、SIR/, () => e6.snr, 'experiments.js params(E6).snr');
claim(8, 'SIR dB', /SIR (−10) dB、withSignal/, () => e6.sir, 'experiments.js params(E6).sir');
claim(8, 'K dB (infinity)', /K→∞（(400) dB）/, () => fromCode('diagS', /\['\(i\) K -> inf[^\]]*?\{ kDb: (\d+) \}/), 'tests/diag_s_sensitivity.js s3 variant (i)');
claim(8, 'K dB (original)', /K = (20) dB 的對照/, () => fromCode('diagS', /\['\(iii\) original K = (\d+) dB'/), 'tests/diag_s_sensitivity.js s3 variant (iii)');
claim(8, 'realisations', /每點 (500) 次實現/, () => { const m = /(\d+) realisations/.exec(s3txt); return m ? +m[1] : null; }, 's3_withsignal.txt first line');
claim(8, 'realisations eigenvalue shares', /佔比為 (200) 次/, () => fromCode('diagS', /Math\.min\(trials, (\d+)\)/), 'tests/diag_s_sensitivity.js: Math.min(trials, 200)');
claim(8, 'speed range start km/h', /速度 (0)–300 km\/h/, () => S3.vs[0], 's3_withsignal.txt first speed column');
claim(8, 'speed range end km/h', /速度 0–(300) km\/h 步長/, () => S3.vs[S3.vs.length - 1], 's3_withsignal.txt last speed column');
claim(8, 'speed step km/h', /步長 (25) km\/h/, () => S3.vs[1] - S3.vs[0], 's3_withsignal.txt speed columns');
claim(8, 'K->inf does not rise', /K→∞ 時不上升/, () => Math.max(...S3.rows.i.slice(1)) <= S3.rows.i[0], 's3_withsignal.txt row (i): no value above v=0', { bool: true });
// --- sentence 9
const b4n = (k, th, v, c) => b4val(k, th, v, c);
claim(9, 'theta1 deg', /在 θ₁ = (0)°、K/, () => 0, 'docs/diagnostics/b4_afc.txt block "K = 20 dB, theta1 = 0"', { fileHasText: [b4txt, /K = 20 dB, theta1 = 0, trainMode = signalFree/] });
claim(9, 'K dB', /K = (20) dB、v/, () => 20, 'b4_afc.txt block header', { fileHasText: [b4txt, /K = 20 dB, theta1 = 0/] });
claim(9, 'v lower km/h', /v = (30)–300 km\/h/, () => Math.min(...Object.keys(b4block('20 dB', 0)).map(Number).filter(v => v > 0)), 'b4_afc.txt speed rows');
claim(9, 'v upper km/h', /v = 30–(300) km\/h/, () => Math.max(...Object.keys(b4block('20 dB', 0)).map(Number)), 'b4_afc.txt speed rows');
claim(9, 'L', /L = (100)（fc/, () => { const m = /L = (\d+), /.exec(b4txt); return m ? +m[1] : null; }, 'b4_afc.txt block header "L = 100"');
claim(9, 'fc GHz', /fc = (5) GHz/, () => CFG('fc') / 1e9, 'core.js CONFIG.fc (tests/diag_b4_afc.js does not set fc)');
claim(9, 'df kHz', /Δf = (15) kHz/, () => CFG('scs') / 1e3, 'core.js CONFIG.scs (diag_b4_afc.js does not set scs)');
claim(9, 'd_min m', /d_min = (30) m/, () => CFG('d_min'), 'core.js CONFIG.d_min (diag_b4_afc.js does not set d_min)');
claim(9, 'tau ms', /τ = (0)、θ₂/, () => fromCode('measure', /x/) === null && 0, 'tests/diag_b4_afc.js: latMs: 0', { fileHasText: [rd('tests/diag_b4_afc.js'), /latMs: 0/] });
claim(9, 'theta2 deg', /θ₂ = (40)°、SNR/, () => 40, 'tests/diag_b4_afc.js: aoaJ: 40', { fileHasText: [rd('tests/diag_b4_afc.js'), /aoaJ: 40/] });
claim(9, 'SNR dB', /SNR (20) dB、SIR/, () => 20, 'tests/diag_b4_afc.js: snr: 20', { fileHasText: [rd('tests/diag_b4_afc.js'), /snr: 20/] });
claim(9, 'SIR dB', /SIR (−10) dB）/, () => -10, 'tests/diag_b4_afc.js: sir: -10', { fileHasText: [rd('tests/diag_b4_afc.js'), /sir: -10/] });
claim(9, 'speed of the collapse km/h', /在 (30) km\/h 起崩潰/, () => Math.min(...Object.keys(b4block('20 dB', 0)).map(Number).filter(v => v > 0)), 'b4_afc.txt first non-zero speed');
claim(9, 'MMSE-P SINR dB at 30 km/h', /崩潰（(−20\.0) dB）/, () => b4n('20 dB', 0, 30, 'mmsep'), 'b4_afc.txt, K = 20 dB, theta1 = 0, v = 30, column (i)');
claim(9, 'genie SINR dB at 300', /可回到 (28\.5) dB（300/, () => b4n('20 dB', 0, 300, 'genie'), 'b4_afc.txt, K = 20 dB, theta1 = 0, v = 300, column (ii)');
claim(9, 'speed of that value km/h', /dB（(300) km\/h；30/, () => 300, 'b4_afc.txt row v = 300', { fileHasText: [b4txt, /\n300 /] });
claim(9, 'genie SINR dB at 30', /30 km\/h 為 (28\.6) dB/, () => b4n('20 dB', 0, 30, 'genie'), 'b4_afc.txt, K = 20 dB, theta1 = 0, v = 30, column (ii)');
claim(9, 'theta1 deg (second case)', /θ₁ = (45)° 時/, () => 45, 'b4_afc.txt block "theta1 = 45"', { fileHasText: [b4txt, /K = 20 dB, theta1 = 45/] });
claim(9, 'speed (second case) km/h', /在 (300) km\/h 只回到/, () => 300, 'b4_afc.txt row v = 300', { fileHasText: [b4txt, /\n300 /] });
claim(9, 'genie SINR dB at 300 (theta1 = 45)', /只回到 (14\.76) dB/, () => b4n('20 dB', 45, 300, 'genie'), 'b4_afc.txt, K = 20 dB, theta1 = 45, v = 300, column (ii)');
claim(9, 'K dB (second case)', /（K = (20) dB），此結論/, () => 20, 'b4_afc.txt block header', { fileHasText: [b4txt, /K = 20 dB, theta1 = 45/] });
claim(9, 'realisations', /每格 (500) 次實現/, () => { const m = /(\d+) trials;/.exec(b4txt); return m ? +m[1] : null; }, 'b4_afc.txt block header "500 trials"');
// --- sentence 10
claim(10, 'L', /L = (20000)，signalFree/, () => fromCode('t17', /L: (20000)/), 'tests/t17_pointing.js T17b: L: 20000');
claim(10, 'delta_theta deg (first)', /δθ = (0)° 與 1°/, () => 0, 'tests/t17_pointing.js T17b: dth in [0, 1, 3, 5]', { fileHasText: [code.t17, /for \(const dth of \[0, 1, 3, 5\]\)/] });
claim(10, 'delta_theta deg (second)', /δθ = 0° 與 (1)°/, () => 1, 'tests/t17_pointing.js T17b: dth in [0, 1, 3, 5]', { fileHasText: [code.t17, /for \(const dth of \[0, 1, 3, 5\]\)/] });
claim(10, 'RMB loss dB', /期望損失 (−0\.0015) dB/, () => suite(/expected finite-L loss at L = 20000, N = 8: .*?= (-?[\d.]+) dB/), SUITE + ': "T17b2 expected finite-L loss"', needSuite('T17', '21 min'));
claim(10, 'z at delta 0', /（z = (−2\.07)、/, () => suite(/T17b2 signalFree, delta = 0 deg:.*z = (-?[\d.]+)/), SUITE + ': "T17b2 signalFree, delta = 0 deg"', needSuite('T17', '21 min'));
claim(10, 'z at delta 1', /、(−1\.78)）。/, () => suite(/T17b2 signalFree, delta = 1 deg:.*z = (-?[\d.]+)/), SUITE + ': "T17b2 signalFree, delta = 1 deg"', needSuite('T17', '21 min'));
claim(10, 'realisations', /各格 (100) 次實現/, () => fromCode('t17', /for \(let t = 0; t < (100); t\+\+\) \{\s*\n\s*const s = mk\(\{ algo: 'SMI', L: 20000/), 'tests/t17_pointing.js T17b loop: t < 100');
claim(10, 'N', /N = (8)、K/, () => TB.t17_pointing.N, 'tests/t17_pointing.js BASE.N');
claim(10, 'K dB', /K = (20) dB、v/, () => TB.t17_pointing.kDb, 'tests/t17_pointing.js BASE.kDb');
claim(10, 'v', /v = (0)、τ/, () => TB.t17_pointing.v, 'tests/t17_pointing.js BASE.v');
claim(10, 'tau ms', /τ = (0)、θ₁/, () => TB.t17_pointing.latMs, 'tests/t17_pointing.js BASE.latMs');
claim(10, 'theta1 deg', /θ₁ = (0)°、θ₂/, () => TB.t17_pointing.aoaT, 'tests/t17_pointing.js BASE.aoaT');
claim(10, 'theta2 deg', /θ₂ = (40)°、SNR/, () => TB.t17_pointing.aoaJ, 'tests/t17_pointing.js BASE.aoaJ');
claim(10, 'SNR dB', /SNR (20) dB、SIR/, () => TB.t17_pointing.snr, 'tests/t17_pointing.js BASE.snr');
claim(10, 'SIR dB', /SIR (−10) dB、QPSK/, () => TB.t17_pointing.sir, 'tests/t17_pointing.js BASE.sir');
claim(10, 'fc GHz', /fc = (5) GHz/, () => CFG('fc') / 1e9, 'core.js CONFIG.fc');
claim(10, 'df kHz', /Δf = (15) kHz/, () => CFG('scs') / 1e3, 'core.js CONFIG.scs');
claim(10, 'd_min m', /d_min = (30) m/, () => CFG('d_min'), 'core.js CONFIG.d_min');
claim(10, 'delta_theta threshold deg', /δθ ≥ (3)° 公式/, () => 3, 'tests/t17_pointing.js T17b: the cells with dth >= 3 keep the old criterion', { fileHasText: [code.t17, /useTheory = tm === 'signalFree' && dth <= 1/] });
// --- sentence 11
claim(11, 'MUSIC/Capon target error at L >= 12 (max |x - 0.2|)', /約 (0\.2)°（受/, () => { const v = t18rows('vs L').filter(r => r.key >= 12); return Math.max(...v.map(r => r.musicT), ...v.map(r => r.capT)); }, 'fresh run of T18 (docs/diagnostics/fresh_T18_output.txt), "vs L" table, rows L = 12 ... 100, MUSIC T and Capon T; value = the largest of the 8 numbers', { within: 0.05, rangeLabel: 'all within 0.2 +- 0.05' });
claim(11, 'MUSIC/Capon target error at separation >= 20 (max)', /約 (0\.2)°（受/, () => { const v = t18rows('vs separation').filter(r => r.key >= 20); return Math.max(...v.map(r => r.musicT), ...v.map(r => r.capT)); }, 'fresh run of T18, "vs separation" table, rows 20 and 40 deg', { within: 0.05, second: true });
claim(11, 'L lower bound', /在 L ≥ (12)（θ₁/, () => 12, 'tests/t18_doa.js T18c L list [4, 8, 12, 24, 48, 100]', { fileHasText: [code.t18, /for \(const L of \[4, 8, 12, 24, 48, 100\]\)/] });
claim(11, 'theta1 deg (L sweep)', /（θ₁ = (0)°、θ₂ = 40°）與/, () => 0, 'tests/t18_doa.js T18c "vs L (theta1 = 0, theta2 = 40)"', { fileHasText: [code.t18, /vs L \(theta1 = 0, theta2 = 40\)/] });
claim(11, 'theta2 deg (L sweep)', /θ₂ = (40)°）與/, () => 40, 'tests/t18_doa.js T18c', { fileHasText: [code.t18, /aoaT: 0, aoaJ: 40, L, snr: 20, kDb: 20/] });
claim(11, 'separation lower bound deg', /源分離 ≥ (20)°/, () => 20, 'tests/t18_doa.js T18c separation list [5, 10, 20, 40]', { fileHasText: [code.t18, /for \(const t2 of \[5, 10, 20, 40\]\)/] });
claim(11, 'theta1 deg (separation sweep)', /（θ₁ = (0)°、L = 100）/, () => 0, 'tests/t18_doa.js T18c "vs separation (theta1 = 0, L = 100)"', { fileHasText: [code.t18, /vs separation \(theta1 = 0, L = 100\)/] });
claim(11, 'L (separation sweep)', /θ₁ = 0°、L = (100)）/, () => 100, 'tests/t18_doa.js T18c', { fileHasText: [code.t18, /vs separation \(theta1 = 0, L = 100\)/] });
claim(11, 'grid deg', /受 (0\.5)° 網格/, () => fromCode('core', /\(-90 \+ (0\.5) \* i\) \* Math\.PI \/ 180/), 'core.js DOA grid: (-90 + 0.5 * i) deg');
claim(11, 'realisations', /T18c，各格 (500) 次實現/, () => fromCode('t18', /trials = (\d+), seed = 1818/), 'tests/t18_doa.js default trials');
claim(11, 'N', /N = (8)、K/, () => TB.t18_doa.N, 'tests/t18_doa.js BASE.N');
claim(11, 'K dB', /K = (20) dB、SNR/, () => 20, 'tests/t18_doa.js T18c: kDb: 20', { fileHasText: [code.t18, /snr: 20, kDb: 20/] });
claim(11, 'SNR dB', /SNR (20) dB、SIR/, () => 20, 'tests/t18_doa.js T18c: snr: 20', { fileHasText: [code.t18, /snr: 20, kDb: 20/] });
claim(11, 'SIR dB', /SIR (−10) dB、v/, () => TB.t18_doa.sir, 'tests/t18_doa.js BASE.sir');
claim(11, 'v', /v = (0)、τ/, () => TB.t18_doa.v, 'tests/t18_doa.js BASE.v');
claim(11, 'tau ms', /τ = (0)、withSignal/, () => TB.t18_doa.latMs, 'tests/t18_doa.js BASE.latMs');
claim(11, 'fc GHz', /fc = (5) GHz/, () => CFG('fc') / 1e9, 'core.js CONFIG.fc');
claim(11, 'df kHz', /Δf = (15) kHz/, () => CFG('scs') / 1e3, 'core.js CONFIG.scs');
claim(11, 'd_min m', /d_min = (30) m/, () => CFG('d_min'), 'core.js CONFIG.d_min');


// --- additional claims for the numbers that appear inside words / labels
claim(1, 'E4 delta_theta point deg', /E4 的 (0)° 點/, () => null, 'data/experiment_numbers.json E4.signalFree has the key "0"', { fileHas: ['measure', /for \(const d of \[(?:0|0, 1), 3, 5\]\)|\[0, 1, 3, 5\]/] });
claim(2, 'delta_theta reference deg', /SMI 與 DL 比 (0)° 低/, () => null, 'data/experiment_numbers.json E3.*.d0', { fileHas: ['measure', /d0: sinr\(\{ algo: a, trainMode: 'withSignal', pointErrDeg: 0 \}/] });
claim(2, 'signalFree delta_theta deg', /signalFree 時 (3)° 只降約/, () => null, 'data/experiment_numbers.json E4.signalFree["3"]', { fileHas: ['measure', /for \(const d of \[0, 1, 3, 5\]\)/] });
claim(5, 'QAM order', /越過 (16)-QAM/, () => parseInt(/QAM(\d+)/.exec(XP('E1', 'mod'))[1], 10), 'experiments.js params(E1).mod = QAM16');
claim(5, 'range start km/h', /前者在 (0)–300 km\/h/, () => fromCode('measure', /for \(let v = (\d+); v <= 300; v \+= 10\)/), 'tests/experiment_measure.js: loop v = 0 ...');
claim(8, 'v km/h (K=inf at 0)', /（v=(0) 為 11\.4/, () => S3.vs[0], 's3_withsignal.txt first speed column');
claim(8, 'v km/h (K=inf at 300)', /11\.4 dB，(300) km\/h 為/, () => S3.vs[S3.vs.length - 1], 's3_withsignal.txt last speed column');
claim(8, 'K dB (second statement)', /6\.7 dB），K=(20) dB 時從/, () => fromCode('diagS', /\['\(iii\) original K = (\d+) dB'/), 'tests/diag_s_sensitivity.js s3 variant (iii)');
claim(8, 'tau of the T22c comparison ms', /T22c（τ = (0) 對照）/, () => null, 'tests/t22_mobility_pointing.js: the mobility run with latMs: 0', { fileHasText: [rd('tests/t22_mobility_pointing.js'), /latMs: 0, pointingMode: 'mobility'/] });
claim(9, 'speed of the genie value at 30 km/h', /；(30) km\/h 為 28\.6/, () => null, 'b4_afc.txt row v = 30', { fileHasText: [b4txt, /\n30 /] });
for (const sn of Object.keys(SENT).map(Number)) if (/v = (0) 時/.test(SENT[sn])) claim(sn, 'v = 0: d_min does not matter', /v = (0) 時/, () => null, 'core.js trackAngle: "if (s === 0 || vms === 0) return th" (no angle change at v = 0, so d_min has no effect)', { fileHas: ['core', /vms === 0\) return th/] });

// ------------------------------------------------------------------ evaluation
const NUMRE = /(?<![A-Za-z_Ͱ-Ͽ\d.])[−-]?\d+(?:\.\d+)?/g;
const rows = []; const covered = {};
for (const c of claims) {
    const sent = SENT[c.s], re = new RegExp(c.re.source, c.re.flags.replace('d', '') + 'd'), m = re.exec(sent);
    const row = { s: c.s, label: c.label, src: c.src };
    if (!m) { row.rep = '(pattern not found)'; row.act = '-'; row.res = 'manual'; row.note = 'the report text no longer matches the pattern of this check; look at the sentence by hand'; rows.push(row); continue; }
    if (m[1] !== undefined) { const [a, b] = m.indices[1]; (covered[c.s] = covered[c.s] || []).push([a, b]); }
    row.rep = m[1] !== undefined ? m[1].replace(/−/g, '-') : '(text)';
    let act;
    try {
        if (c.o.fileHas) act = c.o.fileHas[1].test(code[c.o.fileHas[0]]) ? 'found' : 'missing';
        else if (c.o.fileHasText) act = c.o.fileHasText[1].test(c.o.fileHasText[0]) ? 'found' : 'missing';
        else if (c.o.codeRe) { const mm = c.o.codeRe[0].exec(code[c.o.codeRe[1]]); act = mm ? num(mm[1]) : 'missing'; }
        else act = c.get();
    } catch (e) { act = 'error: ' + e.message; }
    if (act === undefined) { row.act = '-'; row.res = 'manual'; row.note = c.o.manualHint || 'needs --suite'; }
    else if (act === null || (typeof act === 'string' && /^(missing|error)/.test(act))) { row.act = String(act); row.res = 'mismatch'; row.note = 'the source value could not be found or computed'; }
    else if (c.o.bool) { row.act = String(act); row.res = act === true ? 'match' : 'mismatch'; }
    else if (act === 'found') { row.act = 'found in ' + c.src; row.res = 'match'; }
    else if (c.o.list) { const rep = m[1].split(',').map(x => +x.trim()), ok = rep.length === act.length && rep.every((x, i) => x === act[i]); row.rep = '{' + rep.join(', ') + '}'; row.act = '{' + act.join(', ') + '}'; row.res = ok ? 'match' : 'mismatch'; }
    else {
        const rv = num(m[1]), d = dec(m[1]);
        if (c.o.atLeast) { row.act = (+act).toFixed(2) + ' (>= ' + m[1] + ')'; row.res = act >= rv ? 'match' : 'mismatch'; }
        else if (c.o.within) { row.act = 'largest ' + (+act).toFixed(2); row.res = Math.abs(act - rv) <= c.o.within + 1e-9 ? 'match' : 'mismatch'; row.note = c.o.rangeLabel || ''; }
        else { const tol = (c.o.tolAbs !== undefined ? c.o.tolAbs : 0.5 * Math.pow(10, -d)) + 1e-9; row.act = Math.abs(act) < 1e-9 ? '0' : (+act).toFixed(Math.max(d, 2)); row.res = Math.abs(act - rv) <= tol ? 'match' : 'mismatch'; if (row.res === 'mismatch') row.note = 'difference ' + (act - rv).toFixed(4); }
    }
    rows.push(row);
}
// numbers of the sentences that no check covers
for (const s of Object.keys(SENT).map(Number)) {
    const sent = SENT[s], cov = covered[s] || [];
    let m; NUMRE.lastIndex = 0;
    while ((m = NUMRE.exec(sent))) {
        const a = m.index, b = a + m[0].length;
        if (cov.some(([x, y]) => a >= x && b <= y)) continue;
        // skip the sentence number itself, identifiers like fc/E4/T17b2 (excluded by the look-behind) and the digits of file names
        if (a < 3 && /^\d+\. /.test(sent)) continue;
        const ctx = sent.slice(Math.max(0, a - 12), b + 10).replace(/\n/g, ' ');
        rows.push({ s, label: 'not mapped', src: '-', rep: m[0].replace(/−/g, '-'), act: '-', res: 'manual', note: 'no automatic source; context: "' + ctx + '"' });
    }
}
rows.sort((p, q) => p.s - q.s);

// ------------------------------------------------------------------ output
const out = [];
out.push(`check_report_numbers: item 10 of docs/FINAL_REPORT.md; archived suite output: ${suitePath ? path.relative(ROOT, suitePath) + ' (complete run of 2026-10-06, 198/198; used for T17 only)' : 'not used (--no-archive): numbers that exist only there are listed as manual'}`);
out.push('tolerance: half a unit of the last digit written in the report; exact for integers; "bool" = a statement of the report that must hold for the source data');
out.push('');
const W = [3, 44, 14, 24, 9];
out.push(['s'.padEnd(W[0]), 'value'.padEnd(W[1]), 'report'.padEnd(W[2]), 'source'.padEnd(W[3]), 'result'.padEnd(W[4]), 'source / note'].join(' | '));
for (const r of rows) out.push([String(r.s).padEnd(W[0]), r.label.slice(0, W[1]).padEnd(W[1]), String(r.rep).slice(0, W[2]).padEnd(W[2]), String(r.act).slice(0, W[3]).padEnd(W[3]), ({ match: 'match', mismatch: 'MISMATCH', manual: 'manual' })[r.res].padEnd(W[4]), (r.res === 'manual' || r.res === 'mismatch' ? (r.note ? r.note + ' | ' : '') : '') + r.src].join(' | '));
const cnt = k => rows.filter(r => r.res === k).length;
out.push('');
out.push(`total ${rows.length}: match ${cnt('match')}, mismatch ${cnt('mismatch')}, manual ${cnt('manual')}`);
const man = rows.filter(r => r.res === 'manual' && r.label !== 'not mapped'), hints = [...new Set(man.map(r => r.note))];
if (man.length) { out.push('manual because of the missing suite output: ' + man.length + ' (rerun hints: ' + hints.join('; ') + ')'); }
const sRows = rows.filter(r => r.src.startsWith('terminal output of the suite'));
if (sRows.length) { out.push('items that rest on the ARCHIVED terminal output of the complete run of 2026-10-06 (198/198), file docs/diagnostics/suite_output_2026-10-06_198of198.txt, NOT re-executed: ' + sRows.length + ' (T17b2, T17c, T17d); to recompute them: rerun T17 (node tests/run_all.js, about 21 min); with --no-archive they are listed as manual. T11b and T18c are read from fresh runs: docs/diagnostics/fresh_T11_output.txt, fresh_T18_output.txt.'); }
console.log(out.join('\n'));
if (opt('--out')) fs.writeFileSync(opt('--out'), out.join('\n') + '\n');
