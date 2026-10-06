'use strict';
/**
 * Deployment check in a real (headless) browser.   node tests/browser/deploy_check.js <url> [--shot <png>] [--json <file>]
 * <url> is a file:// URL of the local index.html or the live address. Needs playwright-core and Microsoft Edge / Chrome (BROWSER_CHANNEL), PW_DIR for the module path.
 * For Traditional Chinese and for English it loads the page fresh (the language is chosen through localStorage, as the button does), and checks:
 *   - no console error and no page error (also while the experiments are applied),
 *   - the defaults: model = unified, trainMode = signalFree, pointingMode = manual,
 *   - the experiments E0 ... E6 can be applied: every parameter of Experiments.params(id) is on Sys afterwards (the apply button of the page),
 *   - the CSV text of the page has csv_schema_version = 2 (Exporter.buildCsv),
 *   - ?model=legacy loads the legacy model without errors.
 * Read-only: nothing is written except the optional screenshot / json report.
 */
const fs = require('fs'), path = require('path');
const root = path.join(__dirname, '..', '..');
let chromium;
try { ({ chromium } = require('playwright-core')); } catch (e) { ({ chromium } = require(require.resolve('playwright-core', { paths: [process.cwd(), root, process.env.PW_DIR || root] }))); }
const argv = process.argv.slice(2), url = argv.find(a => !a.startsWith('--') && argv[argv.indexOf(a) - 1] !== '--shot' && argv[argv.indexOf(a) - 1] !== '--json');
const opt = n => { const i = argv.indexOf(n); return i >= 0 ? argv[i + 1] : null; };
if (!url) { console.error('usage: node tests/browser/deploy_check.js <url> [--shot png] [--json file]'); process.exit(2); }
const EXPS = ['E0', 'E1', 'E2', 'E3', 'E4', 'E5', 'E6'];
const results = []; let fail = 0;
const check = (name, ok, detail = '') => { results.push({ name, ok, detail }); if (!ok) fail++; console.log(`${ok ? 'PASS' : 'FAIL'}  ${name}${detail ? '   ' + detail : ''}`); };

(async () => {
    const browser = await chromium.launch({ channel: process.env.BROWSER_CHANNEL || 'msedge', headless: true });
    for (const lang of ['zh', 'en']) {
        const ctx = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
        await ctx.addInitScript(l => { try { localStorage.setItem('lang', l); } catch (e) { /* ignore */ } }, lang);
        const page = await ctx.newPage(), errors = [];
        page.on('console', m => { if (m.type() === 'error') errors.push('console: ' + m.text()); });
        page.on('pageerror', e => errors.push('pageerror: ' + e.message));
        page.on('requestfailed', r => errors.push('request failed: ' + r.url()));
        page.on('response', r => { if (r.status() >= 400) errors.push('HTTP ' + r.status() + ': ' + r.url()); });
        await page.goto(url, { waitUntil: 'load' });
        await page.waitForTimeout(2500);
        const L = `[${lang}]`;
        const cur = await page.evaluate(() => Lang.cur);
        check(`${L} language`, cur === lang, `Lang.cur = ${cur}`);
        const def = await page.evaluate(() => ({ model: Sys.model, train: Sys.trainMode, pm: Sys.pointingMode, algo: Sys.algo, title: document.title }));
        check(`${L} defaults model = unified, trainMode = signalFree, pointingMode = manual`, def.model === 'unified' && def.train === 'signalFree' && def.pm === 'manual', JSON.stringify(def));
        for (const id of EXPS) {
            const r = await page.evaluate(async id => {
                const e = Experiments.byId(id); ExpUI.select(id);
                document.getElementById('btn-exp-apply').click();
                await new Promise(res => setTimeout(res, 300));
                const p = Experiments.params(id, e.defaultVariant), bad = [];
                for (const k of Object.keys(p)) { const want = p[k], got = Sys[k]; if (got === undefined) continue; if (typeof want === 'number' ? Math.abs(got - want) > 1e-9 * Math.max(1, Math.abs(want)) : got !== want) bad.push(`${k}: ${got} (wanted ${want})`); }
                return { bad, n: Object.keys(p).length };
            }, id);
            check(`${L} experiment ${id} can be applied`, r.bad.length === 0, r.bad.length ? r.bad.join('; ') : `${r.n} parameters set`);
        }
        // E0: L = 1000 is applied beyond the slider range; the slider keeps working; the next preset restores its own L; the live computation stays fast
        const e0 = await page.evaluate(async () => {
            const wait = ms => new Promise(r => setTimeout(r, ms)), out = {};
            ExpUI.select('E0'); document.getElementById('btn-exp-apply').click(); await wait(400);
            out.L = Sys.L; out.numL = document.getElementById('num-L').value; out.win = document.getElementById('m-win').innerText;
            const t0 = performance.now(); for (let i = 0; i < 5; i++) { Sys.snaps = []; Sys.computeMath(); } out.msPerCompute1000 = (performance.now() - t0) / 5;
            const rng = document.getElementById('rng-L'); rng.value = 150; rng.dispatchEvent(new Event('input', { bubbles: true })); await wait(300);
            out.afterSlider = Sys.L; out.numAfterSlider = document.getElementById('num-L').value;
            ExpUI.select('E1'); document.getElementById('btn-exp-apply').click(); await wait(300); out.afterE1 = Sys.L;
            ExpUI.select('E0'); document.getElementById('btn-exp-apply').click(); await wait(300); out.again = Sys.L;
            ExpUI.select('E6'); document.getElementById('btn-exp-apply').click(); await wait(300); out.afterE6 = Sys.L;
            const t1 = performance.now(); for (let i = 0; i < 5; i++) { Sys.snaps = []; Sys.computeMath(); } out.msPerCompute100 = (performance.now() - t1) / 5;
            return out;
        });
        check(`${L} E0: Sys.L = 1000 and the L read-out shows 1000`, e0.L === 1000 && String(e0.numL) === '1000' && /1000/.test(e0.win), `Sys.L ${e0.L}, input ${e0.numL}, window "${e0.win}"`);
        check(`${L} E0: the L slider still works afterwards`, e0.afterSlider === 150 && String(e0.numAfterSlider) === '150', `Sys.L ${e0.afterSlider}, input ${e0.numAfterSlider}`);
        check(`${L} E1 and E6 restore their own L (100); E0 again gives 1000`, e0.afterE1 === 100 && e0.afterE6 === 100 && e0.again === 1000, `E1 ${e0.afterE1}, E0 again ${e0.again}, E6 ${e0.afterE6}`);
        check(`${L} L = 1000: one computation is fast (< 100 ms)`, e0.msPerCompute1000 < 100, `${e0.msPerCompute1000.toFixed(1)} ms per computation at L = 1000 against ${e0.msPerCompute100.toFixed(1)} ms at L = 100`);
        const csv = await page.evaluate(() => Exporter.buildCsv());
        const m = /(^|\n)meta,csv_schema_version,(\d+),/.exec(csv);
        check(`${L} CSV schema version`, !!m && m[2] === '2', m ? 'csv_schema_version = ' + m[2] : 'not found');
        check(`${L} no console / page errors`, errors.length === 0, errors.slice(0, 5).join(' | '));
        if (lang === 'en' && opt('--shot')) { await page.evaluate(() => { ExpUI.select('E0'); document.getElementById('btn-exp-run').click(); }); await page.waitForTimeout(5000); await page.screenshot({ path: opt('--shot'), fullPage: false }); }
        await ctx.close();
    }
    // ?model=legacy
    { const ctx = await browser.newContext({ viewport: { width: 1600, height: 1000 } }), page = await ctx.newPage(), errors = [];
      page.on('console', m => { if (m.type() === 'error') errors.push(m.text()); }); page.on('pageerror', e => errors.push(e.message));
      await page.goto(url + (url.includes('?') ? '&' : '?') + 'model=legacy', { waitUntil: 'load' }); await page.waitForTimeout(2500);
      const m = await page.evaluate(() => ({ model: Sys.model, tag: (document.getElementById('model-tag') || {}).textContent }));
      check('?model=legacy selects the legacy model', m.model === 'legacy', JSON.stringify(m)); check('?model=legacy: no console / page errors', errors.length === 0, errors.slice(0, 3).join(' | ')); await ctx.close(); }
    await browser.close();
    console.log(`\n${results.length - fail}/${results.length} checks passed, ${fail} failed   (${url})`);
    if (opt('--json')) fs.writeFileSync(opt('--json'), JSON.stringify({ url, when: new Date().toISOString(), results }, null, 1));
    process.exit(fail ? 1 : 0);
})().catch(e => { console.error(e); process.exit(2); });
