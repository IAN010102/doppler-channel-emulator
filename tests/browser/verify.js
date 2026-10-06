'use strict';
/**
 * Real-browser verification of the page (Commit 25).  node tests/browser/verify.js [--no-shots]
 * Needs playwright-core (npm install --no-save playwright-core) and an installed Microsoft Edge / Chrome (channel 'msedge', or set BROWSER_CHANNEL=chrome).
 * For Traditional Chinese and for English it
 *   - loads the page, switches model / training data / modulation / fc / numerology / seed / pointing error, looks at the DOA panel,
 *   - applies and runs the five experiments E1-E5 (E1 with both (fc, df) pairs), runs the velocity sweep until it finishes, exports the CSV,
 *   - collects console errors, page errors, clipped text (scrollWidth > clientWidth in label-like elements), overlapping label/value boxes, blank canvases,
 *   - repeats the default screen at viewport widths 1280 and 1024,
 *   - writes screenshots to docs/screenshots/ and a report to docs/screenshots/report.json.
 */
const http = require('http'), fs = require('fs'), path = require('path');
const root = path.join(__dirname, '..', '..');
let chromium;
try { ({ chromium } = require('playwright-core')); } catch (e) { ({ chromium } = require(require.resolve('playwright-core', { paths: [process.cwd(), root, process.env.PW_DIR || root] }))); }
const OUT = path.join(root, 'docs', 'screenshots'); fs.mkdirSync(OUT, { recursive: true });
const noShots = process.argv.includes('--no-shots');

const MIME = { '.html': 'text/html; charset=utf-8', '.js': 'text/javascript; charset=utf-8', '.json': 'application/json', '.md': 'text/plain; charset=utf-8', '.csv': 'text/csv' };
const server = http.createServer((req, res) => {
    const p = path.join(root, decodeURIComponent(req.url.split('?')[0]).replace(/^\/+/, '') || 'index.html');
    if (!p.startsWith(root) || !fs.existsSync(p) || fs.statSync(p).isDirectory()) { res.writeHead(404); res.end('not found'); return; }
    res.writeHead(200, { 'Content-Type': MIME[path.extname(p)] || 'application/octet-stream' }); fs.createReadStream(p).pipe(res);
});

const report = { generated: new Date().toISOString(), runs: {} };
const issues = [], seenUntr = new Set();
function note(lang, where, msg) { issues.push({ lang, where, msg }); console.log(`  ISSUE [${lang}] ${where}: ${msg}`); }

async function layoutChecks(page, lang, where) {
    if (lang === 'en') {                                    // English mode: no Chinese text may be left on the page
        const left = await page.evaluate(() => {
            const out = new Set(), re = /[一-鿿　-〿＀-￯]/, tw = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT); let n;
            while ((n = tw.nextNode())) { const p = n.parentNode; if (p && (p.id === 'btn-lang' || p.tagName === 'SCRIPT' || p.tagName === 'STYLE')) continue; if (re.test(n.nodeValue)) out.add(n.nodeValue.trim().slice(0, 90)); }
            document.querySelectorAll('[title]').forEach(el => { if (!el.dataset.tip && re.test(el.getAttribute('title'))) out.add('[title] ' + el.getAttribute('title').slice(0, 90)); });
            return [...out];
        });
        left.forEach(x => { if (!seenUntr.has(x)) { seenUntr.add(x); note(lang, where, 'untranslated: ' + x); } });
    }
    const r = await page.evaluate(() => {
        const vis = el => el.offsetParent !== null && el.getBoundingClientRect().width > 0;
        const clipped = [], overlap = [];
        document.querySelectorAll('.kv .k, .lab .zh, .lab .en, .btn, .p-head, .g-title, .s-top, .tx-info .k, select, .chip, .exp-btns .btn').forEach(el => {
            if (!vis(el)) return;
            const cs = getComputedStyle(el), isSelect = el.tagName === 'SELECT';
            if (el.scrollWidth > el.clientWidth + 1 && (cs.overflow === 'hidden' || cs.overflowX === 'hidden' || cs.textOverflow === 'ellipsis' || isSelect)) {
                if (isSelect) { const o = el.options[el.selectedIndex]; if (o && o.text.length < 24) return; }
                clipped.push((el.id ? '#' + el.id + ' ' : '') + (el.className || el.tagName) + ': "' + (el.innerText || el.textContent || '').trim().slice(0, 60) + '" (' + el.scrollWidth + ' > ' + el.clientWidth + ')');
            }
        });
        document.querySelectorAll('.kv').forEach(kv => {
            if (!vis(kv)) return; const k = kv.querySelector('.k'), v = kv.querySelector('.v'); if (!k || !v) return;
            const a = k.getBoundingClientRect(), b = v.getBoundingClientRect();
            if (a.right > b.left + 1 && a.top < b.bottom - 1 && a.bottom > b.top + 1) overlap.push((k.innerText || '').trim().slice(0, 40) + ' | ' + (v.innerText || '').trim().slice(0, 20));
        });
        const pageOverflowX = document.documentElement.scrollWidth > document.documentElement.clientWidth + 1;
        return { clipped, overlap, pageOverflowX };
    });
    r.clipped.forEach(x => note(lang, where, 'clipped text: ' + x)); r.overlap.forEach(x => note(lang, where, 'label/value overlap: ' + x));
    if (r.pageOverflowX) note(lang, where, 'horizontal page scroll');
    return r;
}
async function blankCanvases(page, ids) {
    return page.evaluate(ids => ids.filter(id => {
        const c = document.getElementById(id); if (!c || !c.width) return true;
        const d = c.getContext('2d').getImageData(0, 0, c.width, c.height).data; let first = [d[0], d[1], d[2]], diff = 0;
        for (let i = 0; i < d.length; i += 4 * 17) if (d[i] !== first[0] || d[i + 1] !== first[1] || d[i + 2] !== first[2]) { diff++; if (diff > 5) return false; }
        return true;
    }), ids);
}

(async () => {
    await new Promise(r => server.listen(0, r)); const port = server.address().port, base = `http://127.0.0.1:${port}/index.html`;
    const browser = await chromium.launch({ channel: process.env.BROWSER_CHANNEL || 'msedge', headless: true });
    for (const lang of (process.env.LANGS || 'zh,en').split(',')) {
        console.log(`\n=== ${lang === 'zh' ? '繁體中文' : 'English'} ===`);
        const ctx = await browser.newContext({ viewport: { width: 1500, height: 2400 }, locale: lang === 'zh' ? 'zh-TW' : 'en-US', acceptDownloads: true, deviceScaleFactor: 1 });
        await ctx.addInitScript(l => { try { localStorage.setItem('lang', l); } catch (e) {} }, lang);
        const page = await ctx.newPage(), errors = []; page.setDefaultTimeout(20000);
        page.on('console', m => { if (m.type() === 'error' && !/favicon|404/.test(m.text() + (m.location().url || ''))) errors.push('console.error: ' + m.text()); });
        page.on('pageerror', e => errors.push('pageerror: ' + e.message));
        const shot = async name => { if (!noShots) await page.screenshot({ path: path.join(OUT, `${name}-${lang}.png`), fullPage: false }); };
        const settle = ms => page.waitForTimeout(ms || 700);
        const R = report.runs[lang] = { steps: [], blank: [], csv: null };

        await page.goto(base + '?seed=12345'); await settle(2000);
        const htmlLang = await page.evaluate(() => document.documentElement.lang);
        R.steps.push(`html lang = ${htmlLang}`); if ((lang === 'zh') !== /zh/.test(htmlLang)) note(lang, 'start', 'document language is ' + htmlLang);
        await shot('00-default'); await layoutChecks(page, lang, 'default 1500');
        const b0 = await blankCanvases(page, ['cvs-beam', 'cvs-ofdm', 'cvs-const', 'cvs-matrix', 'cvs-doa']); if (b0.length) note(lang, 'default', 'blank canvas: ' + b0.join(', ')); R.blank.push(...b0);

        // ---- controls
        const sel = async (id, v) => { if (process.env.TRACE) console.log('  select', id, v); await page.selectOption('#' + id, v); await settle(500); if (lang === 'en') await layoutChecks(page, lang, `${id}=${v}`); };
        await sel('sel-model', 'legacy'); await shot('01-legacy'); await sel('sel-model', 'unified');
        await sel('sel-train', 'withSignal'); await shot('02-withsignal'); await sel('sel-train', 'signalFree');
        for (const m of ['QPSK', 'QAM64', 'QAM16']) await sel('sel-mod', m);
        await sel('sel-fc', '28'); await sel('sel-scs', '120'); await shot('03-28ghz-120khz'); await sel('sel-fc', '5'); await sel('sel-scs', '15');
        await page.fill('#num-seed', '777'); await page.press('#num-seed', 'Enter'); await page.dispatchEvent('#num-seed', 'change'); await settle(600);
        await page.fill('#num-pointErrDeg', '3'); await page.dispatchEvent('#num-pointErrDeg', 'change'); await settle(800); await shot('04-pointing-3deg');
        await page.fill('#num-pointErrDeg', '0'); await page.dispatchEvent('#num-pointErrDeg', 'change');
        for (const a of ['FOURIER', 'MMSEP', 'MMSE', 'SMI', 'DL', 'BEAMSPACE']) await sel('sel-algo', a);
        await sel('sel-algo', 'SMI'); await sel('sel-train', 'withSignal'); await settle(800); await shot('05-doa-withsignal'); await sel('sel-train', 'signalFree'); await shot('06-doa-signalfree');
        R.steps.push('controls: model, training data, modulation, fc, numerology, seed, pointing error, algorithms, DOA panel');

        // ---- experiments
        const exps = [['E0', ''], ['E1', '28'], ['E1', '5'], ['E2', ''], ['E3', ''], ['E4', ''], ['E5', '']];
        for (const [id, variant] of exps) {
            await page.click(`#exp-btns .btn[data-exp="${id}"]`); await settle(300);
            if (variant) await page.selectOption('#sel-exp-variant', variant);
            await page.click('#btn-exp-run');
            if (id === 'E1') await page.waitForFunction(() => typeof Bench !== 'undefined' && !Bench.running && Bench.result && Bench.result.v.length === Bench.result.ici.length, null, { timeout: 600000 });
            else await settle(id === 'E4' || id === 'E5' ? 6000 : (id === 'E0' ? 4000 : 1500));
            await settle(800); await shot(`exp-${id.toLowerCase()}${variant ? '-' + variant : ''}`);
            await layoutChecks(page, lang, `experiment ${id}${variant ? ' ' + variant : ''}`);
            if (id === 'E4' || id === 'E5') { const b = await blankCanvases(page, ['cvs-exp']); if (b.length) note(lang, 'experiment ' + id, 'blank result canvas'); R.blank.push(...b.map(x => id + ':' + x)); }
            if (id === 'E1') { const b = await blankCanvases(page, ['cvs-bench']); if (b.length) note(lang, 'sweep', 'blank benchmark canvas'); R.blank.push(...b); }
            R.steps.push(`experiment ${id}${variant ? ' (' + variant + ')' : ''} applied and run`);
        }
        // ---- the standard sweep with the default parameters
        await page.click('#exp-btns .btn[data-exp="E1"]'); await page.click('#btn-exp-apply'); await settle(500);
        await page.click('#btn-sweep');
        await page.waitForFunction(() => !Bench.running && Bench.result && Bench.result.v.length === Bench.result.ici.length, null, { timeout: 600000 }); await settle(800);
        await shot('07-sweep'); R.steps.push('velocity sweep finished');

        // ---- CSV
        const [dl] = await Promise.all([page.waitForEvent('download', { timeout: 20000 }), page.click('#btn-csv')]);
        const csvPath = path.join(OUT, `export-${lang}.csv`); await dl.saveAs(csvPath); const csv = fs.readFileSync(csvPath, 'utf8');
        const need = ['meta,csv_schema_version,2', 'param,model,', 'param,train_mode,', 'param,seed,', 'param,delta_theta_deg,', 'param,fc_GHz,', 'param,algorithm_legacy,', 'result,SINR_opt_dB,', 'sweep_info,train=', 'sweep,velocity_kmh,point_seed'];
        const miss = need.filter(k => !csv.includes(k)); R.csv = { bytes: csv.length, missing: miss }; if (miss.length) note(lang, 'csv', 'missing: ' + miss.join(' | '));
        fs.unlinkSync(csvPath);   // the CSV itself is not kept (it contains a timestamp); only the check result is
        R.steps.push('csv exported: ' + (miss.length ? 'MISSING ' + miss.join(',') : 'schema version 2 and all required keys present'));

        // ---- narrow viewports
        for (const w of [1280, 1024]) {
            await page.setViewportSize({ width: w, height: 2600 }); await settle(1200);
            await shot(`08-width-${w}`); await layoutChecks(page, lang, `width ${w}`);
            const b = await blankCanvases(page, ['cvs-beam', 'cvs-ofdm', 'cvs-const', 'cvs-matrix', 'cvs-doa']); if (b.length) note(lang, `width ${w}`, 'blank canvas: ' + b.join(', '));
        }
        R.errors = errors; errors.forEach(e => note(lang, 'console', e));
        await ctx.close();
    }
    await browser.close(); server.close();
    report.issues = issues; fs.writeFileSync(path.join(OUT, 'report.json'), JSON.stringify(report, null, 1));
    console.log(`\n${issues.length} issue(s). Report: docs/screenshots/report.json`);
    process.exit(issues.length ? 1 : 0);
})().catch(e => { console.error(e); try { server.close(); } catch (x) {} process.exit(2); });
