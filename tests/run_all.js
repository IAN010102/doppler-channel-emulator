'use strict';
/**
 * node tests/run_all.js [--quick]
 *
 * T0  legacy equivalence (core.js vs the pre-refactor page, bit for bit) and UI smoke test
 * T1  regression unified vs legacy (T1a check, T1b informational)     T2  single-path invariance     T3  aging (signalFree / MPDR)
 * T4  diffuse-path ICI (H1)     T5  time-domain OFDM chain     T6  gamma_rel sweep (informational)     T7  straight-track geometry
 * Informational items are listed separately and never count as PASS/FAIL.
 *
 * Every check prints PASS/FAIL with its value and tolerance. Tolerances are fixed in the test files and are never
 * adjusted to make a check pass; failures carry the reason. --quick lowers the trial counts (for development only).
 */
const { spawnSync } = require('child_process');
const path = require('path');
const quick = process.argv.includes('--quick');

const results = [];
function banner(s) { console.log('\n' + '='.repeat(100) + '\n' + s + '\n' + '='.repeat(100)); }

function script(name, label) {
    banner(label);
    const t0 = Date.now();
    const r = spawnSync(process.execPath, [path.join(__dirname, name)], { encoding: 'utf8', maxBuffer: 1 << 26 });
    const lines = (r.stdout || '').trim().split('\n');
    console.log(lines.slice(-3).join('\n') + (r.stderr ? '\n' + r.stderr.trim().split('\n').slice(0, 5).join('\n') : ''));
    results.push({ id: label.split(' ')[0], title: label, checks: [{ name: label, value: `exit ${r.status}`, tol: 'exit 0', pass: r.status === 0, note: '' }], ms: Date.now() - t0 });
}

(async () => {
    script('legacy_equivalence.js', 'T0a core.js reproduces the pre-refactor page bit for bit (seeded)');
    script('ui_smoke.js', 'T0b UI smoke test (stubbed DOM)');
    const opts = quick ? { trials: 200 } : {};
    for (const f of ['t1_regression.js', 't2_single_path.js', 't3_aging.js', 't4_diffuse_ici.js', 't5_ofdm_time_domain.js', 't6_gamma_rel.js', 't7_geometry.js']) {
        const mod = require('./' + f);
        banner(`${mod.id}  ${mod.title}`);
        const t0 = Date.now();
        const r = await mod.run(/^t[457]/.test(f) ? {} : opts);
        r.ms = Date.now() - t0; results.push(r);
    }

    banner('SUMMARY');
    let fails = 0, total = 0;
    for (const r of results) {
        console.log(`\n${r.id}  ${r.title}   (${(r.ms / 1000).toFixed(1)} s)`);
        for (const c of r.checks) {
            total++; if (!c.pass) fails++;
            console.log(`  ${c.pass ? 'PASS' : 'FAIL'}  ${c.name}   value: ${c.value}   tolerance: ${c.tol}`);
            if (!c.pass && c.note) console.log(`        reason: ${c.note}`);
        }
    }
    console.log('\nINFORMATIONAL (not counted)');
    for (const r of results) for (const c of (r.info || [])) console.log(`  ${c.name}   ${c.value}${c.note ? '   (' + c.note + ')' : ''}`);
    console.log(`\n${total - fails}/${total} checks passed, ${fails} failed${quick ? '   (--quick: reduced trial counts)' : ''}`);
    process.exit(fails ? 1 : 0);
})();
