// Motion Analysis: evidence accumulation (synthetic, labelled), benchmark
// comparison against baselines, and exploratory search on the current field.
import { api, waitJob, h, s, fmt, chart, waveColor, DET, detColor } from './util.js';
import { state } from './observatory.js';

const METHOD = {
  link: { label: 'Detect + link (≥3 × 5σ)', color: '#6b7685', kind: 'baseline' },
  best_det: { label: 'Stack, best detector', color: '#8a7f6a', kind: 'baseline' },
  stack: { label: 'Stack, constant flux (KBMOD)', color: '#b3bcc8', kind: 'baseline' },
  perdet: { label: 'Per-detector stacks (χ²₆)', color: '#7aa7c9', kind: 'baseline' },
  subspace2: { label: 'Spectral subspace K=2', color: '#9fd3b8', kind: 'parallax' },
  bank: { label: 'Spectral template bank', color: '#8cb8ff', kind: 'parallax' },
  template: { label: 'True SED (oracle)', color: '#e7c27d', kind: 'reference' },
};

function stampCanvas(arr, lim = 3.5) {
  const n = arr.length, c = h('canvas', { width: n, height: n }), ctx = c.getContext('2d'), im = ctx.createImageData(n, n);
  for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) {
    const v = arr[y][x], t = Math.max(0, Math.min(1, (v + lim * 0.4) / (lim * 1.4))), g = Math.round(255 * Math.sqrt(t)), o = 4 * (y * n + x);
    im.data[o] = g; im.data[o + 1] = g; im.data[o + 2] = g + 4; im.data[o + 3] = 255;
  }
  ctx.putImageData(im, 0, 0); return c;
}

export async function renderMotion(el) {
  el.replaceChildren(
    h('h1', {}, 'Motion Analysis'),
    h('p', { class: 'lede' }, 'A faint moving source can sit below the noise in every single exposure. PARALLAX tests trajectory hypotheses directly against all exposures at once. Because each SPHEREx exposure samples a different wavelength, the source spectrum is modelled together with the motion.'),
  );
  // ------------------------------------------------------------ evidence demo
  const demo = h('section', { class: 'panel' });
  el.append(demo);
  const magSel = h('select', { class: 'field' }, ...[19.5, 20, 20.5, 21, 21.5].map(m => h('option', { value: m, selected: m === 20.5 }, `AB ${m} at 1 µm`)));
  const sedSel = h('select', { class: 'field' }, h('option', { value: 'reflected' }, 'Reflected sunlight'), h('option', { value: 'thermal' }, 'Reflected + thermal'), h('option', { value: 'browndwarf' }, 'Cold brown dwarf proxy'));
  const seedIn = h('input', { class: 'field mono', type: 'number', value: 11, min: 0, max: 9999 });
  const runBtn = h('button', { class: 'btn primary' }, 'Simulate');
  const body = h('div');
  demo.append(h('div', { class: 'row', style: 'justify-content:space-between;margin-bottom:6px' }, h('h2', { style: 'margin:0' }, 'Evidence accumulation'), h('span', { class: 'tag synth mono', style: 'font-size:10px' }, 'SIMULATED DATA')),
    h('p', { class: 'sub' }, 'Synthetic SPHEREx-like visit: the real burst cadence measured from DR1 headers (2–4 exposures 2.1 min apart, one detector pair per burst, about 20 bursts over 20 days), undersampled PSF, per-pixel Gaussian noise, 2% masked pixels.'),
    h('div', { class: 'toolbar' }, h('label', { class: 'ctl' }, h('span', {}, 'Brightness'), magSel), h('label', { class: 'ctl' }, h('span', {}, 'Spectrum'), sedSel), h('label', { class: 'ctl', style: 'min-width:90px' }, h('span', {}, 'Seed'), seedIn), runBtn), body);
  const run = async () => {
    body.replaceChildren(h('p', { class: 'note' }, 'Simulating and evaluating…'));
    const d = await api(`/api/motion/demo?mag=${magSel.value}&sed=${sedSel.value}&seed=${seedIn.value}`);
    const stamps = h('div', { class: 'stamps' });
    d.stamps.forEach(st => { const c = stampCanvas(st.data); c.title = `${DET[st.det].name} λ=${st.wave.toFixed(3)} µm, frame S/N ${st.snr_frame.toFixed(2)}`; c.style.outline = `1px solid ${waveColor(st.wave, 0.8)}`; stamps.append(c); });
    const cum = s('svg', { class: 'chart', viewBox: '0 0 560 220' });
    const n = d.t.length, idx = [...Array(n).keys()].map(i => i + 1);
    const stk = h('div', { class: 'grid2', style: 'gap:10px' },
      h('figure', { style: 'margin:0' }, (() => { const c = stampCanvas(d.stack_moving, 6); c.style.cssText = 'width:100%;image-rendering:pixelated;border-radius:3px'; return c; })(), h('figcaption', { class: 'note' }, 'Matched stack along the true trajectory')),
      h('figure', { style: 'margin:0' }, (() => { const c = stampCanvas(d.stack_static, 6); c.style.cssText = 'width:100%;image-rendering:pixelated;border-radius:3px'; return c; })(), h('figcaption', { class: 'note' }, 'Same frames, stationary hypothesis')));
    body.replaceChildren(
      h('div', { class: 'grid2' },
        h('div', {}, h('div', { class: 'note', style: 'margin-bottom:6px' }, `${n} frames, centred on the predicted position (time order, outlined by wavelength)`), stamps),
        h('div', {}, cum, h('div', { class: 'legend' }, h('span', {}, h('i', { style: 'background:#b3bcc8' }), 'constant-flux stack'), h('span', {}, h('i', { style: 'background:#8cb8ff' }), 'template bank (max)'), h('span', {}, h('i', { style: 'background:#5a6472' }), 'single-frame S/N')),
          h('div', { class: 'grid3', style: 'margin-top:12px' },
            h('div', { class: 'metric' }, h('span', {}, 'max single-frame S/N'), h('b', { class: 'big' }, fmt(d.final.max_single, 1))),
            h('div', { class: 'metric' }, h('span', {}, 'stack S/N (all frames)'), h('b', { class: 'big' }, fmt(d.final.stack, 1))),
            h('div', { class: 'metric' }, h('span', {}, 'template bank S/N'), h('b', { class: 'big' }, fmt(d.final.bank, 1)))),
          h('p', { class: 'note' }, 'These are single-trial statistics for the true trajectory. A blind search over all trajectories needs the higher, empirically calibrated threshold shown below (dashed line). Stacking does not create photons: it uses the photons in every frame under one physical hypothesis.'))),
      h('div', { style: 'margin-top:14px;max-width:420px' }, stk));
    const thr = window.__motionThr || null;
    chart(cum, { W: 560, H: 220, series: [
      { x: idx, y: d.single, color: '#5a6472', r: 1.6 },
      { x: idx, y: d.cum_stack, color: '#b3bcc8', type: 'line', w: 1.4 },
      { x: idx, y: d.cum_bank, color: '#8cb8ff', type: 'line', w: 1.8 }],
      xlab: 'frames combined (time order)', ylab: 'S/N', hline: thr });
  };
  runBtn.addEventListener('click', () => run().catch(e => body.replaceChildren(h('p', { class: 'note' }, e.message))));

  // ------------------------------------------------------------ benchmark
  const bench = h('section', { class: 'panel' }); el.append(bench);
  bench.append(h('div', { class: 'row', style: 'justify-content:space-between' }, h('h2', { style: 'margin:0' }, 'Controlled comparison'), h('span', { class: 'tag synth mono', style: 'font-size:10px' }, 'SYNTHETIC INJECTION-RECOVERY')),
    h('p', { class: 'sub' }, 'Same simulated fields, same trajectory grid, same per-field false-alarm probability for every method. Thresholds come from noise-only fields searched in full, so the look-elsewhere effect is included.'));
  const benchBody = h('div'); bench.append(benchBody);
  try {
    const rep = await api('/api/reports');
    const r = rep.motion_benchmark;
    if (!r) throw new Error('benchmark not run yet');
    window.__motionThr = r.thresholds?.bank;
    const mags = r.config.mags, svg = s('svg', { class: 'chart', viewBox: '0 0 620 260' });
    const order = ['link', 'best_det', 'stack', 'perdet', 'subspace2', 'bank', 'template'].filter(k => r.summary[k]);
    const series = order.flatMap(k => [{ x: mags, y: r.summary[k].completeness, color: METHOD[k].color, type: 'line', w: k === 'bank' ? 2.2 : 1.3, dash: k === 'template' ? '4 3' : null },
                                       { x: mags, y: r.summary[k].completeness, color: METHOD[k].color, r: 1.8 }]);
    const tbl = h('table', { class: 't' }, h('tr', {}, h('th', {}, 'Method'), h('th', { class: 'num' }, 'Threshold'), h('th', { class: 'num' }, 'm₅₀ (all)'), h('th', { class: 'num' }, 'Δm₅₀ vs stack'), ...Object.keys(r.summary.stack.by_class).map(c => h('th', { class: 'num' }, c.split(' (')[0]))));
    const base = r.summary.stack.m50;
    for (const k of order) {
      const v = r.summary[k];
      tbl.append(h('tr', {}, h('td', {}, h('i', { style: `display:inline-block;width:8px;height:8px;border-radius:2px;margin-right:6px;background:${METHOD[k].color}` }), METHOD[k].label),
        h('td', { class: 'num' }, v.threshold != null ? fmt(v.threshold, 2) : '—'), h('td', { class: 'num' }, fmt(v.m50, 2)),
        h('td', { class: 'num' }, isFinite(v.m50) && isFinite(base) ? (v.m50 - base >= 0 ? '+' : '') + fmt(v.m50 - base, 2) : '—'),
        ...Object.values(v.by_class).map(b => h('td', { class: 'num' }, fmt(b.m50, 2)))));
    }
    const th = r.theory_spectral_mismatch;
    benchBody.append(h('div', { class: 'grid2' }, h('div', {}, svg, h('div', { class: 'legend' }, ...order.map(k => h('span', {}, h('i', { style: `background:${METHOD[k].color}` }), METHOD[k].label)))), h('div', {}, h('div', { class: 'scroll' }, tbl),
      th ? h('p', { class: 'note', style: 'margin-top:10px' }, 'Theory: the constant-flux statistic loses a factor R = Σsφ / √(Σφ Σs²φ) in S/N against a matched SED, so the expected gain is Δm₅₀ = 2.5 log₁₀(1/R): ',
        ...Object.entries(th).map(([c, v]) => h('span', { class: 'mono' }, `${c.split(' (')[0]} ${fmt(v.expected_dm50_matched_minus_constant, 2)} mag; `))) : null,
      h('p', { class: 'note' }, `${r.config.n_null} null fields, ${r.config.per_mag} injections per magnitude, ${r.config.size_px}-px fields, |v| ≤ ${r.config.vmax_arcsec_day}″/day, FAP ${r.config.fap_per_field} per field. Cold brown-dwarf magnitudes are referenced to the SED normalisation at 4.6 µm. ${r.provenance === 'synthetic' ? 'All data synthetic.' : ''}`))));
    chart(svg, { W: 620, H: 260, series, xlab: 'injected AB magnitude', ylab: 'completeness', yr: [0, 1.04], hline: 0.5 });
  } catch (e) { benchBody.append(h('p', { class: 'note' }, `Benchmark unavailable: ${e.message}`)); }

  // ------------------------------------------------------------ real field
  const real = h('section', { class: 'panel' }); el.append(real);
  const out = h('div');
  const vmaxSel = h('select', { class: 'field' }, ...[1, 2, 4].map(v => h('option', { value: v, selected: v === 2 }, `|v| ≤ ${v}″/day`)));
  const go = h('button', { class: 'btn' }, 'Search current field');
  real.append(h('div', { class: 'row', style: 'justify-content:space-between' }, h('h2', { style: 'margin:0' }, 'Search a real SPHEREx field'), h('span', { class: 'tag real mono', style: 'font-size:10px' }, 'REAL DATA · EXPLORATORY')),
    h('p', { class: 'sub' }, 'Uses the field loaded in the Observatory. It takes the densest 25-day visit, fits static sources jointly with a free amplitude per frame, removes them, then runs the template-bank trajectory search on the residuals.'),
    h('div', { class: 'toolbar' }, h('label', { class: 'ctl' }, h('span', {}, 'Velocity range'), vmaxSel), go), out);
  go.addEventListener('click', async () => {
    const f = state.field; if (!f) { out.replaceChildren(h('p', { class: 'note' }, 'Load a field in the Observatory first.')); return; }
    out.replaceChildren(h('p', { class: 'note' }, 'Searching…'));
    try {
      const { job } = await api(`/api/fields/${f.id}/search?vmax=${vmaxSel.value}`, { method: 'POST' });
      const r = await waitJob(job, j => out.replaceChildren(h('p', { class: 'note mono' }, `${j.message} ${Math.round(100 * j.progress)}%`)));
      const thr = window.__motionThr;
      const tbl = h('table', { class: 't' }, h('tr', {}, h('th', {}, '#'), h('th', { class: 'num' }, 'Δα* ["]'), h('th', { class: 'num' }, 'Δδ ["]'), h('th', { class: 'num' }, 'v ["/day]'), h('th', { class: 'num' }, 'z (bank)'), h('th', {}, 'Status')));
      r.candidates.forEach((c, i) => tbl.append(h('tr', {}, h('td', {}, i + 1), h('td', { class: 'num' }, fmt(c.xi0, 1)), h('td', { class: 'num' }, fmt(c.eta0, 1)),
        h('td', { class: 'num' }, `${fmt(c.vxi, 2)}, ${fmt(c.veta, 2)}`), h('td', { class: 'num' }, fmt(c.z_bank, 2)),
        h('td', {}, h('span', { class: 'status-pill st-exp' }, thr && c.z_bank > thr && c.above_scrambled_null !== false ? 'Candidate · needs vetting' : 'Below calibrated threshold')))));
      out.replaceChildren(h('p', { class: 'note' }, `${r.n_frames} frames from the visit starting MJD ${fmt(r.visit_start_mjd, 2)} · ${r.static_sources} static sources modelled · ${r.n_trials.toLocaleString()} trajectories.`), h('div', { class: 'scroll' }, tbl),
        h('p', { class: 'note callout warn', style: 'margin-top:10px' }, r.significance_note, ' Real fields contain residuals from imperfect PSFs and saturated stars, so the synthetic thresholds are a lower bound on what a real detection needs.'));
    } catch (e) { out.replaceChildren(h('p', { class: 'note' }, e.message)); }
  });
  run().catch(e => body.replaceChildren(h('p', { class: 'note' }, e.message)));
}
