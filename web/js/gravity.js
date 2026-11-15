// Gravity Analysis: point-lens vs intrinsic-variability inference on
// multi-wavelength photometry; feasibility and degeneracy, stated plainly.
import { api, h, s, fmt, chart, waveColor } from './util.js';
import { state } from './observatory.js';

function magn(u) { return (u * u + 2) / (u * Math.sqrt(u * u + 4)); }

export async function renderGravity(el) {
  el.replaceChildren(h('h1', {}, 'Gravity Analysis'),
    h('p', { class: 'lede' }, 'An unseen mass passing in front of a background star brightens it in the same way at every wavelength. Unresolved neighbours dilute that brightening differently at each wavelength. PARALLAX fits lens, source spectrum and blend spectrum together against every measurement, and tests the fit against ordinary variability.'));
  const ctl = {
    u0: h('input', { type: 'range', min: 0.05, max: 1.5, step: 0.05, value: 0.3 }),
    tE: h('input', { type: 'range', min: 1, max: 80, step: 1, value: 15 }),
    blend: h('input', { type: 'range', min: 0, max: 0.9, step: 0.05, value: 0.5 }),
    sigma: h('select', { class: 'field' }, ...[0.005, 0.01, 0.02, 0.05].map(v => h('option', { value: v, selected: v === 0.02 }, `${v * 100}% per point`))),
    cadence: h('select', { class: 'field' }, h('option', { value: 'deep' }, 'Deep-field-like (dense)'), h('option', { value: 'wide' }, 'Wide survey (2 visits)')),
    lensed: h('select', { class: 'field' }, h('option', { value: 'true' }, 'Inject lens'), h('option', { value: 'false' }, 'No lens (null)')),
    seed: h('input', { class: 'field mono', type: 'number', value: 7 }),
  };
  const lbl = {};
  const mk = (k, title) => { lbl[k] = h('span', { class: 'mono' }); return h('label', { class: 'ctl' }, h('span', {}, title, ' ', lbl[k]), ctl[k]); };
  const run = h('button', { class: 'btn primary' }, 'Simulate & fit');
  const realBtn = h('button', { class: 'btn' }, 'Fit Observatory position');
  const top = h('section', { class: 'panel' },
    h('div', { class: 'row', style: 'justify-content:space-between' }, h('h2', { style: 'margin:0' }, 'Hypothesis test'), h('span', { id: 'gProv' })),
    h('div', { class: 'toolbar', style: 'margin-top:10px' }, mk('u0', 'u₀'), mk('tE', 't_E [d]'), mk('blend', 'blend @2µm'),
      h('label', { class: 'ctl' }, h('span', {}, 'Noise'), ctl.sigma), h('label', { class: 'ctl' }, h('span', {}, 'Cadence'), ctl.cadence),
      h('label', { class: 'ctl' }, h('span', {}, 'Truth'), ctl.lensed), h('label', { class: 'ctl', style: 'min-width:80px' }, h('span', {}, 'Seed'), ctl.seed), run, realBtn));
  const lcSvg = s('svg', { class: 'chart', viewBox: '0 0 760 280' });
  const chromSvg = s('svg', { class: 'chart', viewBox: '0 0 360 220' });
  const verdict = h('div'), table = h('div'), mass = h('div');
  top.append(lcSvg, h('div', { class: 'legend' }, h('span', {}, h('i', { style: 'background:linear-gradient(90deg,#7aa7ff,#a8d678,#ecc96e,#e8687d)' }), 'measurements coloured by λ (0.75→5 µm)'),
    h('span', {}, h('i', { style: 'background:#e7c27d' }), 'point lens + blend'), h('span', {}, h('i', { style: 'background:#c49ad9' }), 'intrinsic variability'), h('span', {}, h('i', { style: 'background:#5a6472' }), 'achromatic A(t), no blend')));
  el.append(top, h('div', { class: 'grid3', style: 'margin-top:18px' },
    h('section', { class: 'panel' }, h('h2', {}, 'Model comparison'), verdict, table),
    h('section', { class: 'panel' }, h('h2', {}, 'Chromatic dilution'), h('p', { class: 'sub' }, 'Brightening near peak vs wavelength. Lensing is achromatic, so the wavelength dependence measures the blend spectrum.'), chromSvg),
    h('section', { class: 'panel' }, h('h2', {}, 'What t_E can say about mass'), mass)));
  const feas = h('section', { class: 'panel', style: 'margin-top:18px' }); el.append(feas);
  const bench = h('section', { class: 'panel', style: 'margin-top:18px' }); el.append(bench);

  for (const k of ['u0', 'tE', 'blend']) { const f = () => (lbl[k].textContent = ctl[k].value); ctl[k].addEventListener('input', f); f(); }

  const draw = (d, isReal) => {
    $prov(isReal);
    const t0 = Math.floor(Math.min(...d.t)), tt = d.t.map(t => t - t0);
    const base = d.baseline_flux || d.H0.model_flux;
    const ratio = d.flux.map((f, i) => f / base[i]), rerr = d.err.map((e, i) => e / Math.abs(base[i]));
    const r1 = d.H1.model_flux.map((f, i) => f / base[i]), r2 = d.H2.model_flux.map((f, i) => f / base[i]);
    const p = d.H1.params, grid = [...Array(400).keys()].map(i => Math.min(...tt) + i / 399 * (Math.max(...tt) - Math.min(...tt)));
    const A = grid.map(t => magn(Math.sqrt(p.u0 ** 2 + ((t + t0 - p.t0) / p.tE) ** 2)));
    chart(lcSvg, { W: 760, H: 280, series: [
      { x: grid, y: A, color: '#5a6472', type: 'line', w: 1, dash: '3 3' },
      { x: tt, y: ratio, e: rerr, colors: d.wave.map(w => waveColor(w, 0.9)), r: 1.9 },
      { x: tt, y: r2, color: '#c49ad9', r: 1.1, opacity: 0.8 },
      { x: tt, y: r1, color: '#e7c27d', r: 1.1, opacity: 0.9 }],
      xlab: `days since MJD ${t0}`, ylab: 'flux / fitted unlensed flux (source + blend)', hline: 1 });
    const ok = d.verdict.startsWith('Microlensing'), alt = d.verdict.startsWith('Alternative') || d.verdict.startsWith('Model inadequate');
    verdict.replaceChildren(h('p', { class: 'verdict', style: `color:${ok ? 'var(--amber)' : alt ? '#c49ad9' : 'var(--text-2)'}` }, d.verdict),
      h('p', { class: 'note' }, `ΔBIC (lens − constant) ${fmt(d.dBIC_H1_vs_H0, 1)} · ΔBIC (lens − variability) ${fmt(d.dBIC_H1_vs_H2, 1)} · n = ${d.n}. Decision rule: model adequacy first (reduced χ² ≤ 3), then |ΔBIC| > 10 with a physically valid lens (0 < f_source ≤ 1, t_E ≤ 3× span). The rate of false lens verdicts is measured in the benchmark below.${d.chi2r ? ` Reduced χ²: H0 ${fmt(d.chi2r.H0, 1)}, H1 ${fmt(d.chi2r.H1, 1)}, H2 ${fmt(d.chi2r.H2, 1)}.` : ''}${isReal ? ' Real photometry includes a 1% calibration floor.' : ''}`));
    table.replaceChildren(h('table', { class: 't' }, h('tr', {}, h('th', {}, 'Model'), h('th', { class: 'num' }, 'χ²'), h('th', { class: 'num' }, 'k')),
      ...[d.H0, d.H1, d.H2].map(m => h('tr', {}, h('td', {}, m.model), h('td', { class: 'num' }, fmt(m.chi2, 1)), h('td', { class: 'num' }, m.n_par)))),
      h('dl', { class: 'kv', style: 'margin-top:10px' }, h('dt', {}, 't₀'), h('dd', {}, `MJD ${fmt(p.t0, 2)}`), h('dt', {}, 'u₀'), h('dd', {}, fmt(p.u0, 3)), h('dt', {}, 't_E'), h('dd', {}, `${fmt(p.tE, 2)} d`),
        h('dt', {}, 'A_max'), h('dd', {}, fmt(p.Amax, 3)), h('dt', {}, 'f_source @2µm'), h('dd', {}, fmt(p.source_frac_2um, 3)),
        ...(d.meta && d.meta.tE ? [h('dt', {}, 'truth'), h('dd', {}, `t₀ ${fmt(d.meta.t0, 1)} u₀ ${d.meta.u0} t_E ${d.meta.tE} blend ${d.meta.blend_frac_2um}${d.meta.lensed ? '' : ' (no lens)'}`)] : [])));
    // chromatic dilution near peak
    const near = d.t.map((t, i) => i).filter(i => Math.abs(d.t[i] - p.t0) < 0.4 * p.tE);
    chart(chromSvg, { W: 360, H: 220, xlog: true, series: [
      { x: near.map(i => d.wave[i]), y: near.map(i => ratio[i] - 1), e: near.map(i => rerr[i]), colors: near.map(i => waveColor(d.wave[i])), r: 2 },
      { x: near.map(i => d.wave[i]), y: near.map(i => r1[i] - 1), color: '#e7c27d', r: 1.2 }], xlab: 'λ [µm]', ylab: 'fractional brightening', hline: 0 });
    if (d.mass) {
      const q = d.mass.M_quantiles_Msun;
      mass.replaceChildren(h('p', { class: 'note' }, `With only t_E = ${fmt(p.tE, 1)} d, lens mass is degenerate with distance and relative proper motion. Under the stated Galactic prior the 90% interval spans`),
        h('p', { class: 'big', style: 'font-size:20px' }, `${fmt(q['5%'], 3)} – ${fmt(q['95%'], 2)} M☉`),
        h('p', { class: 'note' }, `median ${fmt(q['50%'], 2)} M☉. ${d.mass.assumptions}. Parallax or finite-source effects would be needed to narrow this. A single light curve cannot establish that the lens is a black hole.`));
    } else mass.replaceChildren(h('p', { class: 'note' }, 'Mass inference needs a lens-favoured fit.'));
  };
  const $prov = isReal => { const t = el.querySelector('#gProv'); t.replaceChildren(h('span', { class: `tag ${isReal ? 'real' : 'synth'} mono`, style: 'font-size:10px' }, isReal ? 'REAL SPHEREx PHOTOMETRY' : 'SIMULATED EVENT')); };
  const go = async () => {
    verdict.replaceChildren(h('p', { class: 'note' }, 'Fitting…'));
    try {
      const d = await api(`/api/gravity/demo?u0=${ctl.u0.value}&tE=${ctl.tE.value}&blend=${ctl.blend.value}&cadence=${ctl.cadence.value}&lensed=${ctl.lensed.value}&seed=${ctl.seed.value}&sigma=${ctl.sigma.value}`);
      draw(d, false);
    } catch (e) { verdict.replaceChildren(h('p', { class: 'note' }, e.message)); }
  };
  run.addEventListener('click', go);
  realBtn.addEventListener('click', async () => {
    const f = state.field; if (!f) { verdict.replaceChildren(h('p', { class: 'note' }, 'Load a field in the Observatory and click a source first.')); return; }
    const sel = state.viewer.sel || { xi: 0, eta: 0 };
    verdict.replaceChildren(h('p', { class: 'note' }, 'Fitting real photometry…'));
    try {
      const follow = !!f.track; let xi = sel.xi, eta = sel.eta;
      if (follow) { xi -= f.track.xi[state.viewer.i]; eta -= f.track.eta[state.viewer.i]; }
      const d = await api(`/api/fields/${f.id}/gravity?xi=${xi}&eta=${eta}&follow=${follow}`); d.mass = null; draw(d, true);
    } catch (e) { verdict.replaceChildren(h('p', { class: 'note' }, e.message)); }
  });
  go();

  // feasibility + benchmark reports
  try {
    const rep = await api('/api/reports');
    const fz = rep.gravity_feasibility;
    const rng = (a, d = 2) => a ? `${fmt(a[0], d)}–${fmt(a[1], d)}` : '—';
    if (fz) {
      feas.append(h('div', { class: 'row', style: 'justify-content:space-between' }, h('h2', { style: 'margin:0' }, 'Can SPHEREx do this?'), h('span', { class: 'tag proc mono', style: 'font-size:10px' }, 'MODEL ESTIMATE')),
        h('p', { class: 'sub' }, fz.summary),
        h('div', { class: 'scroll' }, h('table', { class: 't' }, h('tr', {}, h('th', {}, 'Scenario'), h('th', { class: 'num' }, 'Stars N'), h('th', { class: 'num' }, 'τ'), h('th', { class: 'num' }, '⟨t_E⟩'), h('th', { class: 'num' }, 'Events / yr (range)'), h('th', { class: 'num' }, 'ε (E3)'), h('th', { class: 'num' }, 'Detectable / yr (range)')),
          ...fz.scenarios.map(sc => h('tr', {}, h('td', {}, sc.name), h('td', { class: 'num' }, sc.n_stars.toExponential(0)), h('td', { class: 'num' }, sc.tau.toExponential(0)), h('td', { class: 'num' }, `${sc.tE_days} d`),
            h('td', { class: 'num' }, `${fmt(sc.events_per_yr, 2)} (${rng(sc.events_range)})`), h('td', { class: 'num' }, fmt(sc.efficiency, 2)), h('td', { class: 'num' }, `${fmt(sc.detectable_per_yr, 2)} (${rng(sc.detectable_range)})`))))),
        h('p', { class: 'note', style: 'margin-top:8px' }, fz.caveats));
    } else feas.append(h('p', { class: 'note' }, 'Feasibility report not generated yet.'));
    const gbs = ['deep', 'wide'].map(c => rep[`gravity_benchmark_${c}`]).filter(Boolean);
    if (gbs.length) bench.append(h('div', { class: 'row', style: 'justify-content:space-between' }, h('h2', { style: 'margin:0' }, 'Joint chromatic fit vs grey and single-band fits'), h('span', { class: 'tag synth mono', style: 'font-size:10px' }, 'SYNTHETIC · E3')));
    for (const gb of gbs) {
      const cl = gb.classification;
      bench.append(h('h3', { style: 'margin-top:14px' }, gb.cadence === 'deep' ? 'Deep-field-like cadence' : 'Wide-survey cadence (two visits)'), h('p', { class: 'sub' }, gb.description),
        h('div', { class: 'scroll' }, h('table', { class: 't' }, h('tr', {}, h('th', {}, 'Method'), h('th', { class: 'num' }, `Δχ² threshold (FAP ${gb.fap})`), ...gb.u0_values.map(u => h('th', { class: 'num' }, `ε u₀=${u}`))),
          ...Object.values(gb.methods).map(m => h('tr', {}, h('td', {}, m.label), h('td', { class: 'num' }, fmt(m.threshold, 1)), ...m.efficiency.map(e => h('td', { class: 'num' }, fmt(e, 2))))))),
        h('p', { class: 'note', style: 'margin-top:8px' }, `${gb.conclusion} Verdicts use the same guarded rule as the app. Lens injections classified lens / ambiguous / variability / none / inadequate: ${cl.lens_injections.lens} / ${cl.lens_injections.ambiguous} / ${cl.lens_injections.var} / ${cl.lens_injections.none} / ${cl.lens_injections.inadequate ?? 0}. Symmetric variables: ${cl.variable_injections.lens} / ${cl.variable_injections.ambiguous} / ${cl.variable_injections.var} / ${cl.variable_injections.none} / ${cl.variable_injections.inadequate ?? 0}.`));
    }
  } catch (_) {}
}
