// Results: what was measured, at which evidence level, against which baseline.
import { api, h, s, fmt, chart } from './util.js';

const LEVELS = {
  impl: ['Implemented', 'st-impl'], test: ['Tested', 'st-test'], syn: ['Validated · synthetic', 'st-syn'],
  real: ['Validated · real', 'st-real'], imp: ['Improvement demonstrated', 'st-imp'], exp: ['Exploratory', 'st-exp'],
  neg: ['Negative result', 'st-exp'],
};
const pill = k => h('span', { class: `status-pill ${LEVELS[k][1]}` }, LEVELS[k][0]);

function matrix(rep) {
  const mb = rep.motion_benchmark, rv = rep.real_validation;
  const bd = mb && mb.summary.bank.by_class['cold brown dwarf proxy'].m50 - mb.summary.stack.by_class['cold brown dwarf proxy'].m50;
  const rows = [
    ['SPHEREx L2 access (AWS mirror, HTTP range reads, RICE flags, ePSF, WCS-WAVE)', 'real', 'Cutouts of real L2 images with per-pixel wavelength, flags and variance; IRSA SIA used when reachable.'],
    ['Coverage discovery without IRSA (solar-elongation prior + header index)', 'real', 'Finds all exposures covering a position by exact WCS containment.'],
    ['Time-ordered viewer, epoch comparison, difference, co-add', 'real', 'Observatory view; every frame labelled with λ, detector, UTC, observation ID.'],
    ['Proper-motion recovery of a known star from SPHEREx centroids', rv && rv.barnard ? 'real' : 'test', rv && rv.barnard ? `Barnard's Star: motion detected from ${rv.barnard.fit.n} SPHEREx centroids; each visit within 1.2″ of the Gaia-predicted position. Formal errors miss ~0.7″ per-visit systematics (details in E2 below).` : 'Pipeline implemented; real-data run pending.'],
    ['Joint spectro-orbital matched filter (template bank)', mb ? 'imp' : 'syn', mb ? `Synthetic, null-calibrated: m50 gain over constant-flux stacking ${fmt(mb.summary.bank.m50 - mb.summary.stack.m50, 2)} mag overall, ${fmt(bd, 2)} mag for cold brown-dwarf SEDs.` : ''],
    ['Exact branch-and-bound over velocity space', 'neg', 'Exactness verified (70/70 vs exhaustive); bound too loose to prune at useful thresholds. No speed-up claimed.'],
    ['Search on real SPHEREx fields', 'exp', rv && rv.barnard_search ? `Runs end-to-end on real frames; limited by bright-star residuals (top statistic ${fmt(rv.barnard_search.candidates[0]?.z_bank, 0)} vs ${fmt(rv.barnard_search.scrambled_null_max, 0)} on time-scrambled frames). No candidate.` : 'Runs end-to-end; outputs are unvetted candidates.'],
    ['Gravity mode: PSPL + chromatic blend, model comparison', rep.gravity_benchmark_deep ? 'syn' : 'test', rep.gravity_benchmark_deep ? (() => { const c = rep.gravity_benchmark_deep.classification, b = rep.gravity_benchmark_deep.classification_bic_only; return `No detection gain over a grey multi-band fit; better timescale recovery; strong gain over single-band. Deep cadence: ${c.variable_injections.lens}/100 symmetric variables and ${c.null_injections.lens}/300 nulls given a lens verdict${b ? ` (BIC alone: ${b.variable_injections.lens}/100)` : ''}; ${c.lens_injections.lens} true lenses so labelled (E3 below).`; })() : 'Nesting and recovery unit-tested.'],
    ['Microlensing yield for SPHEREx', 'exp', 'Order-of-magnitude interval from literature optical depths and measured efficiency; not a detection.'],
    ['Real microlensing / new moving-object discovery', 'exp', 'None claimed.'],
  ];
  return h('div', { class: 'scroll' }, h('table', { class: 't matrix' }, h('tr', {}, h('th', {}, 'Component'), h('th', {}, 'Evidence'), h('th', {}, 'Basis')),
    ...rows.map(([a, k, c]) => h('tr', {}, h('td', { style: 'white-space:normal' }, a), h('td', {}, pill(k)), h('td', {}, c)))));
}

export async function renderResults(el) {
  el.replaceChildren(h('h1', {}, 'Results'), h('p', { class: 'lede' }, 'Every number on this page comes from a script in the repository with fixed seeds. Synthetic experiments are marked as such; nothing here is a discovery.'));
  let rep = {};
  try { rep = await api('/api/reports'); } catch (e) { el.append(h('p', { class: 'note' }, e.message)); return; }
  el.append(h('section', { class: 'panel' }, h('h2', {}, 'Evidence matrix'), matrix(rep)));

  const mb = rep.motion_benchmark;
  if (mb) {
    const svg = s('svg', { class: 'chart', viewBox: '0 0 760 300' });
    const col = { stack: '#7d8794', perdet: '#5f8fbf', bank: '#e7c27d', template: '#86cfa6', link: '#c4746e', best_det: '#4a525e', subspace2: '#b59ad9' };
    const lab = { stack: 'Constant-flux stack (KBMOD-like)', perdet: 'Per-detector χ²', bank: 'Spectral template bank (PARALLAX)', template: 'True-SED oracle', link: 'Detect-then-link (≥3 exp. at 5σ)', best_det: 'Best single detector', subspace2: 'Spectral subspace K=2' };
    const keys = ['link', 'best_det', 'stack', 'perdet', 'subspace2', 'bank', 'template'];
    const mags = mb.config.mags;
    chart(svg, { W: 760, H: 300, xlab: 'AB magnitude at 1 µm', ylab: 'completeness', yr: [0, 1.02],
      series: keys.flatMap(k => [{ x: mags, y: mb.summary[k].completeness, color: col[k], type: 'line', w: k === 'bank' ? 2.2 : 1.3, dash: k === 'template' ? '4 3' : null }, { x: mags, y: mb.summary[k].completeness, color: col[k], r: k === 'bank' ? 2.4 : 1.6 }]), hline: 0.5 });
    const cls = Object.keys(mb.summary.stack.by_class);
    el.append(h('section', { class: 'panel', style: 'margin-top:18px' },
      h('div', { class: 'row', style: 'justify-content:space-between' }, h('h2', { style: 'margin:0' }, 'E1 · Faint moving-object completeness'), h('span', { class: 'tag synth mono', style: 'font-size:10px' }, 'SYNTHETIC')),
      h('p', { class: 'sub' }, `${mb.config.per_mag} injections per magnitude, ${mb.config.n_null} null fields, real SPHEREx burst cadence, FAP ${mb.config.fap_per_field} per field from the null maximum of each statistic (look-elsewhere included).`),
      svg, h('div', { class: 'legend' }, ...keys.map(k => h('span', {}, h('i', { style: `background:${col[k]}` }), lab[k]))),
      h('div', { class: 'scroll', style: 'margin-top:12px' }, h('table', { class: 't' },
        h('tr', {}, h('th', {}, 'Method'), h('th', { class: 'num' }, 'Threshold'), h('th', { class: 'num' }, 'm50 all'), ...cls.map(c => h('th', { class: 'num' }, `m50 ${c.split(' (')[0]}`))),
        ...keys.map(k => h('tr', {}, h('td', {}, lab[k]), h('td', { class: 'num' }, fmt(mb.summary[k].threshold, 2)), h('td', { class: 'num' }, fmt(mb.summary[k].m50, 2)), ...cls.map(c => h('td', { class: 'num' }, fmt(mb.summary[k].by_class[c].m50, 2))))))),
      h('p', { class: 'note', style: 'margin-top:10px' }, 'Expected gain from the SNR ratio R = Σsφ / √(Σφ·Σs²φ) between the constant-flux and matched statistics: ',
        ...Object.entries(mb.theory_spectral_mismatch).map(([c, t]) => `${c.split(' (')[0]} ${fmt(t.expected_dm50_matched_minus_constant, 2)} mag; `),
        'measured gains follow this ordering. Solar-like reflectance spectra gain little because a flat spectrum is already nearly matched; the large gain is for objects whose spectrum differs strongly from flat, such as cold brown dwarfs.')));
  }

  for (const c of ['deep', 'wide']) {
    const gb = rep[`gravity_benchmark_${c}`]; if (!gb) continue;
    const svg = s('svg', { class: 'chart', viewBox: '0 0 560 240' });
    const col = { joint: '#e7c27d', grey: '#7d8794', single: '#5f8fbf' };
    chart(svg, { W: 560, H: 240, xlab: 'impact parameter u₀', ylab: 'detection efficiency', yr: [0, 1.02],
      series: Object.entries(gb.methods).flatMap(([k, m]) => [{ x: gb.u0_values, y: m.efficiency, color: col[k], type: 'line', w: k === 'joint' ? 2 : 1.3 }, { x: gb.u0_values, y: m.efficiency, color: col[k], r: 2 }]) });
    el.append(h('section', { class: 'panel', style: 'margin-top:18px' },
      h('div', { class: 'row', style: 'justify-content:space-between' }, h('h2', { style: 'margin:0' }, `E3 · Gravity detection, ${c} cadence`), h('span', { class: 'tag synth mono', style: 'font-size:10px' }, 'SYNTHETIC')),
      h('p', { class: 'sub' }, gb.description), svg,
      h('div', { class: 'legend' }, ...Object.entries(gb.methods).map(([k, m]) => h('span', {}, h('i', { style: `background:${col[k]}` }), `${m.label} (Δχ² ≥ ${fmt(m.threshold, 1)})`))),
      h('p', { class: 'note', style: 'margin-top:8px' }, gb.conclusion)));
  }

  const rv = rep.real_validation;
  if (rv) {
    el.append(h('section', { class: 'panel', style: 'margin-top:18px' },
      h('div', { class: 'row', style: 'justify-content:space-between' }, h('h2', { style: 'margin:0' }, 'E2 · Real SPHEREx data checks'), h('span', { class: 'tag real mono', style: 'font-size:10px' }, 'REAL DATA')),
      ...Object.entries(rv).filter(([k]) => typeof rv[k] === 'object').map(([k, v]) => h('div', { style: 'margin-top:10px' }, h('h3', {}, v.title || k), h('p', { class: 'note' }, v.summary),
        v.table ? h('div', { class: 'scroll' }, h('table', { class: 't' }, h('tr', {}, ...v.table.columns.map(c => h('th', {}, c))), ...v.table.rows.map(r => h('tr', {}, ...r.map(x => h('td', { class: typeof x === 'number' ? 'num' : '' }, typeof x === 'number' ? fmt(x, 3) : x)))))) : ''))));
  }
}
