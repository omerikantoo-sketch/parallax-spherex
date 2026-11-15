// Observatory workspace: field loading, timeline, playback, inspector, spectra.
import { api, waitJob, $, $$, h, s, fmt, sexRA, sexDec, isoShort, mjdToDate, waveColor, DET, detColor, chart } from './util.js';
import { Viewer } from './viewer.js';

export const state = { field: null, viewer: null, dets: new Set([1, 2, 3, 4, 5, 6]), visible: [], playing: false, timer: null, spec: null, specView: 'wave', target: null };

export function initObservatory(targets) {
  const v = new Viewer($('#stage'), $('#imgCanvas'), $('#ovCanvas'), $('#readout'));
  state.viewer = v;
  v.onSelect = sel => {
    // follow the target's proper motion only when the click is on the target itself
    const f = state.field;
    if (f && f.track) {
      const d = Math.hypot(sel.xi - f.track.xi[v.i], sel.eta - f.track.eta[v.i]);
      $('#followPm').checked = d < 12;
    }
    loadSpectrum(sel);
  };

  // target quick-picks
  const list = $('#targetList');
  for (const t of targets) {
    list.append(h('button', { type: 'button', title: `${t.kind}\n${t.note}`, onclick: () => openField(t.ra, t.dec, t.slug) },
      h('span', {}, t.name), h('span', { class: 'mono' }, `${t.ra.toFixed(3)} ${t.dec >= 0 ? '+' : ''}${t.dec.toFixed(3)}`)));
  }
  // detector chips
  const chips = $('#detChips');
  for (const d of [1, 2, 3, 4, 5, 6]) {
    const b = h('button', { class: 'chip on', type: 'button', 'aria-pressed': 'true', title: `${DET[d].lo}–${DET[d].hi} µm` },
      h('i', { style: `background:${detColor(d)}` }), `${DET[d].name} ${DET[d].lo}–${DET[d].hi}`);
    b.addEventListener('click', () => {
      state.dets.has(d) ? state.dets.delete(d) : state.dets.add(d);
      b.classList.toggle('on', state.dets.has(d)); b.setAttribute('aria-pressed', state.dets.has(d));
      applyFilter();
    });
    chips.append(b);
  }
  $('#coordForm').addEventListener('submit', async e => {
    e.preventDefault();
    const q = $('#coordInput').value.trim(); if (!q) return;
    try { const r = await api(`/api/resolve?q=${encodeURIComponent(q)}`); openField(r.ra, r.dec, r.target); }
    catch (err) { showEmpty('Could not read that position', err.message); }
  });
  // display controls
  const seg = (id, key, attr) => $$(`#${id} button`).forEach(b => b.addEventListener('click', () => {
    $$(`#${id} button`).forEach(x => x.classList.toggle('on', x === b)); v.opts[key] = b.dataset[attr]; v._coadd = null; v._global = null; v.draw(); updateBadge();
  }));
  seg('modeSeg', 'mode', 'm'); seg('scaleSeg', 'scale', 's');
  $('#stretchSel').addEventListener('change', e => { v.opts.stretch = e.target.value; v.draw(); });
  $('#loSl').addEventListener('input', e => { v.opts.lo = +e.target.value; v._global = null; v.draw(); });
  $('#hiSl').addEventListener('input', e => { v.opts.hi = +e.target.value; v._global = null; v.draw(); });
  const tog = (id, key) => $(id).addEventListener('click', e => { v.opts[key] = !v.opts[key]; e.currentTarget.classList.toggle('on', v.opts[key]); v.draw(); });
  tog('#ovTrack', 'track'); tog('#ovFlags', 'flags'); tog('#ovGrid', 'compass');
  $('#resetView').addEventListener('click', () => { v.fit(); v.draw(); });
  $('#smoothBtn').addEventListener('click', e => { v.smooth = !v.smooth; e.currentTarget.classList.toggle('on', v.smooth); v.draw(); });
  $('#playBtn').addEventListener('click', togglePlay);
  $('#prevBtn').addEventListener('click', () => step(-1));
  $('#nextBtn').addEventListener('click', () => step(1));
  $('#fpsSel').addEventListener('change', () => { if (state.playing) { togglePlay(); togglePlay(); } });
  $$('#specSeg button').forEach(b => b.addEventListener('click', () => {
    $$('#specSeg button').forEach(x => x.classList.toggle('on', x === b)); state.specView = b.dataset.v; drawSpectrum();
  }));
  $('#followPm').addEventListener('change', () => v.sel && loadSpectrum(v.sel));
  $('#astroBtn').addEventListener('click', runAstrometry);
  window.addEventListener('keydown', e => {
    if (!location.hash.startsWith('#/observatory') && location.hash) return;
    if (e.target.matches('input, select, textarea')) return;
    if (e.key === 'ArrowRight') { step(1); e.preventDefault(); }
    else if (e.key === 'ArrowLeft') { step(-1); e.preventDefault(); }
    else if (e.key === ' ') { togglePlay(); e.preventDefault(); }
  });
  window.addEventListener('resize', drawTimeline);
}

function showEmpty(title, msg, progress = null) {
  const el = $('#stageEmpty'); el.hidden = false;
  el.querySelector('h2').textContent = title; el.querySelector('p').textContent = msg;
  const bar = $('#loadBar'); bar.hidden = progress == null; if (progress != null) bar.firstElementChild.style.width = `${Math.round(progress * 100)}%`;
}

let openGen = 0;
export async function openField(ra, dec, target = null, size = null) {
  const gen = ++openGen;
  size = size || +$('#fovSel').value; const release = $('#relSel').value;
  if (size) $('#fovSel').value = String(size);
  $('#coordInput').value = `${ra.toFixed(5)} ${dec >= 0 ? '+' : ''}${dec.toFixed(5)}`;
  const q = `ra=${ra.toFixed(5)}&dec=${dec.toFixed(5)}&fov=${size}${target ? `&target=${target}` : ''}`;
  history.replaceState(null, '', `#/observatory?${q}`); $('#permalink').textContent = `#/observatory?${q}`;
  stopPlay(); state.target = target;
  showEmpty('Locating SPHEREx exposures', `${sexRA(ra)} ${sexDec(dec)} · searching ${release.toUpperCase()}`, 0.01);
  $('#loadMsg').textContent = '';
  try {
    const r = await api(`/api/fields?ra=${ra}&dec=${dec}&size=${size}&release=${release}${target ? `&target=${target}` : ''}`, { method: 'POST' });
    if (r.job) await waitJob(r.job, j => { if (gen === openGen) showEmpty('Retrieving observations', j.message || 'working…', j.progress); }, () => gen !== openGen);
    if (gen !== openGen) return;   // superseded by a newer request
    await loadField(r.field);
  } catch (err) { if (gen === openGen) showEmpty('Observations unavailable', err.message); }
}

export async function loadField(fid) {
  const meta = await api(`/api/fields/${fid}`);
  if (!meta.frames || !meta.frames.length) { showEmpty('No coverage', meta.message || 'No SPHEREx exposures found for this position in the selected release.'); return; }
  showEmpty('Loading image stack', `${meta.frames.length} exposures`, 0.95);
  const buf = await (await fetch(`/api/fields/${fid}/stack.f32`)).arrayBuffer();
  const fbuf = await (await fetch(`/api/fields/${fid}/flags.i32`)).arrayBuffer();
  state.field = meta;
  state.viewer.setField(meta, new Float32Array(buf), new Int32Array(fbuf));
  $('#stageEmpty').hidden = true;
  $('#astroBtn').disabled = !meta.track; $('#astroOut').replaceChildren();
  state.spec = null; drawSpectrum(); $('#specPos').textContent = 'click the image';
  applyFilter();
  // default selection: the field centre / target
  const sel = meta.track ? { xi: meta.track.xi[0], eta: meta.track.eta[0] } : { xi: 0, eta: 0 };
  $('#followPm').checked = !!meta.track; state.viewer.sel = sel; loadSpectrum(sel);
}

function applyFilter() {
  const f = state.field; if (!f) return;
  state.visible = f.frames.map((r, i) => i).filter(i => state.dets.has(f.frames[i].detector));
  $('#detCount').textContent = `${state.visible.length} / ${f.frames.length} frames`;
  const v = state.viewer;
  if (!state.visible.includes(v.i)) v.i = state.visible[0] ?? 0;
  if (v.j != null && !state.visible.includes(v.j)) v.j = null;
  v.setVisible(state.visible); updateFrame(); drawTimeline();
}

function step(d) {
  const vis = state.visible; if (!vis.length) return;
  const k = vis.indexOf(state.viewer.i); state.viewer.i = vis[(k + d + vis.length) % vis.length]; updateFrame();
}
function togglePlay() { state.playing ? stopPlay() : startPlay(); }
function startPlay() {
  state.playing = true; $('#playBtn').textContent = 'Pause';
  const fps = +$('#fpsSel').value; state.timer = setInterval(() => step(1), 1000 / fps);
}
function stopPlay() { state.playing = false; clearInterval(state.timer); $('#playBtn').textContent = 'Play'; }

function updateFrame() {
  const f = state.field; if (!f) return;
  const v = state.viewer; v.draw();
  const r = f.frames[v.i];
  $('#whenLbl').textContent = `${isoShort(r.mjd)} · ${DET[r.detector].name} · λ ${r.wave_um.toFixed(3)} µm · frame ${state.visible.indexOf(v.i) + 1}/${state.visible.length}`;
  const kv = $('#frameKv'); kv.replaceChildren();
  const rows = [['Source', h('span', {}, 'SPHEREx ', f.release.toUpperCase(), ' · ', h('a', { href: r.url, target: '_blank', rel: 'noopener' }, 'MEF'))],
    ['Obs ID', r.obs_id], ['Detector', `${DET[r.detector].name} (${DET[r.detector].lo}–${DET[r.detector].hi} µm)`],
    ['λ centre', `${r.wave_um.toFixed(4)} µm`], ['Bandwidth', `${r.band_um.toFixed(4)} µm (R≈${Math.round(r.wave_um / r.band_um)})`],
    ['Mid-exposure', `${r.utc} UTC`], ['MJD (UTC)', r.mjd.toFixed(6)], ['Exposure', `${r.exptime_s.toFixed(1)} s`],
    ['Zodi model', `${r.zodi_mjysr.toFixed(4)} MJy/sr (subtracted)`], ['σ / pixel', `${r.sigma_ujy.toFixed(1)} µJy`],
    ['Flagged', `${(100 * r.flag_frac).toFixed(1)}% of cutout`], ['Wavelength', r.wave_source], ['Discovery', f.discovery]];
  for (const [k, val] of rows) kv.append(h('dt', {}, k), h('dd', {}, val));
  updateBadge(); drawTimelineCursor();
}

const MODE_NOTES = {
  diff: 'Current frame minus a flux-matched comparison (auto: same detector, another epoch, closest λ; Shift-click a tick to choose). Static stars never cancel perfectly: the PSF is undersampled (FWHM ≈ 1.2 px) and the two exposures sample different wavelengths. Moving or variable sources appear as isolated or paired residuals.',
  coadd: 'Mean of all visible frames. It mixes wavelengths, so it shows depth, not colour. A source moving faster than ~1 px per visit is smeared.',
};
function updateBadge() {
  const f = state.field, v = state.viewer, b = $('#badge');
  const mn = $('#modeNote'); if (mn) { mn.textContent = MODE_NOTES[v.opts.mode] || ''; mn.hidden = !MODE_NOTES[v.opts.mode]; }
  if (!f) return;
  const r = f.frames[v.i]; b.replaceChildren();
  b.append(h('span', { class: 'tag real' }, 'REAL'), h('span', {}, `SPHEREx ${f.release.toUpperCase()} · ${DET[r.detector].name} · ${r.wave_um.toFixed(3)} µm · ${r.utc.replace('T', ' ')} UTC`));
  if (v.opts.mode === 'diff') { const j = v.cmpIndex; if (j != null) b.append(h('span', { class: 'tag proc' }, 'DIFFERENCE'), h('span', {}, `minus ${(v.diffScale ?? 1).toFixed(2)} × ${f.frames[j].utc.slice(0, 10)} · ${f.frames[j].wave_um.toFixed(3)} µm (flux-matched)`)); }
  if (v.opts.mode === 'coadd') b.append(h('span', { class: 'tag proc' }, 'CO-ADD'), h('span', {}, `mean of ${state.visible.length} frames · mixed λ`));
  b.append(h('span', { class: 'tag proc', title: 'Bilinear resampling to a north-up grid, for display only' }, 'RESAMPLED'));
}

// ------------------------------------------------------------------ timeline
let tlScale = null;
function drawTimeline() {
  const svg = $('#tlSvg'), f = state.field; svg.replaceChildren(); if (!f) return;
  const W = svg.clientWidth || 800, H = 74; svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const fr = f.frames, t = fr.map(r => r.mjd);
  // compress long gaps (> 40 d, i.e. between survey visits) into fixed-width breaks so visits stay readable; label the gap.
  const order = t.map((v, i) => [v, i]).sort((a, b) => a[0] - b[0]);
  const segs = []; let cur = [order[0][0], order[0][0]];
  for (const [v] of order) { if (v - cur[1] > 40) { segs.push(cur); cur = [v, v]; } else cur[1] = v; }
  segs.push(cur);
  const gapW = 34, padL = 8, padR = 44, usable = W - padL - padR - gapW * (segs.length - 1);
  const tot = segs.reduce((a, s2) => a + Math.max(s2[1] - s2[0], 0.5), 0);
  let x = padL; const map = [];
  segs.forEach((sg, k) => { const w = usable * Math.max(sg[1] - sg[0], 0.5) / tot; map.push({ t0: sg[0], t1: sg[1], x0: x, w }); x += w + gapW; });
  tlScale = mjd => { for (const m of map) if (mjd >= m.t0 - 1e-6 && mjd <= m.t1 + 1e-6) return m.x0 + (m.t1 > m.t0 ? (mjd - m.t0) / (m.t1 - m.t0) * m.w : m.w / 2); return null; };
  const axis = s('g', { class: 'tl-axis' });
  const placed = [];   // avoid overlapping axis labels (approx. 6.2 px per mono character)
  const label = (x, y, txt, anchor = 'start') => {
    const w = txt.length * 6.2, x0 = anchor === 'start' ? x : anchor === 'end' ? x - w : x - w / 2;
    if (x0 < 0 || x0 + w > W - padR + 40 || placed.some(p => p.y === y && x0 < p.x1 + 6 && x0 + w > p.x0 - 6)) return;
    placed.push({ y, x0, x1: x0 + w }); axis.append(s('text', { x, y, 'text-anchor': anchor }, txt));
  };
  map.forEach((m, k) => {
    axis.append(s('line', { x1: m.x0, x2: m.x0 + m.w, y1: 52, y2: 52 }));
    const d0 = mjdToDate(m.t0), d1 = mjdToDate(m.t1);
    label(m.x0, 66, d0.toISOString().slice(0, 10));
    if (m.w > 150) label(m.x0 + m.w, 66, d1.toISOString().slice(5, 10), 'end');
    if (k < map.length - 1) {
      const gx = m.x0 + m.w + gapW / 2, days = map[k + 1].t0 - m.t1;
      axis.append(s('path', { d: `M${gx - 4},48 l3,8 M${gx + 1},48 l3,8`, stroke: '#4a5361', fill: 'none' }));
      label(gx, 42, `${Math.round(days)} d`, 'middle');
    }
  });
  svg.append(axis);
  // detector lanes (two rows: D1-3 SWIR, D4-6 MWIR)
  const g = s('g');
  for (let i = 0; i < fr.length; i++) {
    const r = fr[i], xx = tlScale(r.mjd); if (xx == null) continue;
    const on = state.dets.has(r.detector), y = r.detector <= 3 ? 8 : 26;
    const el = s('rect', { x: xx - 1, y, width: 2.2, height: 14, rx: 0.6, fill: waveColor(r.wave_um, on ? 0.95 : 0.15), 'data-i': i, style: 'cursor:pointer' });
    el.append(s('title', {}, `${r.utc} · ${DET[r.detector].name} · ${r.wave_um.toFixed(3)} µm`));
    g.append(el);
  }
  g.addEventListener('click', e => {
    const i = e.target.getAttribute('data-i'); if (i == null) return;
    const k = +i; if (!state.dets.has(f.frames[k].detector)) return;
    if (e.shiftKey) { state.viewer.j = k; state.viewer.opts.mode = 'diff'; $$('#modeSeg button').forEach(b => b.classList.toggle('on', b.dataset.m === 'diff')); }
    else state.viewer.i = k;
    updateFrame();
  });
  svg.append(g, s('g', { id: 'tlCursor' }));
  svg.append(s('text', { x: W - 6, y: 18, 'text-anchor': 'end', fill: '#5a6472', 'font-size': 9, 'font-family': 'JetBrains Mono' }, 'D1–D3'),
    s('text', { x: W - 6, y: 36, 'text-anchor': 'end', fill: '#5a6472', 'font-size': 9, 'font-family': 'JetBrains Mono' }, 'D4–D6'));
  drawTimelineCursor();
}
function drawTimelineCursor() {
  const g = $('#tlCursor'), f = state.field; if (!g || !f || !tlScale) return; g.replaceChildren();
  const v = state.viewer, x = tlScale(f.frames[v.i].mjd);
  if (x != null) g.append(s('path', { d: `M${x},3 V46`, stroke: '#dfe4ec', 'stroke-width': 1 }), s('path', { d: `M${x - 4},1 h8 l-4,5z`, fill: '#dfe4ec' }));
  if (v.opts.mode === 'diff' && v.cmpIndex != null) { const xj = tlScale(f.frames[v.cmpIndex].mjd); if (xj != null) g.append(s('path', { d: `M${xj},3 V46`, stroke: '#8cb8ff', 'stroke-dasharray': '2 2' })); }
}

// ------------------------------------------------------------------ spectrum
async function loadSpectrum(sel) {
  const f = state.field; if (!f) return;
  const follow = $('#followPm').checked && !!f.track;
  let xi = sel.xi, eta = sel.eta;
  if (follow) { xi -= f.track.xi[state.viewer.i]; eta -= f.track.eta[state.viewer.i]; }
  $('#specPos').textContent = 'measuring…';
  try {
    state.spec = await api(`/api/fields/${f.id}/spectrum?xi=${xi}&eta=${eta}&follow=${follow}`);
    $('#specPos').textContent = `${sexRA(state.spec.ra)} ${sexDec(state.spec.dec)}`;
  } catch (e) { state.spec = null; $('#specPos').textContent = e.message; }
  drawSpectrum();
}
function drawSpectrum() {
  const svg = $('#specSvg'), sp = state.spec, lg = $('#specLegend');
  lg.replaceChildren(); $('#specNote').textContent = '';
  if (!sp || !sp.points.length) {
    svg.replaceChildren();
    if (sp) $('#specNote').textContent = 'No usable measurement at this position: the PSF core is flagged (saturated, nonlinear or missing) or off the detector in every exposure. Very bright stars such as Barnard’s Star saturate in most SPHEREx bands. Click a fainter source.';
    return;
  }
  const p = sp.points.filter(q => state.dets.has(q.det));
  const byWave = state.specView === 'wave';
  chart(svg, {
    W: 312, H: 190, xlog: byWave,
    series: [{ x: p.map(q => byWave ? q.wave : q.t - Math.floor(p[0]?.t || 0)), y: p.map(q => q.flux), e: p.map(q => q.err), colors: p.map(q => waveColor(q.wave)), r: 2 }],
    xlab: byWave ? 'λ [µm]' : `days since MJD ${Math.floor(p[0]?.t || 0)}`, ylab: 'F_ν [µJy]', hline: 0,
  });
  const dets = [...new Set(p.map(q => q.det))].sort();
  for (const d of dets) lg.append(h('span', {}, h('i', { style: `background:${detColor(d)}` }), DET[d].name));
  const snr = p.map(q => q.flux / q.err), med = snr.sort((a, b) => a - b)[Math.floor(snr.length / 2)] || 0;
  $('#specNote').textContent = `${p.length} measurements · median S/N ${fmt(med, 1)} · ${sp.method}. Each point is one exposure; SPHEREx samples a different wavelength in every exposure.`;
}

// ---------------------------------------------------------------- astrometry
async function runAstrometry() {
  const f = state.field, out = $('#astroOut'); if (!f) return;
  out.replaceChildren(h('p', { class: 'note' }, 'Fitting PSF centroids…'));
  try {
    const r = await api(`/api/fields/${f.id}/astrometry`);
    if (r.error) { out.replaceChildren(h('p', { class: 'note' }, r.error)); return; }
    const fit = r.fit, ex = r.expected || {};
    const tbl = h('table', { class: 't' },
      h('tr', {}, h('th', {}, ''), h('th', { class: 'num' }, 'SPHEREx'), h('th', { class: 'num' }, 'Catalogue')),
      h('tr', {}, h('td', {}, 'μα* [mas/yr]'), h('td', { class: 'num' }, `${fmt(fit.pmra_masyr, 0)} ± ${fmt(fit.e_pmra_masyr, 0)}`), h('td', { class: 'num' }, fmt(ex.pmra_masyr, 0))),
      h('tr', {}, h('td', {}, 'μδ [mas/yr]'), h('td', { class: 'num' }, `${fmt(fit.pmdec_masyr, 0)} ± ${fmt(fit.e_pmdec_masyr, 0)}`), h('td', { class: 'num' }, fmt(ex.pmdec_masyr, 0))),
      h('tr', {}, h('td', {}, 'N centroids'), h('td', { class: 'num' }, fit.n), h('td', { class: 'num' }, '')),
      h('tr', {}, h('td', {}, 'χ²/dof'), h('td', { class: 'num' }, `${fmt(fit.chi2 / fit.dof, 2)}`), h('td', { class: 'num' }, '')));
    const svg = s('svg', { class: 'chart', viewBox: '0 0 312 150' });
    const t0 = Math.floor(r.points[0].t);
    out.replaceChildren(tbl, svg, h('p', { class: 'note' }, `Offsets of the measured centroid from the field centre (east, north) vs time. ${r.method}. Catalogue: ${ex.source || '—'}.`),
      fit.chi2 / fit.dof > 2 ? h('p', { class: 'note callout warn' }, `χ²/dof = ${fmt(fit.chi2 / fit.dof, 1)}: the scatter exceeds the centroid errors, which are already scaled by √(χ²/dof). Per-visit systematics (saturated cores, absolute WCS between visits) dominate. Treat the uncertainties as lower limits; see Results → E2.`) : '');
    chart(svg, { W: 312, H: 150, series: [
      { x: r.points.map(p => p.t - t0), y: r.points.map(p => p.eta), e: r.points.map(p => p.seta), color: '#8cb8ff', r: 1.8 },
      { x: r.points.map(p => p.t - t0), y: r.points.map(p => p.xi), e: r.points.map(p => p.sxi), color: '#e7c27d', r: 1.8 }],
      xlab: `days since MJD ${t0}`, ylab: 'offset ["]  · north (blue) east (amber)' });
  } catch (e) { out.replaceChildren(h('p', { class: 'note' }, e.message)); }
}
