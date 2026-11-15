// Shared helpers: API access, formatting, wavelength colours, SVG building.

export async function api(path, opts = {}) {
  const r = await fetch(path, opts);
  if (!r.ok) {
    let msg = r.statusText;
    try { const j = await r.json(); msg = j.detail || JSON.stringify(j); } catch (_) {}
    throw new Error(`${r.status} ${msg}`);
  }
  const ct = r.headers.get('content-type') || '';
  return ct.includes('json') ? r.json() : r;
}

export async function waitJob(id, onProgress, cancelled = () => false) {
  for (;;) {
    if (cancelled()) return null;
    const j = await api(`/api/jobs/${id}`);
    onProgress && onProgress(j);
    if (j.status === 'done') return j.result;
    if (j.status === 'error') throw new Error(j.error || 'job failed');
    await new Promise(r => setTimeout(r, 700));
  }
}

export const $ = (s, el = document) => el.querySelector(s);
export const $$ = (s, el = document) => [...el.querySelectorAll(s)];
export function h(tag, attrs = {}, ...kids) {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (k === 'class') el.className = v;
    else if (k.startsWith('on')) el.addEventListener(k.slice(2), v);
    else if (v !== false && v != null) el.setAttribute(k, v === true ? '' : v);
  }
  for (const c of kids.flat()) if (c != null) el.append(c.nodeType ? c : document.createTextNode(c));
  return el;
}

// ---------------------------------------------------------------- formatting
export const fmt = (v, d = 2) => (v == null || !isFinite(v)) ? '—' : Number(v).toFixed(d);
export function sexRA(deg) {
  const h = ((deg % 360) + 360) % 360 / 15; const H = Math.floor(h); const m = (h - H) * 60; const M = Math.floor(m);
  const s = (m - M) * 60; return `${String(H).padStart(2, '0')}h${String(M).padStart(2, '0')}m${s.toFixed(2).padStart(5, '0')}s`;
}
export function sexDec(deg) {
  const sg = deg < 0 ? '−' : '+'; const a = Math.abs(deg); const D = Math.floor(a); const m = (a - D) * 60; const M = Math.floor(m);
  const s = (m - M) * 60; return `${sg}${String(D).padStart(2, '0')}°${String(M).padStart(2, '0')}′${s.toFixed(1).padStart(4, '0')}″`;
}
export function mjdToDate(mjd) { return new Date((mjd - 40587) * 86400000); }
export function isoShort(mjd) { return mjdToDate(mjd).toISOString().slice(0, 16).replace('T', ' ') + ' UTC'; }

// gnomonic inverse: tangent offsets (arcsec) -> RA/Dec (deg)
export function tan2sky(xi, eta, ra0, dec0) {
  const r = Math.PI / 180, x = xi / 3600 * r, y = eta / 3600 * r, a0 = ra0 * r, d0 = dec0 * r;
  const rho = Math.hypot(x, y), c = Math.atan(rho);
  if (rho === 0) return [ra0, dec0];
  const dec = Math.asin(Math.cos(c) * Math.sin(d0) + y * Math.sin(c) * Math.cos(d0) / rho);
  const ra = a0 + Math.atan2(x * Math.sin(c), rho * Math.cos(d0) * Math.cos(c) - y * Math.sin(d0) * Math.sin(c));
  return [((ra / r) % 360 + 360) % 360, dec / r];
}

// ------------------------------------------------------------- wavelength map
const WSTOPS = [[0.75, [122, 167, 255]], [1.25, [110, 205, 214]], [1.9, [168, 214, 120]], [2.8, [236, 205, 110]],
                [3.9, [242, 152, 98]], [5.0, [232, 104, 125]]];
export function waveColor(w, a = 1) {
  const x = Math.log(Math.min(5, Math.max(0.75, w)));
  for (let i = 0; i < WSTOPS.length - 1; i++) {
    const [w0, c0] = WSTOPS[i], [w1, c1] = WSTOPS[i + 1];
    if (x <= Math.log(w1)) {
      const f = (x - Math.log(w0)) / (Math.log(w1) - Math.log(w0));
      const c = c0.map((v, k) => Math.round(v + f * (c1[k] - v)));
      return `rgba(${c[0]},${c[1]},${c[2]},${a})`;
    }
  }
  return `rgba(232,104,125,${a})`;
}
export const DET = {
  1: { name: 'D1', lo: 0.75, hi: 1.09 }, 2: { name: 'D2', lo: 1.10, hi: 1.62 }, 3: { name: 'D3', lo: 1.63, hi: 2.41 },
  4: { name: 'D4', lo: 2.42, hi: 3.82 }, 5: { name: 'D5', lo: 3.83, hi: 4.41 }, 6: { name: 'D6', lo: 4.42, hi: 5.00 },
};
export const detColor = d => waveColor(Math.sqrt(DET[d].lo * DET[d].hi));

// ---------------------------------------------------------------- SVG charts
const NS = 'http://www.w3.org/2000/svg';
export function s(tag, attrs = {}, ...kids) {
  const el = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) if (v != null) el.setAttribute(k, v);
  for (const c of kids.flat()) if (c != null) el.append(c.nodeType ? c : document.createTextNode(c));
  return el;
}
export function niceTicks(lo, hi, n = 5) {
  if (!isFinite(lo) || !isFinite(hi) || lo === hi) return [lo];
  const span = hi - lo, step0 = span / n, mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(st => span / st <= n) || 10 * mag;
  const out = []; for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9 * span; v += step) out.push(+v.toPrecision(10));
  return out;
}
export function logTicks(lo, hi) {
  const c = [0.75, 1, 1.5, 2, 3, 4, 5]; return c.filter(v => v >= lo * 0.999 && v <= hi * 1.001);
}

/** Generic scatter/line chart into an <svg>. series: [{x:[],y:[],e:[],color|colors,type:'pts'|'line',w,dash,label}] */
export function chart(svg, { series, xlab = '', ylab = '', xlog = false, ylog = false, xr, yr, W, H, pad = [34, 12, 22, 42],
                              xfmt = v => +v.toPrecision(3), yfmt = v => +v.toPrecision(3), hline = null, vlines = [] }) {
  const vb = svg.viewBox.baseVal; W = W || vb.width || 400; H = H || vb.height || 200;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`); svg.replaceChildren();
  const [pt, pr, pb, pl] = pad;
  const xs = series.flatMap(s => s.x).filter(isFinite), ys = series.flatMap(s => s.y.map((v, i) => [v - (s.e?.[i] || 0), v + (s.e?.[i] || 0)]).flat()).filter(isFinite);
  let [x0, x1] = xr || [Math.min(...xs), Math.max(...xs)]; let [y0, y1] = yr || [Math.min(...ys), Math.max(...ys)];
  if (!xs.length) { x0 = 0; x1 = 1; } if (!ys.length) { y0 = 0; y1 = 1; }
  if (!xr && !xlog) { const d = (x1 - x0) || 1; x0 -= d * 0.03; x1 += d * 0.03; }
  if (!yr && !ylog) { const d = (y1 - y0) || Math.abs(y0) || 1; y0 -= d * 0.08; y1 += d * 0.08; }
  const tx = xlog ? v => Math.log(v) : v => v, ty = ylog ? v => Math.log(Math.max(v, 1e-12)) : v => v;
  const X = v => pl + (tx(v) - tx(x0)) / (tx(x1) - tx(x0)) * (W - pl - pr);
  const Y = v => H - pb - (ty(v) - ty(y0)) / (ty(y1) - ty(y0)) * (H - pt - pb);
  const g = s('g', { class: 'grid' }), ax = s('g', { class: 'ax' });
  const xt = xlog ? logTicks(x0, x1) : niceTicks(x0, x1, 5), yt = niceTicks(y0, y1, 4);
  for (const v of yt) { g.append(s('line', { x1: pl, x2: W - pr, y1: Y(v), y2: Y(v) })); ax.append(s('text', { x: pl - 5, y: Y(v) + 3, 'text-anchor': 'end' }, yfmt(v))); }
  for (const v of xt) { g.append(s('line', { y1: pt, y2: H - pb, x1: X(v), x2: X(v) })); ax.append(s('text', { x: X(v), y: H - pb + 13, 'text-anchor': 'middle' }, xfmt(v))); }
  ax.append(s('path', { d: `M${pl},${pt}V${H - pb}H${W - pr}`, fill: 'none' }));
  if (xlab) ax.append(s('text', { x: W - pr, y: H - 3, 'text-anchor': 'end' }, xlab));
  if (ylab) ax.append(s('text', { x: 4, y: pt - 4 }, ylab));
  svg.append(g, ax);
  if (hline != null && hline >= y0 && hline <= y1) svg.append(s('line', { x1: pl, x2: W - pr, y1: Y(hline), y2: Y(hline), stroke: '#48515e', 'stroke-dasharray': '3 3' }));
  for (const vl of vlines) if (vl.x >= x0 && vl.x <= x1) svg.append(s('line', { x1: X(vl.x), x2: X(vl.x), y1: pt, y2: H - pb, stroke: vl.color || '#48515e', 'stroke-dasharray': vl.dash || '2 3' }));
  const clipId = 'c' + Math.random().toString(36).slice(2);
  svg.append(s('clipPath', { id: clipId }, s('rect', { x: pl, y: pt, width: W - pl - pr, height: H - pt - pb })));
  const body = s('g', { 'clip-path': `url(#${clipId})` }); svg.append(body);
  for (const se of series) {
    if (se.type === 'line') {
      const pts = se.x.map((x, i) => [x, se.y[i]]).filter(p => isFinite(p[0]) && isFinite(p[1]));
      if (!pts.length) continue;
      body.append(s('path', { d: 'M' + pts.map(p => `${X(p[0]).toFixed(1)},${Y(p[1]).toFixed(1)}`).join('L'), fill: 'none',
        stroke: se.color, 'stroke-width': se.w || 1.4, 'stroke-dasharray': se.dash || null, opacity: se.opacity ?? 1 }));
    } else {
      se.x.forEach((x, i) => {
        const y = se.y[i]; if (!isFinite(x) || !isFinite(y)) return;
        const c = se.colors ? se.colors[i] : se.color, cx = X(x), cy = Y(y);
        if (se.e && isFinite(se.e[i])) body.append(s('line', { x1: cx, x2: cx, y1: Y(y - se.e[i]), y2: Y(y + se.e[i]), stroke: c, 'stroke-width': 1, opacity: 0.55 }));
        body.append(s('circle', { cx, cy, r: se.r || 2.2, fill: c, opacity: se.opacity ?? 0.95 }));
      });
    }
  }
  return { X, Y };
}
