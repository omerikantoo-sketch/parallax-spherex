// Canvas sky viewer: per-frame stretch, difference / co-add modes, pan & zoom,
// overlays (compass, scale bar, predicted track, flags, selection) and a
// cursor readout in RA/Dec and calibrated surface brightness.
import { tan2sky, sexRA, sexDec, fmt } from './util.js';

const REJECT = (1 << 0) | (1 << 1) | (1 << 2) | (1 << 6) | (1 << 9) | (1 << 10) | (1 << 11) | (1 << 15) | (1 << 17) | (1 << 19);

function percentiles(arr, ps) {
  const v = [];
  const step = Math.max(1, Math.floor(arr.length / 40000));
  for (let i = 0; i < arr.length; i += step) { const x = arr[i]; if (Number.isFinite(x)) v.push(x); }
  if (!v.length) return ps.map(() => 0);
  v.sort((a, b) => a - b);
  return ps.map(p => v[Math.min(v.length - 1, Math.max(0, Math.round(p / 100 * (v.length - 1))))]);
}

export class Viewer {
  constructor(stage, imgCanvas, ovCanvas, readout) {
    Object.assign(this, { stage, imgCanvas, ovCanvas, readout });
    this.ctx = imgCanvas.getContext('2d'); this.octx = ovCanvas.getContext('2d');
    this.off = document.createElement('canvas'); this.offCtx = this.off.getContext('2d');
    this.zoom = 1; this.panX = 0; this.panY = 0; this.smooth = false;
    this.opts = { mode: 'frame', stretch: 'asinh', lo: 1, hi: 99.6, scale: 'frame', track: true, flags: false, compass: true };
    this.sel = null; this.onSelect = null; this.field = null; this.visible = [];
    this.i = 0; this.j = null; this._limits = null;
    new ResizeObserver(() => this.resize()).observe(stage);
    this._bindMouse();
  }

  setField(meta, stack, flags) {
    this.field = meta; this.stack = stack; this.flagStack = flags; this.S = meta.size; this.n = meta.frames.length;
    this.off.width = this.off.height = this.S; this.sel = null; this._coadd = null; this._global = null;
    this.fit(); this.resize();
  }
  frameData(k) { return this.stack.subarray(k * this.S * this.S, (k + 1) * this.S * this.S); }
  setVisible(list) { this.visible = list; this._coadd = null; this.draw(); }

  fit() { this.zoom = 1; this.panX = 0; this.panY = 0; }
  resize() {
    const r = this.stage.getBoundingClientRect(), dpr = window.devicePixelRatio || 1;
    for (const c of [this.imgCanvas, this.ovCanvas]) { c.width = Math.round(r.width * dpr); c.height = Math.round(r.height * dpr); }
    this.dpr = dpr; this.W = r.width; this.H = r.height; this.draw();
  }
  // screen <-> display-pixel transform
  get base() { return Math.min(this.W, this.H) * 0.94 / (this.S || 1); }
  get scale() { return this.base * this.zoom; }
  toScreen(px, py) { const sc = this.scale; return [this.W / 2 + this.panX + (px - this.S / 2) * sc, this.H / 2 + this.panY + (py - this.S / 2) * sc]; }
  toPix(sx, sy) { const sc = this.scale; return [(sx - this.W / 2 - this.panX) / sc + this.S / 2, (sy - this.H / 2 - this.panY) / sc + this.S / 2]; }
  pixToTan(px, py) { const c = (this.S - 1) / 2, p = this.field.pixscale_arcsec; return [-(px - 0.5 - c) * p, (c - (py - 0.5)) * p]; }
  tanToPix(xi, eta) { const c = (this.S - 1) / 2, p = this.field.pixscale_arcsec; return [c - xi / p + 0.5, c - eta / p + 0.5]; }

  _composite() {
    const S2 = this.S * this.S, m = this.opts.mode;
    if (m === 'frame') return { data: this.frameData(this.i), diverging: false };
    if (m === 'diff') {
      const j = this.j ?? this._autoCompare();
      const a = this.frameData(this.i), b = this.frameData(j), out = new Float32Array(S2);
      // flux-matched difference: the two exposures sample different wavelengths, so
      // static stars differ by their spectral slope.  Scale b by the median a/b over
      // the brightest common pixels (robust), so static sources cancel to first order.
      const both = []; for (let k = 0; k < S2; k++) if (Number.isFinite(a[k]) && Number.isFinite(b[k])) both.push(b[k]);
      both.sort((x, y) => x - y); const thr = both[Math.floor(both.length * 0.97)] ?? Infinity;
      const r = []; for (let k = 0; k < S2; k++) if (b[k] >= thr && Number.isFinite(a[k]) && b[k] > 0) r.push(a[k] / b[k]);
      r.sort((x, y) => x - y); const sc = r.length >= 5 ? r[Math.floor(r.length / 2)] : 1;
      this.diffScale = sc;
      for (let k = 0; k < S2; k++) out[k] = a[k] - sc * b[k];
      return { data: out, diverging: true, j };
    }
    if (!this._coadd) {
      const sum = new Float64Array(S2), cnt = new Uint16Array(S2);
      for (const k of this.visible) { const a = this.frameData(k); for (let q = 0; q < S2; q++) if (Number.isFinite(a[q])) { sum[q] += a[q]; cnt[q]++; } }
      const out = new Float32Array(S2); for (let q = 0; q < S2; q++) out[q] = cnt[q] ? sum[q] / cnt[q] : NaN;
      this._coadd = out;
    }
    return { data: this._coadd, diverging: false };
  }
  _autoCompare() {
    // nearest visible frame of the same detector from a different epoch (>1 d away), else previous visible
    const fr = this.field.frames, me = fr[this.i];
    let best = null, bd = Infinity;
    for (const k of this.visible) {
      if (k === this.i || fr[k].detector !== me.detector) continue;
      const dt = Math.abs(fr[k].mjd - me.mjd); if (dt < 1) continue;
      const dw = Math.abs(fr[k].wave_um - me.wave_um);
      const score = dw * 50 + Math.abs(Math.log(dt));
      if (score < bd) { bd = score; best = k; }
    }
    return best ?? (this.visible[Math.max(0, this.visible.indexOf(this.i) - 1)] ?? this.i);
  }

  draw() {
    const ctx = this.ctx; if (!ctx) return;
    ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.clearRect(0, 0, this.imgCanvas.width, this.imgCanvas.height);
    if (!this.field || !this.n) { this.drawOverlay(); return; }
    const { data, diverging, j } = this._composite(); this.cmpIndex = j;
    let lo, hi;
    if (diverging) { const [p] = percentiles(data.map(Math.abs), [this.opts.hi]); lo = -p; hi = p; }
    else if (this.opts.scale === 'global' && this.opts.mode === 'frame') {
      if (!this._global) this._global = percentiles(this.stack, [this.opts.lo, this.opts.hi]);
      [lo, hi] = this._global;
    } else [lo, hi] = percentiles(data, [this.opts.lo, this.opts.hi]);
    this._limits = [lo, hi];
    const img = this.offCtx.createImageData(this.S, this.S), px = img.data, st = this.opts.stretch;
    const span = (hi - lo) || 1, a = 10;
    const f = st === 'linear' ? x => x : st === 'sqrt' ? x => Math.sqrt(x) : st === 'log' ? x => Math.log10(1 + 999 * x) / 3 : x => Math.asinh(a * x) / Math.asinh(a);
    for (let k = 0; k < data.length; k++) {
      const v = data[k], o = 4 * k;
      if (!Number.isFinite(v)) { px[o] = 26; px[o + 1] = 14; px[o + 2] = 18; px[o + 3] = 255; continue; }
      if (diverging) {
        const t = Math.max(-1, Math.min(1, v / hi)), u = Math.asinh(4 * Math.abs(t)) / Math.asinh(4);
        if (t >= 0) { px[o] = 20 + 211 * u; px[o + 1] = 20 + 174 * u; px[o + 2] = 22 + 103 * u; }
        else { px[o] = 20 + 120 * u; px[o + 1] = 20 + 164 * u; px[o + 2] = 22 + 233 * u; }
      } else {
        const t = Math.max(0, Math.min(1, (v - lo) / span)), g = Math.round(255 * f(t));
        px[o] = g; px[o + 1] = g; px[o + 2] = Math.min(255, g + 4);
      }
      px[o + 3] = 255;
    }
    if (this.opts.flags && this.opts.mode === 'frame' && this.flagStack) {
      const fl = this.flagStack.subarray(this.i * this.S * this.S, (this.i + 1) * this.S * this.S);
      for (let k = 0; k < fl.length; k++) if (fl[k] & REJECT) { const o = 4 * k; px[o] = Math.max(px[o], 200); px[o + 1] *= 0.45; px[o + 2] *= 0.45; }
    }
    this.offCtx.putImageData(img, 0, 0);
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    ctx.imageSmoothingEnabled = this.smooth;
    const [x0, y0] = this.toScreen(0, 0), sc = this.scale;
    ctx.drawImage(this.off, x0, y0, this.S * sc, this.S * sc);
    ctx.strokeStyle = '#1f252d'; ctx.lineWidth = 1; ctx.strokeRect(x0 - 0.5, y0 - 0.5, this.S * sc + 1, this.S * sc + 1);
    this.drawOverlay();
  }

  drawOverlay() {
    const o = this.octx; o.setTransform(1, 0, 0, 1, 0, 0); o.clearRect(0, 0, this.ovCanvas.width, this.ovCanvas.height);
    if (!this.field || !this.n) return;
    o.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    const tr = this.field.track;
    if (this.opts.track && tr) {
      o.lineWidth = 1;
      for (const k of this.visible) {
        const [px, py] = this.tanToPix(tr.xi[k], tr.eta[k]); const [sx, sy] = this.toScreen(px, py);
        o.fillStyle = 'rgba(231,194,125,0.35)'; o.beginPath(); o.arc(sx, sy, 1.6, 0, 7); o.fill();
      }
      const [px, py] = this.tanToPix(tr.xi[this.i], tr.eta[this.i]); const [sx, sy] = this.toScreen(px, py);
      const r = Math.max(9, this.scale * 1.6);
      o.strokeStyle = 'rgba(231,194,125,0.95)'; o.beginPath();
      for (const [a0, a1] of [[0.15, 1.42], [1.72, 2.99], [3.29, 4.56], [4.86, 6.13]]) { o.moveTo(sx + r * Math.cos(a0), sy + r * Math.sin(a0)); o.arc(sx, sy, r, a0, a1); }
      o.stroke();
      o.fillStyle = 'rgba(231,194,125,0.9)'; o.font = '10px "JetBrains Mono", monospace';
      o.fillText(`${tr.name} · predicted`, sx + r + 5, sy - r + 2);
    }
    if (this.sel) {
      const [px, py] = this.tanToPix(this.sel.xi, this.sel.eta); const [sx, sy] = this.toScreen(px, py);
      o.strokeStyle = 'rgba(140,184,255,0.95)'; o.lineWidth = 1; const g = 5, L = 12;
      o.beginPath(); o.moveTo(sx - L, sy); o.lineTo(sx - g, sy); o.moveTo(sx + g, sy); o.lineTo(sx + L, sy);
      o.moveTo(sx, sy - L); o.lineTo(sx, sy - g); o.moveTo(sx, sy + g); o.lineTo(sx, sy + L); o.stroke();
    }
    if (this.opts.compass) {
      const x = this.W - 54, y = this.H - 58; o.strokeStyle = '#8590a0'; o.fillStyle = '#8590a0'; o.lineWidth = 1;
      o.font = '10px "JetBrains Mono", monospace';
      o.beginPath(); o.moveTo(x, y); o.lineTo(x, y - 26); o.moveTo(x, y); o.lineTo(x - 26, y); o.stroke();
      o.fillText('N', x - 3, y - 30); o.fillText('E', x - 37, y + 3);
      const arcmin = 60 / this.field.pixscale_arcsec * this.scale;
      const bx = this.W - 24 - arcmin, by = this.H - 20;
      o.beginPath(); o.moveTo(bx, by); o.lineTo(bx + arcmin, by); o.moveTo(bx, by - 3); o.lineTo(bx, by + 3); o.moveTo(bx + arcmin, by - 3); o.lineTo(bx + arcmin, by + 3); o.stroke();
      o.fillText("1′", bx + arcmin / 2 - 5, by - 6);
    }
  }

  _bindMouse() {
    let drag = null;
    const st = this.ovCanvas;
    st.style.cursor = 'crosshair';
    st.addEventListener('wheel', e => {
      if (!this.field) return; e.preventDefault();
      const r = st.getBoundingClientRect(), mx = e.clientX - r.left, my = e.clientY - r.top;
      const [px, py] = this.toPix(mx, my);
      this.zoom = Math.min(24, Math.max(0.5, this.zoom * Math.exp(-e.deltaY * 0.0015)));
      const [sx, sy] = this.toScreen(px, py); this.panX += mx - sx; this.panY += my - sy; this.draw();
    }, { passive: false });
    st.addEventListener('pointerdown', e => { drag = { x: e.clientX, y: e.clientY, px: this.panX, py: this.panY, moved: false }; st.setPointerCapture(e.pointerId); });
    st.addEventListener('pointermove', e => {
      const r = st.getBoundingClientRect();
      if (drag) {
        const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
        if (Math.hypot(dx, dy) > 3) drag.moved = true;
        if (drag.moved) { this.panX = drag.px + dx; this.panY = drag.py + dy; this.draw(); }
      }
      this._readout(e.clientX - r.left, e.clientY - r.top);
    });
    st.addEventListener('pointerup', e => {
      const r = st.getBoundingClientRect();
      if (drag && !drag.moved && this.field) {
        const [px, py] = this.toPix(e.clientX - r.left, e.clientY - r.top);
        if (px >= 0 && py >= 0 && px < this.S && py < this.S) {
          const [xi, eta] = this.pixToTan(px, py); this.sel = { xi, eta }; this.drawOverlay(); this.onSelect && this.onSelect(this.sel);
        }
      }
      drag = null;
    });
    st.addEventListener('pointerleave', () => { this.readout.hidden = true; });
  }
  _readout(mx, my) {
    if (!this.field || !this.n) return;
    const [px, py] = this.toPix(mx, my);
    if (px < 0 || py < 0 || px >= this.S || py >= this.S) { this.readout.hidden = true; return; }
    const [xi, eta] = this.pixToTan(px, py), [ra, dec] = tan2sky(xi, eta, this.field.ra, this.field.dec);
    const k = Math.floor(py) * this.S + Math.floor(px);
    const { data, diverging } = this._composite();
    const v = data[k];
    const fl = this.flagStack ? this.flagStack[this.i * this.S * this.S + k] : 0;
    this.readout.hidden = false;
    this.readout.textContent = `${sexRA(ra)}  ${sexDec(dec)}\nα ${ra.toFixed(5)}°  δ ${dec.toFixed(5)}°   Δ ${fmt(xi, 1)}″, ${fmt(eta, 1)}″\n${diverging ? 'Δ' : ''}I ${Number.isFinite(v) ? v.toFixed(4) : 'masked'} MJy/sr${fl ? `   flags 0x${fl.toString(16)}` : ''}`;
    this.readout.style.whiteSpace = 'pre';
  }
}
