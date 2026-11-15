// PARALLAX front-end entry: router, archive status, page initialisation.
import { api, $, $$, h } from './util.js';
import { initObservatory, openField, state } from './observatory.js';
import { renderMotion } from './motion.js';
import { renderGravity } from './gravity.js';
import { renderResults } from './results.js';
import { renderAbout } from './about.js';

const VIEWS = { observatory: 'obs', motion: 'motion', gravity: 'gravity', results: 'results', about: 'about' };
const rendered = new Set();

async function refreshStatus() {
  const el = $('#status');
  try {
    const s = await api('/api/status');
    const item = (label, v) => h('span', { title: v }, h('i', { class: `dot ${v === 'online' ? 'ok' : 'bad'}` }), `${label} ${v === 'online' ? 'online' : 'offline'}`);
    el.replaceChildren(item('IRSA SIA', s.irsa_sia), item('AWS mirror', s.aws_mirror), h('span', { class: 'mono faint' }, `v${s.version}`));
  } catch (_) { el.replaceChildren(h('span', {}, h('i', { class: 'dot bad' }), 'API offline')); }
}

function parseHash() {
  const [path, query] = (location.hash.replace(/^#\/?/, '') || 'observatory').split('?');
  return { view: VIEWS[path] ? path : 'observatory', params: new URLSearchParams(query || '') };
}

function route() {
  const { view, params } = parseHash();
  const id = VIEWS[view];
  $$('.view').forEach(v => v.classList.toggle('on', v.id === id));
  $$('.nav a').forEach(a => a.classList.toggle('on', a.dataset.view === id));
  document.title = view === 'observatory' ? 'PARALLAX Observatory' : `PARALLAX · ${view[0].toUpperCase() + view.slice(1)}`;
  if (view === 'observatory') {
    state.viewer && state.viewer.resize();
    const ra = parseFloat(params.get('ra')), dec = parseFloat(params.get('dec'));
    if (isFinite(ra) && isFinite(dec)) {
      const f = state.field;
      if (!f || Math.abs(f.ra - ra) > 1e-4 || Math.abs(f.dec - dec) > 1e-4) openField(ra, dec, params.get('target'), +params.get('fov') || null);
    }
  }
  if (!rendered.has(view)) {
    rendered.add(view);
    if (view === 'motion') renderMotion($('#motionInner'));
    if (view === 'gravity') renderGravity($('#gravityInner'));
    if (view === 'results') renderResults($('#resultsInner'));
    if (view === 'about') renderAbout($('#aboutInner'));
  }
}

async function main() {
  refreshStatus(); setInterval(refreshStatus, 120000);
  let targets = [];
  try { targets = await api('/api/targets'); } catch (_) {}
  initObservatory(targets);
  window.addEventListener('hashchange', route);
  route();
}
main();
