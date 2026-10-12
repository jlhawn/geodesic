import * as THREE from '../js/three.module.js';
import { Grid } from '../js/grid.module.js';
import { Globe } from './globe.module.js';
import { slider, choice, buttons, legend, caption, readout, paletteVersion } from './runtime.module.js';
import { loadFrames } from './frames.module.js';

export function playbackGlobe(root, { url, N, center, sequences, fields, framesPerSecond = 8, note = '', when = null, globalMean = false, arrows = null }) {
  const controls = root.querySelector('.controls');
  const globe = new Globe(root, { height: 460, center, lighting: { directional: 1.5, ambient: 1.9 } });
  const cells = [...new Grid(N)], polygons = cells.map((c) => ({ center: c.centerVertex, vertices: c.vertices }));
  const areas = globalMean ? Float64Array.from(cells, (c) => c.area) : null, totalArea = areas ? areas.reduce((sum, a) => sum + a, 0) : 0;
  const paint = globe.dynamicCells(polygons);
  const drawArrows = arrows ? arrowLayer(globe, cells, arrows) : null;
  const top = caption(root, note);
  const key = legend(root, []);
  let data = null, sequence = sequences[0].name, field = 0, position = 0, playing = true, values = null, shown = null, scrub = null;
  const decoded = new Map(), scratch = {};
  const fieldIndex = () => data.header.fields.findIndex((f) => f.name === fields[field].name);
  const times = () => data.times(sequence);
  const rate = () => sequences.find((s) => s.name === sequence).framesPerSecond ?? framesPerSecond;

  function load() {
    const frame = Math.min(times().length - 1, Math.floor(position));
    const id = `${sequence}/${frame}/${field}/${paletteVersion}`;
    if (id === shown) return;
    const spec = fields[field];
    decoded.clear();
    const read = (name) => {
      if (!decoded.has(name)) decoded.set(name, scratch[name] = data.values(sequence, frame, data.header.fields.findIndex((f) => f.name === name), scratch[name]));
      return decoded.get(name);
    };
    values = spec.compute ? spec.compute(read, values ?? new Float32Array(data.C)) : data.values(sequence, frame, fieldIndex(), values ?? undefined);
    paint((i, rgb) => spec.color(values[i], rgb));
    drawArrows?.set(spec.arrow ? (i, out) => spec.arrow(read, i, out) : null);
    shown = id;
    let lo = Infinity, hi = -Infinity, sum = 0;
    for (let i = 0; i < values.length; i++) { const v = values[i]; lo = Math.min(lo, v); hi = Math.max(hi, v); if (areas) sum += v * areas[i]; }
    const t = times()[frame], days = t / 86400;
    const number = spec.format ?? ((v, digits) => v.toFixed(digits));
    const range = [spec.label, `${number(lo, spec.digits ?? 0)} to ${number(hi, spec.digits ?? 0)}\u00a0${spec.unit}`];
    out.set([...(when ? when.readout(t, sequence) : [['day', days < 10 ? days.toFixed(2) : days.toFixed(1)]]), range, ...(areas && spec.mean !== false ? [['global mean', `${number(sum / totalArea, spec.meanDigits ?? 1)}\u00a0${spec.unit}`]] : [])]);
    if (scrub) scrub.value = frame;
  }

  const render = globe.render.bind(globe);
  globe.render = () => { if (data) load(); drawArrows?.face(); render(); };
  globe.fig.step = (dt) => {
    if (!data || !playing) return;
    position += dt * rate();
    if (position >= times().length) position = 0;
  };

  if (sequences.length > 1) choice(controls, { label: 'Show', options: sequences.map((s) => [s.label, s.name]), value: sequence, onChange: (v) => { sequence = v; position = 0; if (data) rebuild(); }, span: true });
  choice(controls, { label: 'Color shows', options: fields.map((f, k) => [f.label, String(k)]), value: '0', onChange: (v) => { field = Number(v); setKey(); shown = null; globe.fig.render(); }, span: true });
  const holder = document.createElement('div');
  holder.style.gridColumn = '1 / -1';
  controls.append(holder);
  const [toggle] = buttons(controls, [['Pause', (b) => { playing = !playing; globe.fig.play(playing); b.textContent = playing ? 'Pause' : 'Play'; }]]);
  const out = readout(controls);

  function setKey() { key.set(fields[field].legend); }
  function rebuild() {
    holder.replaceChildren();
    scrub = slider(holder, { label: 'Time', min: 0, max: times().length - 1, step: 1, value: 0, format: (v) => (when ? when.slider(times()[v], sequence) : `day ${(times()[v] / 86400).toFixed(times()[times().length - 1] < 864000 * 2 ? 2 : 1)}`), onInput: (v) => { position = v; playing = false; globe.fig.play(false); toggle.textContent = 'Play'; globe.fig.render(); } });
    shown = null;
    globe.fig.render();
  }

  setKey();
  loadFrames(url).then((loaded) => { data = loaded; rebuild(); globe.fig.play(playing); });
  return { globe, setNote: (text) => { top.textContent = text; } };
}

function arrowLayer(globe, cells, { coarse = 8, opacity = 0.85, color = 0xffffff, limb = 0.2 }) {
  const picks = [];
  for (const cell of new Grid(coarse, { relax: 0 })) {
    const c = cell.centerVertex;
    let best = 0, closest = -2;
    cells.forEach((other, i) => { const d = c.dot(other.centerVertex); if (d > closest) { closest = d; best = i; } });
    picks.push(best);
  }
  const n = picks.length, centers = new Float32Array(3 * n), basis = new Float32Array(6 * n), vectors = new Float32Array(3 * n), pair = [0, 0];
  const shownCenters = new Float32Array(3 * n), shownVectors = new Float32Array(3 * n), turn = new THREE.Quaternion(), view = new THREE.Vector3(), last = new THREE.Vector3();
  let changed = true;
  picks.forEach((i, k) => {
    const { x, y, z } = cells[i].centerVertex, rho = Math.hypot(x, y) || 1e-12;
    centers.set([x, y, z], 3 * k);
    basis.set([-y / rho, x / rho, 0, -z * x / rho, -z * y / rho, rho], 6 * k);
  });
  const draw = globe.dynamicArrows(n, { color, opacity }), lines = globe.layers.at(-1);
  lines.visible = false;
  return {
    set(arrowAt) {
      lines.visible = Boolean(arrowAt);
      if (!arrowAt) return;
      for (let k = 0; k < n; k++) {
        arrowAt(picks[k], pair);
        for (let c = 0; c < 3; c++) vectors[3 * k + c] = pair[0] * basis[6 * k + c] + pair[1] * basis[6 * k + 3 + c];
      }
      changed = true;
    },
    face() {
      if (!lines.visible) return;
      view.set(0, 0, 1).applyQuaternion(turn.copy(globe.spin).multiply(globe.base).invert());
      if (!changed && view.equals(last)) return;
      last.copy(view);
      changed = false;
      for (let k = 0; k < n; k++) {
        const front = centers[3 * k] * view.x + centers[3 * k + 1] * view.y + centers[3 * k + 2] * view.z >= limb;
        for (let c = 0; c < 3; c++) { shownCenters[3 * k + c] = front ? centers[3 * k + c] : 0; shownVectors[3 * k + c] = front ? vectors[3 * k + c] : 0; }
      }
      draw(shownCenters, shownVectors, 1);
    },
  };
}
