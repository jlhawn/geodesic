import { Grid } from '../js/grid.module.js';
import { Globe } from './globe.module.js';
import { slider, choice, buttons, legend, caption, readout, paletteVersion } from './runtime.module.js';
import { loadFrames } from './frames.module.js';

export function playbackGlobe(root, { url, N, center, sequences, fields, framesPerSecond = 8, note = '' }) {
  const controls = root.querySelector('.controls');
  const globe = new Globe(root, { height: 460, center, lighting: { directional: 1.5, ambient: 1.9 } });
  const polygons = [...new Grid(N)].map((c) => ({ center: c.centerVertex, vertices: c.vertices }));
  const paint = globe.dynamicCells(polygons);
  const top = caption(root, note);
  const key = legend(root, []);
  let data = null, sequence = sequences[0].name, field = 0, position = 0, playing = true, values = null, shown = null, scrub = null;
  const fieldIndex = () => data.header.fields.findIndex((f) => f.name === fields[field].name);
  const times = () => data.times(sequence);

  function load() {
    const frame = Math.min(times().length - 1, Math.floor(position));
    const id = `${sequence}/${frame}/${field}/${paletteVersion}`;
    if (id === shown) return;
    values = data.values(sequence, frame, fieldIndex(), values ?? undefined);
    const spec = fields[field];
    paint((i, rgb) => spec.color(values[i], rgb));
    shown = id;
    let lo = Infinity, hi = -Infinity;
    for (const v of values) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
    const t = times()[frame], days = t / 86400;
    out.set([['day', days < 10 ? days.toFixed(2) : days.toFixed(1)], [spec.label, `${lo.toFixed(spec.digits ?? 0)} to ${hi.toFixed(spec.digits ?? 0)} ${spec.unit}`]]);
    if (scrub) scrub.value = frame;
  }

  const render = globe.render.bind(globe);
  globe.render = () => { if (data) load(); render(); };
  globe.fig.step = (dt) => {
    if (!data || !playing) return;
    position += dt * framesPerSecond;
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
    scrub = slider(holder, { label: 'Time', min: 0, max: times().length - 1, step: 1, value: 0, format: (v) => `day ${(times()[v] / 86400).toFixed(times()[times().length - 1] < 864000 * 2 ? 2 : 1)}`, onInput: (v) => { position = v; playing = false; globe.fig.play(false); toggle.textContent = 'Play'; globe.fig.render(); } });
    shown = null;
    globe.fig.render();
  }

  setKey();
  loadFrames(url).then((loaded) => { data = loaded; rebuild(); globe.fig.play(true); });
  return { globe, setNote: (text) => { top.textContent = text; } };
}
