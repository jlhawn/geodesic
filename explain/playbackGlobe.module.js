import { Grid } from '../js/grid.module.js';
import { Globe } from './globe.module.js';
import { slider, choice, buttons, legend, caption, readout, paletteVersion } from './runtime.module.js';
import { loadFrames } from './frames.module.js';

export function playbackGlobe(root, { url, N, center, sequences, fields, framesPerSecond = 8, note = '', when = null, globalMean = false }) {
  const controls = root.querySelector('.controls');
  const globe = new Globe(root, { height: 460, center, lighting: { directional: 1.5, ambient: 1.9 } });
  const cells = [...new Grid(N)], polygons = cells.map((c) => ({ center: c.centerVertex, vertices: c.vertices }));
  const areas = globalMean ? Float64Array.from(cells, (c) => c.area) : null, totalArea = areas ? areas.reduce((sum, a) => sum + a, 0) : 0;
  const paint = globe.dynamicCells(polygons);
  const top = caption(root, note);
  const key = legend(root, []);
  let data = null, sequence = sequences[0].name, field = 0, position = 0, playing = true, values = null, shown = null, scrub = null;
  const fieldIndex = () => data.header.fields.findIndex((f) => f.name === fields[field].name);
  const times = () => data.times(sequence);
  const rate = () => sequences.find((s) => s.name === sequence).framesPerSecond ?? framesPerSecond;

  function load() {
    const frame = Math.min(times().length - 1, Math.floor(position));
    const id = `${sequence}/${frame}/${field}/${paletteVersion}`;
    if (id === shown) return;
    values = data.values(sequence, frame, fieldIndex(), values ?? undefined);
    const spec = fields[field];
    paint((i, rgb) => spec.color(values[i], rgb));
    shown = id;
    let lo = Infinity, hi = -Infinity, sum = 0;
    for (let i = 0; i < values.length; i++) { const v = values[i]; lo = Math.min(lo, v); hi = Math.max(hi, v); if (areas) sum += v * areas[i]; }
    const t = times()[frame], days = t / 86400;
    const range = [spec.label, `${lo.toFixed(spec.digits ?? 0)} to ${hi.toFixed(spec.digits ?? 0)}\u00a0${spec.unit}`];
    out.set([...(when ? when.readout(t, sequence) : [['day', days < 10 ? days.toFixed(2) : days.toFixed(1)]]), range, ...(areas ? [['global mean', `${(sum / totalArea).toFixed(spec.meanDigits ?? 1)}\u00a0${spec.unit}`]] : [])]);
    if (scrub) scrub.value = frame;
  }

  const render = globe.render.bind(globe);
  globe.render = () => { if (data) load(); render(); };
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
