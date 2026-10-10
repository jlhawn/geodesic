import * as THREE from '../../js/three.module.js';
import { Grid } from '../../js/grid.module.js';
import { Globe } from '../globe.module.js';
import { slider, choice, buttons, legend, caption, readout, rampRGB, sequentialRGB, DARK_NEUTRAL, paletteVersion } from '../runtime.module.js';
import { loadFrames, loadTracks } from '../frames.module.js';
import { heightOf } from '../physics.module.js';
import { linear } from '../scene3d.module.js';

const EXAG = 30, A = 6.371e6, HOURS_PER_SECOND = 8;

export function mountStorm3d(root) {
  const controls = root.querySelector('.controls');
  const globe = new Globe(root, { height: 480, center: { lat: 48 * Math.PI / 180, lon: 165 * Math.PI / 180 }, lighting: { directional: 1.5, ambient: 1.9 } });
  const polygons = [...new Grid(16)].map((c) => ({ center: c.centerVertex, vertices: c.vertices }));
  const paint = globe.dynamicCells(polygons);
  caption(root, 'Computed in advance by the model’s layered core, 27 layers on the N = 16 grid: 192 parcels released on day 6 at two heights, followed hour by hour to day 11. Heights are drawn 30 times larger. Drag to turn the globe.');
  const pressure = ['ramp', 'surface pressure, 970 to 1030 hPa', 'cool', 'warm', 'rgb(52, 55, 62)'];
  const key = legend(root, [pressure, ['ramp', 'each parcel’s path, colored by its height, 0 to 5 km', 'rgb(52, 55, 62)', 'warm']]);
  let frames = null, tracks = null, hour = 0, playing = true, colorBy = 'height', group = 'all', shown = null, lines = [], heads = null, scrub = null;

  const point = (lat, lon, p) => { const r = 1 + EXAG * heightOf(p * 100) / A, c = Math.cos(lat * Math.PI / 180); return [r * c * Math.cos(lon * Math.PI / 180), r * c * Math.sin(lon * Math.PI / 180), r * Math.sin(lat * Math.PI / 180)]; };

  function build() {
    const { header, data } = tracks, P = header.parcels, H = header.times.length;
    for (let q = 0; q < P; q++) {
      const positions = new Float32Array(3 * H), colors = new Float32Array(3 * H);
      for (let t = 0; t < H; t++) { const o = (t * P + q) * 4; positions.set(point(data[o], data[o + 1], data[o + 2]), 3 * t); }
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
      geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
      const line = new THREE.Line(geometry, new THREE.LineBasicMaterial({ vertexColors: true }));
      globe.group.add(line);
      lines.push(line);
    }
    const headGeometry = new THREE.BufferGeometry();
    headGeometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(3 * P), 3));
    heads = new THREE.Points(headGeometry, new THREE.PointsMaterial({ color: 0xffffff, size: 4, sizeAttenuation: false }));
    globe.group.add(heads);
  }

  function recolor() {
    const { header, data } = tracks, P = header.parcels, H = header.times.length, rgb = [0, 0, 0];
    lines.forEach((line, q) => {
      const colors = line.geometry.attributes.color;
      for (let t = 0; t < H; t++) {
        const o = (t * P + q) * 4;
        if (colorBy === 'height') sequentialRGB(heightOf(data[o + 2] * 100) / 5000, rgb); else rampRGB((data[o + 3] - 270) / 50, rgb);
        colors.setXYZ(t, linear(rgb[0]), linear(rgb[1]), linear(rgb[2]));
      }
      colors.needsUpdate = true;
    });
  }

  let values = null, recolored = -1;
  function show() {
    if (!frames || !tracks) return;
    const { header, data } = tracks, P = header.parcels, H = header.times.length, t = Math.min(H - 1, Math.floor(hour));
    if (recolored !== `${paletteVersion}/${colorBy}`) { recolor(); recolored = `${paletteVersion}/${colorBy}`; shown = null; }
    const time = header.times[t], fTimes = frames.times('wave'), frame = fTimes.findIndex((s) => s >= time - 1);
    const id = `${t}/${group}/${paletteVersion}/${colorBy}`;
    if (id === shown) return;
    shown = id;
    values = frames.values('wave', Math.max(0, frame), 0, values ?? undefined);
    paint((i, rgb) => rampRGB(0.5 + (values[i] - 1000) / 60, rgb, DARK_NEUTRAL));
    const positions = heads.geometry.attributes.position;
    let rising = 0, risers = 0;
    lines.forEach((line, q) => {
      const startP = data[q * 4 + 2], visible = group === 'all' || (group === 'low' ? startP > 800 : startP <= 800);
      line.visible = visible;
      line.geometry.setDrawRange(0, t + 1);
      const o = (t * P + q) * 4;
      positions.setXYZ(q, ...(visible ? point(data[o], data[o + 1], data[o + 2]) : [0, 0, 0]));
      if (visible) { const climb = startP - data[o + 2]; rising = Math.max(rising, climb); if (climb > 100) risers++; }
    });
    positions.needsUpdate = true;
    out.set([['day', (time / 86400).toFixed(2)], ['highest climb so far', `${rising.toFixed(0)} hPa`], ['parcels risen more than 100 hPa', `${risers}`]]);
    if (scrub) scrub.value = t;
  }

  const render = globe.render.bind(globe);
  globe.render = () => { show(); render(); };
  globe.fig.step = (dt) => {
    if (!tracks || !playing) return false;
    const before = Math.floor(hour);
    hour += dt * HOURS_PER_SECOND;
    if (hour >= tracks.header.times.length) hour = 0;
    return Math.floor(hour) !== before;
  };

  choice(controls, { label: 'Color the paths by', options: [['height', 'height'], ['potential temperature', 'theta']], value: colorBy, onChange: (v) => { colorBy = v; key.set([pressure, v === 'height' ? ['ramp', 'each parcel’s path, colored by its height, 0 to 5 km', 'rgb(52, 55, 62)', 'warm'] : ['ramp', 'each parcel’s path, colored by its θ, 270 to 320 K', 'cool', 'warm', 'neutral']]); globe.fig.render(); }, span: true });
  choice(controls, { label: 'Parcels', options: [['all', 'all'], ['released near 900 hPa', 'low'], ['released near 700 hPa', 'mid']], value: group, onChange: (v) => { group = v; shown = null; globe.fig.render(); }, span: true });
  const holder = document.createElement('div'); holder.style.gridColumn = '1 / -1'; controls.append(holder);
  const [toggle] = buttons(controls, [['Pause', (b) => { playing = !playing; globe.fig.play(playing); b.textContent = playing ? 'Pause' : 'Play'; }]]);
  const out = readout(controls);

  Promise.all([loadFrames(new URL('../data/jw06_N16.bin', import.meta.url).href), loadTracks(new URL('../data/tracks_N16.bin', import.meta.url).href)]).then(([f, t]) => {
    frames = f; tracks = t; build();
    scrub = slider(holder, { label: 'Time', min: 0, max: t.header.times.length - 1, step: 1, value: 0, format: (v) => `day ${(t.header.times[v] / 86400).toFixed(2)}`, onInput: (v) => { hour = v; playing = false; globe.fig.play(false); toggle.textContent = 'Play'; globe.fig.render(); } });
    globe.fig.render();
    globe.fig.play(playing);
  });
}
