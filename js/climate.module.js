import { Grid } from "./grid.module.js";
import { initUnifiedViewer } from "./unifiedViewer.module.js";
import { createWindParticles } from "./windParticles.module.js";

const FIELDS = {
  ps: { label: 'Surface pressure (hPa)', scale: 0.01, min: 960, max: 1060, center: 1013, palette: 'diverging' },
  ts: { label: 'Surface temperature (K)', scale: 1, min: 240, max: 310, palette: 'sequential' },
  wind: { label: 'Surface wind speed (m/s)', scale: 1, min: 0, max: 25, palette: 'sequential', arrows: 'surface' },
  jet: { label: 'Wind speed at 250 hPa (m/s)', scale: 1, min: 0, max: 60, palette: 'sequential', arrows: 'jet' },
  jetZonal: { label: 'Zonal wind at 250 hPa (m/s)', scale: 1, min: -40, max: 40, center: 0, palette: 'diverging', arrows: 'jet' },
  none: { label: 'No field' },
};

const WINDS = {
  surface: { vector: 'windVector', referenceSpeed: 15, label: 'surface wind' },
  jet: { vector: 'jetVector', referenceSpeed: 40, label: 'wind at 250 hPa' },
};

/*
 * Palettes are kept mid-dark so that the white wind trails read on top
 * of them: the diverging one runs blue → neutral gray → red, the
 * sequential one indigo → blue → teal → olive → red.
 */
const PALETTES = {
  diverging: [[0.13, 0.30, 0.75], [0.35, 0.42, 0.62], [0.50, 0.50, 0.50], [0.68, 0.38, 0.32], [0.80, 0.20, 0.15]],
  sequential: [[0.16, 0.10, 0.45], [0.15, 0.35, 0.65], [0.15, 0.60, 0.45], [0.55, 0.70, 0.20], [0.85, 0.30, 0.15]],
};

function color(t, palette, out, at) {
  const stops = PALETTES[palette];
  const x = Math.max(0, Math.min(1, t)) * (stops.length - 1);
  const k = Math.min(stops.length - 2, Math.floor(x)), f = x - k;
  const a = stops[k], b = stops[k + 1];
  out[at] = Math.round(255 * (a[0] + f * (b[0] - a[0])));
  out[at + 1] = Math.round(255 * (a[1] + f * (b[1] - a[1])));
  out[at + 2] = Math.round(255 * (a[2] + f * (b[2] - a[2])));
}

// The viewer takes vertex colors in linear light and encodes to sRGB on
// output, so palette bytes (sRGB) are linearized before upload.
const LINEAR = Uint8Array.from({ length: 256 }, (_, v) => {
  const s = v / 255;
  return Math.round(255 * (s <= 0.04045 ? s / 12.92 : ((s + 0.055) / 1.055) ** 2.4));
});

function paintCells(rgb, cells, palette, value) {
  for (let i = 0; i < cells; i++) {
    color(value(i), palette, rgb, 3 * i);
    rgb[3 * i] = LINEAR[rgb[3 * i]]; rgb[3 * i + 1] = LINEAR[rgb[3 * i + 1]]; rgb[3 * i + 2] = LINEAR[rgb[3 * i + 2]];
  }
}

function renderLegend(spec, min, max) {
  document.getElementById('legend').style.display = spec.palette ? 'block' : 'none';
  if (!spec.palette) return;
  const canvas = document.querySelector('#legend canvas');
  const context = canvas.getContext('2d');
  const image = context.createImageData(canvas.width, 1);
  const rgb = new Uint8Array(3);
  for (let x = 0; x < canvas.width; x++) {
    color(x / (canvas.width - 1), spec.palette, rgb, 0);
    image.data.set([rgb[0], rgb[1], rgb[2], 255], 4 * x);
  }
  context.putImageData(image, 0, 0);
  document.getElementById('legendTitle').textContent = spec.label;
  const ticks = document.querySelectorAll('#legendTicks span');
  const digits = max - min < 10 ? 1 : 0;
  ticks[0].textContent = `≤ ${min.toFixed(digits)}`;
  ticks[1].textContent = ((min + max) / 2).toFixed(digits);
  ticks[2].textContent = `≥ ${max.toFixed(digits)}`;
}

/*
 * The wind layer draws the wind at the chosen level as particle trails
 * or arrows; choosing a wind field switches the level to that field's.
 * Arrows are strided so that no more than a few thousand are drawn.
 */
function createWindLayer(container, viewer, grid) {
  const arrows = viewer.addArrowLayer();
  const particles = createWindParticles(container, viewer, grid);
  const mode = document.getElementById('windLayer');
  const level = document.getElementById('windLevel');
  const note = document.getElementById('legendNote');
  return function paint(frame, spec) {
    if (spec.arrows) level.value = spec.arrows;
    const kind = WINDS[level.value];
    const vectors = frame[kind.vector];
    const shown = vectors ? mode.value : 'none';
    arrows.setVisible(shown === 'arrows');
    particles.setVisible(shown === 'particles');
    if (shown === 'arrows') {
      arrows.update(vectors, { referenceSpeed: kind.referenceSpeed, stride: Math.ceil(grid.size / 4000) });
      note.textContent = `Arrows point downwind (${kind.label}); full length at ${kind.referenceSpeed} m/s`;
    } else if (shown === 'particles') {
      particles.setField(vectors, kind.referenceSpeed);
      note.textContent = `Particles drift downwind (${kind.label}); brighter trails are faster`;
    } else {
      note.textContent = vectors ? '' : 'No wind vectors in this snapshot';
    }
  };
}

function addOption(select, value, label) {
  const option = document.createElement('option');
  option.value = value;
  option.textContent = label;
  select.appendChild(option);
}

/*
 * Shows one saved snapshot (the JSON the emergence runs write) instead
 * of running the model: ?snapshot=<url>. Sibling days of the same run
 * are reachable with the prev/next buttons.
 */
export async function showSnapshot(url) {
  const response = await fetch(url);
  const snapshot = await response.json();
  const grid = new Grid(snapshot.N);
  const cells = grid.size;
  const rgb = new Uint8Array(3 * cells);
  const viewer = initUnifiedViewer(document.getElementById('globe'), grid, { backgroundColor: 0x151515, dynamicColors: true, getColor: () => ({ r: 0.25, g: 0.25, b: 0.25 }) });
  const readout = document.getElementById('readout');
  const select = document.getElementById('field');
  addOption(select, 'jetZonal', 'Zonal wind at 250 hPa');
  const paintWind = createWindLayer(document.getElementById('globe'), viewer, grid);
  const pauseButton = document.getElementById('pause');
  pauseButton.textContent = 'Previous day';
  const nextButton = document.createElement('button');
  nextButton.textContent = 'Next day';
  pauseButton.after(nextButton);
  function paint() {
    const spec = FIELDS[select.value];
    if (!spec.palette) {
      rgb.fill(LINEAR[40]);
      viewer.updateColors(rgb);
      renderLegend(spec);
      paintWind(snapshot, spec);
      readout.textContent = `snapshot day ${snapshot.day}, N=${snapshot.N}`;
      return;
    }
    const values = snapshot[select.value];
    const sorted = Float64Array.from(values, (v) => v * spec.scale).sort();
    const lo = sorted[Math.floor(0.02 * cells)], hi = sorted[Math.ceil(0.98 * cells) - 1];
    let min = lo, max = hi;
    if (spec.palette === 'diverging') { const half = Math.max(Math.abs(lo - spec.center), Math.abs(hi - spec.center)); min = spec.center - half; max = spec.center + half; }
    paintCells(rgb, cells, spec.palette, (i) => (values[i] * spec.scale - min) / (max - min));
    viewer.updateColors(rgb);
    renderLegend(spec, min, max);
    paintWind(snapshot, spec);
    readout.textContent = `${spec.label} — snapshot day ${snapshot.day}, N=${snapshot.N} · field range ${sorted[0].toFixed(1)} to ${sorted[cells - 1].toFixed(1)}`;
  }
  const step = (delta) => {
    const match = url.match(/day(\d{3})\.json$/);
    if (!match) return;
    const day = Math.max(0, Number(match[1]) + delta);
    location.search = `?snapshot=${encodeURIComponent(url.replace(/day\d{3}\.json$/, `day${String(day).padStart(3, '0')}.json`))}`;
  };
  pauseButton.addEventListener('click', () => step(-10));
  nextButton.addEventListener('click', () => step(10));
  select.addEventListener('change', paint);
  document.getElementById('windLayer').addEventListener('change', paint);
  document.getElementById('windLevel').addEventListener('change', paint);
  paint();
}

export default function runClimate(N = 16) {
  const grid = new Grid(N);
  const cells = grid.size;
  const rgb = new Uint8Array(3 * cells).fill(60);
  const viewer = initUnifiedViewer(document.getElementById('globe'), grid, {
    backgroundColor: 0x151515,
    dynamicColors: true,
    getColor: () => ({ r: 0.25, g: 0.25, b: 0.25 }),
  });

  const readout = document.getElementById('readout');
  const select = document.getElementById('field');
  const pauseButton = document.getElementById('pause');
  const paintWind = createWindLayer(document.getElementById('globe'), viewer, grid);
  let latest = null, paused = false;

  function paint() {
    if (!latest) return;
    const spec = FIELDS[select.value];
    const values = latest[select.value];
    if (spec.palette) paintCells(rgb, cells, spec.palette, (i) => (values[i] * spec.scale - spec.min) / (spec.max - spec.min));
    else rgb.fill(LINEAR[40]);
    viewer.updateColors(rgb);
    renderLegend(spec, spec.min, spec.max);
    paintWind(latest, spec);
    const d = latest.diagnostics;
    readout.textContent = `${spec.label} — day ${latest.day.toFixed(2)} · ps ${(d.piMin / 100).toFixed(1)}–${(d.piMax / 100).toFixed(1)} hPa · Ts ${d.meanSurfaceT.toFixed(2)} K · max wind ${d.maxWind.toFixed(1)} m/s · solar ${d.absorbedSolar.toFixed(0)} / OLR ${d.outgoingLongwave.toFixed(0)} W/m²`;
  }

  const worker = new Worker(new URL('./model.worker.js', import.meta.url), { type: 'module' });
  worker.onmessage = (event) => {
    const message = event.data;
    if (message.type === 'ready') readout.textContent = `model ready: ${message.cells} cells × ${message.layers} layers, dt ${message.dt} s`;
    if (message.type === 'frame') { latest = message; paint(); }
  };
  worker.onerror = (error) => { readout.textContent = `worker error: ${error.message}`; };
  worker.postMessage({ type: 'start', N });

  select.addEventListener('change', paint);
  document.getElementById('windLayer').addEventListener('change', paint);
  document.getElementById('windLevel').addEventListener('change', paint);
  pauseButton.addEventListener('click', () => {
    paused = !paused;
    worker.postMessage({ type: paused ? 'pause' : 'resume' });
    pauseButton.textContent = paused ? 'Resume' : 'Pause';
  });
}
