import { Grid } from "./grid.module.js";
import { initUnifiedViewer } from "./unifiedViewer.module.js";

const FIELDS = {
  ps: { label: 'Surface pressure (hPa)', scale: 0.01, min: 960, max: 1060, palette: 'diverging' },
  ts: { label: 'Surface temperature (K)', scale: 1, min: 240, max: 310, palette: 'sequential' },
  wind: { label: 'Surface wind speed (m/s)', scale: 1, min: 0, max: 25, palette: 'sequential' },
};

function color(t, palette, out, at) {
  const x = Math.max(0, Math.min(1, t));
  let r, g, b;
  if (palette === 'diverging') {
    if (x < 0.5) { const s = x / 0.5; r = 0.15 + 0.8 * s; g = 0.3 + 0.65 * s; b = 0.85 + 0.1 * s; }
    else { const s = (x - 0.5) / 0.5; r = 0.95; g = 0.95 - 0.7 * s; b = 0.95 - 0.8 * s; }
  } else {
    r = 0.1 + 0.85 * x; g = 0.1 + 0.7 * Math.sin(Math.PI * x); b = 0.9 - 0.8 * x;
  }
  out[at] = Math.round(255 * r); out[at + 1] = Math.round(255 * g); out[at + 2] = Math.round(255 * b);
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
  let latest = null, paused = false;

  function paint() {
    if (!latest) return;
    const spec = FIELDS[select.value];
    const values = latest[select.value];
    for (let i = 0; i < cells; i++) {
      color((values[i] * spec.scale - spec.min) / (spec.max - spec.min), spec.palette, rgb, 3 * i);
    }
    viewer.updateColors(rgb);
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
  pauseButton.addEventListener('click', () => {
    paused = !paused;
    worker.postMessage({ type: paused ? 'pause' : 'resume' });
    pauseButton.textContent = paused ? 'Resume' : 'Pause';
  });
}
