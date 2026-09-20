import { Grid } from "./grid.module.js";
import { initUnifiedViewer } from "./unifiedViewer.module.js";
import { createWindParticles } from "./windParticles.module.js";
import { seasonPhrase } from "./levels.module.js";

const WIND_MAX = { surface: 25, 1000: 30, 850: 40, 700: 40, 500: 50, 250: 70, 70: 100, 10: 150 };
const TEMP_RANGE = { surface: [240, 310], 1000: [240, 310], 850: [230, 300], 700: [220, 290], 500: [210, 280], 250: [190, 250], 70: [180, 240], 10: [200, 280] };
const REFERENCE_SPEED = { surface: 15, 1000: 20, 850: 25, 700: 25, 500: 30, 250: 40, 70: 50, 10: 60 };

const OVERLAYS = {
  wind: { label: 'Wind', unit: 'm/s', kind: 'sequential', field: 'speed', scale: 1, range: (level) => [0, WIND_MAX[level]] },
  temp: { label: 'Temp', unit: 'K', kind: 'sequential', field: 'temperature', scale: 1, range: (level) => TEMP_RANGE[level] },
  rh: { label: 'RH', unit: '%', kind: 'sequential', field: 'humidity', scale: 100, range: () => [0, 100] },
  precip: { label: 'Precip', unit: 'mm/day', kind: 'sequential', field: 'precipitation', scale: 1, range: () => [0, 30] },
  tpw: { label: 'TPW', unit: 'kg/m²', kind: 'sequential', field: 'water', scale: 1, range: () => [0, 60] },
  tcw: { label: 'TCW', unit: 'g/m²', kind: 'sequential', field: 'cloud', scale: 1000, range: () => [0, 500] },
  ice: { label: 'Ice', unit: 'm', kind: 'sequential', field: 'ice', scale: 1, range: () => [0, 3] },
  albedo: { label: 'Albedo', unit: '', kind: 'sequential', field: 'albedo', scale: 1, range: () => [0, 0.8] },
  mslp: { label: 'MSLP', unit: 'hPa', kind: 'diverging', field: 'ps', scale: 0.01, range: () => [960, 1060] },
  none: { label: 'None' },
};
const OVERLAY_NAMES = { wind: 'Wind speed', temp: 'Temperature', rh: 'Relative humidity', precip: 'Precipitation', tpw: 'Total precipitable water', tcw: 'Total cloud water', ice: 'Sea ice thickness', albedo: 'Surface albedo', mslp: 'Mean sea level pressure' };

/*
 * Palettes as sRGB stops. The sequential ones are perceptually uniform
 * colormaps that stay distinguishable under the common colour-vision
 * deficiencies; the diverging ones pair hues that do the same, around a
 * neutral gray so the white wind trails read on top.
 */
const PALETTES = {
  sequential: {
    viridis: [[0.267, 0.005, 0.329], [0.283, 0.141, 0.458], [0.254, 0.265, 0.530], [0.207, 0.372, 0.553], [0.164, 0.471, 0.558], [0.128, 0.567, 0.551], [0.135, 0.659, 0.518], [0.267, 0.749, 0.441], [0.478, 0.821, 0.318], [0.741, 0.873, 0.150], [0.993, 0.906, 0.144]],
    cividis: [[0.000, 0.135, 0.304], [0.127, 0.196, 0.416], [0.256, 0.263, 0.437], [0.367, 0.336, 0.446], [0.472, 0.409, 0.458], [0.578, 0.484, 0.463], [0.690, 0.562, 0.450], [0.807, 0.645, 0.411], [0.926, 0.734, 0.339], [0.994, 0.831, 0.243]],
    inferno: [[0.001, 0.000, 0.014], [0.087, 0.036, 0.209], [0.258, 0.039, 0.406], [0.416, 0.090, 0.433], [0.578, 0.148, 0.404], [0.735, 0.216, 0.330], [0.865, 0.317, 0.226], [0.954, 0.462, 0.109], [0.988, 0.645, 0.040], [0.965, 0.844, 0.146], [0.988, 0.998, 0.645]],
    magma: [[0.001, 0.000, 0.014], [0.098, 0.062, 0.259], [0.269, 0.060, 0.478], [0.446, 0.122, 0.507], [0.617, 0.183, 0.499], [0.792, 0.253, 0.446], [0.933, 0.372, 0.375], [0.987, 0.541, 0.383], [0.996, 0.719, 0.518], [0.987, 0.898, 0.729]],
    dusk: [[0.16, 0.10, 0.45], [0.15, 0.35, 0.65], [0.15, 0.60, 0.45], [0.55, 0.70, 0.20], [0.85, 0.30, 0.15]],
  },
  diverging: {
    'blue-gray-red': [[0.13, 0.30, 0.75], [0.35, 0.42, 0.62], [0.50, 0.50, 0.50], [0.68, 0.38, 0.32], [0.80, 0.20, 0.15]],
    'purple-gray-orange': [[0.33, 0.15, 0.53], [0.50, 0.45, 0.67], [0.50, 0.50, 0.50], [0.88, 0.51, 0.08], [0.70, 0.35, 0.02]],
    'teal-gray-brown': [[0.00, 0.40, 0.37], [0.35, 0.64, 0.60], [0.50, 0.50, 0.50], [0.75, 0.55, 0.30], [0.55, 0.32, 0.04]],
  },
};

const DEFAULTS = { overlay: 'wind', level: 'surface', animate: 'particles', isobars: 'off', isobarStep: 4, heightStep: 60, graticule: '15', projection: 'sphere', palettes: { sequential: 'viridis', diverging: 'blue-gray-red' }, panel: 'open' };

/*
 * The contour row draws isobars of surface pressure at the surface and
 * height contours of the pressure surface at any other level.
 */
const ISOLINES = {
  surface: { label: 'Isobars', unit: 'hPa', setting: 'isobarStep', steps: [2, 4, 5, 10, 25], field: (frame) => Float32Array.from(frame.ps, (p) => p / 100) },
  level: { label: 'Height lines', unit: 'm', setting: 'heightStep', steps: [20, 30, 60, 120, 240], field: (frame) => frame.height },
};
const isolinesFor = (level) => (level === 'surface' ? ISOLINES.surface : ISOLINES.level);

function loadSettings() {
  try { return { ...DEFAULTS, ...JSON.parse(localStorage.getItem('climate.settings') || '{}'), palettes: { ...DEFAULTS.palettes, ...JSON.parse(localStorage.getItem('climate.settings') || '{}').palettes } }; } catch { return { ...DEFAULTS }; }
}
function saveSettings(settings) {
  try { localStorage.setItem('climate.settings', JSON.stringify(settings)); } catch { /* storage unavailable */ }
}

function color(t, stops, out, at) {
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

function renderScale(stops, min, max, unit) {
  const canvas = document.getElementById('scale');
  const context = canvas.getContext('2d');
  const image = context.createImageData(canvas.width, 1);
  const rgb = new Uint8Array(3);
  for (let x = 0; x < canvas.width; x++) {
    color(x / (canvas.width - 1), stops, rgb, 0);
    image.data.set([rgb[0], rgb[1], rgb[2], 255], 4 * x);
  }
  context.putImageData(image, 0, 0);
  const ticks = document.querySelectorAll('#scaleTicks span');
  ticks[0].textContent = `≤ ${min} ${unit}`;
  ticks[1].textContent = `${(min + max) / 2}`;
  ticks[2].textContent = `≥ ${max} ${unit}`;
}

function levelLabel(level) { return level === 'surface' ? 'Surface' : `${level} hPa`; }

function formatDate(time) {
  const day = Math.floor(time / 86400), seconds = time - 86400 * day;
  const hh = String(Math.floor(seconds / 3600)).padStart(2, '0'), mm = String(Math.floor((seconds % 3600) / 60)).padStart(2, '0');
  return `Day ${day} ${hh}:${mm} · ${seasonPhrase(time)}`;
}

/*
 * The neighbours of a saved state among the run's files, for the ‹ ›
 * controls: files of the same tag listed by the server's directory
 * index, ordered by day.
 */
async function siblingStates(from) {
  const match = from.match(/^(.*\/)?([^/]+)_state_day(\d+)\.json$/);
  if (!match) return null;
  const [, directory = '', tag, day] = match;
  const html = await (await fetch(directory || './')).text();
  const days = [...html.matchAll(new RegExp(`href="${tag}_state_day(\\d+)\\.json"`, 'g'))].map((m) => Number(m[1])).sort((a, b) => a - b);
  const at = days.indexOf(Number(day));
  const url = (d) => `${directory}${tag}_state_day${String(d).padStart(3, '0')}.json`;
  return { prev: at > 0 ? url(days[at - 1]) : null, next: at >= 0 && at < days.length - 1 ? url(days[at + 1]) : null };
}

export default function runClimate({ N = null, from = null, workers = 1, paused = false } = {}) {
  const settings = loadSettings();
  const panel = document.getElementById('panel');
  const note = document.getElementById('note');
  const status = document.getElementById('status');
  let latest = null, grid = null, viewer = null, particles = null, arrows = null, isobars = null, graticule = null, rgb = null, running = !paused, siblings = null;
  const clock = [];
  function simulatedHoursPerMinute() {
    if (clock.length < 2) return null;
    const first = clock[0], last = clock[clock.length - 1];
    if (last.wall - first.wall < 2000) return null;
    return (last.time - first.time) / 3600 / ((last.wall - first.wall) / 60000);
  }

  const worker = new Worker(new URL('./model.worker.js', import.meta.url), { type: 'module' });

  function setup(cells) {
    grid = new Grid(cells);
    rgb = new Uint8Array(3 * grid.size);
    viewer = initUnifiedViewer(document.getElementById('globe'), grid, { backgroundColor: 0x151515, dynamicColors: true, controls: false, getColor: () => ({ r: 0.25, g: 0.25, b: 0.25 }) });
    arrows = viewer.addArrowLayer();
    isobars = viewer.addContourLayer();
    graticule = viewer.addGraticuleLayer();
    particles = createWindParticles(document.getElementById('globe'), viewer, grid);
    viewer.setProjection(settings.projection);
  }

  function paintOverlay() {
    const overlay = OVERLAYS[settings.overlay];
    const scaleRow = document.querySelector('.scaleRow');
    if (!overlay.field) {
      rgb.fill(LINEAR[40]);
      viewer.updateColors(rgb);
      scaleRow.style.visibility = 'hidden';
      document.getElementById('data').textContent = `Wind @ ${levelLabel(settings.level)} · no overlay`;
      return;
    }
    const stops = PALETTES[overlay.kind][settings.palettes[overlay.kind]] ?? Object.values(PALETTES[overlay.kind])[0];
    const [min, max] = overlay.range(settings.level);
    const values = latest[overlay.field];
    for (let i = 0; i < grid.size; i++) {
      color((values[i] * overlay.scale - min) / (max - min), stops, rgb, 3 * i);
      rgb[3 * i] = LINEAR[rgb[3 * i]]; rgb[3 * i + 1] = LINEAR[rgb[3 * i + 1]]; rgb[3 * i + 2] = LINEAR[rgb[3 * i + 2]];
    }
    viewer.updateColors(rgb);
    scaleRow.style.visibility = 'visible';
    renderScale(stops, min, max, overlay.unit);
    const columnField = ['ps', 'precipitation', 'water', 'cloud', 'ice', 'albedo'].includes(overlay.field);
    document.getElementById('data').textContent = `${OVERLAY_NAMES[settings.overlay]} @ ${columnField ? 'Surface' : levelLabel(settings.level)} · wind @ ${levelLabel(settings.level)}`;
  }

  function paintWind() {
    const reference = REFERENCE_SPEED[settings.level];
    arrows.setVisible(settings.animate === 'arrows');
    particles.setVisible(settings.animate === 'particles');
    if (settings.animate === 'arrows') arrows.update(latest.vector, { referenceSpeed: reference, stride: Math.ceil(grid.size / 4000) });
    if (settings.animate === 'particles') particles.setField(latest.vector, reference);
    note.textContent = settings.animate === 'particles' ? `trails brighten toward ${reference} m/s` : settings.animate === 'arrows' ? `full arrow at ${reference} m/s` : '';
  }

  function paintGraticule() {
    if (!graticule) return;
    graticule.setVisible(settings.graticule !== 'off');
    if (settings.graticule !== 'off') graticule.update(Number(settings.graticule));
  }

  function paintIsobars() {
    isobars.setVisible(settings.isobars === 'on');
    if (settings.isobars !== 'on') return;
    const isolines = isolinesFor(settings.level);
    isobars.setColor(settings.overlay === 'none' ? 0xffffff : 0x000000);
    isobars.update(isolines.field(latest), Number(settings[isolines.setting]));
  }

  function render() {
    for (const group of panel.querySelectorAll('.options[data-setting]')) {
      for (const button of group.querySelectorAll('button[data-value]')) button.classList.toggle('selected', button.dataset.value === String(settings[group.dataset.setting]));
    }
    const isolines = isolinesFor(settings.level);
    document.getElementById('isolineLabel').textContent = isolines.label;
    document.getElementById('isolineUnit').textContent = isolines.unit;
    const stepSelect = document.getElementById('isolineStep');
    stepSelect.replaceChildren(...isolines.steps.map((step) => { const option = document.createElement('option'); option.value = String(step); option.textContent = String(step); return option; }));
    stepSelect.value = String(settings[isolines.setting]);
    const paletteSelect = document.getElementById('palette');
    const overlay = OVERLAYS[settings.overlay];
    paletteSelect.style.display = overlay.kind ? '' : 'none';
    if (overlay.kind) {
      paletteSelect.replaceChildren(...Object.keys(PALETTES[overlay.kind]).map((name) => { const option = document.createElement('option'); option.value = name; option.textContent = name; return option; }));
      paletteSelect.value = settings.palettes[overlay.kind];
    }
    document.querySelector('[data-control="play"]').textContent = running ? '❚❚' : '▶';
    panel.classList.toggle('hidden', settings.panel !== 'open');
    if (!latest) return;
    paintOverlay();
    paintWind();
    paintIsobars();
    paintGraticule();
    const d = latest.diagnostics;
    document.getElementById('date').textContent = formatDate(latest.time);
    const rate = simulatedHoursPerMinute();
    document.getElementById('rate').textContent = !running ? 'paused' : rate === null ? 'measuring…' : `${rate.toFixed(1)} simulated hours per minute`;
    status.textContent = `ps ${(d.piMin / 100).toFixed(0)}–${(d.piMax / 100).toFixed(0)} hPa · Ts ${d.meanSurfaceT.toFixed(1)} K · solar ${d.absorbedSolar.toFixed(0)} / OLR ${d.outgoingLongwave.toFixed(0)} W/m² · LH ${d.latentHeat.toFixed(0)} SH ${d.sensibleHeat.toFixed(0)} · rain ${(d.precipitation * 86400).toFixed(2)} mm/d · TPW ${d.columnWater.toFixed(1)} · TCW ${(1000 * d.columnCloud).toFixed(0)} g/m² · ice ${(100 * d.iceFraction).toFixed(0)}% · albedo ${d.planetaryAlbedo.toFixed(2)} · N=${latest.N ?? ''} ${latest.workers > 1 ? `· ${latest.workers} workers` : ''}`;
  }

  function update(changes) {
    Object.assign(settings, changes);
    saveSettings(settings);
    if ('level' in changes) worker.postMessage({ type: 'level', level: settings.level });
    if ('projection' in changes && viewer) viewer.setProjection(settings.projection);
    render();
  }

  let ready = null;
  worker.onmessage = (event) => {
    const message = event.data;
    if (message.type === 'status') document.getElementById('date').textContent = message.text;
    if (message.type === 'ready') {
      ready = message;
      setup(message.N);
      document.getElementById('date').textContent = `model ready: ${message.cells} cells × ${message.layers} layers, dt ${message.dt} s, ${message.workers > 1 ? `${message.workers} workers` : 'one thread'}`;
    }
    if (message.type === 'frame') {
      latest = { ...message, N: ready.N, workers: ready.workers };
      clock.push({ wall: performance.now(), time: message.time });
      while (clock.length > 2 && clock[clock.length - 1].wall - clock[0].wall > 30000) clock.shift();
      render();
    }
  };
  worker.onerror = (error) => { document.getElementById('date').textContent = `worker error: ${error.message}`; };
  worker.postMessage({ type: 'start', N, from: from ? new URL(from, location.href).href : null, workers: crossOriginIsolated ? workers : 1, paused, level: settings.level });

  for (const group of panel.querySelectorAll('.options[data-setting]')) {
    group.addEventListener('click', (event) => {
      const button = event.target.closest('button[data-value]');
      if (button) update({ [group.dataset.setting]: button.dataset.value });
    });
  }
  document.getElementById('isolineStep').addEventListener('change', (event) => update({ [isolinesFor(settings.level).setting]: Number(event.target.value) }));
  document.getElementById('palette').addEventListener('change', (event) => {
    const kind = OVERLAYS[settings.overlay].kind;
    update({ palettes: { ...settings.palettes, [kind]: event.target.value } });
  });
  document.querySelector('[data-control="play"]').addEventListener('click', () => {
    running = !running;
    clock.length = 0;
    worker.postMessage({ type: running ? 'resume' : 'pause' });
    render();
  });
  const step = async (direction) => {
    siblings ??= from ? await siblingStates(from) : {};
    const target = siblings[direction];
    if (!target) { note.textContent = `no ${direction === 'prev' ? 'earlier' : 'later'} saved day`; return; }
    const params = new URLSearchParams(location.search);
    params.set('snapshot', target);
    location.search = params.toString();
  };
  document.querySelector('[data-control="prev"]').addEventListener('click', () => step('prev'));
  document.querySelector('[data-control="next"]').addEventListener('click', () => step('next'));
  document.getElementById('menu').addEventListener('click', () => update({ panel: settings.panel === 'open' ? 'closed' : 'open' }));
  render();
}
