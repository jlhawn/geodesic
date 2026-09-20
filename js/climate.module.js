import { Grid } from "./grid.module.js";
import { initUnifiedViewer } from "./unifiedViewer.module.js";
import { createWindParticles } from "./windParticles.module.js";
import { seasonPhrase } from "./levels.module.js";
import { sunDirection, DAY, YEAR } from "./physics/radiation.module.js";
import { createDisplayClock } from "./displayClock.module.js";
import { listSnapshots, saveSnapshot, getSnapshot, renameSnapshot, deleteSnapshot, cloneSnapshot } from "./snapshots.module.js";

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
  swdown: { label: 'SW↓', unit: 'W/m²', kind: 'sequential', field: 'shortwave', scale: 1, range: () => [0, 1200] },
  olr: { label: 'OLR', unit: 'W/m²', kind: 'sequential', field: 'longwave', scale: 1, range: () => [100, 320] },
  mslp: { label: 'MSLP', unit: 'hPa', kind: 'diverging', field: 'ps', scale: 0.01, range: () => [960, 1060] },
  none: { label: 'None' },
};
const OVERLAY_NAMES = { wind: 'Wind speed', temp: 'Temperature', rh: 'Relative humidity', precip: 'Precipitation', tpw: 'Precipitable water', tcw: 'Cloud water', ice: 'Sea ice thickness', albedo: 'Surface albedo', swdown: 'Sunlight reaching the surface', olr: 'Outgoing longwave at the top', mslp: 'Sea-level pressure' };

/*
 * The cloud view: open water is ocean blue, ice whitens with thickness,
 * and cloud is white composited on top with an opacity that rises with
 * the column's cloud water, 1 − exp(−TCW / CLOUD_OPACITY_SCALE), so
 * clear sky is transparent and 40 g/m² is two-thirds opaque.
 */
const OCEAN_COLOR = [0.05, 0.22, 0.45], ICE_COLOR = [0.85, 0.90, 0.95], CLOUD_COLOR = [1, 1, 1], CLOUD_OPACITY_SCALE = 40;
const cloudOpacity = (grams) => 1 - Math.exp(-Math.max(0, grams) / CLOUD_OPACITY_SCALE);

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

const DEFAULTS = { view: 'data', overlay: 'wind', level: 'surface', animate: 'particles', isobars: 'off', isobarStep: 5, heightStep: 60, graticule: '15', projection: 'sphere', palettes: { sequential: 'viridis', diverging: 'blue-gray-red' }, panel: 'open' };

/*
 * The contour row draws isobars of surface pressure at the surface and
 * height contours of the pressure surface at any other level.
 */
const ISOLINES = {
  surface: { label: 'Isobars', unit: 'hPa', setting: 'isobarStep', steps: [1, 2.5, 5, 10], field: (frame) => Float32Array.from(frame.ps, (p) => p / 100) },
  level: { label: 'Height lines', unit: 'm', setting: 'heightStep', steps: [20, 30, 60, 120, 240], field: (frame) => frame.height },
};
const isolinesFor = (level) => (level === 'surface' ? ISOLINES.surface : ISOLINES.level);

function loadSettings() {
  try {
    const stored = JSON.parse(localStorage.getItem('climate.settings') || '{}');
    const settings = { ...DEFAULTS, ...stored, palettes: { ...DEFAULTS.palettes, ...stored.palettes } };
    if (settings.overlay === 'clouds') { settings.overlay = 'none'; settings.view = 'space'; }
    if (!OVERLAYS[settings.overlay]) settings.overlay = DEFAULTS.overlay;
    return settings;
  } catch { return { ...DEFAULTS }; }
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
/*
 * The saved runs the server lists in its runs/ directory index, as
 * built-in snapshots the page can download into its own store.
 */
async function builtinSnapshots() {
  try {
    const html = await (await fetch('runs/')).text();
    const files = [...new Set([...html.matchAll(/href="([^"]+_state_day\d+\.json)"/g)].map((m) => m[1]))].sort();
    return files.map((file) => ({ file, url: new URL(`runs/${file}`, location.href).href, name: file.replace(/\.json$/, '') }));
  } catch { return []; }
}

const HEIGHT_OVERLAYS = new Set(['wind', 'temp', 'rh', 'none']);

const VIEW_NOTES = [
  ['Mode', 'Data paints the chosen overlay on an evenly lit globe. Satellite renders the planet as it would look from space: ocean, ice and cloud lit by the sun in its true direction for the model date and time, a dark ambient on the night side, and the stars turning behind it once a sidereal day.'],
  ['Wind animation', 'Particles trace the wind at the chosen height as fading trails, brighter where it blows faster; Vectors draw one arrow per cell; None hides the motion.'],
  ['Height', 'The pressure level shown by the wind, temperature and humidity views and followed by the animation: Sfc is the lowest layer, about 60 m up; the others are hPa. Column views hide it and use the surface wind.'],
  ['Wind speed', 'Speed at the chosen height.'],
  ['Temperature', 'Air temperature at the chosen height.'],
  ['Relative humidity', 'At the chosen height.'],
  ['Sea-level pressure', 'Surface pressure; with no terrain it is the sea-level pressure.'],
  ['Precipitation', 'Rain rate over the last frame, from convection and from cloud that rained out.'],
  ['Precipitable water', 'All the vapour in the column, as the depth of rain it would make.'],
  ['Cloud water', 'All the condensed water in the column.'],
  ['Sea ice', 'Sea-ice thickness.'],
  ['Albedo', 'The surface albedo for diffuse light: 0.06 over water, rising to 0.5 over half a metre of ice.'],
  ['Surface sunlight', 'Shortwave reaching the surface, direct and diffuse, before the surface reflects its share.'],
  ['Outgoing longwave', 'Infrared leaving the top of the atmosphere: low over cold cloud tops and the poles, high over clear warm regions.'],
  ['Isobars / Height lines', 'Contours of surface pressure at the surface, of geopotential height on a pressure level, at the chosen interval.'],
  ['Graticule', 'Parallels and meridians at the chosen spacing; the meridians stop at the outermost parallel.'],
  ['Projection', 'The orthographic globe, or the Equal Earth map; both can be dragged to any orientation.'],
  ['Snapshots', 'Save the paused state in this browser, restore it later, or download one of the runs saved on the server.'],
];

export default function runClimate({ N = null, from = null, workers = 1, engine = 'cpu', paused = false } = {}) {
  const settings = loadSettings();
  const panel = document.getElementById('panel');
  const activeLevel = () => (settings.view === 'space' || HEIGHT_OVERLAYS.has(settings.overlay) ? settings.level : 'surface');
  let latest = null, grid = null, viewer = null, particles = null, arrows = null, isobars = null, graticule = null, rgb = null, running = !paused;
  const clock = [];
  function simulatedHoursPerMinute() {
    if (clock.length < 2) return null;
    const first = clock[0], last = clock[clock.length - 1];
    if (last.wall - first.wall < 2000) return null;
    return (last.time - first.time) / 3600 / ((last.wall - first.wall) / 60000);
  }

  const display = createDisplayClock();
  function recentFrames(count = 5) {
    if (clock.length < 2) return { rate: 0, interval: 500 };
    const from = Math.max(0, clock.length - count);
    const first = clock[from], last = clock[clock.length - 1];
    if (last.wall <= first.wall) return { rate: 0, interval: 500 };
    return { rate: (last.time - first.time) / (last.wall - first.wall), interval: (last.wall - first.wall) / (clock.length - 1 - from) };
  }
  function tick(now) {
    requestAnimationFrame(tick);
    if (!latest) return;
    const time = display.advance(now, latest.time, { ...recentFrames(), running });
    if (running) document.getElementById('date').textContent = formatDate(time);
    if (settings.view === 'space' && viewer) viewer.setSpace({ enabled: true, sun: sunDirection(time), sidereal: 2 * Math.PI * time * (1 / DAY + 1 / YEAR) });
  }
  requestAnimationFrame(tick);

  const worker = new Worker(new URL('./model.worker.js', import.meta.url), { type: 'module' });

  function setup(cells) {
    grid = new Grid(cells);
    rgb = new Uint8Array(3 * grid.size);
    viewer = initUnifiedViewer(document.getElementById('globe'), grid, { backgroundColor: 0x151515, dynamicColors: true, controls: false, getColor: () => ({ r: 0.25, g: 0.25, b: 0.25 }) });
    arrows = viewer.addArrowLayer({ opacity: 0.5 });
    isobars = viewer.addContourLayer({ opacity: 0.25 });
    graticule = viewer.addGraticuleLayer({ opacity: 0.25 });
    particles = createWindParticles(document.getElementById('globe'), viewer, grid);
    viewer.setProjection(settings.projection);
  }

  function paintSatellite() {
    const cloud = latest.cloud, ice = latest.ice;
    for (let i = 0; i < grid.size; i++) {
      const frozen = Math.min(1, ice[i] / 0.5), opacity = cloudOpacity(cloud[i] * 1000);
      for (let j = 0; j < 3; j++) {
        const base = OCEAN_COLOR[j] + frozen * (ICE_COLOR[j] - OCEAN_COLOR[j]);
        rgb[3 * i + j] = LINEAR[Math.round(255 * (base + opacity * (CLOUD_COLOR[j] - base)))];
      }
    }
    viewer.updateColors(rgb);
    document.querySelector('.scaleRow').classList.add('hidden');
    document.getElementById('data').textContent = `Satellite view · wind @ ${levelLabel(activeLevel())}`;
  }

  function paintOverlay() {
    if (settings.view === 'space') { paintSatellite(); return; }
    const overlay = OVERLAYS[settings.overlay];
    const scaleRow = document.querySelector('.scaleRow');
    if (!overlay.field) {
      rgb.fill(LINEAR[40]);
      viewer.updateColors(rgb);
      scaleRow.classList.add('hidden');
      document.getElementById('data').textContent = `Wind @ ${levelLabel(activeLevel())} · no overlay`;
      return;
    }
    const [min, max] = overlay.range(activeLevel());
    const values = latest[overlay.field];
    const stops = PALETTES[overlay.kind][settings.palettes[overlay.kind]] ?? Object.values(PALETTES[overlay.kind])[0];
    for (let i = 0; i < grid.size; i++) {
      color((values[i] * overlay.scale - min) / (max - min), stops, rgb, 3 * i);
      rgb[3 * i] = LINEAR[rgb[3 * i]]; rgb[3 * i + 1] = LINEAR[rgb[3 * i + 1]]; rgb[3 * i + 2] = LINEAR[rgb[3 * i + 2]];
    }
    viewer.updateColors(rgb);
    scaleRow.classList.remove('hidden');
    renderScale(stops, min, max, overlay.unit);
    const columnField = ['ps', 'precipitation', 'water', 'cloud', 'ice', 'albedo', 'shortwave', 'longwave'].includes(overlay.field);
    document.getElementById('data').textContent = columnField ? `${OVERLAY_NAMES[settings.overlay]} · wind @ ${levelLabel(activeLevel())}` : `${OVERLAY_NAMES[settings.overlay]} @ ${levelLabel(activeLevel())}`;
  }

  function paintWind() {
    const reference = REFERENCE_SPEED[activeLevel()];
    arrows.setVisible(settings.animate === 'arrows');
    particles.setVisible(settings.animate === 'particles');
    if (settings.animate === 'arrows') arrows.update(latest.vector, { referenceSpeed: reference });
    if (settings.animate === 'particles') particles.setField(latest.vector, reference);
    const note = settings.animate === 'particles' ? `trails brighten toward ${reference} m/s` : settings.animate === 'arrows' ? `full arrow at ${reference} m/s` : '';
    if (note) document.getElementById('data').textContent += ` · ${note}`;
  }

  function paintGraticule() {
    if (!graticule) return;
    graticule.setVisible(settings.graticule !== 'off');
    if (settings.graticule !== 'off') graticule.update(Number(settings.graticule));
  }

  function paintIsobars() {
    isobars.setVisible(settings.isobars === 'on');
    if (settings.isobars !== 'on') return;
    const isolines = isolinesFor(activeLevel());
    isobars.update(isolines.field(latest), Number(settings[isolines.setting]));
  }

  function fillSelect(select, values, value) {
    if ([...select.options].map((option) => option.value).join('\n') !== values.join('\n')) {
      select.replaceChildren(...values.map((v) => { const option = document.createElement('option'); option.value = v; option.textContent = v; return option; }));
    }
    if (select.value !== value) select.value = value;
  }

  function fillSegments(group, values, labels) {
    if ([...group.children].map((button) => button.dataset.value).join('\n') !== values.join('\n')) {
      group.replaceChildren(...values.map((v, k) => { const button = document.createElement('button'); button.dataset.value = v; button.textContent = labels[k]; return button; }));
    }
  }

  function isolineChoice() {
    const isolines = isolinesFor(activeLevel());
    return settings.isobars === 'off' ? 'off' : String(settings[isolines.setting]);
  }

  function render() {
    const isolines = isolinesFor(activeLevel());
    if (!isolines.steps.includes(settings[isolines.setting])) {
      settings[isolines.setting] = isolines.steps.reduce((a, b) => (Math.abs(b - settings[isolines.setting]) < Math.abs(a - settings[isolines.setting]) ? b : a));
      saveSettings(settings);
    }
    fillSegments(document.getElementById('isolineSegments'), ['off', ...isolines.steps.map(String)], ['Off', ...isolines.steps.map(String)]);
    for (const group of panel.querySelectorAll('.options[data-setting]')) {
      const current = group.dataset.setting === 'isolines' ? isolineChoice() : String(settings[group.dataset.setting]);
      for (const button of group.querySelectorAll('button[data-value]')) button.classList.toggle('selected', button.dataset.value === current);
    }
    const space = settings.view === 'space';
    const heights = space || HEIGHT_OVERLAYS.has(settings.overlay);
    document.getElementById('heightLabel').classList.toggle('hidden', !heights);
    document.getElementById('heightOptions').classList.toggle('hidden', !heights);
    document.getElementById('overlayLabel').classList.toggle('hidden', space);
    document.getElementById('overlayOptions').classList.toggle('hidden', space);
    document.getElementById('isolineLabel').textContent = isolines.label;
    document.getElementById('isolineUnit').textContent = isolines.unit;
    const paletteSelect = document.getElementById('palette');
    const overlay = OVERLAYS[settings.overlay];
    paletteSelect.style.display = overlay.kind && PALETTES[overlay.kind] ? '' : 'none';
    if (overlay.kind && PALETTES[overlay.kind]) fillSelect(paletteSelect, Object.keys(PALETTES[overlay.kind]), settings.palettes[overlay.kind]);
    document.querySelector('[data-control="play"]').textContent = running ? '❚❚' : '▶';
    panel.classList.toggle('hidden', settings.panel !== 'open');
    if (viewer) viewer.setSpace({ enabled: space });
    if (!latest) return;
    paintOverlay();
    paintWind();
    paintIsobars();
    paintGraticule();
    const d = latest.diagnostics;
    document.getElementById('date').textContent = formatDate(latest.time);
    const rate = simulatedHoursPerMinute();
    document.getElementById('rate').textContent = !running ? 'paused' : rate === null ? 'measuring…' : `${rate.toFixed(1)} simulated hours per minute`;
    if (!document.getElementById('modelModal').classList.contains('hidden')) renderModelDetails();
  }

  function renderModelDetails() {
    if (!latest) return;
    const d = latest.diagnostics;
    const rows = [
      ['Surface pressure', `<b>${(d.piMin / 100).toFixed(0)}–${(d.piMax / 100).toFixed(0)} hPa</b> — the lowest and highest on the globe right now.`],
      ['Surface temperature', `<b>${d.meanSurfaceT.toFixed(1)} K</b> — area-weighted global mean of the skin temperature.`],
      ['Absorbed solar', `<b>${d.absorbedSolar.toFixed(0)} W/m²</b> — global mean sunlight absorbed by atmosphere and surface.`],
      ['Outgoing longwave', `<b>${d.outgoingLongwave.toFixed(0)} W/m²</b> — infrared leaving the top; absorbed solar minus this is the planet's energy imbalance, <b>${(d.absorbedSolar - d.outgoingLongwave).toFixed(0)} W/m²</b>.`],
      ['Latent heat', `<b>${d.latentHeat.toFixed(0)} W/m²</b> — heat leaving the surface as evaporation.`],
      ['Sensible heat', `<b>${d.sensibleHeat.toFixed(0)} W/m²</b> — heat conducted from the surface into the air.`],
      ['Precipitation', `<b>${(d.precipitation * 86400).toFixed(2)} mm/day</b> — global mean rain rate over the last frame.`],
      ['Precipitable water', `<b>${d.columnWater.toFixed(1)} kg/m²</b> — all the vapour in a column, global mean; the same number in mm of rain.`],
      ['Cloud water', `<b>${(1000 * d.columnCloud).toFixed(0)} g/m²</b> — condensed water in a column, global mean.`],
      ['Sea ice', `<b>${(100 * d.iceFraction).toFixed(0)}%</b> of the area${d.iceThickness ? `, <b>${d.iceThickness.toFixed(2)} m</b> thick on average` : ''}.`],
      ['Planetary albedo', `<b>${d.planetaryAlbedo.toFixed(2)}</b> — the fraction of sunlight reflected back to space by clouds, ice and water.`],
      ['Resolution', `<b>N=${latest.N}</b> — icosahedral grid, cells about <b>${(7720 / latest.N).toFixed(0)} km</b> across, ${ready ? `${ready.cells.toLocaleString()} cells × ${ready.layers} layers` : ''}.`],
      ['Time step', ready ? `<b>${ready.dt} s</b> per step.` : ''],
      ['Engine', latest.engine === 'gpu' ? '<b>GPU</b> — every kernel runs on the graphics processor through WebGPU in single precision.' : `<b>${latest.workers > 1 ? `${latest.workers} worker threads` : 'one thread'}</b> — the CPU engine in double precision.`],
    ];
    if (d.oceanUpperDepth !== undefined) rows.push(['Ocean', `upper layer <b>${d.oceanUpperDepth.toFixed(0)} m</b> deep on average, currents to <b>${d.oceanSpeed.toFixed(2)} m/s</b>, thermocline <b>${d.oceanThermoclineT.toFixed(1)} K</b>.`]);
    document.getElementById('modelDetails').innerHTML = rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
  }

  function update(changes) {
    const before = activeLevel();
    if ('isolines' in changes) {
      const { isolines, ...rest } = changes;
      changes = isolines === 'off' ? { ...rest, isobars: 'off' } : { ...rest, isobars: 'on', [isolinesFor(before).setting]: Number(isolines) };
    }
    Object.assign(settings, changes);
    saveSettings(settings);
    if (activeLevel() !== before) worker.postMessage({ type: 'level', level: activeLevel() });
    if ('projection' in changes && viewer) viewer.setProjection(settings.projection);
    render();
  }

  let ready = null;
  worker.onmessage = (event) => {
    const message = event.data;
    if (message.type === 'status') {
      document.getElementById('date').textContent = message.text;
      const progress = document.getElementById('progress');
      if (message.fraction !== null && message.fraction !== undefined) {
        progress.classList.add('visible');
        document.getElementById('progressBar').style.width = `${Math.round(100 * message.fraction)}%`;
        document.getElementById('progressText').textContent = message.text;
      }
    }
    if (message.type === 'ready') {
      ready = message;
      setup(message.N);
      document.getElementById('date').textContent = `model ready: ${message.cells} cells × ${message.layers} layers, dt ${message.dt} s, ${message.workers > 1 ? `${message.workers} workers` : 'one thread'}`;
    }
    if (message.type === 'snapshotData') storeSnapshot(message);
    if (message.type === 'frame') {
      document.getElementById('progress').classList.remove('visible');
      latest = { ...message, N: ready.N, workers: ready.workers };
      clock.push({ wall: performance.now(), time: message.time });
      while (clock.length > 2 && clock[clock.length - 1].wall - clock[0].wall > 30000) clock.shift();
      render();
    }
  };
  worker.onerror = (error) => { document.getElementById('date').textContent = `worker error: ${error.message}`; };
  worker.postMessage({ type: 'start', N, from: from ? new URL(from, location.href).href : null, workers: crossOriginIsolated ? workers : 1, engine, paused, level: activeLevel() });

  for (const group of panel.querySelectorAll('.options[data-setting]')) {
    group.addEventListener('click', (event) => {
      const button = event.target.closest('button[data-value]');
      if (button) update({ [group.dataset.setting]: button.dataset.value });
    });
  }
  document.getElementById('palette').addEventListener('change', (event) => {
    const kind = OVERLAYS[settings.overlay].kind;
    update({ palettes: { ...settings.palettes, [kind]: event.target.value } });
  });
  const localList = document.getElementById('localList'), builtinList = document.getElementById('builtinList');
  const item = (name, meta, tag, buttons) => {
    const li = document.createElement('li');
    const label = document.createElement('span'); label.className = 'name'; label.textContent = name; li.append(label);
    if (tag) { const t = document.createElement('span'); t.className = 'tag'; t.textContent = tag; li.append(t); }
    const m = document.createElement('span'); m.className = 'meta'; m.textContent = meta; li.append(m);
    for (const [text, handler] of buttons) { const b = document.createElement('button'); b.className = 'flat'; b.textContent = text; b.addEventListener('click', handler); li.append(b); }
    return li;
  };
  const toSnapshot = (saved, url = null) => {
    const arrays = Object.fromEntries(['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].filter((k) => saved[k]).map((k) => [k, Float64Array.from(saved[k]).buffer]));
    const ocean = saved.ocean ? Object.fromEntries(Object.entries(saved.ocean).map(([k, v]) => [k, Float64Array.from(v).buffer])) : null;
    const bytes = Object.values(arrays).reduce((n, b) => n + b.byteLength, 0) + (ocean ? Object.values(ocean).reduce((n, b) => n + b.byteLength, 0) : 0);
    return { meta: { N: saved.N, K: saved.K, day: saved.day, time: saved.time, bytes, source: url }, data: { arrays, ocean } };
  };
  async function refreshSnapshots() {
    const list = await listSnapshots();
    localList.replaceChildren(...list.map((meta) => item(meta.name, `day ${meta.day.toFixed(0)} · N=${meta.N} · ${(meta.bytes / 1048576).toFixed(0)} MB`, meta.source ? 'built in' : '', [
      ['Restore', () => restoreSnapshot(meta.id)],
      ['Rename', async () => { const name = prompt('Snapshot name', meta.name); if (name && name !== meta.name) { await renameSnapshot(meta.id, name); refreshSnapshots(); } }],
      ['Clone', async () => { const name = prompt('Name for the copy', `${meta.name} (copy)`); if (name) { await cloneSnapshot(meta.id, name); refreshSnapshots(); } }],
      ['Delete', async () => { if (confirm(`Delete "${meta.name}"?`)) { await deleteSnapshot(meta.id); refreshSnapshots(); } }],
    ])));
    document.getElementById('localEmpty').style.display = list.length ? 'none' : '';
    for (const button of document.querySelectorAll('[data-snapshot="latest"]')) button.disabled = list.length === 0;
    const defaultUrl = from ? new URL(from, location.href).href : null;
    const builtin = (await builtinSnapshots()).sort((a, b) => (b.url === defaultUrl) - (a.url === defaultUrl));
    builtinList.replaceChildren(...builtin.map((entry) => {
      const local = list.find((meta) => meta.source === entry.url);
      const day = entry.file.match(/_state_day(\d+)/)?.[1];
      return item(entry.name, `${day ? `day ${Number(day)}` : ''}${entry.url === defaultUrl ? ' · the page default' : ''}`, local ? 'downloaded' : '', [
        [local ? 'Restore' : 'Download and restore', async () => { const id = local ? local.id : await download(entry); if (id) restoreSnapshot(id); }],
        ...(local ? [] : [['Download', async () => { await download(entry); }]]),
      ]);
    }));
    document.getElementById('builtinNote').style.display = builtin.length ? '' : 'none';
    return list;
  }
  async function download(entry) {
    try {
      document.getElementById('date').textContent = `downloading ${entry.name}…`;
      const saved = await (await fetch(entry.url)).json();
      const { meta, data } = toSnapshot(saved, entry.url);
      const id = await saveSnapshot({ name: entry.name, created: Date.now(), ...meta }, data);
      await refreshSnapshots();
      render();
      return id;
    } catch (error) { document.getElementById('date').textContent = `download failed: ${error.message}`; return null; }
  }
  async function storeSnapshot(message) {
    const { meta, data } = toSnapshot({ N: message.N, K: message.K, day: message.day, time: message.time, ...Object.fromEntries(Object.entries(message.arrays).map(([k, b]) => [k, new Float64Array(b)])), ocean: message.ocean ? Object.fromEntries(Object.entries(message.ocean).map(([k, b]) => [k, new Float64Array(b)])) : null });
    const name = `Day ${Math.floor(message.day)} · ${new Date().toLocaleString()}`;
    await saveSnapshot({ name, created: Date.now(), ...meta }, data);
    await refreshSnapshots();
    document.getElementById('date').textContent = `saved "${name}"`;
  }
  async function restoreSnapshot(id) {
    const { meta, data } = await getSnapshot(id);
    if (!meta) return;
    running = false;
    clock.length = 0;
    closeModals();
    render();
    const transfer = [...Object.values(data.arrays), ...(data.ocean ? Object.values(data.ocean) : [])];
    worker.postMessage({ type: 'restore', snapshot: { N: meta.N, K: meta.K, day: meta.day, time: meta.time, arrays: data.arrays, ocean: data.ocean } }, transfer);
  }
  for (const button of document.querySelectorAll('[data-snapshot]')) {
    button.addEventListener('click', async () => {
      const action = button.dataset.snapshot;
      if (action === 'save') { running = false; clock.length = 0; render(); worker.postMessage({ type: 'pause' }); worker.postMessage({ type: 'snapshot' }); }
      if (action === 'latest') { const list = await listSnapshots(); if (list.length) await restoreSnapshot(list[0].id); }
    });
  }
  const closeModals = () => { for (const modal of document.querySelectorAll('.modal')) modal.classList.add('hidden'); };
  for (const tab of document.querySelectorAll('.tabs button[data-tab]')) {
    tab.addEventListener('click', () => {
      for (const button of tab.parentElement.children) button.classList.toggle('selected', button === tab);
      for (const section of document.querySelectorAll('.tab[data-tab]')) section.classList.toggle('hidden', section.dataset.tab !== tab.dataset.tab);
    });
  }
  document.getElementById('modelButton').addEventListener('click', () => { renderModelDetails(); document.getElementById('modelModal').classList.remove('hidden'); });
  document.getElementById('snapshotsButton').addEventListener('click', () => { refreshSnapshots(); document.getElementById('snapshotModal').classList.remove('hidden'); });
  for (const button of document.querySelectorAll('[data-close]')) button.addEventListener('click', closeModals);
  for (const modal of document.querySelectorAll('.modal')) modal.addEventListener('click', (event) => { if (event.target === modal) closeModals(); });
  window.addEventListener('keydown', (event) => { if (event.key === 'Escape') closeModals(); });
  document.getElementById('viewDetails').innerHTML = VIEW_NOTES.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');

  document.querySelector('[data-control="play"]').addEventListener('click', () => {
    running = !running;
    clock.length = 0;
    worker.postMessage({ type: running ? 'resume' : 'pause' });
    render();
  });
  document.getElementById('menu').addEventListener('click', () => update({ panel: settings.panel === 'open' ? 'closed' : 'open' }));
  const bottom = document.getElementById('bottom');
  new ResizeObserver(() => { panel.style.bottom = `${bottom.offsetHeight + 20}px`; }).observe(bottom);
  render();
}
