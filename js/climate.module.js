import { Grid } from "./grid.module.js";
import { initUnifiedViewer } from "./unifiedViewer.module.js";
import { createWindParticles } from "./windParticles.module.js";
import { fetchState, stateName, decodeState } from './stateFile.module.js';
import { seasonPhrase } from "./levels.module.js";
import { sunDirection, DAY, YEAR } from "./physics/radiation.module.js";
import { createDisplayClock } from "./displayClock.module.js";
import { listSnapshots, saveSnapshot, getSnapshot, renameSnapshot, deleteSnapshot, cloneSnapshot } from "./snapshots.module.js";

const WIND_MAX = { surface: 25, 1000: 30, 850: 40, 700: 40, 500: 50, 250: 70, 70: 100, 10: 150 };
const TEMP_RANGE = { surface: [-35, 35], 1000: [-35, 35], 850: [-45, 25], 700: [-55, 15], 500: [-65, 5], 250: [-85, -25], 70: [-95, -35], 10: [-75, 5] };
const CELSIUS = -273.15;
const REFERENCE_SPEED = { surface: 15, 1000: 20, 850: 25, 700: 25, 500: 30, 250: 40, 70: 50, 10: 60 };

const OVERLAYS = {
  wind: { label: 'Wind speed', short: 'WIND', unit: 'm/s', kind: 'sequential', field: 'speed', scale: 1, range: (level) => [0, WIND_MAX[level]] },
  temp: { label: 'Temperature', short: 'TEMP', unit: '°C', kind: 'sequential', field: 'temperature', scale: 1, offset: CELSIUS, range: (level) => TEMP_RANGE[level] },
  rh: { label: 'Relative humidity', short: 'RH', unit: '%', kind: 'sequential', field: 'humidity', scale: 100, range: () => [0, 100] },
  mi: { label: 'Misery index', short: 'MI', unit: '°C', kind: 'sequential', derive: (frame) => deriveField(frame, miseryIndex), point: miseryIndex, scale: 1, range: () => [-40, 45] },
  wbt: { label: 'Wet-bulb temperature', short: 'WBT', unit: '°C', kind: 'sequential', derive: (frame) => deriveField(frame, wetBulb), point: wetBulb, scale: 1, range: () => [-40, 35] },
  dp: { label: 'Dew point', short: 'DP', unit: '°C', kind: 'sequential', derive: (frame) => deriveField(frame, dewPoint), point: dewPoint, scale: 1, range: () => [-40, 30] },
  rain: { label: 'Recent rain', short: 'RAIN', unit: 'mm', kind: 'sequential', field: 'rain', scale: 1, range: () => [0, 20] },
  tpw: { label: 'Total precipitable water', short: 'TPW', unit: 'kg/m²', kind: 'sequential', field: 'water', scale: 1, range: () => [0, 60] },
  tcw: { label: 'Total cloud water', short: 'TCW', unit: 'g/m²', kind: 'sequential', field: 'cloud', scale: 1000, range: () => [0, 500] },
  cloudcover: { label: 'Cloud cover', short: 'CC', unit: 'g/m²', kind: 'clouds', field: 'cloud', scale: 1000, range: () => [0, 100] },
  albedo: { label: 'Surface albedo', short: 'ALB', unit: '', kind: 'sequential', field: 'albedo', scale: 1, range: () => [0, 0.8] },
  swdown: { label: 'Surface sunlight', short: 'SSI', unit: 'W/m²', kind: 'sequential', field: 'shortwave', scale: 1, range: () => [0, 1200] },
  olr: { label: 'Outgoing longwave radiation', short: 'OLR', unit: 'W/m²', kind: 'sequential', field: 'longwave', scale: 1, range: () => [100, 320] },
  ice: { label: 'Sea ice thickness', short: 'ICE', unit: 'm', kind: 'sequential', field: 'ice', scale: 1, decimals: 2, range: () => [0, 3] },
  mslp: { label: 'Sea-level pressure', short: 'MSLP', unit: 'hPa', kind: 'diverging', field: 'mslp', scale: 0.01, range: () => [960, 1060] },
  ps: { label: 'Surface pressure', short: 'PS', unit: 'hPa', kind: 'sequential', field: 'ps', scale: 0.01, range: () => [500, 1050] },
  soil: { label: 'Soil water', short: 'SOIL', unit: 'kg/m²', kind: 'sequential', field: 'soil', scale: 1, range: () => [0, 150] },
  snow: { label: 'Snow', short: 'SNOW', unit: 'kg/m²', kind: 'sequential', field: 'snow', scale: 1, range: () => [0, 100] },
  elevation: { label: 'Elevation', short: 'ELEV', unit: 'm', kind: 'diverging', field: 'elevation', scale: 1, range: () => [-4000, 4000] },
  sst: { label: 'Sea surface temperature', short: 'SST', unit: '°C', kind: 'sequential', field: 'sst', scale: 1, offset: CELSIUS, range: () => [-2, 32] },
  current: { label: 'Current speed', short: 'CUR', unit: 'm/s', kind: 'sequential', field: 'current', scale: 1, range: () => [0, 1] },
  layer: { label: 'Mixed layer depth', short: 'MLD', unit: 'm', kind: 'sequential', field: 'layerDepth', scale: 1, range: () => [10, 300] },
  thermocline: { label: 'Thermocline depth', short: 'THD', unit: 'm', kind: 'sequential', field: 'thermocline', scale: 1, range: () => [0, 1200] },
  sss: { label: 'Sea surface salinity', short: 'SSS', unit: 'psu', kind: 'sequential', field: 'sss', scale: 1, range: () => [32, 38] },
  ssh: { label: 'Sea surface height', short: 'SSH', unit: 'm', kind: 'sequential', field: 'ssh', scale: 1, decimals: 2, range: () => [-1.5, 1.5] },
  none: { label: 'None', short: 'None' },
};
const MODE_OVERLAYS = {
  atmosphere: [['none'], ['wind', 'temp', 'rh'], ['mi', 'wbt', 'dp'], ['rain', 'tpw', 'tcw', 'cloudcover'], ['albedo', 'swdown', 'olr'], ['mslp', 'ps', 'elevation'], ['ice', 'snow', 'soil']],
  ocean: [['none'], ['sst', 'current', 'layer', 'thermocline', 'sss', 'ssh']],
};
const MODE_DEFAULT_OVERLAY = { atmosphere: 'wind', ocean: 'sst' };

/*
 * Fields derived on the page from the frame's temperature, relative
 * humidity and wind at the chosen height. Dew point by the Magnus
 * formula; wet-bulb by Stull's fit; the misery index is the NWS heat
 * index above 26.7 °C, the wind chill below 10 °C, and the air
 * temperature between.
 */
function dewPoint(t, rh) {
  const r = Math.max(1e-3, Math.min(1, rh)), a = 17.625, b = 243.04, g = Math.log(r) + a * t / (b + t);
  return b * g / (a - g);
}
function wetBulb(t, rh) {
  const p = 100 * Math.max(0.05, Math.min(0.99, rh));
  return t * Math.atan(0.151977 * Math.sqrt(p + 8.313659)) + Math.atan(t + p) - Math.atan(p - 1.676331) + 0.00391838 * Math.pow(p, 1.5) * Math.atan(0.023101 * p) - 4.686035;
}
function heatIndex(t, rh) {
  const T = t * 9 / 5 + 32, R = 100 * Math.max(0, Math.min(1, rh));
  let hi = -42.379 + 2.04901523 * T + 10.14333127 * R - 0.22475541 * T * R - 6.83783e-3 * T * T - 5.481717e-2 * R * R + 1.22874e-3 * T * T * R + 8.5282e-4 * T * R * R - 1.99e-6 * T * T * R * R;
  if (R < 13 && T <= 112) hi -= ((13 - R) / 4) * Math.sqrt((17 - Math.abs(T - 95)) / 17);
  else if (R > 85 && T <= 87) hi += ((R - 85) / 10) * ((87 - T) / 5);
  return (hi - 32) * 5 / 9;
}
function windChill(t, v) {
  const k = 3.6 * v;
  if (k < 4.8) return t;
  const p = Math.pow(k, 0.16);
  return 13.12 + 0.6215 * t - 11.37 * p + 0.3965 * t * p;
}
function miseryIndex(t, rh, v) {
  if (t >= 26.7) return Math.max(t, heatIndex(t, rh));
  if (t <= 10) return Math.min(t, windChill(t, v));
  return t;
}
function deriveField(frame, fn) {
  const T = frame.temperature, H = frame.humidity, V = frame.speed, out = new Float32Array(T.length);
  for (let i = 0; i < T.length; i++) out[i] = fn(T[i] + CELSIUS, H[i], V[i]);
  return out;
}
const OVERLAY_NAMES = Object.fromEntries(Object.entries(OVERLAYS).map(([key, overlay]) => [key, overlay.label]));

/*
 * The cloud view: open water is ocean blue, ice whitens with thickness,
 * and cloud is white composited on top with an opacity that rises with
 * the column's cloud water, 1 − exp(−TCW / CLOUD_OPACITY_SCALE), so
 * clear sky is transparent and 40 g/m² is two-thirds opaque.
 */
const OCEAN_COLOR = [0.05, 0.22, 0.45], ICE_COLOR = [0.85, 0.90, 0.95], CLOUD_COLOR = [1, 1, 1], CLOUD_OPACITY_SCALE = 40;
const DRY_LAND = [0.45, 0.36, 0.22], WET_LAND = [0.16, 0.30, 0.12], SNOW_COLOR = [0.9, 0.92, 0.95];
const COVER_BASE = [0.22, 0.22, 0.22];
const cloudOpacity = (grams) => 1 - Math.exp(-Math.max(0, grams) / CLOUD_OPACITY_SCALE);
const COVER_STOPS = Array.from({ length: 11 }, (_, k) => { const a = cloudOpacity(10 * k); return COVER_BASE.map((c) => c + a * (1 - c)); });

/*
 * Palettes as sRGB stops. The sequential ones are perceptually uniform
 * colormaps that stay distinguishable under the common colour-vision
 * deficiencies; the diverging ones pair hues that do the same, around a
 * neutral gray so the white wind trails read on top.
 */
const PALETTES = {
  viridis: [[0.267, 0.005, 0.329], [0.283, 0.141, 0.458], [0.254, 0.265, 0.530], [0.207, 0.372, 0.553], [0.164, 0.471, 0.558], [0.128, 0.567, 0.551], [0.135, 0.659, 0.518], [0.267, 0.749, 0.441], [0.478, 0.821, 0.318], [0.741, 0.873, 0.150], [0.993, 0.906, 0.144]],
  cividis: [[0.000, 0.135, 0.304], [0.127, 0.196, 0.416], [0.256, 0.263, 0.437], [0.367, 0.336, 0.446], [0.472, 0.409, 0.458], [0.578, 0.484, 0.463], [0.690, 0.562, 0.450], [0.807, 0.645, 0.411], [0.926, 0.734, 0.339], [0.994, 0.831, 0.243]],
  inferno: [[0.001, 0.000, 0.014], [0.087, 0.036, 0.209], [0.258, 0.039, 0.406], [0.416, 0.090, 0.433], [0.578, 0.148, 0.404], [0.735, 0.216, 0.330], [0.865, 0.317, 0.226], [0.954, 0.462, 0.109], [0.988, 0.645, 0.040], [0.965, 0.844, 0.146], [0.988, 0.998, 0.645]],
  magma: [[0.001, 0.000, 0.014], [0.098, 0.062, 0.259], [0.269, 0.060, 0.478], [0.446, 0.122, 0.507], [0.617, 0.183, 0.499], [0.792, 0.253, 0.446], [0.933, 0.372, 0.375], [0.987, 0.541, 0.383], [0.996, 0.719, 0.518], [0.987, 0.898, 0.729]],
  dusk: [[0.16, 0.10, 0.45], [0.15, 0.35, 0.65], [0.15, 0.60, 0.45], [0.55, 0.70, 0.20], [0.85, 0.30, 0.15]],
  'blue-gray-red': [[0.13, 0.30, 0.75], [0.35, 0.42, 0.62], [0.50, 0.50, 0.50], [0.68, 0.38, 0.32], [0.80, 0.20, 0.15]],
  'purple-gray-orange': [[0.33, 0.15, 0.53], [0.50, 0.45, 0.67], [0.50, 0.50, 0.50], [0.88, 0.51, 0.08], [0.70, 0.35, 0.02]],
  'teal-gray-brown': [[0.00, 0.40, 0.37], [0.35, 0.64, 0.60], [0.50, 0.50, 0.50], [0.75, 0.55, 0.30], [0.55, 0.32, 0.04]],
};

const DEFAULTS = { view: 'atmosphere', overlay: 'wind', level: 'surface', animate: 'particles', isobars: 'off', isobarStep: 5, heightStep: 60, graticule: '15', projection: 'sphere', palette: 'viridis', panel: 'open' };

/*
 * The contour row draws isobars of surface pressure at the surface and
 * height contours of the pressure surface at any other level.
 */
const ISOLINES = {
  surface: { label: 'Isobars', unit: 'hPa', setting: 'isobarStep', steps: [1, 2.5, 5, 10], field: (frame) => Float32Array.from(frame.mslp ?? frame.ps, (p) => p / 100) },
  level: { label: 'Height lines', unit: 'm', setting: 'heightStep', steps: [20, 30, 60, 120, 240], field: (frame) => frame.height },
};
const isolinesFor = (level) => (level === 'surface' ? ISOLINES.surface : ISOLINES.level);

function loadSettings() {
  try {
    const stored = JSON.parse(localStorage.getItem('climate.settings') || '{}');
    const { palettes: oldPalettes, ...rest } = stored;
    const settings = { ...DEFAULTS, ...rest };
    if (!PALETTES[settings.palette]) settings.palette = (oldPalettes && PALETTES[oldPalettes.sequential]) ? oldPalettes.sequential : DEFAULTS.palette;
    if (settings.overlay === 'clouds') { settings.overlay = 'none'; settings.view = 'space'; }
    reconcile(settings);
    return settings;
  } catch { return { ...DEFAULTS }; }
}
/*
 * Settings named in the page's query string override the stored ones;
 * a value is accepted when the panel offers it, and naming something
 * only the data view shows implies that view.
 */
function applyOverrides(settings, overrides) {
  if (!('view' in overrides)) {
    if ('overlay' in overrides) settings.view = modeOf(OVERLAY_ALIASES[overrides.overlay] ?? overrides.overlay);
    else if (['level', 'isobars', 'isobarStep', 'heightStep'].some((key) => key in overrides)) settings.view = 'atmosphere';
    else if ('animate' in overrides && settings.view === 'space') settings.view = modeOf(settings.overlay);
  }
  for (const key of Object.keys(DEFAULTS)) {
    if (!(key in overrides)) continue;
    let value = overrides[key];
    if (typeof DEFAULTS[key] === 'number') { if (Number(value) > 0) settings[key] = Number(value); continue; }
    if (key === 'overlay') value = OVERLAY_ALIASES[value] ?? value;
    const known = key === 'palette' ? PALETTES[value] : key === 'overlay' ? OVERLAYS[value] : key === 'view' ? ['space', 'atmosphere', 'ocean', 'data'].includes(value) : key === 'panel' ? ['open', 'closed'].includes(value) : key === 'isobars' ? ['on', 'off'].includes(value) : document.querySelector(`[data-setting="${key}"] [data-value="${CSS.escape(value)}"]`);
    if (known) settings[key] = value;
  }
  reconcile(settings);
}
/*
 * Every mode shows only its own overlays: an overlay names its mode,
 * the old data view is whichever mode its overlay belongs to, and a
 * mode switch that leaves an overlay behind takes the mode's default.
 */
const OVERLAY_ALIASES = { precip: 'rain', p3h: 'rain' };
const OCEAN_OVERLAYS = new Set(MODE_OVERLAYS.ocean.flat().filter((key) => key !== 'none'));
function modeOf(overlay) { return OCEAN_OVERLAYS.has(overlay) ? 'ocean' : 'atmosphere'; }
function reconcile(settings) {
  settings.overlay = OVERLAY_ALIASES[settings.overlay] ?? settings.overlay;
  if (!OVERLAYS[settings.overlay]) settings.overlay = DEFAULTS.overlay;
  if (!['space', 'atmosphere', 'ocean'].includes(settings.view)) settings.view = modeOf(settings.overlay);
  if (settings.view !== 'space' && !MODE_OVERLAYS[settings.view].flat().includes(settings.overlay)) settings.overlay = MODE_DEFAULT_OVERLAY[settings.view];
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
 * built-in snapshots the page can download into its own store; a
 * server without directory listings still offers the page default.
 */
async function builtinSnapshots(fallback = null) {
  const entry = (file) => ({ file, url: new URL(file, location.href).href, name: stateName(file.replace(/.*\//, '')) });
  try {
    const html = await (await fetch('runs/')).text();
    const files = [...new Set([...html.matchAll(/href="([^"]+_state_day\d+(?:\.json(?:\.gz)?|\.parts\.json))"/g)].map((m) => decodeURIComponent(m[1])))].sort();
    const newest = new Map(files.map((file) => [stateName(file), file]));
    if (newest.size) return [...newest.values()].map((file) => entry(`runs/${file}`));
  } catch {}
  return fallback ? [entry(fallback)] : [];
}

const HEIGHT_OVERLAYS = new Set(['wind', 'temp', 'rh', 'mi', 'wbt', 'dp', 'none']);
const CURRENT_REFERENCE = 0.2;

const VIEW_NOTES = [
  ['Mode', 'Atmosphere and Ocean paint the chosen overlay on an evenly lit globe, each with its own overlays: the wind or the current is what the animation follows, and only Atmosphere offers isobars and height lines. Satellite renders the planet as it would look from space: ocean, ice and cloud lit by the sun in its true direction for the model date and time, a dark ambient on the night side, and the stars turning behind it once a sidereal day.'],
  ['Wind animation', 'Particles trace the wind at the chosen height as fading trails, brighter where it blows faster; Vectors draw one arrow per cell; None hides the motion.'],
  ['Height', 'The pressure level shown by the wind, temperature and humidity views and followed by the animation: Sfc is the lowest layer, about 60 m up; the others are hPa. Column views hide it and use the surface wind.'],
  ['Wind speed', 'Speed at the chosen height.'],
  ['Temperature', 'Air temperature at the chosen height.'],
  ['Relative humidity', 'At the chosen height.'],
  ['Surface pressure', 'The pressure at the ground itself, about 1000 hPa at the coast and 550 hPa on the Tibetan plateau; the weather signal is the small variation on top of the elevation.'],
  ['Cloud cover', 'Cloud as white over grey with the opacity the Satellite view uses, from the column\'s cloud water.'],
  ['Sea surface temperature', 'The temperature of the ocean\'s wind-driven upper layer, the freezing point under ice; grey over land.'],
  ['Current speed', 'The upper layer\'s current, up to a metre a second in the boundary currents. With any ocean view selected, Particles and Vectors trace the current instead of the wind.'],
  ['Mixed layer depth', 'The thickness of the surface mixed layer: deep where winter cooling and wind stirring reach down, shallow under summer warming and along upwelling coasts.'],
  ['Thermocline depth', 'The depth of the base of the second interior layer, deep in the subtropical gyres where the wind piles warm water up and shallow at the equator and toward the poles.'],
  ['Sea surface height', 'The free surface: high over the subtropical gyres and low around the poles, with the currents flowing along its contours.'],
  ['Thermocline temperature', 'The layer below the upper one, which entrains into it where the upper layer thins.'],
  ['Sea-level pressure', 'Surface pressure reduced to sea level through a standard-lapse-rate column below the terrain; the isobars use it too. Where a pressure level lies below the ground, wind, temperature and humidity show the lowest layer of that column, and only the height is extrapolated hydrostatically so its contours stay a pressure field.'],
  ['Soil water', 'The land bucket: up to 150 kg/m² of soil water; evaporation slows as it dries and rain beyond its capacity runs off.'],
  ['Snow', 'Snow on land in water equivalent; it falls when the lowest air is below freezing and melts into the bucket.'],
  ['Elevation', 'The mean elevation of each cell from ETOPO 2022; negative under the sea.'],
  ['Coastlines', 'The mesh edges between land and ocean cells, drawn in the Atmosphere and Ocean modes.'],
  ['Misery index', 'How the air feels: the heat index where it is warmer than 26.7 °C, the wind chill where it is colder than 10 °C, the air temperature between.'],
  ['Wet-bulb temperature', 'The coolest a wet surface can get by evaporation at the chosen height; above about 35 °C the body can no longer shed heat.'],
  ['Dew point', 'The temperature the air would have to cool to for its vapour to condense; close to the air temperature means humid air.'],
  ['Recent rain', 'Rain from convection and from cloud that rained out, as an exponentially weighted accumulation with a three-hour memory: steady rain settles at its three-hour total and a shower fades over the hours after it.'],
  ['Precipitable water', 'All the vapour in the column, as the depth of rain it would make.'],
  ['Cloud water', 'All the condensed water in the column.'],
  ['Sea ice', 'Sea-ice thickness.'],
  ['Albedo', 'The surface albedo for diffuse light: 0.06 over water, rising to 0.5 over half a metre of ice.'],
  ['Surface sunlight', 'Shortwave reaching the surface, direct and diffuse, before the surface reflects its share.'],
  ['Outgoing longwave', 'Infrared leaving the top of the atmosphere: low over cold cloud tops and the poles, high over clear warm regions.'],
  ['Isobars / Height lines', 'Contours of surface pressure at the surface, of geopotential height on a pressure level, at the chosen interval.'],
  ['Graticule', 'Parallels and meridians at the chosen spacing; the meridians stop at the outermost parallel.'],
  ['Projection', 'The orthographic globe, or the Equal Earth map; both can be dragged to any orientation.'],
  ['Snapshots', 'Save the paused state in this browser, restore it later, download one of the runs saved on the server, or import and export snapshot files to share them.'],
];

export default function runClimate({ N = null, from = null, workers = 1, engine = 'cpu', paused = false, land = true, topography = null, terrain = true, settings: overrides = {}, view = null } = {}) {
  const settings = loadSettings();
  applyOverrides(settings, overrides);
  const panel = document.getElementById('panel');
  const activeLevel = () => (settings.view === 'atmosphere' && HEIGHT_OVERLAYS.has(settings.overlay) ? settings.level : 'surface');
  const shownLevel = () => latest?.level ?? activeLevel();
  let latest = null, grid = null, viewer = null, particles = null, arrows = null, isobars = null, graticule = null, coast = null, highlight = null, rgb = null, running = !paused, animatedSource = null, seaCells = null;
  let cells = null, centres = null, selected = -1, hoverTip = null;
  let geographyFields = {}, hasLand = false;
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
  /*
   * The address bar always holds a link to the current view: every
   * setting, the globe's orientation and whether the model is paused,
   * rewritten in place a moment after the last change. The page reads
   * the query string only when it starts.
   */
  let urlVersion = -1, urlRunning = null, urlTimer = 0;
  function reflectUrl() {
    const url = new URL(location.href);
    for (const key of Object.keys(DEFAULTS)) url.searchParams.set(key, String(settings[key]));
    if (viewer) {
      const orientation = viewer.view();
      const tenth = (x) => String(Math.round(10 * x) / 10 || 0);
      url.searchParams.set('lat', tenth(orientation.lat));
      url.searchParams.set('lon', tenth(orientation.lon));
      url.searchParams.set('zoom', String(Math.round(orientation.zoom)));
      if (Math.abs(orientation.roll) >= 0.05) url.searchParams.set('roll', tenth(orientation.roll)); else url.searchParams.delete('roll');
    }
    if (running) url.searchParams.delete('paused'); else url.searchParams.set('paused', '');
    history.replaceState(null, '', url);
  }
  function scheduleUrl() { clearTimeout(urlTimer); urlTimer = setTimeout(reflectUrl, 300); }

  function tick(now) {
    requestAnimationFrame(tick);
    if (viewer && (viewer.viewVersion() !== urlVersion || running !== urlRunning)) { urlVersion = viewer.viewVersion(); urlRunning = running; scheduleUrl(); }
    if (!latest) return;
    const time = display.advance(now, latest.time, { ...recentFrames(), running });
    if (running) document.getElementById('date').textContent = formatDate(time);
    if (settings.view === 'space' && viewer) viewer.setSpace({ enabled: true, sun: sunDirection(time), sidereal: 2 * Math.PI * time * (1 / DAY + 1 / YEAR) });
  }
  requestAnimationFrame(tick);

  const worker = new Worker(new URL('./model.worker.js', import.meta.url), { type: 'module' });

  function setup(size) {
    grid = new Grid(size);
    rgb = new Uint8Array(3 * grid.size);
    viewer = initUnifiedViewer(document.getElementById('globe'), grid, { backgroundColor: 0x151515, dynamicColors: true, controls: false, getColor: () => ({ r: 0.25, g: 0.25, b: 0.25 }) });
    arrows = viewer.addArrowLayer({ opacity: 0.5 });
    isobars = viewer.addContourLayer({ opacity: 0.25 });
    graticule = viewer.addGraticuleLayer({ opacity: 0.25 });
    coast = viewer.addSegmentLayer({ opacity: 0.6 });
    highlight = viewer.addSegmentLayer({ color: 0xffe8a0, opacity: 1 });
    cells = [...grid];
    centres = new Float32Array(3 * cells.length);
    for (const cell of cells) { const c = cell.centerVertex, r = Math.hypot(c.x, c.y, c.z); centres[3 * cell.index] = c.x / r; centres[3 * cell.index + 1] = c.y / r; centres[3 * cell.index + 2] = c.z / r; }
    selected = -1;
    particles = createWindParticles(document.getElementById('globe'), viewer, grid);
    viewer.setProjection(settings.projection);
    if (view) viewer.setView(view);
  }
  function teardown() {
    if (particles) particles.dispose();
    if (viewer) viewer.dispose();
    particles = viewer = arrows = isobars = graticule = coast = highlight = null;
    cells = centres = null; selected = -1;
  }

  function paintSatellite() {
    const cloud = latest.cloud, ice = latest.ice, land = latest.land, soil = latest.soil, snow = latest.snow;
    for (let i = 0; i < grid.size; i++) {
      const opacity = cloudOpacity(cloud[i] * 1000);
      const onLand = land && land[i];
      const frozen = onLand ? Math.min(1, snow[i] / 20) : Math.min(1, ice[i] / 0.5);
      const wet = onLand ? Math.min(1, soil[i] / 150) : 0;
      for (let j = 0; j < 3; j++) {
        const ground = onLand ? DRY_LAND[j] + wet * (WET_LAND[j] - DRY_LAND[j]) : OCEAN_COLOR[j];
        const white = onLand ? SNOW_COLOR[j] : ICE_COLOR[j];
        const base = ground + frozen * (white - ground);
        rgb[3 * i + j] = LINEAR[Math.round(255 * (base + opacity * (CLOUD_COLOR[j] - base)))];
      }
    }
    viewer.updateColors(rgb);
    document.querySelector('.scaleRow').classList.add('hidden');
    document.getElementById('data').textContent = 'Satellite view';
  }

  function paintOverlay() {
    if (settings.view === 'space') { paintSatellite(); return; }
    const overlay = OVERLAYS[settings.overlay];
    const scaleRow = document.querySelector('.scaleRow');
    if (!overlay.field && !overlay.derive) {
      rgb.fill(LINEAR[40]);
      viewer.updateColors(rgb);
      scaleRow.classList.add('hidden');
      document.getElementById('data').textContent = settings.view === 'ocean' ? 'Surface current · no overlay' : `Wind @ ${levelLabel(shownLevel())} · no overlay`;
      return;
    }
    const [min, max] = overlay.range(shownLevel());
    const values = overlay.derive ? overlay.derive(latest) : latest[overlay.field];
    if (overlay.kind === 'clouds') {
      for (let i = 0; i < grid.size; i++) {
        const opacity = cloudOpacity(values[i] * overlay.scale);
        for (let j = 0; j < 3; j++) rgb[3 * i + j] = LINEAR[Math.round(255 * (COVER_BASE[j] + opacity * (1 - COVER_BASE[j])))];
      }
      viewer.updateColors(rgb);
      scaleRow.classList.remove('hidden');
      renderScale(COVER_STOPS, min, max, overlay.unit);
      document.getElementById('data').textContent = `${OVERLAY_NAMES[settings.overlay]} · wind @ ${levelLabel(shownLevel())}`;
      return;
    }
    const stops = PALETTES[settings.palette] ?? PALETTES.viridis;
    for (let i = 0; i < grid.size; i++) {
      if (Number.isNaN(values[i])) { rgb[3 * i] = rgb[3 * i + 1] = rgb[3 * i + 2] = LINEAR[70]; continue; }
      color((values[i] * overlay.scale + (overlay.offset || 0) - min) / (max - min), stops, rgb, 3 * i);
      rgb[3 * i] = LINEAR[rgb[3 * i]]; rgb[3 * i + 1] = LINEAR[rgb[3 * i + 1]]; rgb[3 * i + 2] = LINEAR[rgb[3 * i + 2]];
    }
    viewer.updateColors(rgb);
    scaleRow.classList.remove('hidden');
    renderScale(stops, min, max, overlay.unit);
    const columnField = ['ps', 'mslp', 'rain', 'water', 'cloud', 'ice', 'albedo', 'shortwave', 'longwave', 'soil', 'snow', 'elevation'].includes(overlay.field);
    document.getElementById('data').textContent = settings.view === 'ocean' ? `${OVERLAY_NAMES[settings.overlay]} · surface current` : columnField ? `${OVERLAY_NAMES[settings.overlay]} · wind @ ${levelLabel(shownLevel())}` : `${OVERLAY_NAMES[settings.overlay]} @ ${levelLabel(shownLevel())}`;
  }

  function paintWind() {
    const ocean = settings.view === 'ocean' && latest.currentVector && latest.currentVector.length > 0;
    const reference = ocean ? CURRENT_REFERENCE : REFERENCE_SPEED[shownLevel()];
    const field = ocean ? latest.currentVector : latest.vector;
    const animate = settings.view === 'space' ? 'none' : settings.animate;
    arrows.setVisible(animate === 'arrows');
    particles.setVisible(animate === 'particles');
    if (animate === 'arrows') arrows.update(field, { referenceSpeed: reference });
    const source = ocean ? 'current' : `wind ${shownLevel()}`;
    if (animate === 'particles' && source !== animatedSource) particles.reset();
    animatedSource = source;
    if (animate === 'particles') particles.setField(field, reference, ocean ? seaCells : null);
    const note = animate === 'particles' ? `full pace at ${reference} m/s` : animate === 'arrows' ? `full arrow at ${reference} m/s` : '';
    if (note) document.getElementById('data').textContent += ` · ${note}`;
  }

  function paintGraticule() {
    if (!graticule) return;
    graticule.setVisible(settings.graticule !== 'off');
    if (settings.graticule !== 'off') graticule.update(Number(settings.graticule));
  }

  function paintIsobars() {
    const on = settings.isobars === 'on' && settings.view !== 'space';
    isobars.setVisible(on);
    if (!on) return;
    const isolines = isolinesFor(shownLevel());
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
    const mode = settings.view === 'ocean' ? 'ocean' : 'atmosphere', overlayBox = document.getElementById('overlayOptions');
    if (overlayBox.dataset.mode !== mode) {
      overlayBox.dataset.mode = mode;
      overlayBox.replaceChildren(...MODE_OVERLAYS[mode].map((row) => {
        const segment = document.createElement('div');
        segment.className = 'segmented';
        for (const key of row) { const button = document.createElement('button'); button.dataset.value = key; button.textContent = OVERLAYS[key].short; button.dataset.tip = OVERLAYS[key].label; segment.append(button); }
        return segment;
      }));
    }
    for (const group of panel.querySelectorAll('.options[data-setting]')) {
      const current = group.dataset.setting === 'isolines' ? isolineChoice() : String(settings[group.dataset.setting]);
      for (const button of group.querySelectorAll('button[data-value]')) button.classList.toggle('selected', button.dataset.value === current);
    }
    const space = settings.view === 'space';
    const heights = settings.view === 'atmosphere' && HEIGHT_OVERLAYS.has(settings.overlay);
    document.getElementById('heightLabel').classList.toggle('hidden', !heights);
    document.getElementById('heightOptions').classList.toggle('hidden', !heights);
    for (const id of ['overlayLabel', 'overlayOptions', 'animateLabel', 'animateOptions']) document.getElementById(id).classList.toggle('hidden', space);
    for (const id of ['isolineLabel', 'isolineOptions']) document.getElementById(id).classList.toggle('hidden', settings.view !== 'atmosphere');
    document.getElementById('animateLabel').textContent = settings.view === 'ocean' ? 'Current animation' : 'Wind animation';
    document.getElementById('isolineLabel').textContent = isolines.label;
    document.getElementById('isolineUnit').textContent = isolines.unit;
    const paletteSelect = document.getElementById('palette');
    const overlay = OVERLAYS[settings.overlay];
    const palettes = !!overlay.kind && overlay.kind !== 'clouds';
    paletteSelect.style.display = palettes ? '' : 'none';
    if (palettes) fillSelect(paletteSelect, Object.keys(PALETTES), settings.palette);
    document.querySelector('[data-control="play"]').textContent = running ? '❚❚' : '▶';
    panel.classList.toggle('hidden', settings.panel !== 'open');
    if (viewer) viewer.setSpace({ enabled: space });
    if (!latest) return;
    paintOverlay();
    paintWind();
    paintIsobars();
    paintGraticule();
    refreshTip();
    if (coast) coast.setVisible(hasLand && settings.view !== 'space');
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
      ['Surface temperature', `<b>${(d.meanSurfaceT + CELSIUS).toFixed(1)} °C</b> — area-weighted global mean of the skin temperature.`],
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
    if (d.oceanUpperDepth !== undefined) rows.push(['Ocean', `mixed layer <b>${d.oceanUpperDepth.toFixed(0)} m</b> deep on average, currents to <b>${d.oceanSpeed.toFixed(2)} m/s</b>${d.oceanTransport !== undefined ? `, the strongest transport <b>${d.oceanTransport.toFixed(0)} Sv</b>` : ''}${d.oceanThermoclineDepth !== undefined ? `, thermocline <b>${d.oceanThermoclineDepth.toFixed(0)} m</b>` : ''}.`]);
    if (d.landFraction !== undefined) rows.push(['Land', `<b>${(100 * d.landFraction).toFixed(0)}%</b> of the area${ready && ready.terrain ? ' with terrain' : ', flat'}; surface <b>${(d.landMeanT + CELSIUS).toFixed(1)} °C</b>, soil water <b>${d.soilWater.toFixed(0)} kg/m²</b>, snow on <b>${(100 * d.snowFraction).toFixed(0)}%</b> of it.`]);
    document.getElementById('modelDetails').innerHTML = rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
  }

  function update(changes) {
    const before = activeLevel();
    if ('isolines' in changes) {
      const { isolines, ...rest } = changes;
      changes = isolines === 'off' ? { ...rest, isobars: 'off' } : { ...rest, isobars: 'on', [isolinesFor(before).setting]: Number(isolines) };
    }
    Object.assign(settings, changes);
    reconcile(settings);
    saveSettings(settings);
    scheduleUrl();
    if (activeLevel() !== before) worker.postMessage({ type: 'level', level: activeLevel() });
    if ('projection' in changes && viewer) viewer.setProjection(settings.projection);
    render();
  }

  /*
   * The rain field is an exponentially weighted accumulation of each
   * frame's rain with a three-hour memory, S ← S·e^(−Δt/3h) + rate·Δt,
   * which for steady rain settles at the three-hour total and needs no
   * history; a frame from an earlier time than the last one starts over.
   */
  const RAIN_MEMORY = 3 * 3600;
  const rain = { time: null, total: null };
  function accumulateRain(message) {
    const rate = message.precipitation;
    if (rain.total === null || rain.total.length !== rate.length || message.time < rain.time) { rain.total = new Float32Array(rate.length); rain.time = message.time; }
    const span = message.time - rain.time, keep = Math.exp(-span / RAIN_MEMORY), days = span / 86400;
    for (let i = 0; i < rate.length; i++) rain.total[i] = rain.total[i] * keep + rate[i] * days;
    rain.time = message.time;
    return rain.total;
  }

  /*
   * The tip box in the bottom-right corner shows the full name of the
   * panel button under the pointer, or, when a cell has been clicked,
   * that cell's value of the current overlay and its position; the
   * highlighted cell's boundary follows the globe.
   */
  const DECIMALS = { '°C': 1, '%': 0, 'm/s': 2, mm: 1, 'kg/m²': 1, 'g/m²': 0, 'W/m²': 0, hPa: 1, m: 0, psu: 2, '': 2 };
  function valueAt(overlay, i) {
    if (overlay.point) return overlay.point(latest.temperature[i] + CELSIUS, latest.humidity[i], latest.speed[i]);
    if (!overlay.field) return null;
    const field = latest[overlay.field];
    return field && field.length > i ? field[i] * overlay.scale + (overlay.offset || 0) : NaN;
  }
  function positionText(i) {
    const lat = Math.asin(Math.max(-1, Math.min(1, centres[3 * i + 2]))) * 180 / Math.PI, lon = Math.atan2(centres[3 * i + 1], centres[3 * i]) * 180 / Math.PI;
    const place = latest && latest.land ? (latest.land[i] ? 'land' : 'sea') : null;
    return `${Math.abs(lat).toFixed(1)}°${lat < 0 ? 'S' : 'N'} ${Math.abs(lon).toFixed(1)}°${lon < 0 ? 'W' : 'E'}${place ? ` · ${place}` : ''}`;
  }
  function readoutHtml() {
    if (selected < 0 || !latest) return null;
    const overlay = settings.view === 'space' ? null : OVERLAYS[settings.overlay];
    const value = overlay ? valueAt(overlay, selected) : null;
    const text = value === null ? 'no overlay' : Number.isNaN(value) ? `${overlay.label}: —` : `${overlay.label}: ${value.toFixed(overlay.decimals ?? DECIMALS[overlay.unit] ?? 1)} ${overlay.unit}`.trim();
    return `<div class="value">${text}</div><div>${positionText(selected)}</div>`;
  }
  function refreshTip() {
    const tip = document.getElementById('tip');
    const html = hoverTip !== null ? `<div class="value">${hoverTip}</div>` : readoutHtml();
    tip.classList.toggle('hidden', html === null);
    if (html !== null && tip.innerHTML !== html) tip.innerHTML = html;
  }
  function select(i) {
    selected = i;
    if (!highlight) return;
    if (i < 0) { highlight.setVisible(false); refreshTip(); return; }
    const verts = cells[i].vertices || [], positions = new Float32Array(6 * verts.length), which = new Int32Array(verts.length);
    for (let k = 0; k < verts.length; k++) {
      const a = verts[k], b = verts[(k + 1) % verts.length], ra = Math.hypot(a.x, a.y, a.z), rb = Math.hypot(b.x, b.y, b.z);
      positions.set([a.x / ra, a.y / ra, a.z / ra, b.x / rb, b.y / rb, b.z / rb], 6 * k);
      which[k] = i;
    }
    highlight.set(positions, which);
    highlight.setVisible(true);
    refreshTip();
  }
  function nearestCell(point) {
    let best = -1, bestDot = -Infinity;
    for (let i = 0; i < cells.length; i++) {
      const d = point[0] * centres[3 * i] + point[1] * centres[3 * i + 1] + point[2] * centres[3 * i + 2];
      if (d > bestDot) { bestDot = d; best = i; }
    }
    return best;
  }
  {
    const globe = document.getElementById('globe');
    let press = null;
    globe.addEventListener('mousedown', (event) => { if (event.button === 0) press = { x: event.clientX, y: event.clientY, at: performance.now() }; });
    globe.addEventListener('mouseup', (event) => {
      if (!press || event.button !== 0) return;
      const moved = Math.hypot(event.clientX - press.x, event.clientY - press.y), held = performance.now() - press.at;
      press = null;
      if (moved > 4 || held > 500 || !viewer || !cells) return;
      const rect = globe.getBoundingClientRect(), point = viewer.unprojectPoint(event.clientX - rect.left, event.clientY - rect.top, [0, 0, 0]);
      select(point ? nearestCell(point) : -1);
    });
    window.addEventListener('keydown', (event) => { if (event.key === 'Escape' && selected >= 0) select(-1); });
    panel.addEventListener('mouseover', (event) => { const button = event.target.closest('button[data-tip]'); if (button) { hoverTip = button.dataset.tip; refreshTip(); } });
    panel.addEventListener('mouseout', (event) => { const button = event.target.closest('button[data-tip]'); if (button && hoverTip !== null) { hoverTip = null; refreshTip(); } });
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
      rain.total = null;
      const rebuild = !viewer || !ready || ready.cells !== message.cells;
      ready = message;
      if (rebuild) { teardown(); setup(message.N); latest = null; }
      animatedSource = null;
      hasLand = !!message.land;
      geographyFields = hasLand ? { land: message.land, landFraction: message.landFraction, elevation: message.elevation } : {};
      seaCells = hasLand ? Uint8Array.from(message.land, (l) => 1 - l) : null;
      if (hasLand) coast.set(message.coast, message.coastCells);
      document.getElementById('date').textContent = `model ready: ${message.cells} cells × ${message.layers} layers, dt ${message.dt} s, ${message.workers > 1 ? `${message.workers} workers` : 'one thread'}`;
    }
    if (message.type === 'snapshotData') storeSnapshot(message);
    if (message.type === 'frame') {
      document.getElementById('progress').classList.remove('visible');
      latest = { ...message, ...geographyFields, N: ready.N, workers: ready.workers, rain: accumulateRain(message) };
      clock.push({ wall: performance.now(), time: message.time });
      while (clock.length > 2 && clock[clock.length - 1].wall - clock[0].wall > 30000) clock.shift();
      render();
    }
  };
  worker.onerror = (error) => { document.getElementById('date').textContent = `worker error: ${error.message}`; };
  worker.postMessage({ type: 'start', N, from: from ? new URL(from, location.href).href : null, workers: crossOriginIsolated ? workers : 1, engine, paused, level: activeLevel(), land, terrain, topography: topography ? new URL(topography, location.href).href : null });

  for (const group of panel.querySelectorAll('.options[data-setting]')) {
    group.addEventListener('click', (event) => {
      const button = event.target.closest('button[data-value]');
      if (button) update({ [group.dataset.setting]: button.dataset.value });
    });
  }
  document.getElementById('palette').addEventListener('change', (event) => update({ palette: event.target.value }));
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
    const land = saved.land ? Object.fromEntries(Object.entries(saved.land).map(([k, v]) => [k, Float64Array.from(v).buffer])) : null;
    const size = (o) => (o ? Object.values(o).reduce((n, b) => n + b.byteLength, 0) : 0);
    const bytes = size(arrays) + size(ocean) + size(land);
    return { meta: { N: saved.N, K: saved.K, day: saved.day, time: saved.time, terrain: !!saved.terrain, bytes, source: url }, data: { arrays, ocean, land } };
  };
  async function refreshSnapshots() {
    const list = await listSnapshots();
    localList.replaceChildren(...list.map((meta) => item(meta.name, `day ${meta.day.toFixed(0)} · N=${meta.N} · ${(meta.bytes / 1048576).toFixed(0)} MB`, meta.source ? 'built in' : '', [
      ['Restore', () => restoreSnapshot(meta.id)],
      ['Rename', async () => { const name = prompt('Snapshot name', meta.name); if (name && name !== meta.name) { await renameSnapshot(meta.id, name); refreshSnapshots(); } }],
      ['Clone', async () => { const name = prompt('Name for the copy', `${meta.name} (copy)`); if (name) { await cloneSnapshot(meta.id, name); refreshSnapshots(); } }],
      ['Export', () => exportSnapshot(meta)],
      ['Delete', async () => { if (confirm(`Delete "${meta.name}"?`)) { await deleteSnapshot(meta.id); refreshSnapshots(); } }],
    ])));
    document.getElementById('localEmpty').style.display = list.length ? 'none' : '';
    for (const button of document.querySelectorAll('[data-snapshot="latest"]')) button.disabled = list.length === 0;
    const defaultUrl = from ? new URL(from, location.href).href : null;
    const builtin = (await builtinSnapshots(from)).sort((a, b) => (b.url === defaultUrl) - (a.url === defaultUrl));
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
      const saved = await fetchState(entry.url);
      const { meta, data } = toSnapshot(saved, entry.url);
      const id = await saveSnapshot({ name: entry.name, created: Date.now(), ...meta }, data);
      await refreshSnapshots();
      render();
      return id;
    } catch (error) { document.getElementById('date').textContent = `download failed: ${error.message}`; return null; }
  }
  async function storeSnapshot(message) {
    const { meta, data } = toSnapshot({ N: message.N, K: message.K, day: message.day, time: message.time, terrain: message.terrain, ...Object.fromEntries(Object.entries(message.arrays).map(([k, b]) => [k, new Float64Array(b)])), ocean: message.ocean ? Object.fromEntries(Object.entries(message.ocean).map(([k, b]) => [k, new Float64Array(b)])) : null, land: message.land ? Object.fromEntries(Object.entries(message.land).map(([k, b]) => [k, new Float64Array(b)])) : null });
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
    scheduleUrl();
    render();
    const transfer = [...Object.values(data.arrays), ...(data.ocean ? Object.values(data.ocean) : []), ...(data.land ? Object.values(data.land) : [])];
    worker.postMessage({ type: 'restore', snapshot: { N: meta.N, K: meta.K, day: meta.day, time: meta.time, terrain: !!meta.terrain, arrays: data.arrays, ocean: data.ocean, land: data.land ?? null } }, transfer);
  }
  const slug = (name) => name.replace(/[^A-Za-z0-9]+/g, '_').replace(/^_+|_+$/g, '') || 'snapshot';
  async function exportSnapshot(meta) {
    const { data } = await getSnapshot(meta.id);
    if (!data) return;
    document.getElementById('date').textContent = `exporting "${meta.name}"…`;
    const saved = { N: meta.N, K: meta.K, day: meta.day, time: meta.time, terrain: !!meta.terrain };
    for (const [key, buffer] of Object.entries(data.arrays)) saved[key] = Array.from(new Float64Array(buffer));
    if (data.ocean) saved.ocean = Object.fromEntries(Object.entries(data.ocean).map(([k, b]) => [k, Array.from(new Float64Array(b))]));
    if (data.land) saved.land = Object.fromEntries(Object.entries(data.land).map(([k, b]) => [k, Array.from(new Float64Array(b))]));
    const gz = await new Response(new Blob([JSON.stringify(saved)]).stream().pipeThrough(new CompressionStream('gzip'))).blob();
    const link = document.createElement('a');
    link.href = URL.createObjectURL(gz);
    link.download = `${slug(meta.name)}_state_day${Math.round(meta.day)}.json.gz`;
    link.click();
    setTimeout(() => URL.revokeObjectURL(link.href), 60000);
    document.getElementById('date').textContent = `exported ${link.download} (${(gz.size / 1048576).toFixed(0)} MB)`;
  }
  async function importSnapshots(files) {
    for (const file of files) {
      try {
        document.getElementById('date').textContent = `importing ${file.name}…`;
        const saved = await decodeState(new Uint8Array(await file.arrayBuffer()));
        if (!saved || !saved.N || !saved.pi) throw new Error('not a saved state');
        const { meta, data } = toSnapshot(saved);
        await saveSnapshot({ name: stateName(file.name), created: Date.now(), ...meta }, data);
        document.getElementById('date').textContent = `imported ${file.name}`;
      } catch (error) { document.getElementById('date').textContent = `import failed: ${error.message}`; }
    }
    await refreshSnapshots();
  }
  document.getElementById('snapshotFile').addEventListener('change', async (event) => { const files = [...event.target.files]; event.target.value = ''; await importSnapshots(files); });
  for (const button of document.querySelectorAll('[data-snapshot]')) {
    button.addEventListener('click', async () => {
      const action = button.dataset.snapshot;
      if (action === 'import') document.getElementById('snapshotFile').click();
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
    scheduleUrl();
    render();
  });
  document.getElementById('menu').addEventListener('click', () => update({ panel: settings.panel === 'open' ? 'closed' : 'open' }));
  const bottom = document.getElementById('bottom');
  new ResizeObserver(() => { panel.style.bottom = `${bottom.offsetHeight + 20}px`; }).observe(bottom);
  render();
}
