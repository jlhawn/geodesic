/*
 * The page and the model worker talk as client and server. The page
 * subscribes to what it draws,
 *   { type: 'subscribe', subscription: { level, depth, fields, diagnostics } },
 * with level 'surface' (the lowest layer) or a pressure in hPa, depth
 * 'surface' (the mixed layer) or a depth in metres for the ocean
 * fields that vary with depth, fields a list of names from FIELDS, and
 * diagnostics true while it shows the global summary. Every frame the
 * worker posts carries exactly that:
 *   { type: 'frame', frame, time, day, level, depth, engine, pause, fields: { name: Float32Array }, diagnostics }
 * with diagnostics null unless subscribed and pause the worker's idle
 * time after each GPU step in ms. A frame built before the latest
 * subscription arrived may lack a field the page now wants. The page
 * also reports { type: 'pace', late } once a second while it runs.
 */
export const FIELDS = {
  temperature: 'air temperature at the level, K',
  height: 'geopotential height of the level, m',
  humidity: 'relative humidity at the level, fraction',
  speed: 'wind speed at the level, m/s',
  wind: 'wind at the level, m/s, three components per cell',
  dewPoint: 'dew point at the level, K',
  wetBulb: 'wet-bulb temperature at the level, K',
  misery: 'heat index above 26.7 °C, wind chill below 10 °C, air temperature between, at the level, K',
  vertical: 'vertical velocity at the level, m/s, positive upward, averaged with the neighbouring cells and over a two-hour memory',
  ps: 'surface pressure, Pa',
  mslp: 'sea-level pressure, Pa',
  water: 'precipitable water, kg/m²',
  cloud: 'column cloud water: the resolved condensate, the cumulus cloud\'s condensate times its cover and the stratocumulus deck\'s water times its cover, kg/m²',
  cloudLow: 'the resolved condensate in the layers below 800 hPa, kg/m²',
  cloudMid: 'the resolved condensate in the layers between 800 and 500 hPa, kg/m²',
  cloudHigh: 'the resolved condensate in the layers above 500 hPa, kg/m²',
  cloudCumulus: 'the cumulus cloud\'s condensate times its cover, kg/m²',
  cloudDeck: 'the stratocumulus deck\'s water times its cover, kg/m²',
  cloudTop: 'height above sea level where light from above first meets the column\'s cloud, averaged over where it does, 0 where the column is less than 1% opaque, m',
  cloudBase: 'height above sea level where light from below first meets the column\'s cloud, averaged over where it does, 0 where the column is less than 1% opaque, m',
  rain: 'recent rain, mm, with a three-hour exponential memory',
  ice: 'sea-ice thickness over the part of the cell the ice covers, m',
  concentration: 'sea-ice concentration, the fraction of the cell the ice covers, 0 to 1',
  albedo: 'surface albedo for diffuse light',
  shortwave: 'sunlight reaching the surface, W/m²',
  longwave: 'outgoing longwave radiation, W/m²',
  soil: 'soil water, kg/m²',
  snow: 'snow water, kg/m²',
  vegetation: 'vegetation cover, 0 (bare ground) to 1 (dense forest)',
  sst: 'sea temperature at the depth, K',
  sss: 'sea surface salinity, psu',
  layerDepth: 'mixed layer depth, m',
  thermocline: 'thermocline depth, m',
  ssh: 'sea surface height, m',
  current: 'current speed at the depth, m/s',
  currents: 'current at the depth, m/s, three components per cell',
  upwelling: 'vertical velocity at the depth, m/s, positive upward',
};

export const LEVEL_FIELDS = new Set(['temperature', 'height', 'humidity', 'speed', 'wind', 'dewPoint', 'wetBulb', 'misery', 'vertical']);
export const OCEAN_FIELDS = new Set(['sst', 'sss', 'layerDepth', 'thermocline', 'ssh', 'current', 'currents', 'upwelling']);
export const DEPTH_FIELDS = new Set(['sst', 'current', 'currents', 'upwelling']);

// A layer's group goes by the pressure at its middle, Pa.
export const CLOUD_LOW_PRESSURE = 80000, CLOUD_HIGH_PRESSURE = 50000;
export const CLOUD_TYPES = ['cloudLow', 'cloudMid', 'cloudHigh', 'cloudCumulus', 'cloudDeck'];
// The cloud overlays' legend ranges, g/m², the top near the 95th percentile of the cells holding each type.
export const CLOUD_RANGES = { cloud: 100, cloudLow: 200, cloudMid: 500, cloudHigh: 400, cloudCumulus: 40, cloudDeck: 150 };

// The cloud water path, kg/m², over which a column's opacity rises by a factor e, and the opacity below which it counts as clear.
export const CLOUD_OPACITY_PATH = 0.040, CLOUD_SEEN = 0.01;

/*
 * Whether a frame still comes from a finite model. A model gone non-finite
 * fills every prognostic field within a few steps, so each field is read
 * on a stride of FINITE_STRIDE and counts as gone when most of those
 * samples are not numbers; the ocean's fields are left out, as they carry
 * NaN over land by design, and the mean surface temperature, when the
 * diagnostics carry it, is read whole.
 */
export const FINITE_STRIDE = 61;
export function frameIsFinite(fields, diagnostics = null) {
  if (diagnostics && diagnostics.meanSurfaceT !== undefined && !Number.isFinite(diagnostics.meanSurfaceT)) return false;
  for (const [name, values] of Object.entries(fields)) {
    if (OCEAN_FIELDS.has(name) || !values || typeof values.length !== 'number') continue;
    let bad = 0, seen = 0;
    for (let i = 0; i < values.length; i += FINITE_STRIDE) { seen++; if (!Number.isFinite(values[i])) bad++; }
    if (seen && 2 * bad > seen) return false;
  }
  return true;
}

/*
 * Where light first meets a column's cloud, for layers k = 0 (top) to K-1
 * (bottom) with mid heights z[k] above sea level over ground at 'ground':
 * each layer spans halfway to its neighbours (the top layer as far above
 * its middle as halfway down, the bottom one down to the ground) and holds
 * the cloud water path[k], kg/m², and the deck's path joins the lowest
 * layer reaching deckTop (the top layer if none does). With the opacity
 * 1 - exp(-P / CLOUD_OPACITY_PATH) of a path P, light from above stops in a
 * layer with the layer's own opacity times the clearness of all above it;
 * out[0] is that layer's upper bound averaged with those weights, which
 * sum to the column's opacity, and out[1] the same for light from below
 * and lower bounds. Both are 0 where the column is less than CLOUD_SEEN
 * opaque.
 */
export function visibleCloudHeights(K, z, ground, path, deck, deckTop, out) {
  let total = deck;
  for (let k = 0; k < K; k++) total += path[k];
  const opacity = 1 - Math.exp(-total / CLOUD_OPACITY_PATH);
  out[0] = 0; out[1] = 0;
  if (opacity < CLOUD_SEEN) return out;
  let below = 0, fromBelow = 1, left = deck, top = 0, base = 0;
  for (let k = K - 1; k >= 0; k--) {
    const lower = k === K - 1 ? ground : 0.5 * (z[k] + z[k + 1]);
    const upper = k === 0 ? z[0] + 0.5 * (z[0] - z[1]) : 0.5 * (z[k] + z[k - 1]);
    let p = path[k];
    if (left > 0 && (upper >= deckTop || k === 0)) { p += left; left = 0; }
    if (p > 0) {
      const through = Math.exp(-p / CLOUD_OPACITY_PATH);
      base += fromBelow * (1 - through) * lower;
      fromBelow *= through;
      top += Math.exp(-Math.max(0, total - below - p) / CLOUD_OPACITY_PATH) * (1 - through) * upper;
      below += p;
    }
  }
  out[0] = top / opacity; out[1] = base / opacity;
  return out;
}

export const RAIN_MEMORY = 3 * 3600;
export const VERTICAL_MEMORY = 2 * 3600;
