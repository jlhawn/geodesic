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
  cloud: 'column cloud water, kg/m²',
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

export const RAIN_MEMORY = 3 * 3600;
export const VERTICAL_MEMORY = 2 * 3600;
