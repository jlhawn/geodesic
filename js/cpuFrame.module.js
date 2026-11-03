import { cellVector } from './dynamics/operators.module.js';
import { levelFields, dewPoint, wetBulb, miseryIndex, verticalVelocity, smoothCells } from './levels.module.js';
import { LEVEL_FIELDS, OCEAN_FIELDS, CLOUD_TYPES, VERTICAL_MEMORY, SURFACE_UPWELLING_DEPTH } from './frames.module.js';
import { LATENT_HEAT } from './physics/moist.module.js';

const FREEZING = 273.15;

/*
 * The fields of frames.module.js for the CPU engines, from the model's
 * own arrays, only those named: the cell-center winds of a layer are
 * reconstructed only when that layer bounds the level somewhere. `step`
 * is the seconds of the model's last step, which the precipitation rate
 * is over. A frame's `memories` are the three-hour rain, `rain`, and
 * the vertical velocity's, `vertical` ({ values, level, time }, which
 * this folds the level's into); a caller that passes none gets none of
 * the rain, vertical and upwelling fields, and changes nothing of the
 * model's.
 */
export function cpuFrameFields(model, { level = 'surface', depth = 'surface', fields = [] } = {}, { step = 0, memories = null } = {}) {
  const { mesh, core, state } = model;
  const C = mesh.nCells, E = mesh.nEdges, [pi, theta, u, , q, , iceField] = state, time = model.time;
  const want = new Set(fields), out = {}, { rain = null, vertical = null } = memories ?? {};
  const onLand = model.geography ? model.geography.land : null;
  const sea = (source) => Float32Array.from({ length: C }, (_, i) => (onLand && onLand[i] ? NaN : source[i]));
  if (fields.some((name) => LEVEL_FIELDS.has(name))) {
    const layerWinds = new Array(core.K).fill(null);
    const layerWind = (k) => layerWinds[k] ??= cellVector(mesh, u.subarray(k * E, (k + 1) * E), new Float64Array(3 * C));
    const f = levelFields(core, pi, theta, layerWind, level, q);
    for (const [name, values] of Object.entries({ temperature: f.temperature, height: f.height, humidity: f.humidity, speed: f.speed, wind: f.vector })) if (want.has(name)) out[name] = values;
    const comfort = (fn) => Float32Array.from(f.temperature, (t, i) => fn(t - FREEZING, f.humidity[i], f.speed[i]) + FREEZING);
    if (want.has('dewPoint')) out.dewPoint = comfort(dewPoint);
    if (want.has('wetBulb')) out.wetBulb = comfort(wetBulb);
    if (want.has('misery')) out.misery = comfort(miseryIndex);
    if (vertical && want.has('vertical')) {
      const smoothed = smoothCells(mesh, verticalVelocity(mesh, core, pi, u, level, theta));
      const keep = vertical.values && vertical.level === level ? Math.exp(-Math.max(0, time - vertical.time) / VERTICAL_MEMORY) : 0;
      if (!vertical.values || vertical.values.length !== C) vertical.values = new Float32Array(C);
      for (let i = 0; i < C; i++) vertical.values[i] = keep * vertical.values[i] + (1 - keep) * smoothed[i];
      vertical.level = level; vertical.time = time;
      out.vertical = Float32Array.from(vertical.values);
    } else if (vertical) vertical.level = null;
  }
  if (want.has('ps')) out.ps = Float32Array.from(pi);
  if (want.has('mslp')) {
    out.mslp = Float32Array.from(pi);
    if (model.surfaceGeopotential) {
      const phis = model.surfaceGeopotential, { R, g, exnerLayer } = core.diagnostics, bottom = (core.K - 1) * C;
      for (let i = 0; i < C; i++) out.mslp[i] = pi[i] * Math.exp(phis[i] / (R * (theta[bottom + i] * exnerLayer[bottom + i] + 0.00325 * phis[i] / g)));
    }
  }
  if (want.has('water')) out.water = Float32Array.from({ length: C }, (_, i) => model.moist.columnWater(pi, q, i));
  if (want.has('cloud')) out.cloud = Float32Array.from({ length: C }, (_, i) => model.cloudWater(i));
  if (CLOUD_TYPES.some((name) => want.has(name))) {
    const parts = Array.from({ length: C }, (_, i) => model.cloudParts(i));
    for (const name of CLOUD_TYPES) if (want.has(name)) out[name] = Float32Array.from(parts, (part) => part[name]);
  }
  if (want.has('cloudTop') || want.has('cloudBase')) {
    const heights = model.cloudHeights();
    if (want.has('cloudTop')) out.cloudTop = heights.top;
    if (want.has('cloudBase')) out.cloudBase = heights.base;
  }
  if (rain && want.has('rain')) out.rain = Float32Array.from(rain);
  if (want.has('ice')) out.ice = Float32Array.from(iceField);
  if (want.has('concentration')) out.concentration = Float32Array.from(model.seaIce.concentration);
  if (want.has('albedo')) out.albedo = Float32Array.from(iceField, (h, i) => (onLand && onLand[i] ? model.land.albedo(i) : model.seaIce.albedo(h, null, model.seaIce.snow[i], model.seaIce.cover(i, h), model.state[3][i], model.seaIce.snowAlbedo[i])));
  if (want.has('shortwave')) out.shortwave = Float32Array.from(model.radiation.surfaceShortwave);
  if (want.has('longwave')) out.longwave = Float32Array.from(model.radiation.outgoing);
  if (model.land && want.has('soil')) out.soil = Float32Array.from(model.land.soil, (w, i) => w + model.land.surface[i] + (model.land.deepWater ? model.land.deepWater[i] : 0));
  if (model.land && want.has('snow')) out.snow = Float32Array.from(model.land.snow);
  if (model.land && want.has('vegetation')) out.vegetation = Float32Array.from(model.land.vegetation);
  if (want.has('surfaceT')) out.surfaceT = Float32Array.from(state[3]);
  if (model.moist && model.moist.rain && step > 0 && want.has('precipitation')) out.precipitation = Float32Array.from(model.moist.rain, (r) => r * 86400 / step);
  if (want.has('latentHeat')) out.latentHeat = Float32Array.from(model.radiation.evaporation, (e) => LATENT_HEAT * e);
  if (want.has('sensibleHeat')) out.sensibleHeat = Float32Array.from(model.radiation.sensibleHeat);
  if (model.surface && model.surface.stress && (want.has('zonalStress') || want.has('meridionalStress'))) {
    const stress = cellVector(mesh, model.surface.stress(state, new Float64Array(E))), { latCell, lonCell } = mesh;
    if (want.has('zonalStress')) out.zonalStress = Float32Array.from({ length: C }, (_, i) => -Math.sin(lonCell[i]) * stress[3 * i] + Math.cos(lonCell[i]) * stress[3 * i + 1]);
    if (want.has('meridionalStress')) out.meridionalStress = Float32Array.from({ length: C }, (_, i) => { const s = Math.sin(latCell[i]); return -s * Math.cos(lonCell[i]) * stress[3 * i] - s * Math.sin(lonCell[i]) * stress[3 * i + 1] + Math.cos(latCell[i]) * stress[3 * i + 2]; });
  }
  const ocean = model.oceanFields && fields.some((name) => OCEAN_FIELDS.has(name)) ? model.oceanFields(depth === 'surface' ? 0 : Number(depth), depth === 'surface' ? SURFACE_UPWELLING_DEPTH : Number(depth)) : null;
  if (ocean) {
    for (const [name, values] of Object.entries({ sst: ocean.temperature, sss: ocean.S1, layerDepth: ocean.h1, thermocline: ocean.thermoclineDepth, ssh: ocean.eta, current: ocean.speed, upwelling: ocean.upwelling })) if (want.has(name) && (name !== 'upwelling' || memories)) out[name] = sea(values);
    if (want.has('currents')) {
      const vector = ocean.current;
      if (onLand) for (let i = 0; i < C; i++) if (onLand[i]) vector.fill(0, 3 * i, 3 * i + 3);
      out.currents = Float32Array.from(vector);
    }
  }
  return out;
}
