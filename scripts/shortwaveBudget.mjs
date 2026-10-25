// The model's own CPU shortwave on one saved state, term by term, over one
// day lit at TIMES (8) instants:
//   [RADIATION='{...}'] [EXTRA='{...}'] [PASS=clear] node scripts/shortwaveBudget.mjs <state.bin> [out.json]
// Every column is taken as the saved state holds it (cloud, deck and cumulus
// included; PASS=clear takes them away) and its sunlight split into what the
// column absorbs (asr), the surface (surf) and the atmosphere (atm) take of
// it, what leaves the top (refl), the cloud's own absorption (cloud), the
// gases' (gas), ozone's, aerosol's and the upward pass's, and the shielding
// term of radiation.vaporShielding; closure is the insolation less absorbed
// and reflected. Area-weighted means in W/m2 over the globe, its sea and land,
// 0-20N, 40-60N and each 10-degree band, printed and written to out.json.
// RADIATION (JSON) replaces the radiation options the page's states are spun
// up with and EXTRA (JSON) adds to them, MOIST (JSON) sets the moist physics'.
import { readFileSync, writeFileSync } from 'node:fs';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { savedGeography, checkGeography, readTopography, geographyOptions } from '../js/geography.module.js';
import { savedDeckField, DECK_FIELDS, savedMoistField, MOIST_FIELDS } from '../js/physics/regrid.module.js';
import { withCadence } from '../js/cadence.module.js';

const FILE = process.argv[2], OUT = process.argv[3], TIMES = Number(process.env.TIMES ?? 8), PASS = process.env.PASS ?? 'all';
const RADIATION_OPTIONS = JSON.parse(process.env.RADIATION ?? '{"cloudColumns":true,"cloudSpread":0.75,"varianceCover":"water","buoyancyLength":0.53}');
const MOIST = JSON.parse(process.env.MOIST ?? '{"stratocumulusInversion":7,"cumulusEntrainmentHeight":1,"cumulusAscent":"kinetic"}');
const EXTRA = JSON.parse(process.env.EXTRA ?? '{}');
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const N = saved.N, dt = 1350 * 16 / N, deg = 180 / Math.PI;
const RADIATION = withCadence({ radiation: { clearSkyPass: true, ...RADIATION_OPTIONS, ...EXTRA } }, dt, N).radiation;
const dragEvery = withCadence({ dragEvery: undefined }, dt, N).dragEvery;
const GEOGRAPHY = checkGeography({ ...savedGeography(saved) });
const topography = readTopography(GEOGRAPHY.topography), sheets = geographyOptions(GEOGRAPHY);
const model = createModel(stateGrid(saved, topography).mesh, { topography, geography: sheets.geography, ocean: false, radiation: RADIATION, moist: MOIST, boundaryLayer: {}, surface: {}, land: { ...sheets.land }, dynamics: {}, dragEvery, levels: savedLevels(saved) });
const { mesh, core, state, radiation, boundaryLayer: bl, seaIce, land } = model;
const { K, C, g, geopotential } = core.diagnostics;
const [pi, theta, , surfaceT, q, qc, ice] = state;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
for (const field of Object.keys(DECK_FIELDS)) radiation[field].set(savedDeckField(saved, field, model));
for (const field of Object.keys(MOIST_FIELDS)) model.moist[field].set(savedMoistField(saved, field, model));
if (saved.boundaryDepth) bl.depth.set(saved.boundaryDepth);
if (saved.mixingTop) bl.mixingTop.set(saved.mixingTop);
if (saved.boundaryRegime) bl.regime.set(saved.boundaryRegime);
if (saved.boundaryBuoyancy) bl.buoyancyFlux.set(saved.boundaryBuoyancy);
if (saved.windSpeed) model.surface.windSpeed.set(saved.windSpeed);
for (const field of ['cumulusCover', 'cumulusWater']) if (saved[field] && saved[field].length === model.moist[field].length) model.moist[field].set(saved[field]);
if (saved.subcloudVirtual && saved.subcloudVirtual.length === model.moist.subcloudVirtual.length) model.moist.subcloudVirtual.set(saved.subcloudVirtual);
land.load(saved.land);
for (let i = 0; i < C; i++) core.diagnoseColumn(i, pi, theta, q, qc);
const bottom = (K - 1) * C, landMask = model.geography.land, wind = model.surface.windSpeed;
const day = Math.round(saved.time / 86400);
const KEYS = ['beam', 'asr', 'surf', 'atm', 'refl', 'cloud', 'gas', 'ozone', 'aerosol', 'upward', 'shielded', 'miss'];
const cells = [];
for (let i = 0; i < C; i++) {
  const isLand = !!landMask[i];
  cells.push({ i, la: mesh.latCell[i] * deg, w: mesh.areaCell[i], isLand, h: ice[i], area: isLand ? 0 : seaIce.cover(i, ice[i]), sum: Object.fromEntries(KEYS.map((k) => [k, 0])) });
}
const albedos = (c, mu) => {
  if (c.isLand) { const a = land.albedo(c.i); return [a, a]; }
  const skin = surfaceT[c.i], snowy = seaIce.snowAlbedo[c.i];
  return [seaIce.albedo(c.h, mu, seaIce.snow[c.i], c.area, skin, snowy), seaIce.albedo(c.h, null, seaIce.snow[c.i], c.area, skin, snowy)];
};
const gasOf = (b) => b.vaporSolar + b.oxygenSolar + b.carbonDioxideSolar + b.methaneSolar;
if (PASS === 'clear') radiation.useCumulus(null, null);
for (let h = 0; h < TIMES; h++) {
  radiation.setTime(saved.time + h * 86400 / TIMES);
  for (const c of cells) {
    const { i } = c, mu = radiation.cosZenith(i);
    if (!(mu > 0)) continue;
    const [direct, diffuse] = albedos(c, mu), top = bl.turbulence === 'moist' ? bl.mixingTop[i] - geopotential[bottom + i] / g : 0;
    if (PASS === 'all') radiation.column(i, pi[i], theta, surfaceT[i], wind[i], undefined, undefined, q[bottom + i], q, qc, direct, diffuse, 1, 0, c.isLand ? 0 : 1 - c.area, bl.depth[i] - geopotential[bottom + i] / g, 0, top);
    else radiation.column(i, pi[i], theta, surfaceT[i], wind[i], undefined, undefined, q[bottom + i], q, null, direct, diffuse, 1, 0, 0, 0, 0, 0);
    const b = radiation.budget, s = c.sum;
    s.beam += b.insolation; s.asr += b.absorbedSolar; s.surf += b.surfaceSolar; s.atm += b.atmosphereSolar; s.refl += b.reflectedSolar; s.cloud += b.cloudSolar;
    s.gas += gasOf(b); s.ozone += b.ozoneSolar; s.aerosol += b.aerosolSolar; s.upward += b.upwardGasSolar; s.shielded += b.shieldedGasSolar;
    s.miss += b.insolation - b.absorbedSolar - b.reflectedSolar;
  }
}
const sets = {};
const add = (name, c) => {
  const s = (sets[name] ??= { area: 0, sum: Object.fromEntries(KEYS.map((k) => [k, 0])), cells: 0 });
  s.area += c.w; s.cells++;
  for (const k of KEYS) s.sum[k] += c.w * c.sum[k] / TIMES;
};
for (const c of cells) {
  add('global', c); add(c.isLand ? 'globalLand' : 'globalSea', c);
  add(`band${Math.min(17, Math.floor((c.la + 90) / 10))}`, c);
  if (c.la >= 0 && c.la < 20) add('z0020', c);
  if (c.la >= 40 && c.la < 60) { add('z4060', c); add(c.isLand ? 'z4060Land' : 'z4060Sea', c); }
}
const result = { file: FILE, day, pass: PASS, extra: EXTRA, times: TIMES, sets: {} };
for (const [name, s] of Object.entries(sets)) result.sets[name] = { cells: s.cells, ...Object.fromEntries(KEYS.map((k) => [k, s.sum[k] / s.area])) };
if (OUT) writeFileSync(OUT, JSON.stringify(result));
const f = (x) => x.toFixed(2).padStart(7);
console.log(`${FILE.split('/').slice(-1)[0]}: day ${day}, ${PASS} sky, ${JSON.stringify(EXTRA)}; W/m2`);
for (const name of ['global', 'globalSea', 'globalLand', 'z0020', 'z4060']) {
  const r = result.sets[name];
  console.log(`${name.padEnd(10)} beam ${f(r.beam)} ASR ${f(r.asr)} surface ${f(r.surf)} atmosphere ${f(r.atm)} reflected ${f(r.refl)} | cloud ${f(r.cloud)} gas ${f(r.gas)} ozone ${f(r.ozone)} aerosol ${f(r.aerosol)} upward ${f(r.upward)} shield ${f(r.shielded)} | closure ${r.miss.toExponential(1)}`);
}
