// The atmosphere alone on the GPU from a saved state over a sea surface
// held fixed: the ocean off and the sea a slab of CAP J/m²/K (1e12, so the
// surface does not move), stepped DAYS days as the spin-up steps, a state
// saved every SAVE_DAYS days from day SAVE_FROM on (stateFile.module.js
// format, no ocean; the land and the atmosphere's carried fields as
// scripts/spinup.mjs saves them). SST=own keeps the state's sea surface;
// SST=woa sets every open-sea cell to WOA23's monthly SST
// (data/woa_sst_monthly_1deg.bin) on the state's date. MOIST,
// BOUNDARY_LAYER, SURFACE, RADIATION and GRAVITY_WAVES (JSON, or false)
// pass options as scripts/spinup.mjs does; LEVELS (a grid of SIGMA_GRIDS)
// runs the arm on that grid with the state remapped onto it. The AMIP
// question for one state: what the atmosphere does over the observed sea.
//   SST=woa DAYS=12 SAVE_DAYS=3 node scripts/fixedSst.mjs STATE.bin OUT TAG
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createGpuModel } from '../js/gpu/model.gpu.js';
import { withCadence } from '../js/cadence.module.js';
import { decodeState, encodeState, savedLevels } from '../js/stateFile.module.js';
import { savedDeckField, DECK_FIELDS, savedMoistField, MOIST_FIELDS, savedRadiationField, RADIATION_FIELDS, remapLevels } from '../js/physics/regrid.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { stateGrid, gridHeader } from '../js/stateGrid.module.js';
import { decodeMonthlySst, monthWeights, climatologyAt } from '../js/ocean/climatology.module.js';

const [FILE, OUT, TAG] = process.argv.slice(2);
const DAYS = Number(process.env.DAYS ?? 12), SAVE_DAYS = Number(process.env.SAVE_DAYS ?? 3), SAVE_FROM = Number(process.env.SAVE_FROM ?? 0), SST = process.env.SST ?? 'own', CAP = Number(process.env.CAP ?? 1e12), BATCH = 8;
const env = (name) => JSON.parse(process.env[name] ?? '{}');
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const mesh = stateGrid(saved, topography).mesh, dt = 1350 * 16 / saved.N, perDay = Math.round(86400 / dt);
const levels = process.env.LEVELS ? sigmaInterfaces(process.env.LEVELS) : savedLevels(saved);
const model = await createGpuModel(mesh, { topography, levels, ...withCadence({ ocean: false, radiation: { clearSkyPass: true, ...env('RADIATION') } }, dt, saved.N), ice: { slabHeatCapacity: CAP }, moist: env('MOIST'), boundaryLayer: env('BOUNDARY_LAYER'), surface: env('SURFACE'), gravityWaves: env('GRAVITY_WAVES') });
const { state, core } = model, C = mesh.nCells;
const carried = { ...saved, ...(process.env.LEVELS ? remapLevels(savedLevels(saved), levels, saved, mesh) : {}) };
['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => state[a].set(carried[name]));
model.seaIce.load(state[6], saved.concentration ?? null);
for (const field of Object.keys(DECK_FIELDS)) model.radiation[field].set(savedDeckField(saved, field, model));
for (const field of Object.keys(MOIST_FIELDS)) model.moist[field].set(savedMoistField(saved, field, model));
for (const field of Object.keys(RADIATION_FIELDS)) model.radiation[field].set(savedRadiationField(saved, field, model));
for (const [key, name] of [['boundaryDepth', 'depth'], ['mixingTop', 'mixingTop'], ['boundaryRegime', 'regime'], ['boundaryBuoyancy', 'buoyancyFlux']]) if (saved[key]) model.boundaryLayer[name].set(saved[key]);
if (saved.windSpeed) model.surface.windSpeed.set(saved.windSpeed);
if (saved.evaporation) model.radiation.evaporation.set(saved.evaporation);
if (saved.exchangeHeat && model.exchange && !model.exchange.fixed) model.exchange.heat.set(saved.exchangeHeat);
if (saved.exchangeWind && model.exchange) model.exchange.wind.set(saved.exchangeWind);
for (const field of ['cumulusCover', 'cumulusWater']) if (saved[field] && saved[field].length === model.moist[field].length) model.moist[field].set(saved[field]);
if (saved.subcloudVirtual && saved.subcloudVirtual.length === model.moist.subcloudVirtual.length) model.moist.subcloudVirtual.set(saved.subcloudVirtual); else model.moist.subcloudVirtual.fill(0);
model.time = saved.time;
const land = model.geography.land;
let changed = 0, shift = 0;
if (SST === 'woa') {
  const sst = decodeMonthlySst(readFileSync(new URL('../data/woa_sst_monthly_1deg.bin', import.meta.url)));
  const months = monthWeights(sst, (((saved.day + 78) % 365) + 365) % 365);
  for (let i = 0; i < C; i++) {
    if (land[i] || state[6][i] > 0) continue;
    const t = climatologyAt(sst, months, mesh.latCell[i], mesh.lonCell[i]);
    if (!Number.isFinite(t)) continue;
    shift += t - state[3][i]; state[3][i] = t; changed++;
  }
}
console.log(`${TAG}: day ${saved.day}, SST ${SST}${changed ? ` set on ${changed} open-sea cells, mean shift ${(shift / changed).toFixed(2)} K` : ''}, slab ${CAP.toExponential(0)} J/m²/K, ${DAYS} days`);
model.load({ ocean: false });
model.land.load(saved.land);
mkdirSync(OUT, { recursive: true });
const seaMean = () => { let s = 0, n = 0; for (let i = 0; i < C; i++) if (!land[i]) { s += state[3][i]; n++; } return s / n; };
const t0 = performance.now();
for (let day = 1; day <= DAYS; day++) {
  for (let done = 0; done < perDay;) { const count = Math.min(BATCH, perDay - done); await model.stepBatch(count, dt); done += count; }
  if (day >= SAVE_FROM && (day % SAVE_DAYS === 0 || day === DAYS)) {
    await model.sync();
    const [pi, theta, u, surfaceT, q, qc, ice] = state, { concentration } = model.seaIce, r = model.radiation, b = model.boundaryLayer;
    const header = { N: saved.N, ...gridHeader(mesh), K: core.K, day: saved.day + day, time: model.time, terrain: !!model.surfaceGeopotential, levels: core.levels };
    const name = `${OUT}/${TAG}_day${String(saved.day + day).padStart(4, '0')}.bin`;
    writeFileSync(name, encodeState({ ...header, pi, theta, u, surfaceT, q, qc, ice, concentration, mlmSubsidence: r.mlmSubsidence, mlmHeight: r.mlmHeight, mlmGate: r.mlmGate, convectiveRain: model.moist.convectiveRain, largeScaleRain: model.moist.largeScaleRain, meanAbsorbedSolar: r.meanAbsorbedSolar, meanOutgoingLongwave: r.meanOutgoingLongwave, meanPlanetaryAlbedo: r.meanPlanetaryAlbedo, boundaryDepth: b.depth, mixingTop: b.mixingTop, boundaryRegime: b.regime, boundaryBuoyancy: b.buoyancyFlux, windSpeed: model.surface.windSpeed, evaporation: r.evaporation, ...(model.exchange ? { ...(model.exchange.fixed ? {} : { exchangeHeat: model.exchange.heat }), exchangeWind: model.exchange.wind } : {}), cumulusCover: model.moist.cumulusCover, cumulusWater: model.moist.cumulusWater, subcloudVirtual: model.moist.subcloudVirtual, land: await model.land.serialize() }));
    console.log(`day +${day}: sea-mean surface ${(seaMean() - 273.15).toFixed(3)} °C, ${((performance.now() - t0) / 1000 / day).toFixed(1)} s/day, saved ${name.replace(/.*\//, '')}`);
  }
}
process.exit(0);
