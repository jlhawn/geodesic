// The Arctic (60-90N) surface sunlight of a saved state on the CPU engine, all-sky and cloud-free, and the cloud effect over the sea ice:
//   [OFFSETS=-25,-15,-5] [HOURS=8] node scripts/arcticSunlight.mjs STATE.bin [OBSDIR]
// Each column north of 60N is lit at HOURS instants on each day OFFSETS days from the save (the sun of the 30-day window the save closes),
// the state held as saved, through radiation.apply one cell at a time, once as is and once under a cloud mask that takes the resolved
// cloud, cumulus and deck away. The summary rows are area-weighted means over the ice70 cells: sea cells north of 70N with a
// concentration above 0.15, by the saved state's concentration (the snapshot set) and by the save's 30-day mean concentration (the window set).
//   snapshot cloud effect = all-sky offline sunlight minus cloud-free offline sunlight, the state as saved;
//   window cloud effect   = the save's 30-day mean all-sky surface sunlight minus the same cloud-free sunlight.
// With OBSDIR (the cache scripts/arcticObs.mjs writes) the window cloud effect of ERA5 (msdwswrf - msdwswrfcs, bilinear on its 1-degree
// subsample) and of CERES-EBAF (rsds - rsdscs, nearest cell) at the window set's cells follows, the months weighted by the window's days.
import { readFileSync, existsSync } from 'node:fs';
import { stateGeography } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels, stateTerrain } from '../js/stateFile.module.js';
import { savedDeckField, DECK_FIELDS, savedMoistField, MOIST_FIELDS } from '../js/physics/regrid.module.js';
import { placeCloudIce, placeCumulus } from '../js/physics/moist.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { loadSavedLake, lakeOptionsFor, savedIceOptions } from './savedLand.mjs';
import { refuseUnreadOptions } from './optionEnvironment.mjs';

refuseUnreadOptions(import.meta.url);

const [FILE, OBS] = process.argv.slice(2);
const OFFSETS = (process.env.OFFSETS ?? '-25,-15,-5').split(',').map(Number), HOURS = Number(process.env.HOURS ?? 8), LAT0 = 60, DEG = 180 / Math.PI;
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const ground = stateGeography(saved), { topography } = ground;
const model = createModel(stateGrid(saved, topography).mesh, { lake: lakeOptionsFor(saved), ice: savedIceOptions(saved), terrain: stateTerrain(saved), topography, geography: ground.geography, land: ground.land, levels: savedLevels(saved), ocean: false, radiation: { ...(saved.options?.radiation ?? {}) }, moist: { ...(saved.options?.moist ?? {}) }, boundaryLayer: { ...(saved.options?.boundaryLayer ?? {}) } });
const { mesh, core, state, radiation, boundaryLayer: bl, seaIce, land, geography } = model;
const { K, C } = core.diagnostics, g = core.diagnostics.g;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
for (const field of Object.keys(DECK_FIELDS)) radiation[field].set(savedDeckField(saved, field, model));
for (const field of Object.keys(MOIST_FIELDS)) model.moist[field].set(savedMoistField(saved, field, model));
land.load({ soil: Float64Array.from(saved.land.soil), snow: Float64Array.from(saved.land.snow), ...(saved.land.vegetation ? { vegetation: Float64Array.from(saved.land.vegetation) } : {}), ...(saved.land.surface ? { surface: Float64Array.from(saved.land.surface) } : {}), ...(saved.land.snowAlbedo ? { snowAlbedo: Float64Array.from(saved.land.snowAlbedo) } : {}) }, state[6]);
loadSavedLake(model, saved);
model.time = saved.time;
const zs = Float64Array.from({ length: C }, (_, i) => (model.surfaceGeopotential ? model.surfaceGeopotential[i] / g : 0));
if (saved.boundaryDepth) bl.depth.set(Float64Array.from(saved.boundaryDepth, (z, i) => z + zs[i]));
if (saved.mixingTop) bl.mixingTop.set(Float64Array.from(saved.mixingTop, (z, i) => z + zs[i]));
if (saved.boundaryRegime) bl.regime.set(saved.boundaryRegime);
if (saved.boundaryBuoyancy) bl.buoyancyFlux.set(saved.boundaryBuoyancy);
if (saved.windSpeed) model.surface.windSpeed.set(saved.windSpeed);
placeCumulus(model, saved);
placeCloudIce(model, saved);
const [pi, theta, , surfaceT, q, qc, ice] = state;
for (let i = 0; i < C; i++) core.diagnoseColumn(i, pi, theta, q, qc);
const landMask = geography.land, cells = [];
for (let i = 0; i < C; i++) if (mesh.latCell[i] * DEG >= LAT0) cells.push(i);
const n = cells.length, swAll = new Float64Array(n), swClear = new Float64Array(n), insol = new Float64Array(n), samples = new Float64Array(n);
const fluxT = Float64Array.from(surfaceT), adir = new Float64Array(C), adif = new Float64Array(C), openSea = new Float64Array(C), zeroQc = new Float64Array(K * C);
const out = [null, new Float64Array(K * C), null, null, new Float64Array(K * C)];
const stateAll = [pi, theta, state[2], fluxT, q, qc], stateClear = [pi, theta, state[2], fluxT, q, zeroQc];
const clearMask = { resolved: new Float64Array(K), cumulus: new Float64Array(K), fraction: 0, deck: 0 };
for (const i of cells) {
  if (landMask[i]) { adif[i] = land.albedo(i); openSea[i] = 0; continue; }
  const h = ice[i], cover = seaIce.cover(i, h);
  adif[i] = seaIce.albedo(h, null, seaIce.snow[i], cover, surfaceT[i], seaIce.snowAlbedo[i]);
  openSea[i] = 1 - cover;
  if (h > 0 && cover < 1) fluxT[i] = cover * surfaceT[i] + (1 - cover) * 271.35;
}
for (const off of OFFSETS) for (let h = 0; h < HOURS; h++) {
  radiation.setTime(saved.time + off * 86400 + (h + 0.5) * 86400 / HOURS);
  for (let c = 0; c < n; c++) {
    const i = cells[c], mu = radiation.cosZenith(i);
    samples[c]++;
    if (!(mu > 0)) continue;
    if (landMask[i]) adir[i] = adif[i]; else { const hh = ice[i], cover = seaIce.cover(i, hh); adir[i] = seaIce.albedo(hh, mu, seaIce.snow[i], cover, surfaceT[i], seaIce.snowAlbedo[i]); }
    radiation.useCloudMask(null);
    radiation.apply(stateAll, out, model.surface.windSpeed, null, i, i + 1, adir, adif, null, openSea, bl.depth, 0);
    swAll[c] += radiation.surfaceShortwave[i]; insol[c] += radiation.budget.insolation;
    radiation.useCloudMask(clearMask);
    radiation.apply(stateClear, out, model.surface.windSpeed, null, i, i + 1, adir, adif, null, openSea, bl.depth, 0);
    swClear[c] += radiation.surfaceShortwave[i];
    if (radiation.budget.cloudCover > 0 || radiation.budget.stratus > 0) throw new Error(`cell ${i}: the clear pass keeps cloud ${radiation.budget.cloudCover} deck ${radiation.budget.stratus}`);
  }
}
radiation.useCloudMask(null);

const day = Math.round(saved.time / 86400), file = FILE.split('/').pop();
const meanShortwave = saved.means?.shortwave ?? null, meanConcentration = saved.means?.concentration ?? null;
const sets = {
  snapshot: (i) => !landMask[i] && mesh.latCell[i] * DEG >= 70 && saved.concentration?.[i] > 0.15,
  window: (i) => !landMask[i] && mesh.latCell[i] * DEG >= 70 && meanConcentration?.[i] > 0.15,
};
const average = (test, value) => { let s = 0, w = 0; for (let c = 0; c < n; c++) { const i = cells[c]; if (!test(i)) continue; s += mesh.areaCell[i] * value(c, i); w += mesh.areaCell[i]; } return s / w; };
const result = { file, day };
for (const [name, test] of Object.entries(sets)) {
  const count = cells.filter(test).length, all = average(test, (c) => swAll[c] / samples[c]), clear = average(test, (c) => swClear[c] / samples[c]), insolation = average(test, (c) => insol[c] / samples[c]);
  const window = meanShortwave ? average(test, (c, i) => meanShortwave[i]) : NaN;
  result[name] = { cells: count, insolation, allSky: all, clear, snapshotCloudEffect: all - clear, windowAllSky: window, windowCloudEffect: window - clear };
  console.log(`${file} day ${day}, ice70 by the ${name} concentration (${count} cells): insolation ${insolation.toFixed(1)}; offline all-sky ${all.toFixed(1)}, cloud-free ${clear.toFixed(1)}, snapshot cloud effect ${(all - clear).toFixed(2)}; window mean all-sky ${window.toFixed(1)}, window cloud effect ${(window - clear).toFixed(2)} W/m2`);
}

if (OBS) {
  const f32 = (name) => { const path = `${OBS}/${name}.f32`; if (!existsSync(path)) return null; const b = readFileSync(path); return new Float32Array(b.buffer, b.byteOffset, b.length / 4); };
  const era = { all: f32('era5_msdwswrf'), clear: f32('era5_msdwswrfcs') }, ceres = { all: f32('ceres_rsds'), clear: f32('ceres_rsdscs') };
  const weights = new Array(12).fill(0);
  for (let d = day - 29; d <= day; d++) weights[new Date(Date.UTC(2019, 2, 20) + (((d % 365) + 365) % 365) * 86400e3).getUTCMonth()] += 1 / 30;
  const eraAt = (a, m, la, lo) => {
    const r = Math.min(29.999, Math.max(0, 90 - la)), r0 = Math.floor(r), fr = r - r0, x = ((lo % 360) + 360) % 360, c0 = Math.floor(x) % 360, c1 = (c0 + 1) % 360, fc = x - Math.floor(x), at = (rr, cc) => a[(m * 31 + rr) * 360 + cc];
    return (1 - fr) * ((1 - fc) * at(r0, c0) + fc * at(r0, c1)) + fr * ((1 - fc) * at(r0 + 1, c0) + fc * at(r0 + 1, c1));
  };
  const ceresAt = (a, m, la, lo) => a[(m * 30 + Math.min(29, Math.max(0, Math.floor(la - 60)))) * 360 + Math.floor(((lo % 360) + 360) % 360)];
  const sample = (a, at, i) => { let s = 0; for (let m = 0; m < 12; m++) if (weights[m]) s += weights[m] * at(a, m, mesh.latCell[i] * DEG, mesh.lonCell[i] * DEG); return s; };
  for (const [name, test] of Object.entries(sets)) {
    if (!era.all || !era.clear || !ceres.all || !ceres.clear) { console.log(`${OBS} lacks the ERA5 or CERES files`); break; }
    const eraEffect = average(test, (c, i) => sample(era.all, eraAt, i) - sample(era.clear, eraAt, i)), ceresEffect = average(test, (c, i) => sample(ceres.all, ceresAt, i) - sample(ceres.clear, ceresAt, i));
    result[name].eraCloudEffect = eraEffect; result[name].ceresCloudEffect = ceresEffect;
    console.log(`${file} day ${day}, ice70 by the ${name} concentration: window cloud effect of ERA5 ${eraEffect.toFixed(2)}, of CERES-EBAF ${ceresEffect.toFixed(2)} W/m2 (model window ${result[name].windowCloudEffect.toFixed(2)}, snapshot ${result[name].snapshotCloudEffect.toFixed(2)})`);
  }
}
console.log(JSON.stringify(result));
