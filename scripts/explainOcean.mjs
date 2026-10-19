// Precomputes the moving sea that explain/ocean.html plays back: Part 7's
// wet planet (explainLayers.mjs wet) with its layered ocean let go, over a
// flat 4000 m floor, on the GPU engine (Dawn in Node):
//   node scripts/explainOcean.mjs [SPINUP_DAYS] [RECORD_DAYS] [OUT]
// writing explain/data/ocean_N16.bin (3650 days of spin-up and 365
// recorded by default) and, at every March equinox of the spin-up, the
// whole model state to CHECKPOINT (runs/explainOcean_N16_day{day}.bin by
// default, {day} the run day; only the newest is kept unless
// KEEP_CHECKPOINTS=all). TEST_CHECKPOINT_DAYS=d1,d2 adds checkpoints at
// those spin-up run days. RESUME=<checkpoint> continues from such a state
// with SPINUP_DAYS more days of spin-up, which must end on a March
// equinox, then the record; it refuses a checkpoint written by another
// engine (other js/ contents) unless FORCE_RESUME=1, and stops if the
// loaded state's stores differ from those saved.
// WET_ICE (JSON, Part 7's dark ice by default), OCEAN_OPTS (JSON merged
// into the ocean options) and SURFACE_OPTS (JSON surface options, such as
// warmLayer and coolSkin) set the options, over the checkpoint's on RESUME.
// SKIN_WATCH, SKIN_FIELD_DAYS and SKIN_WATCH_FROM keep scripts/skinWatch.mjs's
// record.
import { writeFileSync, readFileSync, mkdirSync, renameSync, unlinkSync, existsSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { dirname, basename, relative as relativeTo, resolve, isAbsolute } from 'node:path';
import { execFileSync } from 'node:child_process';
import { Grid } from '../js/grid.module.js';
import { stateGrid, gridHeader } from '../js/stateGrid.module.js';
import { createSigmaCore, sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { FREEZING_POINT, MELTING_POINT } from '../js/physics/ice.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { withCadence } from '../js/cadence.module.js';
import { createGpuModel } from '../js/gpu/model.gpu.js';
import { createForcingRecorder } from '../js/gpu/forcing.gpu.js';
import { watchSkin } from './skinWatch.mjs';
import { edgeLayerFluxes } from './oceanEdgeFlux.mjs';
import { readRanges, getDevice } from '../js/gpu/device.module.js';
import { decodeForcing } from '../js/forcing.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';
import { depthFields } from '../js/ocean/layered.module.js';
import { savedDeckField, DECK_FIELDS, savedMoistField, MOIST_FIELDS, savedRadiationField, RADIATION_FIELDS } from '../js/physics/regrid.module.js';

const FLOAT32_TOLERANCE = 2 ** -21, SECTION_MIN_THICKNESS = 0.1, NET_FRACTION = 0.1;
const N = 16, DAY = 86400, STEP = 450, GRID = 'bl36', FRAME_DAYS = 5, SPINUP_FRAMES = 60, BINS = 36, CAP_STEP = 5, UPWELLING_DEPTH = 50, CELSIUS = 273.15;
const SPINUP_DAYS = Number(process.argv[2] ?? 3650), RECORD_DAYS = Number(process.argv[3] ?? 365);
if (![SPINUP_DAYS, RECORD_DAYS].every((days) => Number.isInteger(days) && days >= 0) || SPINUP_DAYS + RECORD_DAYS === 0) { console.error('usage: node scripts/explainOcean.mjs [SPINUP_DAYS] [RECORD_DAYS] [OUT]'); process.exit(1); }
const out = process.argv[4] ?? new URL(`../explain/data/ocean_N${N}.bin`, import.meta.url).pathname;
const DARK_ICE = { iceAlbedo: 0.07, meltingIceAlbedo: 0.07, iceSnowAlbedo: 0.07, slabHeatCapacity: 4.2e7 };
const EQUINOX_DAY_OF_YEAR = 78;
const MONTHS = [['January', 31], ['February', 28], ['March', 31], ['April', 30], ['May', 31], ['June', 30], ['July', 31], ['August', 31], ['September', 30], ['October', 31], ['November', 30], ['December', 31]];
const monthStart = MONTHS.map((_, m) => MONTHS.slice(0, m).reduce((sum, [, days]) => sum + days, 0));
const monthOf = (day) => { const date = (Math.floor(day) + EQUINOX_DAY_OF_YEAR) % 365; let m = 11; while (monthStart[m] > date) m--; return m; };
const monthDays = (m) => { const from = (monthStart[m] - EQUINOX_DAY_OF_YEAR + 365) % 365, to = from + MONTHS[m][1]; return to <= 365 ? [[from, to]] : [[from, 365], [0, to - 365]]; };
const ROOT = new URL('..', import.meta.url).pathname;
const git = (...args) => { try { return execFileSync('git', args, { cwd: ROOT, encoding: 'utf8', maxBuffer: 1 << 28 }); } catch { return null; } };
const lines = (text) => (text ?? '').trimEnd().split('\n').filter(Boolean);
const sha256 = (data) => createHash('sha256').update(data).digest('hex');
const shown = (path) => { if (!path) return null; const inside = relativeTo(ROOT, resolve(path)); return inside && !inside.startsWith('..') && !isAbsolute(inside) ? inside : basename(path); };
const jsDiff = git('diff', 'HEAD', '--', 'js') ?? '';
const jsFiles = lines(git('ls-files', '--cached', '--others', '--exclude-standard', '--', 'js')).filter((path) => existsSync(ROOT + path)).sort();
const jsBlobs = jsFiles.length ? lines(execFileSync('git', ['hash-object', '--stdin-paths'], { cwd: ROOT, input: jsFiles.join('\n'), encoding: 'utf8' })) : [];
const engine = {
  commit: (git('rev-parse', 'HEAD') ?? '').trim() || null, jsSha256: jsFiles.length ? sha256(jsFiles.map((path, n) => `${jsBlobs[n]} ${path}\n`).join('')) : null,
  jsDiffSha256: jsDiff ? sha256(jsDiff) : null, jsDiffStat: lines(git('diff', '--stat', 'HEAD', '--', 'js')),
  jsUntracked: lines(git('ls-files', '--others', '--exclude-standard', '--', 'js')).map((path) => ({ path, sha256: sha256(readFileSync(ROOT + path)) })),
  scriptSha256: sha256(readFileSync(new URL(import.meta.url))),
};
const sameEngine = (saved) => {
  if (saved.jsSha256) return saved.jsSha256 === engine.jsSha256;
  const diff = /^[0-9a-f]{40}$/.test(saved.commit ?? '') ? git('diff', saved.commit, '--', 'js') : null;
  return diff !== null && (diff ? sha256(diff) : null) === saved.jsDiffSha256 && JSON.stringify(saved.jsUntracked ?? []) === JSON.stringify(engine.jsUntracked);
};
const commit = (git('rev-parse', '--short', 'HEAD') ?? '').trim() || null, modified = lines(git('status', '--porcelain', '--', 'js', 'scripts')), startedAt = new Date().toISOString();

const resumed = process.env.RESUME ? await decodeState(new Uint8Array(readFileSync(process.env.RESUME))) : null;
if (resumed && resumed.kind !== 'explainOcean') throw new Error(`${process.env.RESUME} is not a checkpoint of scripts/explainOcean.mjs`);
const savedEngine = resumed && resumed.engineJson ? JSON.parse(resumed.engineJson) : null;
const engineMatches = !!savedEngine && sameEngine(savedEngine), scriptMatches = !!savedEngine && savedEngine.scriptSha256 === engine.scriptSha256;
if (resumed && !engineMatches) {
  const why = !savedEngine ? 'it does not record the engine that wrote it'
    : savedEngine.jsSha256 ? `its js/ hashes to ${savedEngine.jsSha256} (commit ${savedEngine.commit}), this tree's to ${engine.jsSha256} (commit ${engine.commit})`
    : `it was written at ${savedEngine.commit} with js/ diff ${savedEngine.jsDiffSha256 ?? 'none'} and untracked ${JSON.stringify(savedEngine.jsUntracked)}, which this tree's js/ does not reproduce`;
  if (process.env.FORCE_RESUME !== '1') { console.error(`${process.env.RESUME}: ${why}. FORCE_RESUME=1 resumes it anyway.`); process.exit(1); }
  console.error(`resuming across engines (FORCE_RESUME=1): ${why}`);
}
if (resumed && !scriptMatches) console.error(`warning: ${process.env.RESUME} was written by ${savedEngine?.scriptSha256 ? `scripts/explainOcean.mjs ${savedEngine.scriptSha256}, this is ${engine.scriptSha256}` : 'a scripts/explainOcean.mjs that did not record its hash'}; resuming, and the header's resumeCheck.sameScript is false`);
const base = resumed ? JSON.parse(resumed.optionsJson) : { ocean: { everySteps: 4 }, ice: DARK_ICE };
const ICE = process.env.WET_ICE ? JSON.parse(process.env.WET_ICE) : base.ice;
const OCEAN = { ...base.ocean, ...JSON.parse(process.env.OCEAN_OPTS ?? '{}') };
const cadence = withCadence({ ocean: OCEAN }, STEP, N);
const options = { ocean: cadence.ocean, levels: sigmaInterfaces(GRID), radiation: { ...cadence.radiation, clearSkyPass: true }, dragEvery: cadence.dragEvery, ice: ICE, surface: JSON.parse(process.env.SURFACE_OPTS ?? '{}') };
const runStartTime = resumed ? resumed.runStartTime : (Math.ceil(SPINUP_DAYS / 365) * 365 - SPINUP_DAYS) * DAY;
const invocationStart = resumed ? resumed.time : runStartTime, recordStart = invocationStart + SPINUP_DAYS * DAY, end = RECORD_DAYS * DAY;
if (Math.round(recordStart / DAY) % 365 !== 0) { console.error(`the checkpoint is at model day ${Math.round(invocationStart / DAY)}; a record must start on a March equinox, ${365 - (Math.round(invocationStart / DAY) % 365)} days on (or that plus whole years)`); process.exit(1); }
const segments = [...(resumed ? JSON.parse(resumed.segmentsJson).map(({ from, ...s }) => ({ fromTime: from, ...s })) : []), { fromTime: invocationStart, spinupDays: SPINUP_DAYS, recordDays: RECORD_DAYS, ocean: OCEAN, ice: ICE, ...engine, date: startedAt }];
const checkpointPath = (day) => (process.env.CHECKPOINT ?? `${ROOT}runs/explainOcean_N${N}_day{day}.bin`).replaceAll('{day}', String(day).padStart(5, '0'));
const testCheckpointDays = new Set((process.env.TEST_CHECKPOINT_DAYS ?? '').split(',').filter(Boolean).map(Number));
const firstDay = Math.round((invocationStart - runStartTime) / DAY), lastSpinupDay = Math.round((recordStart - runStartTime) / DAY);
if (![...testCheckpointDays].every((d) => Number.isInteger(d) && d > firstDay && d <= lastSpinupDay)) { console.error(`TEST_CHECKPOINT_DAYS lists run days of this spin-up, ${firstDay + 1} to ${lastSpinupDay}`); process.exit(1); }
const keepCheckpoints = process.env.KEEP_CHECKPOINTS === 'all';

const t0 = performance.now(), seconds = () => (performance.now() - t0) / 1000;
const model = await createGpuModel(resumed ? stateGrid(resumed, null).mesh : new Grid(N), options);
const { mesh, core, state, seaIce, gpu } = model, ocean = model.oceanEngine, { device } = gpu;
const [pi, theta, u, surfaceT, q, qc, ice] = state;
const C = mesh.nCells, E = mesh.nEdges, K = core.K, L = ocean.layers, PH = gpu.layout.PH, SL = gpu.layout.S, OD = ocean.layout.OD;
const { cellsOnEdge, dvEdge, dcEdge, areaCell: weight } = mesh;
const rhoCp = ocean.options.density * ocean.options.specificHeat, LATENT = model.moist.latentHeat, FUSION = seaIce.latentHeatFusion;
const airCp = core.diagnostics.cp, gravity = core.diagnostics.g, { dSigma } = core.diagnostics;
const D = ocean.D, cellOcean = ocean.cellOcean, everySteps = ocean.everySteps, perDay = Math.round(DAY / STEP), quarter = perDay / 4;
if (process.env.SKIN_WATCH) watchSkin(model, { prefix: process.env.SKIN_WATCH, sea: (i) => !!cellOcean[i], surface: options.surface, fieldDays: (process.env.SKIN_FIELD_DAYS ?? '').split(',').filter(Boolean), from: Number(process.env.SKIN_WATCH_FROM ?? 0) });

const restartNames = [];
if (resumed) {
  ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => state[a].set(resumed[name]));
  seaIce.load(ice, resumed.concentration);
  for (const field of Object.keys(DECK_FIELDS)) model.radiation[field].set(savedDeckField(resumed, field, model));
  for (const field of Object.keys(MOIST_FIELDS)) model.moist[field].set(savedMoistField(resumed, field, model));
  for (const field of Object.keys(RADIATION_FIELDS)) model.radiation[field].set(savedRadiationField(resumed, field, model));
  model.boundaryLayer.depth.set(resumed.boundaryDepth); model.boundaryLayer.mixingTop.set(resumed.mixingTop);
  model.boundaryLayer.regime.set(resumed.boundaryRegime); model.boundaryLayer.buoyancyFlux.set(resumed.boundaryBuoyancy);
  model.surface.windSpeed.set(resumed.windSpeed); model.radiation.evaporation.set(resumed.evaporation);
  if (model.exchange && resumed.exchangeHeat && !model.exchange.fixed) model.exchange.heat.set(resumed.exchangeHeat);
  if (model.exchange && resumed.exchangeWind) model.exchange.wind.set(resumed.exchangeWind);
  for (const field of ['cumulusCover', 'cumulusWater', 'subcloudVirtual']) model.moist[field].set(resumed[field]);
  if (model.skin) model.skin.load(resumed);
  model.time = resumed.time;
  model.load({ ocean: false });
  model.ocean.load(resumed.ocean, surfaceT, ice);
  if (resumed.eddyV && resumed.eddyV.length === L * E) device.queue.writeBuffer(ocean.buffers.OD, 4 * OD.EDDYV, Float32Array.from(resumed.eddyV));
  device.queue.writeBuffer(gpu.buffers.PH, 4 * PH.SNOW, Float32Array.from(resumed.seaSnow));
  device.queue.writeBuffer(gpu.buffers.PH, 4 * PH.SNOWALB, Float32Array.from(resumed.seaSnowAlbedo));
  const loaded = await model.ocean.serialize({ restart: true }), differing = [];
  restartNames.push(...new Set([...Object.keys(resumed.ocean), ...Object.keys(loaded)].filter((name) => !['T', 'S', 'densities'].includes(name))));
  for (const name of restartNames) {
    const want = resumed.ocean[name], got = loaded[name];
    if ((want == null && (name === 'penetrated' || name === 'remainder') && got && got.every((x) => x === 0)) || (name === 'transmission' && (want == null || got == null))) continue;
    if (want == null || got == null || want.length !== got.length) { differing.push(`${name} (${want == null ? 'missing' : want.length} saved, ${got == null ? 'missing' : got.length} loaded)`); continue; }
    let count = 0;
    for (let n = 0; n < want.length; n++) if (Math.fround(want[n]) !== Math.fround(got[n]) && !(Number.isNaN(want[n]) && Number.isNaN(got[n]))) count++;
    if (count) differing.push(`${name} (${count} of ${want.length})`);
  }
  if (differing.length) { console.error(`${process.env.RESUME}: the ocean did not load as saved: ${differing.join(', ')}; stopping`); process.exit(1); }
} else {
  initializeState(model).forEach((values, s) => state[s].set(values));
  model.time = runStartTime;
  model.load();
}
const recorder = await createForcingRecorder(model);
const look = createSigmaCore(mesh, { levels: core.levels }), { exnerLayer: lookExner } = look.diagnostics;

const area = weight.reduce((sum, x) => sum + x, 0), radius = Math.sqrt(area / (4 * Math.PI)), bandWidth = Math.PI * radius / BINS;
const bin = Int32Array.from(mesh.latCell, (lat) => Math.min(BINS - 1, Math.floor((lat + Math.PI / 2) / Math.PI * BINS)));
const binArea = new Float64Array(BINS), meanFloor = D.reduce((sum, depth, i) => sum + weight[i] * depth, 0) / area;
for (let i = 0; i < C; i++) binArea[bin[i]] += weight[i];
const eastUnit = new Float64Array(3 * C), northUnit = new Float64Array(3 * C);
for (let i = 0; i < C; i++) {
  const x = mesh.xCell[3 * i], y = mesh.xCell[3 * i + 1], z = mesh.xCell[3 * i + 2], rho = Math.hypot(x, y) || 1e-12;
  eastUnit.set([-y / rho, x / rho, 0], 3 * i); northUnit.set([-z * x / rho, -z * y / rho, rho], 3 * i);
}
const project = (vec, i, unit) => vec[3 * i] * unit[3 * i] + vec[3 * i + 1] * unit[3 * i + 1] + vec[3 * i + 2] * unit[3 * i + 2];
const global = (value) => { let sum = 0; for (let i = 0; i < C; i++) sum += weight[i] * value(i); return sum / area; };
const hemisphere = (value, north) => { let sum = 0, total = 0; for (let i = 0; i < C; i++) if ((mesh.latCell[i] >= 0) === north) { sum += weight[i] * value(i); total += weight[i]; } return sum / total; };
const bands = (value) => {
  const sum = new Float64Array(BINS);
  for (let i = 0; i < C; i++) sum[bin[i]] += weight[i] * value(i);
  return Array.from(sum, (v, b) => v / binArea[b]);
};
const CAPS = Array.from({ length: 180 / CAP_STEP - 1 }, (_, c) => -90 + CAP_STEP * (c + 1));
const capEdges = CAPS.map((degrees) => {
  const phi = degrees * Math.PI / 180, list = [];
  for (let e = 0; e < E; e++) { const south = mesh.latCell[cellsOnEdge[2 * e]] < phi; if (south !== mesh.latCell[cellsOnEdge[2 * e + 1]] < phi) list.push([e, south ? 1 : -1]); }
  return list;
});
const capSouth = CAPS.map((degrees) => { const list = []; for (let i = 0; i < C; i++) if (mesh.latCell[i] < degrees * Math.PI / 180) list.push(i); return list; });
const binCells = new Int32Array(BINS);
for (let i = 0; i < C; i++) binCells[bin[i]]++;
const exactRatio = Array.from(binArea, (a, b) => a / (2 * Math.PI * radius ** 2 * (Math.sin(-Math.PI / 2 + (b + 1) * Math.PI / BINS) - Math.sin(-Math.PI / 2 + b * Math.PI / BINS))));
const PARTITION = `bins and caps partition the cells by center latitude: a cell belongs to the 5° bin its center falls in, a cap holds the cells whose centers lie south of its latitude, and the edges whose two cells lie on opposite sides carry the transport across it; binArea (m²) is the summed area of the bin's cells, from ${(100 * (Math.min(...exactRatio) - 1)).toFixed(0)} % to +${(100 * (Math.max(...exactRatio) - 1)).toFixed(0)} % off the exact 5° band's at N = ${N}, so integrate band means with binArea`;
const LABELS = Array.from(ocean.table.densities), FIXED = ocean.table.fixed ?? 0;
const groupOf = (k) => (k < FIXED ? 0 : LABELS[k] < 1026 ? 1 : LABELS[k] < 1026.8 ? 2 : 3);
const GROUPS = ['fixed-depth layers', 'layers lighter than 1026', 'layers 1026 to 1026.8', 'layers from 1026.8'];
const edgeFluxes = edgeLayerFluxes(mesh, { layers: L, fixed: FIXED, options: ocean.options, cellOcean, depth: D, sill: ocean.sillEdge });
const edgeOut = { volume: new Float64Array(L), net: new Float64Array(L), heat: new Float64Array(L), diffusion: new Float64Array(L), donor: new Int32Array(L) };

async function oceanSnapshot() {
  const [o, [eddy], [iceNow, skin], [concentration]] = await Promise.all([
    ocean.download(),
    readRanges(device, ocean.buffers.OD, [{ offset: OD.EDDYV, length: L * E }]),
    readRanges(device, gpu.buffers.S, [{ offset: SL.ICE, length: C }, { offset: SL.TS, length: C }]),
    readRanges(device, gpu.buffers.PH, [{ offset: PH.CONC, length: C }]),
  ]);
  return { ...o, eddy: ocean.options.eddyDiffusivity > 0 ? Float64Array.from(eddy) : new Float64Array(L * E), ice: Float64Array.from(iceNow), surfaceT: Float64Array.from(skin), concentration: Float64Array.from(concentration) };
}
const coverOf = (iceH, concentration) => (iceH > 0 ? (concentration > 0 ? concentration : 1) : 0);

const STORES = ['ocean', 'ice', 'airEnthalpy', 'airKinetic', 'latent', 'vapor', 'condensate', 'salt', 'iceVolume'];
const store = Object.fromEntries(STORES.map((name) => [name, new Float64Array(C)]));
const kineticLayer = new Float64Array(C);
function storeCells(snap, snow, oceanFlux, pending) {
  look.diagnose(pi, theta);
  const remaining = (everySteps - 1 - (model.oceanCounter % everySteps)) * STEP, { h, T, S } = snap;
  for (let i = 0; i < C; i++) {
    let heat = 0, salt = 0;
    for (let k = 1; k < L; k++) { const n = k * C + i; heat += h[n] * (T[n] - CELSIUS); salt += h[n] * S[n]; }
    const iced = ice[i] > 0, cover = seaIce.cover(i, ice[i]);
    store.ocean[i] = rhoCp * (heat + h[i] * ((iced ? FREEZING_POINT : surfaceT[i]) - CELSIUS)) + oceanFlux[i] * remaining + pending[i];
    store.salt[i] = salt + h[i] * S[i];
    store.ice[i] = iced ? seaIce.energy(surfaceT[i], ice[i], snow[i], cover) : 0;
    store.iceVolume[i] = cover * ice[i];
    let enthalpy = 0;
    for (let k = 0; k < K; k++) { const n = k * C + i; enthalpy += airCp * theta[n] * lookExner[n] * pi[i] * dSigma[k] / gravity; }
    store.airEnthalpy[i] = enthalpy;
    store.vapor[i] = model.moist.columnWater(pi, q, i); store.condensate[i] = model.moist.columnWater(pi, qc, i);
    store.latent[i] = LATENT * store.vapor[i];
  }
  store.airKinetic.fill(0);
  for (let k = 0; k < K; k++) {
    kineticLayer.fill(0);
    for (let e = 0; e < E; e++) { const w = 0.25 * dcEdge[e] * dvEdge[e] * u[k * E + e] ** 2; kineticLayer[cellsOnEdge[2 * e]] += w; kineticLayer[cellsOnEdge[2 * e + 1]] += w; }
    for (let i = 0; i < C; i++) store.airKinetic[i] += kineticLayer[i] / weight[i] * pi[i] * dSigma[k] / gravity;
  }
}
const totalsNow = () => Object.fromEntries(STORES.map((name) => [name, global((i) => store[name][i])]));
const energyOf = (t) => t.ocean + t.ice + t.airEnthalpy + t.latent;
const surfaceOf = (t) => t.ocean + t.ice;

const FIELDS = [
  { name: 'seaSurfaceTemperature', unit: '°C', per: 100, offset: 0 }, { name: 'currentEast', unit: 'm/s', per: 1000, offset: 0 }, { name: 'currentNorth', unit: 'm/s', per: 1000, offset: 0 },
  { name: 'mixedLayerDepth', unit: 'm', per: 10, offset: 0 }, { name: 'seaSurfaceHeight', unit: 'm', per: 1000, offset: 0 }, { name: 'upwelling50', unit: 'm/d', per: 100, offset: 0 },
  { name: 'thermoclineTop', unit: 'm', per: 5, offset: 0 }, { name: 'surfaceSalinity', unit: 'psu', per: 1000, offset: 35 }, { name: 'iceVolume', unit: 'm', per: 100, offset: 0 },
  { name: 'iceCover', unit: '1', per: 10000, offset: 0 }, { name: 'netSurfaceHeat', unit: 'W/m²', per: 10, offset: 0 }, { name: 'rain', unit: 'mm/d', per: 100, offset: 0 },
  { name: 'cloudWater', unit: 'g/m²', per: 4, offset: 0 },
];
const F = FIELDS.length, tally = () => ({ nonFinite: Object.fromEntries(FIELDS.map(({ name }) => [name, 0])), clipped: Object.fromEntries(FIELDS.map(({ name }) => [name, 0])) });
const pack = (values, field, counts) => Int16Array.from(values, (v) => {
  if (!Number.isFinite(v)) { counts.nonFinite[field.name]++; return 0; }
  const raw = Math.round((v - field.offset) * field.per);
  if (Math.abs(raw) > 32767) { counts.clipped[field.name]++; return Math.sign(raw) * 32767; }
  return raw;
});
const recent = resumed && resumed.recentCount ? Array.from({ length: resumed.recentCount }, (_, n) => ({ rain: resumed.recentRain.subarray(n * C, (n + 1) * C), net: resumed.recentNet.subarray(n * C, (n + 1) * C) })) : [];
function intervalMean(days, name) {
  const used = recent.slice(-days), mean = new Float64Array(C);
  for (const day of used) for (let i = 0; i < C; i++) mean[i] += day[name][i] / used.length;
  return mean;
}
function frameValues(snap, cloud, sequence) {
  const days = sequence.every / DAY, temperature = (k, i) => snap.T[k * C + i], surface = depthFields(mesh, L, { h: snap.h, u: snap.u, temperature, cellOcean }, 0), deep = depthFields(mesh, L, { h: snap.h, u: snap.u, temperature, cellOcean }, UPWELLING_DEPTH);
  const net = intervalMean(days, 'net'), rain = intervalMean(days, 'rain');
  const values = [
    Float64Array.from({ length: C }, (_, i) => snap.T[i] - CELSIUS),
    Float64Array.from({ length: C }, (_, i) => project(surface.current, i, eastUnit)), Float64Array.from({ length: C }, (_, i) => project(surface.current, i, northUnit)),
    snap.h1, snap.eta, Float64Array.from(deep.upwelling, (w) => w * DAY), snap.thermoclineDepth, snap.S.subarray(0, C),
    Float64Array.from({ length: C }, (_, i) => coverOf(ice[i], seaIce.concentration[i]) * ice[i]), Float64Array.from({ length: C }, (_, i) => coverOf(ice[i], seaIce.concentration[i])),
    net, Float64Array.from(rain, (x) => x * DAY), Float64Array.from(cloud, (x) => 1000 * x),
  ];
  return values.map((v, f) => pack(v, FIELDS[f], sequence.counts));
}
const year = { name: 'year', every: FRAME_DAYS * DAY, times: [], frames: [], counts: tally() };
const spinup = { name: 'spinup', every: DAY, times: resumed ? Array.from(resumed.spinupFrameTimes) : [], frames: [], counts: tally() };
if (resumed) for (let n = 0; n < spinup.times.length; n++) spinup.frames.push(Array.from({ length: F }, (_, f) => Int16Array.from(resumed.spinupFrameData.subarray((n * F + f) * C, (n * F + f + 1) * C))));
if (resumed && resumed.countsJson) { const saved = JSON.parse(resumed.countsJson), counts = saved.spinup ?? saved; for (const { name } of FIELDS) { spinup.counts.nonFinite[name] += counts.nonFinite[name] ?? 0; spinup.counts.clipped[name] += counts.clipped[name] ?? 0; } }

const zonalSums = () => ({ ...Object.fromEntries(['sst', 'sss', 'mixedLayerDepth', 'seaSurfaceHeight', 'uSurface', 'vSurface', 'uDepthMean', 'vDepthMean', 'zonalTransport', 'upwelling', 'thermoclineTop', 'iceCover', 'iceVolume'].map((name) => [name, new Float64Array(BINS)])), samples: 0 });
const sectionSums = () => ({ ...Object.fromEntries(['h', 'hT', 'hS', 'hU', 'hV'].map((name) => [name, new Float64Array(BINS * L)])), samples: 0 });
const overturningSums = () => ({ ...Object.fromEntries(['mass', 'heat', 'northMass', 'northHeat', 'southMass', 'southHeat', 'eddyMass', 'eddyHeat', 'diffusionHeat'].map((name) => [name, new Float64Array(CAPS.length * L)])), net: new Float64Array(CAPS.length), samples: 0, eddySamples: 0 });
const zonalOcean = MONTHS.map(zonalSums), sections = MONTHS.map(sectionSums), overturning = MONTHS.map(overturningSums);
const layerEast = new Float64Array(L * C), layerNorth = new Float64Array(L * C), vec = new Float64Array(3 * C);
function accumulateOcean(snap, m) {
  const { h, u: velocity, T, S, eta, eddy } = snap, z = zonalOcean[m], sec = sections[m], ov = overturning[m];
  for (let k = 0; k < L; k++) {
    cellVector(mesh, velocity.subarray(k * E, (k + 1) * E), vec);
    for (let i = 0; i < C; i++) { layerEast[k * C + i] = project(vec, i, eastUnit); layerNorth[k * C + i] = project(vec, i, northUnit); }
  }
  const up = depthFields(mesh, L, { h, u: velocity, temperature: (k, i) => T[k * C + i], cellOcean }, UPWELLING_DEPTH).upwelling;
  for (let i = 0; i < C; i++) {
    const b = bin[i], w = weight[i];
    let depth = 0, transportEast = 0, transportNorth = 0;
    for (let k = 0; k < L; k++) {
      const n = k * C + i, hk = h[n], at = b * L + k;
      depth += hk; transportEast += hk * layerEast[n]; transportNorth += hk * layerNorth[n];
      sec.h[at] += w * hk; sec.hT[at] += w * hk * (T[n] - CELSIUS); sec.hS[at] += w * hk * S[n]; sec.hU[at] += w * hk * layerEast[n]; sec.hV[at] += w * hk * layerNorth[n];
    }
    const cover = coverOf(snap.ice[i], snap.concentration[i]);
    z.sst[b] += w * (T[i] - CELSIUS); z.sss[b] += w * S[i]; z.mixedLayerDepth[b] += w * snap.h1[i]; z.seaSurfaceHeight[b] += w * eta[i];
    z.uSurface[b] += w * layerEast[i]; z.vSurface[b] += w * layerNorth[i]; z.uDepthMean[b] += w * transportEast / depth; z.vDepthMean[b] += w * transportNorth / depth;
    z.zonalTransport[b] += w * transportEast; z.upwelling[b] += w * up[i] * DAY; z.thermoclineTop[b] += w * snap.thermoclineDepth[i];
    z.iceCover[b] += w * cover; z.iceVolume[b] += w * cover * snap.ice[i];
  }
  edgeFluxes.load({ h, u: velocity, T, eta });
  CAPS.forEach((_, c) => {
    for (const [e, sign] of capEdges[c]) {
      const a = cellsOnEdge[2 * e], bb = cellsOnEdge[2 * e + 1], out = edgeFluxes.edge(e, edgeOut);
      for (let k = 0; k < L; k++) {
        const at = c * L + k, flux = sign * out.volume[k], carried = rhoCp * sign * (out.heat[k] - CELSIUS * out.volume[k]);
        ov.net[c] += sign * out.net[k];
        ov.mass[at] += flux; ov.heat[at] += carried; ov.diffusionHeat[at] += rhoCp * sign * out.diffusion[k];
        if (flux > 0) { ov.northMass[at] += flux; ov.northHeat[at] += carried; } else { ov.southMass[at] += flux; ov.southHeat[at] += carried; }
        const f = eddy ? eddy[k * E + e] : 0;
        if (f !== 0) { const from = f > 0 ? a : bb; ov.eddyMass[at] += sign * f; ov.eddyHeat[at] += rhoCp * sign * f * (T[k * C + from] - CELSIUS); }
      }
    }
  });
  z.samples++; sec.samples++; ov.samples++;
  if (eddy) ov.eddySamples++;
}

const BUDGET = ['absorbedSolar', 'outgoingLongwave', 'clearAbsorbedSolar', 'clearOutgoingLongwave', 'insolation', 'reflectedSolar', 'netSurfaceHeat', 'snowFusion', 'surfaceShortwave', 'sensibleHeat', 'latentHeat', 'evaporation', 'precipitation', 'convectiveRain', 'largeScaleRain', 'snowfall', 'surfaceTemperature', 'seaSurfaceTemperature', 'iceCover'];
const budgetSums = () => ({ ...Object.fromEntries(BUDGET.map((name) => [name, new Float64Array(BINS)])), days: 0 });
const budgets = MONTHS.map(budgetSums), cellDay = Object.fromEntries(BUDGET.map((name) => [name, new Float64Array(C)]));
function fillCellDay(forcing, sunlight, reflected) {
  const r = model.radiation, m = model.moist;
  for (let i = 0; i < C; i++) {
    cellDay.absorbedSolar[i] = r.meanAbsorbedSolar[i]; cellDay.outgoingLongwave[i] = r.meanOutgoingLongwave[i];
    cellDay.insolation[i] = sunlight[i] / perDay; cellDay.reflectedSolar[i] = reflected[i] / perDay; cellDay.snowFusion[i] = FUSION * forcing.snowfall[i];
    cellDay.clearAbsorbedSolar[i] = r.meanAbsorbedSolar[i] - r.meanShortwaveCloudEffect[i]; cellDay.clearOutgoingLongwave[i] = r.meanOutgoingLongwave[i] + r.meanLongwaveCloudEffect[i];
    cellDay.netSurfaceHeat[i] = forcing.netFlux[i]; cellDay.surfaceShortwave[i] = forcing.shortwave[i]; cellDay.sensibleHeat[i] = forcing.sensible[i];
    cellDay.latentHeat[i] = LATENT * forcing.evaporation[i]; cellDay.evaporation[i] = forcing.evaporation[i] * DAY; cellDay.precipitation[i] = forcing.rain[i] * DAY;
    cellDay.convectiveRain[i] = m.convectiveRain[i]; cellDay.largeScaleRain[i] = m.largeScaleRain[i]; cellDay.snowfall[i] = forcing.snowfall[i] * DAY;
    cellDay.surfaceTemperature[i] = forcing.surfaceT[i]; cellDay.seaSurfaceTemperature[i] = forcing.sst[i]; cellDay.iceCover[i] = forcing.concentration[i];
  }
}

const boundaryDays = new Set([0, RECORD_DAYS]);
for (let d = 1; d < RECORD_DAYS; d++) if (monthOf(d) !== monthOf(d - 1)) boundaryDays.add(d);
const books = { boundaries: [], periods: [] };
let period = null;
const bandStores = () => Object.fromEntries(STORES.map((name) => [name, bands((i) => store[name][i])]));

const DAILY = ['surfaceTemperature', 'absorbedSolar', 'outgoingLongwave', 'clearAbsorbedSolar', 'clearOutgoingLongwave', 'planetaryAlbedo', 'precipitation', 'convectiveRain', 'largeScaleRain', 'evaporation', 'snowfall', 'sensibleHeat', 'precipitableWater', 'cloudCondensate', 'cloudWater', 'ice', 'maxWind', 'waterResidual',
  'seaSurfaceTemperature', 'mixedLayerDepth', 'interiorTemperature', 'surfaceSalinity', 'meanSalinity', 'fastestCurrent', 'fastestCurrentLayer', 'fastestMixedLayerCurrent', 'largestTransport', 'clampedEdges', 'largestSeaSurfaceHeight', 'meanSeaLevel',
  'iceNorth', 'iceSouth', 'iceVolumeNorth', 'iceVolumeSouth', 'toaNet', 'surfaceNet', 'snowFusion', 'oceanHeatChange', 'iceEnergyChange', 'airEnergyChange', 'energyResidual', 'surfaceResidual', 'atmosphereResidual', 'nonFinite'];
const daily = Object.fromEntries(DAILY.map((name) => [name, resumed ? Array.from(resumed[`daily_${name}`] ?? resumed.daily_toaNet.map(() => NaN), (v) => (Number.isFinite(v) ? v : null)) : []]));
const run = resumed ? JSON.parse(resumed.runJson) : { toa: 0, surface: 0, snow: 0, evaporated: 0, rained: 0, seconds: 0, start: null, water0: null, initial: null };

const initial = await model.beginFrame({ fields: ['cloud'], diagnostics: true });
await model.sync();
let snap = await oceanSnapshot(), [snow, oceanFlux, pending] = await readRanges(device, gpu.buffers.PH, [{ offset: PH.SNOW, length: C }, { offset: PH.OFLUX, length: C }, { offset: PH.PEN, length: C }]);
if (resumed && !(resumed.eddyV && resumed.eddyV.length === L * E)) snap.eddy = null;
storeCells(snap, snow, oceanFlux, pending);
let totals = totalsNow(), previous = totals, resumeCheck = null;
if (resumed) {
  const saved = resumed.totalsJson ? JSON.parse(resumed.totalsJson) : null;
  if (!saved && process.env.FORCE_RESUME !== '1') { console.error(`${process.env.RESUME} does not record its storage totals; FORCE_RESUME=1 resumes it anyway`); process.exit(1); }
  const relative = saved && Object.fromEntries(STORES.map((name) => [name, Math.abs(totals[name] - saved[name]) / Math.max(Math.abs(totals[name]), Math.abs(saved[name]), 1e-300)]));
  resumeCheck = {
    forced: !engineMatches || !saved, sameEngine: engineMatches, sameScript: scriptMatches,
    stores: relative && Object.fromEntries(STORES.map((name) => [name, +relative[name].toPrecision(3)])), tolerance: FLOAT32_TOLERANCE,
    oceanArrays: `${restartNames.join(', ')} read back equal to the saved float32 values`, eddyV: !!(resumed.eddyV && resumed.eddyV.length === L * E),
  };
  if (saved) {
    const off = STORES.filter((name) => !(relative[name] <= FLOAT32_TOLERANCE));
    console.error(`resumed stores against the checkpoint's (relative): ${STORES.map((name) => `${name} ${relative[name].toExponential(1)}`).join(', ')}`);
    if (off.length) { console.error(`${process.env.RESUME}: ${off.join(', ')} differ from the saved totals by more than float32 rounding (${FLOAT32_TOLERANCE.toExponential(1)} relative); stopping`); process.exit(1); }
  } else console.error(`${process.env.RESUME} records no storage totals; resuming without the store check (FORCE_RESUME=1), and the header's resumeCheck.stores is null`);
}
const totalWater = () => global((i) => store.vapor[i] + store.condensate[i]);
if (!resumed) {
  const d = initial.diagnostics;
  run.start = totals; run.water0 = totalWater();
  run.initial = { meanMixedLayer: d.oceanUpperDepth, interiorTemperature: d.oceanInteriorT - CELSIUS, largestSeaSurfaceHeight: d.oceanSSH, surfaceSalinity: d.oceanSalinity, meanSeaSurfaceTemperature: global((i) => snap.T[i] - CELSIUS), ice: 100 * d.iceFraction };
}
const recordTotals = { start: null, end: null, toa: 0, surface: 0, snow: 0, seconds: 0 };

async function atInstant(cloud, first) {
  const time = model.time, rt = time - recordStart, runDay = Math.round((time - runStartTime) / DAY);
  if (!(first && resumed) && runDay <= SPINUP_FRAMES && rt <= 0) { spinup.times.push(time); spinup.frames.push(frameValues(snap, cloud, spinup)); }
  if (rt >= 0 && rt <= end && Math.round(rt / DAY) % FRAME_DAYS === 0) { year.times.push(time); year.frames.push(frameValues(snap, cloud, year)); }
  const recordDay = Math.round(rt / DAY);
  if (RECORD_DAYS > 0 && rt >= 0 && boundaryDays.has(recordDay)) {
    const boundary = { time: rt, bands: bandStores(), global: totals, capVolume: capSouth.map((cells) => cells.reduce((sum, i) => sum + weight[i] * snap.eta[i], 0)) };
    if (period) { period.to = rt; books.periods.push(finishPeriod(period, books.boundaries[books.boundaries.length - 1], boundary)); }
    books.boundaries.push(boundary);
    period = recordDay < RECORD_DAYS ? { from: rt, days: 0, sums: Object.fromEntries(BUDGET.map((name) => [name, new Float64Array(BINS)])) } : null;
    if (recordDay === 0) recordTotals.start = totals;
    if (recordDay === RECORD_DAYS) recordTotals.end = totals;
  }
}
function finishPeriod(p, from, to) {
  const seconds = to.time - from.time, mean = (name) => Array.from(p.sums[name], (v, b) => v / (binArea[b] * p.days));
  const fluxes = Object.fromEntries(BUDGET.map((name) => [name, mean(name)]));
  const change = Object.fromEntries(STORES.map((name) => [name, from.bands[name].map((x, b) => (to.bands[name][b] - x) / seconds)]));
  const globalMean = (values) => values.reduce((sum, v, b) => sum + v * binArea[b], 0) / area;
  const toaNet = fluxes.absorbedSolar.map((x, b) => x - fluxes.outgoingLongwave[b]);
  const stored = (energyOf(to.global) - energyOf(from.global)) / seconds, surfaceStored = (surfaceOf(to.global) - surfaceOf(from.global)) / seconds;
  const round = (values, digits) => values.map((v) => +v.toFixed(digits));
  return {
    from: p.from, to: p.to, days: p.days,
    fluxes: { ...Object.fromEntries(BUDGET.map((name) => [name, round(fluxes[name], 4)])), toaNet: round(toaNet, 4) },
    storageChange: Object.fromEntries(STORES.map((name) => [name, round(change[name], 4)])),
    global: {
      toaNet: globalMean(toaNet), surfaceNet: globalMean(fluxes.netSurfaceHeat), snowFusion: globalMean(fluxes.snowFusion), stored, surfaceStored, residual: stored - globalMean(toaNet),
      surfaceResidual: surfaceStored - (globalMean(fluxes.netSurfaceHeat) - globalMean(fluxes.snowFusion)),
      atmosphereResidual: stored - surfaceStored - (globalMean(toaNet) - globalMean(fluxes.netSurfaceHeat) + globalMean(fluxes.snowFusion)),
    },
  };
}
await atInstant(resumed && resumed.cloud && resumed.cloud.length === C ? Float64Array.from(resumed.cloud) : initial.fields.cloud, true);

const checkpoints = [];
async function saveCheckpoint(path, cloud) {
  await model.sync();
  const saved = await model.ocean.serialize({ restart: true });
  const [seaSnow, seaSnowAlbedo] = await readRanges(device, gpu.buffers.PH, [{ offset: PH.SNOW, length: C }, { offset: PH.SNOWALB, length: C }]);
  const r = model.radiation, m = model.moist, b = model.boundaryLayer;
  const spinupFrameData = new Float32Array(spinup.frames.length * F * C);
  spinup.frames.forEach((frame, n) => frame.forEach((field, f) => spinupFrameData.set(field, (n * F + f) * C)));
  const kept = recent.slice(-FRAME_DAYS), recentRain = new Float32Array(kept.length * C), recentNet = new Float32Array(kept.length * C);
  kept.forEach((day, n) => { recentRain.set(day.rain, n * C); recentNet.set(day.net, n * C); });
  const dailyArrays = Object.fromEntries(DAILY.map((name) => [`daily_${name}`, Float64Array.from(daily[name], (v) => v ?? NaN)]));
  const bytes = encodeState({
    kind: 'explainOcean', N, ...gridHeader(mesh), K, day: Math.round((model.time - runStartTime) / DAY), time: model.time, terrain: false, levels: core.levels,
    runStartTime, optionsJson: JSON.stringify({ ocean: OCEAN, ice: ICE }), segmentsJson: JSON.stringify(segments.map((s, n) => (n === segments.length - 1 ? { ...s, toTime: model.time } : s))), runJson: JSON.stringify(run), recentCount: kept.length,
    engineJson: JSON.stringify(engine), totalsJson: JSON.stringify(totals), countsJson: JSON.stringify({ spinup: spinup.counts }), eddyV: snap.eddy ? Float32Array.from(snap.eddy) : new Float32Array(0), cloud: Float32Array.from(cloud),
    pi, theta, u, surfaceT, q, qc, ice, concentration: seaIce.concentration, mlmSubsidence: r.mlmSubsidence, mlmHeight: r.mlmHeight, mlmGate: r.mlmGate,
    convectiveRain: m.convectiveRain, largeScaleRain: m.largeScaleRain, meanAbsorbedSolar: r.meanAbsorbedSolar, meanOutgoingLongwave: r.meanOutgoingLongwave, meanPlanetaryAlbedo: r.meanPlanetaryAlbedo,
    meanShortwaveCloudEffect: r.meanShortwaveCloudEffect, meanLongwaveCloudEffect: r.meanLongwaveCloudEffect, boundaryDepth: b.depth, mixingTop: b.mixingTop, boundaryRegime: b.regime, boundaryBuoyancy: b.buoyancyFlux,
    windSpeed: model.surface.windSpeed, evaporation: r.evaporation, ...(model.exchange ? { ...(model.exchange.fixed ? {} : { exchangeHeat: model.exchange.heat }), exchangeWind: model.exchange.wind } : {}),
    cumulusCover: m.cumulusCover, cumulusWater: m.cumulusWater, subcloudVirtual: m.subcloudVirtual, ...(model.skin ? model.skin.serialize() : {}), ocean: saved, seaSnow, seaSnowAlbedo,
    spinupFrameTimes: Float64Array.from(spinup.times), spinupFrameData, recentRain, recentNet, ...dailyArrays,
  }, { f64: ['spinupFrameTimes', ...Object.keys(dailyArrays)] });
  mkdirSync(dirname(path), { recursive: true });
  writeFileSync(`${path}.partial`, bytes);
  renameSync(`${path}.partial`, path);
  console.error(`checkpoint ${path}: model day ${Math.round(model.time / DAY)}, ${(bytes.length / 1048576).toFixed(1)} MB`);
  const runDay = Math.round((model.time - runStartTime) / DAY), previousWritten = checkpoints.at(-1);
  if (previousWritten && previousWritten.path === path) previousWritten.kept = false;
  else if (!keepCheckpoints && previousWritten && existsSync(previousWritten.path)) { unlinkSync(previousWritten.path); previousWritten.kept = false; }
  checkpoints.push({ path, runDay, time: model.time, kept: true });
}

const totalDays = SPINUP_DAYS + RECORD_DAYS, wallStart = performance.now();
let lastDay = performance.now();
for (let day = 1; day <= totalDays; day++) {
  for (let part = 0; part < 4; part++) {
    const rt = model.time - recordStart;
    if (rt >= 0 && rt < end) accumulateOcean(snap, monthOf(rt / DAY));
    for (let n = 1; n <= quarter; n++) { await model.step(STEP); recorder.step(); if (n % 8 === 0) await model.settle(); }
    await model.settle();
    const next = model.time - recordStart;
    if (part === 3 || (next >= 0 && next < end)) snap = await oceanSnapshot();
  }
  const runDay = Math.round((model.time - runStartTime) / DAY);
  const [sunlight, reflected] = await readRanges(device, gpu.buffers.PH, [{ offset: PH.INSSUM, length: C }, { offset: PH.REFLSUM, length: C }]);
  const frame = await model.beginFrame({ fields: ['cloud'], diagnostics: true }), d = frame.diagnostics, cloud = Float64Array.from(frame.fields.cloud);
  const forcing = (await decodeForcing(await recorder.day(runDay))).fields;
  await model.sync();
  [snow, oceanFlux, pending] = await readRanges(device, gpu.buffers.PH, [{ offset: PH.SNOW, length: C }, { offset: PH.OFLUX, length: C }, { offset: PH.PEN, length: C }]);
  storeCells(snap, snow, oceanFlux, pending);
  previous = totals; totals = totalsNow();
  recent.push({ rain: Float32Array.from(forcing.rain), net: Float32Array.from(forcing.netFlux) });
  if (recent.length > FRAME_DAYS) recent.shift();
  fillCellDay(forcing, sunlight, reflected);
  const dayStart = model.time - DAY - recordStart;
  if (dayStart >= 0 && dayStart < end) {
    const sums = budgets[monthOf(dayStart / DAY)];
    for (let i = 0; i < C; i++) for (const name of BUDGET) { sums[name][bin[i]] += weight[i] * cellDay[name][i]; period.sums[name][bin[i]] += weight[i] * cellDay[name][i]; }
    sums.days++; period.days++;
  }

  let bad = 0;
  for (const array of [...state, snap.h, snap.u, snap.T, snap.S, snap.eta]) for (let n = 0; n < array.length; n++) if (!Number.isFinite(array[n])) bad++;
  let fastest = 0, fastestAt = 0;
  for (let n = 0; n < L * E; n++) if (Math.abs(snap.u[n]) > fastest) { fastest = Math.abs(snap.u[n]); fastestAt = n; }
  let volume = 0, salt = 0;
  for (let i = 0; i < C; i++) for (let k = 0; k < L; k++) { volume += weight[i] * snap.h[k * C + i]; salt += weight[i] * snap.h[k * C + i] * snap.S[k * C + i]; }
  const toaNet = d.absorbedSolar - d.outgoingLongwave, surfaceNet = global((i) => forcing.netFlux[i]), snowfall = global((i) => forcing.snowfall[i]);
  run.toa += toaNet * DAY; run.surface += surfaceNet * DAY; run.snow += FUSION * snowfall * DAY; run.seconds += DAY;
  run.evaporated += global((i) => forcing.evaporation[i]) * DAY; run.rained += global((i) => forcing.rain[i]) * DAY;
  if (dayStart >= 0 && dayStart < end) { recordTotals.toa += toaNet * DAY; recordTotals.surface += surfaceNet * DAY; recordTotals.snow += FUSION * snowfall * DAY; recordTotals.seconds += DAY; }
  const means = Object.fromEntries(BUDGET.map((name) => [name, global((i) => cellDay[name][i])]));
  const coverMean = (north) => 100 * hemisphere((i) => coverOf(ice[i], seaIce.concentration[i]), north), volumeMean = (north) => hemisphere((i) => store.iceVolume[i], north);
  const values = {
    surfaceTemperature: means.surfaceTemperature, absorbedSolar: d.absorbedSolar, outgoingLongwave: d.outgoingLongwave, clearAbsorbedSolar: d.clearAbsorbedSolar, clearOutgoingLongwave: d.clearOutgoingLongwave, planetaryAlbedo: d.planetaryAlbedo,
    precipitation: means.precipitation, convectiveRain: means.convectiveRain, largeScaleRain: means.largeScaleRain, evaporation: means.evaporation, snowfall: means.snowfall, sensibleHeat: means.sensibleHeat,
    precipitableWater: totals.vapor, cloudCondensate: 1000 * totals.condensate, cloudWater: 1000 * global((i) => cloud[i]), ice: 100 * means.iceCover, maxWind: d.maxWind,
    waterResidual: totalWater() - run.water0 - run.evaporated + run.rained,
    seaSurfaceTemperature: global((i) => snap.T[i] - CELSIUS), mixedLayerDepth: d.oceanUpperDepth, interiorTemperature: d.oceanInteriorT - CELSIUS, surfaceSalinity: d.oceanSalinity, meanSalinity: salt / volume,
    fastestCurrent: fastest, fastestCurrentLayer: Math.floor(fastestAt / E), fastestMixedLayerCurrent: d.oceanSpeed, largestTransport: d.oceanTransport, clampedEdges: d.oceanLimited, largestSeaSurfaceHeight: d.oceanSSH, meanSeaLevel: volume / area - meanFloor,
    iceNorth: coverMean(true), iceSouth: coverMean(false), iceVolumeNorth: volumeMean(true), iceVolumeSouth: volumeMean(false),
    toaNet, surfaceNet, snowFusion: FUSION * snowfall, oceanHeatChange: (totals.ocean - previous.ocean) / DAY, iceEnergyChange: (totals.ice - previous.ice) / DAY, airEnergyChange: (totals.airEnthalpy + totals.latent - previous.airEnthalpy - previous.latent) / DAY,
    energyResidual: (energyOf(totals) - energyOf(run.start) - run.toa) / run.seconds, surfaceResidual: (surfaceOf(totals) - surfaceOf(run.start) - (run.surface - run.snow)) / run.seconds,
    atmosphereResidual: (energyOf(totals) - surfaceOf(totals) - energyOf(run.start) + surfaceOf(run.start) - (run.toa - run.surface + run.snow)) / run.seconds, nonFinite: bad,
  };
  for (const name of DAILY) daily[name].push(Number.isFinite(values[name]) ? +values[name].toPrecision(7) : null);
  const pace = (performance.now() - lastDay) / 1000;
  lastDay = performance.now();
  console.error(`ocean day ${runDay}${dayStart < 0 ? ' (spin-up)' : ''}: Ts ${values.surfaceTemperature.toFixed(2)} K, SST ${values.seaSurfaceTemperature.toFixed(2)} °C, ASR ${d.absorbedSolar.toFixed(1)} OLR ${d.outgoingLongwave.toFixed(1)} surface ${surfaceNet.toFixed(1)} W/m², albedo ${d.planetaryAlbedo.toFixed(3)}, P ${means.precipitation.toFixed(2)} E ${means.evaporation.toFixed(2)} mm/d, ice N ${values.iceNorth.toFixed(1)} S ${values.iceSouth.toFixed(1)} %, ML ${d.oceanUpperDepth.toFixed(1)} m, interior ${values.interiorTemperature.toFixed(3)} °C, fastest ${fastest.toFixed(2)} m/s (layer ${values.fastestCurrentLayer}; ML ${d.oceanSpeed.toFixed(2)}), transport ${d.oceanTransport.toFixed(0)} Sv, |η| ${d.oceanSSH.toFixed(2)} m, clamped ${d.oceanLimited}, sea level ${(1000 * values.meanSeaLevel).toFixed(3)} mm, ocean ${values.oceanHeatChange.toFixed(1)} W/m², residual ${values.energyResidual.toFixed(2)} (surface ${values.surfaceResidual.toFixed(2)}, air ${values.atmosphereResidual.toFixed(2)}, snow ${values.snowFusion.toFixed(2)}) W/m², ${pace.toFixed(2)} s/day, ${seconds().toFixed(0)} s`);
  if (bad || ![values.surfaceTemperature, d.absorbedSolar, d.outgoingLongwave, d.maxWind, fastest, d.oceanSSH].every(Number.isFinite) || d.maxWind > 400) { console.error(`unstable: ${bad} non-finite values; stopping`); process.exit(1); }
  await atInstant(cloud, false);
  if (model.time <= recordStart + 1 && (Math.round(model.time / DAY) % 365 === 0 || testCheckpointDays.has(runDay))) await saveCheckpoint(checkpointPath(runDay), cloud);
}
const wall = (performance.now() - wallStart) / 1000;

const sumMonths = (list, make) => {
  const total = make();
  for (const sums of list) for (const [name, value] of Object.entries(sums)) { if (typeof value === 'number') total[name] += value; else value.forEach((v, m) => { total[name][m] += v; }); }
  return total;
};
const ZONAL_DIGITS = { sst: 3, sss: 4, mixedLayerDepth: 2, seaSurfaceHeight: 4, uSurface: 5, vSurface: 5, uDepthMean: 5, vDepthMean: 5, upwelling: 4, thermoclineTop: 2, iceCover: 5, iceVolume: 4 };
const finishZonal = (z) => {
  if (!z.samples) return { samples: 0 };
  const mean = (name) => Array.from(z[name], (v, b) => v / (binArea[b] * z.samples));
  return {
    samples: z.samples,
    ...Object.fromEntries(Object.entries(ZONAL_DIGITS).map(([name, digits]) => [name, mean(name).map((v) => +v.toFixed(digits))])),
    zonalTransport: mean('zonalTransport').map((v) => +(v * bandWidth / 1e6).toFixed(2)),
    iceThickness: Array.from(z.iceVolume, (v, b) => (z.iceCover[b] > 0 ? +(v / z.iceCover[b]).toFixed(4) : null)),
  };
};
const finishSection = (s, full) => {
  if (!s.samples) return { samples: 0 };
  const base = new Array(BINS * L);
  for (let b = 0; b < BINS; b++) { let depth = 0; for (let k = 0; k < L; k++) { depth += s.h[b * L + k] / (binArea[b] * s.samples); base[b * L + k] = +depth.toFixed(1); } }
  const weighted = (name, scale, digits) => Array.from(s[name], (v, at) => (s.h[at] / (binArea[Math.floor(at / L)] * s.samples) >= SECTION_MIN_THICKNESS ? +(scale * v / s.h[at]).toFixed(digits) : null));
  return { samples: s.samples, base, T: weighted('hT', 1, 3), S: weighted('hS', 1, 4), u: weighted('hU', 100, 3), ...(full ? { v: weighted('hV', 100, 3) } : {}) };
};
const etaNetOver = (month) => {
  const change = new Float64Array(CAPS.length);
  let seconds = 0;
  for (let n = 1; n < books.boundaries.length; n++) {
    const a = books.boundaries[n - 1], b = books.boundaries[n];
    if (month !== null && monthOf(a.time / DAY) !== month) continue;
    for (let c = 0; c < CAPS.length; c++) change[c] += b.capVolume[c] - a.capVolume[c];
    seconds += b.time - a.time;
  }
  return seconds > 0 ? Array.from(change, (v) => -v / seconds) : null;
};
const finishOverturning = (o, etaNet) => {
  if (!o.samples) return { samples: 0 };
  const mean = (name, scale, count) => Array.from(o[name], (v) => (count ? +(v * scale / count).toPrecision(4) : null));
  const imbalance = etaNet && CAPS.map((_, c) => {
    let sum = 0;
    for (let k = 0; k < L; k++) sum += o.mass[c * L + k] / o.samples + (o.eddySamples ? o.eddyMass[c * L + k] / o.eddySamples : 0);
    return +((sum - etaNet[c]) / 1e6).toPrecision(4);
  });
  return {
    samples: o.samples, eddySamples: o.eddySamples, mass: mean('mass', 1e-6, o.samples), heat: mean('heat', 1e-15, o.samples), eddyMass: mean('eddyMass', 1e-6, o.eddySamples), eddyHeat: mean('eddyHeat', 1e-15, o.eddySamples), diffusionHeat: mean('diffusionHeat', 1e-15, o.samples), net: mean('net', 1e-6, o.samples),
    etaNet: etaNet ? etaNet.map((v) => +(v / 1e6).toPrecision(4)) : null, imbalance,
  };
};
const finishBudget = (sums) => {
  if (!sums.days) return { days: 0 };
  return {
    days: sums.days, ...Object.fromEntries(BUDGET.map((name) => [name, Array.from(sums[name], (v, b) => +(v / (binArea[b] * sums.days)).toFixed(4))])),
    planetaryAlbedo: Array.from(sums.reflectedSolar, (v, b) => (sums.insolation[b] > 0 ? +(v / sums.insolation[b]).toFixed(4) : null)),
  };
};
const recorded = (list, finish) => MONTHS.map(([name], m) => ({ name, days: monthDays(m), ranges: monthDays(m), ...finish(list[m], false) }));
const returnAcross = (o, degrees) => {
  const c = CAPS.indexOf(degrees);
  if (!o.samples || c < 0) return null;
  const sv = (v) => +(v / o.samples / 1e6).toFixed(2), celsius = (heat, mass) => (mass !== 0 ? +(heat / mass / rhoCp).toFixed(2) : null);
  const groups = GROUPS.map((name, g) => {
    let mass = 0, heat = 0, eddy = 0, north = 0, northHeat = 0, south = 0, southHeat = 0;
    for (let k = 0, at = c * L; k < L; k++, at++) if (groupOf(k) === g) { mass += o.mass[at]; heat += o.heat[at]; eddy += o.eddyMass[at]; north += o.northMass[at]; northHeat += o.northHeat[at]; south += o.southMass[at]; southHeat += o.southHeat[at]; }
    return {
      name, mass: sv(mass), temperature: Math.abs(mass) >= NET_FRACTION * (north - south) ? celsius(heat, mass) : null,
      northward: { mass: sv(north), temperature: celsius(northHeat, north) }, southward: { mass: sv(south), temperature: celsius(southHeat, south) },
      eddyMass: o.eddySamples ? +(eddy / o.eddySamples / 1e6).toFixed(2) : null,
    };
  });
  let heat = 0, eddyHeat = 0;
  for (let k = 0; k < L; k++) { heat += o.heat[c * L + k] + o.diffusionHeat[c * L + k]; eddyHeat += o.eddyHeat[c * L + k]; }
  return { latitude: degrees, groups, heat: +((heat / o.samples + (o.eddySamples ? eddyHeat / o.eddySamples : 0)) / 1e15).toFixed(3), net: +(o.net[c] / o.samples / 1e6).toFixed(2) };
};
const annualOverturning = sumMonths(overturning, overturningSums);

const residuals = (start, stop, integral, seconds) => (start && stop && seconds > 0 ? {
  seconds, toaNet: integral.toa / seconds, surfaceNet: integral.surface / seconds, snowFusion: integral.snow / seconds,
  storageChange: Object.fromEntries(STORES.filter((name) => !['vapor', 'condensate', 'salt', 'iceVolume'].includes(name)).map((name) => [name, (stop[name] - start[name]) / seconds])),
  residual: (energyOf(stop) - energyOf(start) - integral.toa) / seconds,
  residualWithKinetic: (energyOf(stop) + stop.airKinetic - energyOf(start) - start.airKinetic - integral.toa) / seconds,
  surfaceResidual: (surfaceOf(stop) - surfaceOf(start) - (integral.surface - integral.snow)) / seconds,
  atmosphereResidual: (stop.airEnthalpy + stop.latent - start.airEnthalpy - start.latent - (integral.toa - integral.surface + integral.snow)) / seconds,
} : null);
const adapter = (await getDevice()).adapter, info = adapter.info ?? {};
const probe = (h, snowCover, T) => +seaIce.albedo(h, null, snowCover, h > 0 ? 1 : 0, T).toFixed(3);
const { geography, climatology, buffers, salinityProfile, layerTable, ...oceanScalars } = ocean.options;
const endTime = model.time;
const relativeSegments = segments.map(({ fromTime, toTime, ...rest }) => ({ from: fromTime - recordStart, to: (toTime ?? endTime) - recordStart, ...rest }));
const bins = { bins: BINS, binArea: Array.from(binArea, (a) => +a.toPrecision(9)), binCells: Array.from(binCells), partition: PARTITION };
const counted = (kind) => Object.fromEntries(FIELDS.map(({ name }) => [name, year.counts[kind][name] + spinup.counts[kind][name]]));
const header = {
  version: 1, case: 'ocean', engine: 'gpu', device: { vendor: info.vendor ?? null, architecture: info.architecture ?? null, description: info.description ?? null, runtime: 'Dawn through the webgpu package in Node' },
  N, K, layers: L, dt: STEP, levels: Array.from(core.levels), grid: GRID,
  fields: FIELDS.map(({ name, unit, per, offset }) => ({ name, unit, scale: 1 / per, offset })),
  sequences: [year, spinup].map((s) => ({ name: s.name, every: s.every, times: s.times.map((t) => t - recordStart) })),
  time: 'seconds from the March equinox (March 20) that begins the recorded year; the model clock reads recordStart + time, and its day 0 and every 365th day are March equinoxes',
  intervals: `netSurfaceHeat and rain in a frame are means over the days since the sequence\'s previous frame (5 for year, 1 for spinup; year\'s first frame takes the spin-up\'s last 5 days), 0 in the run\'s first frame; netSurfaceHeat is the surface flux into the sea and ice and leaves out L_f × snowfall, the heat they give up melting the snow that lands on them (frames hold no snowfall; zonalBudget\'s snowFusion has it by band), so where snow falls it is more than the sea and ice take; every other field is the state at the frame\'s time; cloudWater adds the cumulus and the deck to the resolved condensate and clips at 8191.75 g/m² (clipped counts it); iceVolume is cover × thickness; seaSurfaceTemperature is the top layer\'s, at freezing under ice; mixedLayerDepth is the ocean\'s diagnosed mixed-layer depth (h1: where the density first exceeds the top layer\'s by 0.03 kg/m³); thermoclineTop is the depth at which the density, held from decreasing downward, reaches 1024.0 kg/m³, interpolated between layer centers; upwelling50 is the vertical velocity at 50 m, upward positive`,
  spinupDays: Math.round((recordStart - runStartTime) / DAY), recordDays: RECORD_DAYS, startTime: runStartTime, recordStart,
  options: { ...options, levels: GRID, ocean: { ...options.ocean, layers: L, densities: Array.from(ocean.table.densities) } }, cadence: { radiationEvery: options.radiation.radiationEvery, dragEvery: options.dragEvery, oceanEvery: everySteps },
  segments: relativeSegments, resumedFrom: shown(process.env.RESUME), resumeCheck,
  checkpoint: shown(checkpoints.filter((c) => c.kept).at(-1)?.path), checkpoints: checkpoints.map(({ path, time, ...c }) => ({ path: shown(path), ...c, time: time - recordStart })),
  oceanOptions: Object.fromEntries(Object.entries(oceanScalars).filter(([, v]) => typeof v !== 'function')),
  floor: { kind: 'flat', depth: ocean.options.flatDepth, minimum: Math.min(...D), maximum: Math.max(...D), land: 'none: every cell and edge is sea' },
  initial: 'initializeState(model): open sea at max(T_f + 0.3 K, 30 cos²φ − 2 °C − 4 K·clamp((|φ| − 50°)/20°, 0, 1)), 28 °C at the equator with the 4 K cooling ramping in from 50° to 70°, and at the freezing point T_f under the ice: 1.5 m of ice poleward of 72°N and 0.7 m poleward of 68°S, RH 0.7σ², winds geostrophic outside ±15°; the ocean\'s analytic start on it: every layer at its minimum thickness (the z grid of the hybrid table, the deepest layer to the floor), T = T_d + (T₀ − T_d) e^(−z/λ) averaged over each layer with T_d 2 °C, T₀ that SST and λ = 400 m × (0.7 + 0.6 cos²φ), S the same blend from 34 + 2 exp(−((|φ| − 25°)/20°)²) to 34.7 psu, the free surface level, no currents; the first regrid moves the interfaces below the twelve fixed-depth layers to their isopycnal positions',
  start: run.initial,
  ice: {
    options: ICE, latentHeatFusion: FUSION, freezingPoint: FREEZING_POINT, meltingPoint: MELTING_POINT, note: 'slabHeatCapacity only sets the surface capacity before the ocean writes its own, ρc_p·max(h₀, 1 m), at load',
    diffuseAlbedo: { water: probe(0, 0, 0), coldIce: probe(2, 0, 250), meltingIce: probe(2, 0, MELTING_POINT), freshSnowOnIce: probe(2, 100, 250) },
  },
  constants: { cp: airCp, g: gravity, latentHeat: LATENT, rhoCp, oceanDensity: ocean.options.density, oceanSpecificHeat: ocean.options.specificHeat, radius },
  commit, modified, source: engine, date: new Date().toISOString(), startedAt,
  pace: { wallSeconds: +wall.toFixed(1), modelDays: totalDays, secondsPerDay: +(wall / totalDays).toFixed(3), setupSeconds: +((wallStart - t0) / 1000).toFixed(1) },
  labels: LABELS, fixedLayers: FIXED, minimumThickness: Array.from(ocean.table.minimum),
  zonalMean: { ...bins, from: 0, to: end, units: 'area means of 6-hourly samples: sst °C (the top layer\'s), sss psu, mixedLayerDepth m (diagnosed, h1), seaSurfaceHeight m, uSurface/vSurface the top layer\'s east/north current m/s, uDepthMean/vDepthMean Σh·u/Σh m/s, zonalTransport Sv (the band\'s mean Σh·u_east times its 5° width), upwelling m/d at 50 m, thermoclineTop m (where the density reaches 1024.0 kg/m³), iceCover fraction, iceVolume m (cover × thickness), iceThickness m over the iced part', ...finishZonal(sumMonths(zonalOcean, zonalSums)), months: recorded(zonalOcean, finishZonal) },
  section: { ...bins, labels: LABELS, index: 'bin * layers + layer, bins south to north, layer 0 the top fixed-depth layer', units: `base m (depth of the layer's base below the surface, the cumulative mean thickness), T °C, S psu, u and v cm/s east and north (thickness-weighted means; null where the layer's mean thickness in the bin is under ${SECTION_MIN_THICKNESS} m)`, ...finishSection(sumMonths(sections, sectionSums), true), months: MONTHS.map(([name], m) => ({ name, days: monthDays(m), ranges: monthDays(m), ...finishSection(sections[m], false) })) },
  overturning: {
    caps: CAPS, labels: LABELS, index: 'cap * layers + layer, caps south to north', partition: PARTITION,
    units: `northward Sv and PW across the latitude circle, means of 6-hourly snapshots (samples): mass and heat are the resolved flow, each layer\'s edge thickness by the model\'s rule (layer 0 the mean of its two cells, every other layer the smaller, the column clipped from the top at the sill depth plus mean η, as the ocean\'s edgeColumn, and every layer donor-limited as its thickness flux, scripts/oceanEdgeFlux.mjs), heat at the face value of the ocean\'s tracer scheme in °C, so relative to 0 °C; diffusionHeat is the lateral diffusion of heat across the circle (the background κ in every layer and the shear-scaled κ in the fixed-depth layers), which moves no water; eddyMass and eddyHeat are the Gent–McWilliams transport the last ocean step applied (the engine\'s EDDYV scratch), averaged over eddySamples; net is the same transport without the donor limit; etaNet (one per cap) is the net northward volume transport the change of η south of the cap implies between the period\'s boundaries, exact because freshwater enters the sea as a salt flux and never changes its volume; imbalance (one per cap) is Σ over layers of mass + eddyMass less etaNet, the flow the snapshots miss: the free surface follows the barotropic sub-step\'s averaged transport and the rescale to it moves every layer, and end-of-step snapshots alias the flow; so a summed heat transport depends on its reference, shifting by ρc_p·imbalance per kelvin; returnFlow sums the annual means by group of layers at 10°S, 10°N, 30°S and 30°N: mass is the resolved net northward Sv and temperature its heat over its volume in °C, defined only where the group\'s flow runs mostly one way (null where |mass| is under ${100 * NET_FRACTION} % of the group\'s northward plus southward transport); northward and southward split the resolved transport by the direction it runs at each edge and snapshot, each with its own flow-weighted °C (null where there is none); eddyMass is the eddy Sv (zero summed over all layers); heat is the PW of both, and net the cap\'s net`,
    ...finishOverturning(annualOverturning, etaNetOver(null)), months: MONTHS.map(([name], m) => ({ name, days: monthDays(m), ranges: monthDays(m), ...finishOverturning(overturning[m], etaNetOver(m)) })),
    returnFlow: [-10, 10, -30, 30].map((degrees) => returnAcross(annualOverturning, degrees)),
  },
  zonalBudget: { ...bins, from: 0, to: end, units: 'daily means: W/m² for the radiation (clear-sky pair from the cloud effects; insolation and reflectedSolar at the top of the atmosphere), netSurfaceHeat (into the surface, what the ocean and ice receive through the surface flux), snowFusion (L_f times the snowfall: snow that lands on the sea melts into it or onto its ice at that cost outside netSurfaceHeat, and the air gains it as the snow forms, so the sea and ice take netSurfaceHeat − snowFusion), surfaceShortwave, sensibleHeat (upward) and latentHeat (upward); mm/d for evaporation, precipitation, convectiveRain, largeScaleRain and snowfall (the forcing recorder\'s estimate: precipitation over sea where the lowest layer is below freezing after the step, which can differ from what the step froze); K for surfaceTemperature (skin) and seaSurfaceTemperature (freezing under ice); fraction for iceCover; planetaryAlbedo is reflectedSolar over insolation (null in polar night)', ...finishBudget(sumMonths(budgets, budgetSums)), months: recorded(budgets, finishBudget) },
  heatBooks: {
    ...bins,
    units: 'storage J/m² by band: ocean ρc_p·Σh(T − 273.15 K) over every layer, layer 0 at the surface temperature the ocean takes at its next step (freezing under ice) plus the ocean flux still owed to the ice; ice the sea-ice energy relative to freezing water (skin heat less the latent heat of ice and snow); airEnthalpy c_p·T over the column; airKinetic; latent L·vapor; vapor and condensate kg/m²; salt psu·m; iceVolume m. Period fluxes as in zonalBudget; storageChange per second (W/m² for the energy stores, kg/m²/s for vapor and condensate, psu·m/s for salt, m/s for iceVolume); global: residual = Δ(ocean + ice + airEnthalpy + latent)/Δt − toaNet, surfaceResidual = Δ(ocean + ice)/Δt − (surfaceNet − snowFusion), atmosphereResidual = Δ(airEnthalpy + latent)/Δt − (toaNet − surfaceNet + snowFusion); the sea and ice of a band take netSurfaceHeat − snowFusion',
    boundaries: books.boundaries.map((b) => ({ time: b.time, global: b.global, bands: b.bands })), periods: books.periods,
  },
  books: {
    run: residuals(run.start, totals, run, run.seconds), record: residuals(recordTotals.start, recordTotals.end, recordTotals, recordTotals.seconds),
    units: 'W/m² over the interval: residual = Δ(ocean + ice + airEnthalpy + latent)/Δt − top-of-atmosphere net; surfaceResidual = Δ(ocean + ice)/Δt − (surface net − snowFusion); atmosphereResidual = Δ(airEnthalpy + latent)/Δt − (TOA net − surface net + snowFusion); snowFusion is L_f times the (estimated) snowfall, which the sea and ice give up outside the net surface flux and the air gains, so it cancels in residual',
  },
  nonFinite: counted('nonFinite'), clipped: counted('clipped'),
  daily: { from: runStartTime - recordStart, units: 'one value per model day at its end: Part 7\'s series (precipitableWater and cloudCondensate are the day\'s end, cloudCondensate and cloudWater g/m², ice percent of the globe, from the forcing recorder\'s step means), plus seaSurfaceTemperature °C, mixedLayerDepth m, interiorTemperature °C, surfaceSalinity and meanSalinity psu, fastestCurrent m/s over every layer and edge (fastestCurrentLayer its layer), fastestMixedLayerCurrent m/s, largestTransport Sv (the largest column transport through one edge), clampedEdges, largestSeaSurfaceHeight m, meanSeaLevel m (the global mean of Σh less the floor\'s depth, whose drift is the ocean\'s volume error), ice cover percent and volume m per hemisphere, toaNet, surfaceNet, snowFusion and the storage changes W/m², energyResidual, surfaceResidual and atmosphereResidual W/m² since the run began (defined as in books)', ...daily },
};
const json = Buffer.from(JSON.stringify(header)), pad = (4 - ((8 + json.length) % 4)) % 4;
const blocks = [];
for (const s of [year, spinup]) for (const f of s.frames) for (const field of f) blocks.push(Buffer.from(field.buffer, field.byteOffset, field.byteLength));
const head = Buffer.alloc(8); head.write('EXF1', 0); head.writeUInt32LE(json.length + pad, 4);
mkdirSync(dirname(out), { recursive: true });
writeFileSync(out, Buffer.concat([head, json, Buffer.alloc(pad, 32), ...blocks]));
console.error(`wrote ${out}: year ${year.frames.length} frames, spinup ${spinup.frames.length} frames, ${(json.length / 1024).toFixed(0)} KB header, ${wall.toFixed(1)} s stepping (${(wall / totalDays).toFixed(3)} s per model day), ${seconds().toFixed(1)} s in all`);
model.destroy();
process.exit(0);
