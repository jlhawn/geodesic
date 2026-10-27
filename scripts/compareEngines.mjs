// The CPU engine against the GPU engine at a saved state, field by field.
//
//   node scripts/compareEngines.mjs STATE.bin [STEPS=1,2] [--out FILE] [--cpu-only]
//
// Builds the CPU (f64) and GPU (f32) models of one state, steps both and, after each count in STEPS,
// compares every state array, the ocean (h, u, T, S, eta and the heat and salt contents h·T, h·S), the
// land and every region of the GPU's PH buffer against the CPU array it stands for. The regions are
// matched by the fixed name map PH_MAP below, not by content, so a region that is all zero on one side
// and not on the other cannot pass as a best match of something else. A PH region that is in neither
// PH_MAP nor PH_INPUT is reported first as UNMAPPED, and the run exits 1 when any region is unmapped, is
// missing on the CPU, is all zero on one side only, is all zero on the GPU without being listed in
// MAY_BE_ZERO, holds a non-finite value, differs from its CPU array by more than TOL (env, relative rms,
// default 1e-3), or, for the DECISION_ROWS below, differs in a larger share of cells than the row's bound;
// the state arrays are only reported, with each engine's own rms change since the last report.
//   --cpu-only   builds only the CPU model (no GPU is imported), steps it, and checks that every
//                mapped CPU array exists and is not all zero: the dry run of this tool
//   --gpu-dump P builds only the GPU model, steps it and writes its states at 0 and at each count in STEPS to
//                P.step<n> (each array in its own type): the GPU half, a minute on the GPU at N=128
//   --against P  builds only the CPU model and compares it at each count in STEPS with P.step<n> from a
//                --gpu-dump run of the same state, options and STEPS: the CPU half, which needs no GPU
//   --out FILE   also writes the report to FILE;  env FLAGS sets Dawn toggles for the GPU adapter,
//                env OCEAN a JSON object of ocean options for both models, env RADIATION one of
//                radiation options, env MOIST one of moist options, env SURFACE one of surface options and
//                env LAND one of land options
// Run it with the GPU free (CLAUDE.md): one GPU model at a time.
import { readFileSync, writeFileSync } from 'node:fs';
import { serialize, deserialize } from 'node:v8';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';

const REPO = fileURLToPath(new URL('..', import.meta.url)).replace(/\/$/, '');
const flagValue = (flag) => (process.argv.includes(flag) ? process.argv[process.argv.indexOf(flag) + 1] : null);
const args = process.argv.slice(2).filter((a) => !a.startsWith('--') && !['--out', '--gpu-dump', '--against'].some((flag) => a === flagValue(flag)));
const FILE = args[0], STEPS = (args[1] ?? '1,2').split(',').map(Number), CPU_ONLY = process.argv.includes('--cpu-only');
const GPU_DUMP = flagValue('--gpu-dump'), AGAINST = flagValue('--against'), NO_GPU = CPU_ONLY || !!AGAINST;
const OUT = process.argv.includes('--out') ? process.argv[process.argv.indexOf('--out') + 1] : null, TOL = Number(process.env.TOL ?? 1e-3);
if (!FILE) { console.error('usage: node scripts/compareEngines.mjs STATE.bin [STEPS=1,2] [--out FILE] [--cpu-only]'); process.exit(2); }

// PH region -> the CPU model's array ([tail]: the CPU array is K·C long and the region holds its last layers)
export const PH_MAP = {
  SFLUX: 'radiation.surfaceFlux', OFLUX: 'seaIce.oceanFlux', CAP: 'ocean.capacity', MIX: 'boundaryLayer.mixing', DEPTH: 'boundaryLayer.depth', MIXTOP: 'boundaryLayer.mixingTop',
  RAIN: 'moist.precipitation', STEPRAIN: 'moist.rain', CONV: 'moist.convectivePrecipitation', COND: 'moist.largeScalePrecipitation', FOG: 'moist.fogPrecipitation',
  SH: 'radiation.sensibleHeat', EVAP: 'radiation.evaporation', SWDN: 'radiation.surfaceShortwave', BUOY: 'radiation.surfaceBuoyancy', USTAR: 'boundaryLayer.friction',
  STRESS: 'boundaryLayer.surfaceStress', SDRAG: 'boundaryLayer.surfaceDrag', HEATX: 'exchange.heat', REFX: 'exchange.reference', DRAG: 'exchange.drag', XWIND: 'exchange.wind',
  REGIME: 'boundaryLayer.regime', ENTRAIN: 'boundaryLayer.entrainment', VRAD: 'boundaryLayer.radiativeVelocity', CTCOOL: 'boundaryLayer.cloudTopCooling',
  DECK: 'radiation.stratus', DECKF: 'radiation.stratusFraction', STRAT: 'radiation.stratiform', MLMCOVER: 'radiation.mlmCover', MLMWATER: 'radiation.mlmWater', MLMENT: 'radiation.mlmEntrainment',
  MLMTOP: 'radiation.mlmTop', MLMH: 'radiation.mlmHeight', MLMSUB: 'radiation.mlmSubsidence', MLMGATE: 'radiation.mlmGate', LWH: 'radiation.longwave',
  REFLSUM: 'radiation.summed.reflectedSolar', ABSSUM: 'radiation.summed.absorbedSolar',
  OLR: 'radiation.held.outgoing', OLRSUM: 'radiation.summed.outgoingLongwave', ATMSUM: 'radiation.summed.atmosphereSolar',
  INSSUM: 'radiation.summed.insolation', ABSCLRSUM: 'radiation.summed.clearAbsorbedSolar', OLRCLRSUM: 'radiation.summed.clearOutgoingLongwave',
  ASRMEAN: 'radiation.meanAbsorbedSolar', OLRMEAN: 'radiation.meanOutgoingLongwave', ALBMEAN: 'radiation.meanPlanetaryAlbedo', SWCREMEAN: 'radiation.meanShortwaveCloudEffect', LWCREMEAN: 'radiation.meanLongwaveCloudEffect',
  CUMF: 'moist.cumulusBaseFlux', CUTOP: 'moist.cumulusTop', CUCOVER: 'moist.cumulusCover[tail]', CUWATER: 'moist.cumulusWater[tail]', SUBTV: 'moist.subcloudVirtual', CONVMEAN: 'moist.convectiveRain', CONDMEAN: 'moist.largeScaleRain',
  PEN: 'ocean.penetrated', TRANS: 'ocean.transmission', CONC: 'seaIce.concentration', SNOW: 'seaIce.snow', SNOWALB: 'seaIce.snowAlbedo',
  SURF: 'land.surface', SOIL: 'land.soil', RUNOFF: 'land.runoff', DISCHARGE: 'land.discharge', ICELOSS: 'land.iceLoss', MLMK: 'radiation.mlmLayer', MLMM: 'radiation.mlmShare', DECKK: 'radiation.deckLayer', DECKFROM: 'radiation.deckFrom', VEG: 'land.vegetation', CANOPY: 'land.canopy', SEASONL: 'land.seasonLength', SEASONW: 'land.seasonWarmth', DEMAND: 'land.demandMean', RAINMEAN: 'land.rainMean',
  LITTERM: 'land.litterMean', DECAYM: 'land.decayMean', SOILC: 'land.soilCarbon', SNOWFREEV: 'land.snowFreeCover', SOILT2: 'land.soilT2', SOILT3: 'land.soilT3', SOILT4: 'land.soilT4',
  GWE: 'gravityWaves.east', GWN: 'gravityWaves.north', OLAUNCH: 'orography.launch', OWAVE: 'orography.wave', OSTRESS: 'orography.stress', OBETA: 'orography.beta', OBLOCK: 'orography.blocking', ODIR: 'orography.direction',
  OSTD: 'orography.fields.deviation', OORI: 'orography.fields.orientation', OANI: 'orography.fields.anisotropy', OSLP: 'orography.fields.slope',
  RADT: 'radiation.held.escape', RADTC: 'radiation.held.clearEscape', RADDF: 'radiation.held.share', RADDIR: 'radiation.held.direct', RADREFL: 'radiation.held.reflected', RADABS: 'radiation.held.absorbed',
  RADSWDN: 'radiation.held.down', RADSW: 'radiation.held.shortwave', RADTOA: 'radiation.held.top', RADBACK: 'radiation.held.back', RADOLR: 'radiation.held.outgoing', RADATM: 'radiation.held.atmosphere',
  RADCLROLR: 'radiation.held.clearOutgoing', RADCLRSW: 'radiation.held.clearAbsorbed', RADEMIT: 'radiation.held.emission', RADMU: 'radiation.held.mu',
};
// regions that are inputs or hand-offs with no CPU array of their own (they are checked for non-finite values only)
export const PH_INPUT = {
  TAU: 'relaxation time, an input', ADIF: 'diffuse albedo, a local of model.module.js', LAND: 'surface flags 0/1/2, an input', CUFALL: 'plume-to-rain hand-off', CURESERVE: 'plume-to-rain hand-off', CUFROZEN: 'plume-to-rain hand-off',
  CUDEEP: 'plume-to-rain hand-off', CUBASE: 'plume-to-rain hand-off', CURAIN: 'cumulus part of CONV', LWSFCSUM: 'running sum for the means', CLOUDK: 'cloud layer index, written inside the physics kernel',
  CUSHALLOW: 'plume diagnostic', CUSNOW: 'plume diagnostic', MOMU: 'momentum-layer scratch', MOMK: 'momentum-layer scratch', MOMD: 'momentum-layer scratch', MOMKD: 'momentum-layer scratch', MOMS: 'momentum-layer scratch',
  ABS: 'per-step value whose CPU counterpart is the running sum ABSSUM compares', INS: 'per-step value whose CPU counterpart is the running sum INSSUM compares', REFL: 'per-step value whose CPU counterpart is the running sum REFLSUM compares', ATMSW: 'per-step value whose CPU counterpart is the running sum ATMSUM compares',
  FOGMEAN: 'running mean', FSTRESS: 'form drag, off unless formDrag is on', OFLT: 'orographic filter scratch', TOFD: 'form drag profile', RADP: 'held radiation scratch', LYR: 'layered shortwave hand-off from the radiation kernel to its pass',
};
// mapped regions that may legitimately be all zero in a configuration (env MAY_BE_ZERO_EXTRA adds names, comma separated); the land
// regions are exempt on a planet without land and the radiation.held ones where held radiation is off (no CPU array, no GPU region)
export const MAY_BE_ZERO = new Set(['OFLUX', 'FSTRESS', 'OSTRESS', 'OLAUNCH', 'OWAVE', 'OBETA', 'OBLOCK', 'ODIR', 'GWE', 'GWN', 'RUNOFF', 'DISCHARGE', 'ICELOSS', 'MLMK', 'MLMM', 'SNOW', 'SNOWALB', 'FOG', 'FOGMEAN', ...(process.env.MAY_BE_ZERO_EXTRA ?? '').split(',').filter(Boolean)]);
// mapped regions set by a threshold decision in each column (mixing regime, cumulus trigger, rain, fog, launch): f32 and
// f64 part at a few columns near the threshold, so these are judged by the share of cells that differ by more than TOL of
// the row's largest CPU value and fail past their bound here, max(3 × the worst share measured on
// runs/sixteen/sixteen128_day0730.bin through 17 steps, 1e-3); env SHARE sets one bound for every row, env DECISION_EXTRA
// adds rows (comma separated) under it
export const DECISION_BOUNDS = { MIX: 0.016, DEPTH: 0.048, MIXTOP: 0.048, ENTRAIN: 0.21, REGIME: 0.03, VRAD: 0.15, CTCOOL: 0.14, RAIN: 0.001, STEPRAIN: 0.026, CONV: 0.02, FOG: 0.0033, CUMF: 0.22, CUTOP: 0.072, CUCOVER: 0.051, CUWATER: 0.036, GWE: 0.001, GWN: 0.001, OWAVE: 0.001 };
for (const name of (process.env.DECISION_EXTRA ?? '').split(',').filter(Boolean)) DECISION_BOUNDS[name] ??= 1e-3;
export const DECISION_ROWS = new Set(Object.keys(DECISION_BOUNDS));
const decisionBound = (name) => (process.env.SHARE ? Number(process.env.SHARE) : DECISION_BOUNDS[name]);

const typed = (v) => ArrayBuffer.isView(v) && !(v instanceof DataView);
export function cpuArray(cpu, path) {
  const tail = path.endsWith('[tail]'), parts = path.replace('[tail]', '').split('.');
  let v = cpu;
  for (const key of parts) { if (v == null) return null; try { v = v[key]; } catch { return null; } }
  return typed(v) ? { array: v, tail } : null;
}

if (!NO_GPU) {
  const require = createRequire(`${REPO}/package.json`);
  const { create, globals } = await import(require.resolve('webgpu'));
  Object.assign(globalThis, globals);
  Object.defineProperty(globalThis, 'navigator', { value: { gpu: create((process.env.FLAGS ?? '').split(' ').filter(Boolean)) }, configurable: true });
}
const { createModel } = await import(`${REPO}/js/model.module.js`);
const { createGpuModel } = NO_GPU ? {} : await import(`${REPO}/js/gpu/model.gpu.js`);
const { readBuffer } = NO_GPU ? {} : await import(`${REPO}/js/gpu/device.module.js`);

async function build(file) {
  const fs = await import(`${REPO}/scripts/figures/figureState.mjs`);
  const { withCadence } = await import(`${REPO}/js/cadence.module.js`);
  const { savedLevels } = await import(`${REPO}/js/stateFile.module.js`);
  const { savedDeckField, DECK_FIELDS, savedMoistField, MOIST_FIELDS, savedRadiationField, RADIATION_FIELDS } = await import(`${REPO}/js/physics/regrid.module.js`);
  const { stateGrid } = await import(`${REPO}/js/stateGrid.module.js`);
  const { decodeSubgrid, stateGeography, subgridUrl } = await import(`${REPO}/js/geography.module.js`);
  const saved = await fs.readState(file), sg = stateGeography(saved), topography = sg.topography, dt = fs.dt(saved.N);
  const mesh = stateGrid(saved, topography).mesh;
  const subgridOf = () => { if (process.env.NO_SUBGRID) return undefined; try { const b = readFileSync(subgridUrl(saved.N, sg.record.topography)); return decodeSubgrid(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength)); } catch { return undefined; } };
  const options = { topography, geography: sg.geography, land: { ...sg.land, ...JSON.parse(process.env.LAND ?? '{}') }, levels: savedLevels(saved), ...(process.env.MOIST ? { moist: JSON.parse(process.env.MOIST) } : {}), ...(process.env.SURFACE ? { surface: JSON.parse(process.env.SURFACE) } : {}), ...withCadence({ ocean: JSON.parse(process.env.OCEAN ?? '{}'), radiation: { clearSkyPass: true, ...JSON.parse(process.env.RADIATION ?? '{}') } }, dt, saved.N) };
  const load = (model) => {
    const { state } = model, C = model.mesh.nCells, K = model.core.K;
    ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => state[a].set(saved[name]));
    model.seaIce.load(state[6], saved.concentration ?? null);
    for (const field of Object.keys(DECK_FIELDS)) model.radiation[field].set(savedDeckField(saved, field, model));
    for (const field of Object.keys(MOIST_FIELDS)) model.moist[field].set(savedMoistField(saved, field, model));
    for (const field of Object.keys(RADIATION_FIELDS)) model.radiation[field].set(savedRadiationField(saved, field, model));
    if (saved.boundaryDepth) model.boundaryLayer.depth.set(saved.boundaryDepth);
    if (saved.mixingTop) model.boundaryLayer.mixingTop.set(saved.mixingTop);
    if (saved.boundaryRegime) model.boundaryLayer.regime.set(saved.boundaryRegime);
    if (saved.boundaryBuoyancy) model.boundaryLayer.buoyancyFlux.set(saved.boundaryBuoyancy);
    if (saved.windSpeed) model.surface.windSpeed.set(saved.windSpeed);
    if (saved.evaporation) model.radiation.evaporation.set(saved.evaporation);
    if (saved.exchangeHeat && model.exchange && !model.exchange.fixed) model.exchange.heat.set(saved.exchangeHeat);
    if (saved.exchangeWind && model.exchange) model.exchange.wind.set(saved.exchangeWind);
    if (model.skin) model.skin.load(saved);
    for (const field of ['cumulusCover', 'cumulusWater']) {
      const v = saved[field], t = model.moist[field];
      if (!v) continue;
      if (v.length === t.length) t.set(v); else if (t.length === K * C) { t.fill(0); t.set(v, K * C - v.length); }
    }
    const sv = model.moist.subcloudVirtual;
    if (saved.subcloudVirtual && saved.subcloudVirtual.length === sv.length) sv.set(saved.subcloudVirtual); else sv.fill(0);
    if (!model.gpu && model.surfaceGeopotential) for (let i = 0; i < C; i++) { const zs = model.surfaceGeopotential[i] / 9.80665; if (saved.boundaryDepth) model.boundaryLayer.depth[i] += zs; if (saved.mixingTop) model.boundaryLayer.mixingTop[i] += zs; }
    model.time = saved.time;
    if (model.load) model.load();
    model.ocean.load(saved.ocean, state[3], state[6]);
    if (saved.land) model.land.load(saved.land, state[6], saved.time);
    return model;
  };
  const cpu = GPU_DUMP ? null : load(createModel(mesh, { ...options, subgrid: subgridOf() }));
  const gpu = NO_GPU ? null : load(await createGpuModel(mesh, { ...options, subgrid: subgridOf() }));
  return { cpu, gpu, dt, label: `${file.replace(/.*\//, '')} N=${saved.N}` };
}

const { cpu, gpu, dt, label } = await build(FILE);
const C = (cpu ?? gpu).mesh.nCells;
const lines = [`# ${label}: CPU against GPU, dt ${dt} s (${REPO.replace(/.*\//, '')})`];
let failures = 0;
const fail = (text) => { failures++; lines.push(`**${text}**`); };

function stats(g, c, base = null, n = Math.min(g.length, c.length), gOffset = 0, cOffset = 0, gBase = null) {
  let dd = 0, cc = 0, bb = 0, gg = 0, mc = 0, md = 0, gz = 0, cz = 0, nan = 0;
  for (let x = 0; x < n; x++) {
    const a = g[gOffset + x], b = c[cOffset + x];
    if (!Number.isFinite(a) || !Number.isFinite(b)) { if (Number.isFinite(a) !== Number.isFinite(b)) nan++; continue; }
    const d = a - b;
    dd += d * d; cc += b * b; md = Math.max(md, Math.abs(d)); mc = Math.max(mc, Math.abs(b));
    if (a === 0) gz++; if (b === 0) cz++;
    if (base) { const s = b - base[cOffset + x]; bb += s * s; }
    if (gBase) { const s = a - gBase[gOffset + x]; gg += s * s; }
  }
  let over = 0;
  for (let x = 0; x < n; x++) { const a = g[gOffset + x], b = c[cOffset + x]; if (Number.isFinite(a) && Number.isFinite(b) && Math.abs(a - b) > TOL * mc) over++; }
  return { own: base && gBase ? [Math.sqrt(bb / n), Math.sqrt(gg / n)] : null, share: over / n, rms: cc > 0 ? Math.sqrt(dd / cc) : dd > 0 ? Infinity : 0, max: mc > 0 ? md / mc : md > 0 ? Infinity : 0, step: base ? (bb > 0 ? Math.sqrt(dd / bb) : dd > 0 ? Infinity : 0) : null, gz: gz / n, cz: cz / n, scale: mc, nan, n };
}
const fmt = (v) => (v === null ? '-' : v === Infinity ? 'inf' : v === 0 ? '0' : v.toExponential(1));
const row = (name, s) => `| ${name} | ${fmt(s.rms)} | ${fmt(s.max)} | ${s.step === null ? '-' : fmt(s.step)} | ${(100 * s.gz).toFixed(0)}/${(100 * s.cz).toFixed(0)} | ${s.own ? `${fmt(s.own[0])} / ${fmt(s.own[1])}` : '-'} | ${fmt(s.scale)}${s.nan ? ` (${s.nan} non-finite on one side)` : ''} |`;
const HEAD = '| field | rel rms | rel max | rms / change since last report | zero % gpu/cpu | rms change cpu / gpu | max abs cpu |\n|---|---|---|---|---|---|---|';

function nonFinite(region) {
  let count = 0, first = -1;
  for (let x = 0; x < region.length; x++) if (!Number.isFinite(region[x])) { if (first < 0) first = x; count++; }
  return count ? `NON-FINITE (${count} of ${region.length}, first at layer ${Math.floor(first / C)} cell ${first % C})` : null;
}

function mappedReport(phNames, phArray, step) {
  lines.push(`\n### PH regions against the fixed name map\n\n| PH | CPU array | rel rms | rel max | zero % gpu/cpu | max abs cpu | verdict |\n|---|---|---|---|---|---|---|`);
  const decisions = [];
  const unmapped = phNames.filter((n) => !(n in PH_MAP) && !(n in PH_INPUT));
  if (unmapped.length) fail(`UNMAPPED PH regions (add them to PH_MAP or PH_INPUT): ${unmapped.join(', ')}`);
  for (const [name, region] of phNames.map((n) => [n, phArray(n)])) {
    const path = PH_MAP[name];
    if (!path) { const bad = region ? nonFinite(region) : null; lines.push(`| ${name} | (input: ${PH_INPUT[name] ?? 'unmapped'}) | | | | | ${bad ?? 'ok'} |`); if (bad) failures++; continue; }
    const found = cpuArray(cpu, path);
    if (!found && path.startsWith('radiation.held.') && !cpu.radiation?.held) { lines.push(`| ${name} | ${path} | | | | | n/a (held radiation is off at this N) |`); continue; }
    if (!found) { fail(`MISSING CPU array ${path} for PH ${name}`); lines.push(`| ${name} | ${path} | | | | | MISSING |`); continue; }
    const c = found.tail && region ? found.array.subarray(found.array.length - region.length) : found.array;
    const landless = path.startsWith('land.') && (!cpu.geography || !cpu.geography.land || !Array.from(cpu.geography.land).some(Boolean));
    let cNonzero = false; for (let x = 0; x < c.length; x++) if (c[x] !== 0) { cNonzero = true; break; }
    if (!region) { lines.push(`| ${name} | ${path} | | | | ${cNonzero ? '' : '0'} | ${cNonzero ? 'cpu nonzero' : 'CPU ALL ZERO'} |`); if (!cNonzero && !MAY_BE_ZERO.has(name) && !landless) failures++; continue; }
    const s = stats(region, c, null, Math.min(region.length, c.length));
    let gNonzero = false; for (let x = 0; x < region.length; x++) if (region[x] !== 0) { gNonzero = true; break; }
    let verdict = 'ok';
    if (!gNonzero && cNonzero) verdict = 'GPU ALL ZERO'; else if (gNonzero && !cNonzero) verdict = 'CPU ALL ZERO'; else if (!gNonzero && !cNonzero) verdict = MAY_BE_ZERO.has(name) ? 'zero on both (listed)' : 'ZERO ON BOTH (not exercised)';
    else if (DECISION_ROWS.has(name)) verdict = s.share > decisionBound(name) ? `OFF (share > ${decisionBound(name)})` : 'ok';
    else if (!(s.rms <= TOL)) verdict = `OFF (> ${TOL})`;
    if (verdict !== 'ok' && verdict !== 'zero on both (listed)' && !((MAY_BE_ZERO.has(name) || landless) && verdict.startsWith('ZERO'))) failures++;
    if (DECISION_ROWS.has(name) && gNonzero && cNonzero) decisions.push(`| ${name} | ${path} | ${fmt(s.share)} | ${decisionBound(name)} | ${fmt(s.rms)} | ${fmt(s.max)} | ${fmt(s.scale)} | ${verdict} |`);
    else lines.push(`| ${name} | ${path} | ${fmt(s.rms)} | ${fmt(s.max)} | ${(100 * s.gz).toFixed(0)}/${(100 * s.cz).toFixed(0)} | ${fmt(s.scale)} | ${verdict} |`);
  }
  lines.push(`\n### decision divergence: share of cells\n\n| PH | CPU array | share > TOL | bound | rel rms | rel max | max abs cpu | verdict |\n|---|---|---|---|---|---|---|---|`, ...decisions);
}

const held = (v) => (typed(v) ? v.slice() : Array.isArray(v) && v.every((x) => typeof x === 'number') ? Float64Array.from(v) : v);
const heldAll = (o) => Object.fromEntries(Object.entries(o).map(([k, v]) => [k, held(v)]));
async function gpuSnapshot() {
  await gpu.sync();
  const layout = gpu.gpu.layout.PH, ph = await readBuffer(gpu.gpu.device, gpu.gpu.buffers.PH, gpu.gpu.buffers.PH.size, Float32Array);
  const entries = Object.entries(layout).filter(([k, v]) => k !== 'total' && typeof v === 'number').sort((x, y) => x[1] - y[1]);
  const regions = Object.fromEntries(entries.map(([name, start], e) => [name, ph.slice(start, e + 1 < entries.length ? entries[e + 1][1] : layout.total ?? ph.length)]).filter(([, r]) => r.length >= C));
  return {
    state: gpu.state.map((a) => Float32Array.from(a)), concentration: Float32Array.from(gpu.seaIce.concentration), regions,
    ocean: gpu.ocean ? heldAll(await gpu.ocean.serialize()) : null, land: gpu.land ? heldAll(await gpu.land.serialize()) : null,
  };
}
const gpuAt = async (n) => (AGAINST ? deserialize(readFileSync(`${AGAINST}.step${n}`)) : gpuSnapshot());
let before = cpu ? cpu.state.map((a) => Float64Array.from(a)) : null;
let oceanBefore = cpu && cpu.ocean ? cpu.ocean.serialize() : null;
let gBefore = CPU_ONLY ? null : await gpuAt(0);
if (GPU_DUMP) writeFileSync(`${GPU_DUMP}.step0`, serialize(gBefore));
let done = 0;
for (const target of STEPS) {
  while (done < target) { const t0 = Date.now(); if (cpu) cpu.step(dt); const t1 = Date.now(); if (gpu) { await gpu.step(dt); await gpu.settle(); } done++; console.error(`step ${done}${cpu ? `: cpu ${((t1 - t0) / 1000).toFixed(1)} s` : ''}`); }
  if (GPU_DUMP) { writeFileSync(`${GPU_DUMP}.step${done}`, serialize(await gpuSnapshot())); continue; }
  lines.push(`\n## after ${done} step${done > 1 ? 's' : ''}`);
  if (CPU_ONLY) { mappedReport(Object.keys(PH_MAP), () => null, done); continue; }
  const g = await gpuAt(done);
  lines.push(`\n### state\n\n${HEAD}`);
  ['pi', 'theta', 'u', 'Ts', 'q', 'qc', 'ice'].forEach((name, a) => lines.push(row(name, stats(g.state[a], cpu.state[a], before[a], undefined, 0, 0, gBefore.state[a]))));
  lines.push(row('seaIce.concentration', stats(g.concentration, cpu.seaIce.concentration)));
  if (cpu.ocean && g.ocean) {
    const co = cpu.ocean.serialize(), go = g.ocean, gb = gBefore.ocean;
    lines.push(`\n### ocean\n\n${HEAD}`);
    for (const key of ['h', 'u', 'T', 'S', 'eta']) lines.push(row(key, stats(go[key], co[key], oceanBefore[key], undefined, 0, 0, gb[key])));
    const content = (o, key) => Float64Array.from(o[key], (v, n) => v * o.h[n]);
    for (const key of ['T', 'S']) lines.push(row(`h·${key}`, stats(content(go, key), content(co, key), content(oceanBefore, key), undefined, 0, 0, content(gb, key))));
  }
  if (cpu.land && g.land) {
    const cl = cpu.land.serialize(), gl = g.land;
    lines.push(`\n### land\n\n${HEAD}`);
    for (const key of Object.keys(cl)) if (gl[key] && cl[key] && cl[key].length === gl[key].length && cl[key].length) lines.push(row(key, stats(gl[key], cl[key])));
  }
  mappedReport(Object.keys(g.regions), (n) => {
    let r = g.regions[n];
    if ((n === 'DEPTH' || n === 'MIXTOP') && cpu.surfaceGeopotential) { r = Float32Array.from(r); for (let i = 0; i < C; i++) r[i] += cpu.surfaceGeopotential[i] / 9.80665; }
    return r;
  }, done);
  before = cpu.state.map((x) => Float64Array.from(x));
  if (cpu.ocean) oceanBefore = cpu.ocean.serialize();
  gBefore = g;
}
const text = lines.join('\n');
if (OUT) writeFileSync(OUT, text + '\n');
console.log(text);
console.log(failures ? `\n${failures} failure${failures > 1 ? 's' : ''}` : '\nno failures');
process.exit(failures ? 1 : 0);
