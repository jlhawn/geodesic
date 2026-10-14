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
// MAY_BE_ZERO, or differs from its CPU array by more than TOL (env, relative rms, default 1e-3; the
// state arrays are only reported).
//   --cpu-only   builds only the CPU model (no GPU is imported), steps it, and checks that every
//                mapped CPU array exists and is not all zero: the dry run of this tool
//   --out FILE   also writes the report to FILE;  env FLAGS sets Dawn toggles for the GPU adapter
// Run it with the GPU free (CLAUDE.md): one GPU model at a time.
import { readFileSync, writeFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';

const REPO = fileURLToPath(new URL('..', import.meta.url)).replace(/\/$/, '');
const args = process.argv.slice(2).filter((a) => !a.startsWith('--') && a !== process.argv[process.argv.indexOf('--out') + 1]);
const FILE = args[0], STEPS = (args[1] ?? '1,2').split(',').map(Number), CPU_ONLY = process.argv.includes('--cpu-only');
const OUT = process.argv.includes('--out') ? process.argv[process.argv.indexOf('--out') + 1] : null, TOL = Number(process.env.TOL ?? 1e-3);
if (!FILE) { console.error('usage: node scripts/compareEngines.mjs STATE.bin [STEPS=1,2] [--out FILE] [--cpu-only]'); process.exit(2); }

// PH region -> the CPU model's array ([tail]: the CPU array is K·C long and the region holds its last layers)
export const PH_MAP = {
  SFLUX: 'radiation.surfaceFlux', OFLUX: 'seaIce.oceanFlux', CAP: 'ocean.capacity', MIX: 'boundaryLayer.mixing', DEPTH: 'boundaryLayer.depth', MIXTOP: 'boundaryLayer.mixingTop',
  RAIN: 'moist.precipitation', STEPRAIN: 'moist.precipitation', CONV: 'moist.convectivePrecipitation', COND: 'moist.largeScalePrecipitation', FOG: 'moist.fogPrecipitation',
  SH: 'radiation.sensibleHeat', EVAP: 'radiation.evaporation', SWDN: 'radiation.surfaceShortwave', BUOY: 'radiation.surfaceBuoyancy', USTAR: 'boundaryLayer.friction',
  STRESS: 'boundaryLayer.surfaceStress', SDRAG: 'boundaryLayer.surfaceDrag', HEATX: 'exchange.heat', REFX: 'exchange.reference', DRAG: 'exchange.drag', XWIND: 'exchange.wind',
  REGIME: 'boundaryLayer.regime', ENTRAIN: 'boundaryLayer.entrainment', VRAD: 'boundaryLayer.radiativeVelocity', CTCOOL: 'boundaryLayer.cloudTopCooling',
  DECK: 'radiation.stratus', DECKF: 'radiation.stratusFraction', STRAT: 'radiation.stratiform', MLMCOVER: 'radiation.mlmCover', MLMWATER: 'radiation.mlmWater', MLMENT: 'radiation.mlmEntrainment',
  MLMTOP: 'radiation.mlmTop', MLMH: 'radiation.mlmHeight', MLMSUB: 'radiation.mlmSubsidence', MLMGATE: 'radiation.mlmGate', LWH: 'radiation.longwave',
  REFL: 'radiation.summed.reflectedSolar', REFLSUM: 'radiation.summed.reflectedSolar', ABS: 'radiation.summed.absorbedSolar', ABSSUM: 'radiation.summed.absorbedSolar',
  OLR: 'radiation.held.outgoing', OLRSUM: 'radiation.held.outgoing', ATMSW: 'radiation.summed.atmosphereSolar', ATMSUM: 'radiation.summed.atmosphereSolar',
  INS: 'radiation.summed.insolation', INSSUM: 'radiation.summed.insolation', ABSCLRSUM: 'radiation.summed.clearAbsorbedSolar', OLRCLRSUM: 'radiation.held.clearOutgoing',
  ASRMEAN: 'radiation.meanAbsorbedSolar', OLRMEAN: 'radiation.meanOutgoingLongwave', ALBMEAN: 'radiation.meanPlanetaryAlbedo', SWCREMEAN: 'radiation.meanShortwaveCloudEffect', LWCREMEAN: 'radiation.meanLongwaveCloudEffect',
  CUMF: 'moist.cumulusBaseFlux', CUTOP: 'moist.cumulusTop', CUCOVER: 'moist.cumulusCover[tail]', CUWATER: 'moist.cumulusWater[tail]', SUBTV: 'moist.subcloudVirtual', CONVMEAN: 'moist.convectiveRain', CONDMEAN: 'moist.largeScaleRain',
  PEN: 'ocean.penetrated', TRANS: 'ocean.transmission', CONC: 'seaIce.concentration', SNOW: 'seaIce.snow', SNOWALB: 'seaIce.snowAlbedo',
  SURF: 'land.surface', SOIL: 'land.soil', RUNOFF: 'land.runoff', VEG: 'land.vegetation', CANOPY: 'land.canopy', SEASONL: 'land.seasonLength', SEASONW: 'land.seasonWarmth', DEMAND: 'land.demandMean', RAINMEAN: 'land.rainMean',
  LITTERM: 'land.litterMean', DECAYM: 'land.decayMean', SOILC: 'land.soilCarbon', SNOWFREEV: 'land.snowFreeCover',
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
  FOGMEAN: 'running mean', FSTRESS: 'form drag, off unless formDrag is on', OFLT: 'orographic filter scratch', TOFD: 'form drag profile', RADP: 'held radiation scratch',
};
// mapped regions that may legitimately be all zero in a configuration (env MAY_BE_ZERO_EXTRA adds names, comma separated); the land
// regions are exempt on a planet without land and the radiation.held ones where held radiation is off (no CPU array, no GPU region)
export const MAY_BE_ZERO = new Set(['OFLUX', 'FSTRESS', 'OSTRESS', 'OLAUNCH', 'OWAVE', 'OBETA', 'OBLOCK', 'ODIR', 'GWE', 'GWN', 'RUNOFF', 'SNOW', 'SNOWALB', 'FOG', 'FOGMEAN', ...(process.env.MAY_BE_ZERO_EXTRA ?? '').split(',').filter(Boolean)]);

const typed = (v) => ArrayBuffer.isView(v) && !(v instanceof DataView);
export function cpuArray(cpu, path) {
  const tail = path.endsWith('[tail]'), parts = path.replace('[tail]', '').split('.');
  let v = cpu;
  for (const key of parts) { if (v == null) return null; try { v = v[key]; } catch { return null; } }
  return typed(v) ? { array: v, tail } : null;
}

if (!CPU_ONLY) {
  const require = createRequire(`${REPO}/package.json`);
  const { create, globals } = await import(require.resolve('webgpu'));
  Object.assign(globalThis, globals);
  Object.defineProperty(globalThis, 'navigator', { value: { gpu: create((process.env.FLAGS ?? '').split(' ').filter(Boolean)) }, configurable: true });
}
const { createModel } = await import(`${REPO}/js/model.module.js`);
const { createGpuModel } = CPU_ONLY ? {} : await import(`${REPO}/js/gpu/model.gpu.js`);
const { readBuffer } = CPU_ONLY ? {} : await import(`${REPO}/js/gpu/device.module.js`);

async function build(file) {
  const fs = await import(`${REPO}/scripts/figures/figureState.mjs`);
  const { withCadence } = await import(`${REPO}/js/cadence.module.js`);
  const { savedLevels } = await import(`${REPO}/js/stateFile.module.js`);
  const { savedDeckField, DECK_FIELDS, savedMoistField, MOIST_FIELDS, savedRadiationField, RADIATION_FIELDS } = await import(`${REPO}/js/physics/regrid.module.js`);
  const { stateGrid } = await import(`${REPO}/js/stateGrid.module.js`);
  const { decodeSubgrid } = await import(`${REPO}/js/geography.module.js`);
  const saved = await fs.readState(file), topography = fs.readTopography(), dt = fs.dt(saved.N);
  const mesh = stateGrid(saved, topography).mesh;
  const subgridOf = () => { if (process.env.NO_SUBGRID) return undefined; try { const b = readFileSync(`${REPO}/data/subgrid_N${saved.N}.bin`); return decodeSubgrid(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength)); } catch { return undefined; } };
  const options = { topography, levels: savedLevels(saved), ...withCadence({ ocean: {}, radiation: { clearSkyPass: true } }, dt, saved.N) };
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
    if (saved.land) model.land.load(saved.land);
    return model;
  };
  const cpu = load(createModel(mesh, { ...options, subgrid: subgridOf() }));
  const gpu = CPU_ONLY ? null : load(await createGpuModel(mesh, { ...options, subgrid: subgridOf() }));
  return { cpu, gpu, dt, label: `${file.replace(/.*\//, '')} N=${saved.N}` };
}

const { cpu, gpu, dt, label } = await build(FILE);
const C = cpu.mesh.nCells;
const lines = [`# ${label}: CPU against GPU, dt ${dt} s (${REPO.replace(/.*\//, '')})`];
let failures = 0;
const fail = (text) => { failures++; lines.push(`**${text}**`); };

function stats(g, c, base = null, n = Math.min(g.length, c.length), gOffset = 0, cOffset = 0) {
  let dd = 0, cc = 0, bb = 0, mc = 0, md = 0, gz = 0, cz = 0, nan = 0;
  for (let x = 0; x < n; x++) {
    const a = g[gOffset + x], b = c[cOffset + x];
    if (!Number.isFinite(a) || !Number.isFinite(b)) { if (Number.isFinite(a) !== Number.isFinite(b)) nan++; continue; }
    const d = a - b;
    dd += d * d; cc += b * b; md = Math.max(md, Math.abs(d)); mc = Math.max(mc, Math.abs(b));
    if (a === 0) gz++; if (b === 0) cz++;
    if (base) { const s = b - base[cOffset + x]; bb += s * s; }
  }
  return { rms: cc > 0 ? Math.sqrt(dd / cc) : dd > 0 ? Infinity : 0, max: mc > 0 ? md / mc : md > 0 ? Infinity : 0, step: base ? (bb > 0 ? Math.sqrt(dd / bb) : dd > 0 ? Infinity : 0) : null, gz: gz / n, cz: cz / n, scale: mc, nan, n };
}
const fmt = (v) => (v === null ? '-' : v === Infinity ? 'inf' : v === 0 ? '0' : v.toExponential(1));
const row = (name, s) => `| ${name} | ${fmt(s.rms)} | ${fmt(s.max)} | ${s.step === null ? '-' : fmt(s.step)} | ${(100 * s.gz).toFixed(0)}/${(100 * s.cz).toFixed(0)} | ${fmt(s.scale)}${s.nan ? ` (${s.nan} non-finite on one side)` : ''} |`;
const HEAD = '| field | rel rms | rel max | rms / step change | zero % gpu/cpu | max abs cpu |\n|---|---|---|---|---|---|';

function mappedReport(phNames, phArray, step) {
  lines.push(`\n### PH regions against the fixed name map\n\n| PH | CPU array | rel rms | rel max | zero % gpu/cpu | max abs cpu | verdict |\n|---|---|---|---|---|---|---|`);
  const unmapped = phNames.filter((n) => !(n in PH_MAP) && !(n in PH_INPUT));
  if (unmapped.length) fail(`UNMAPPED PH regions (add them to PH_MAP or PH_INPUT): ${unmapped.join(', ')}`);
  for (const [name, region] of phNames.map((n) => [n, phArray(n)])) {
    const path = PH_MAP[name];
    if (!path) { lines.push(`| ${name} | (input: ${PH_INPUT[name] ?? 'unmapped'}) | | | | | ${region && region.some((x) => !Number.isFinite(x)) ? 'NON-FINITE' : 'ok'} |`); if (region && region.some((x) => !Number.isFinite(x))) failures++; continue; }
    const found = cpuArray(cpu, path);
    if (!found && path.startsWith('radiation.held.') && !cpu.radiation?.held) { lines.push(`| ${name} | ${path} | | | | | n/a (held radiation is off at this N) |`); continue; }
    if (!found) { fail(`MISSING CPU array ${path} for PH ${name}`); lines.push(`| ${name} | ${path} | | | | | MISSING |`); continue; }
    const c = found.tail ? found.array.subarray(found.array.length - (region?.length ?? C)) : found.array;
    const landless = path.startsWith('land.') && (!cpu.geography || !cpu.geography.land || !Array.from(cpu.geography.land).some(Boolean));
    let cNonzero = false; for (let x = 0; x < c.length; x++) if (c[x] !== 0) { cNonzero = true; break; }
    if (!region) { lines.push(`| ${name} | ${path} | | | | ${cNonzero ? '' : '0'} | ${cNonzero ? 'cpu nonzero' : 'CPU ALL ZERO'} |`); if (!cNonzero && !MAY_BE_ZERO.has(name) && !landless) failures++; continue; }
    const s = stats(region, c, null, Math.min(region.length, c.length));
    let gNonzero = false; for (let x = 0; x < region.length; x++) if (region[x] !== 0) { gNonzero = true; break; }
    let verdict = 'ok';
    if (!gNonzero && cNonzero) verdict = 'GPU ALL ZERO'; else if (gNonzero && !cNonzero) verdict = 'CPU ALL ZERO'; else if (!gNonzero && !cNonzero) verdict = MAY_BE_ZERO.has(name) ? 'zero on both (listed)' : 'ZERO ON BOTH (not exercised)';
    else if (!(s.rms <= TOL)) verdict = `OFF (> ${TOL})`;
    if (verdict !== 'ok' && verdict !== 'zero on both (listed)' && !((MAY_BE_ZERO.has(name) || landless) && verdict.startsWith('ZERO'))) failures++;
    lines.push(`| ${name} | ${path} | ${fmt(s.rms)} | ${fmt(s.max)} | ${(100 * s.gz).toFixed(0)}/${(100 * s.cz).toFixed(0)} | ${fmt(s.scale)} | ${verdict} |`);
  }
}

const before = cpu.state.map((a) => Float64Array.from(a));
const oceanBefore = cpu.ocean ? cpu.ocean.serialize() : null;
let done = 0;
for (const target of STEPS) {
  while (done < target) { const t0 = Date.now(); cpu.step(dt); const t1 = Date.now(); if (gpu) { await gpu.step(dt); await gpu.settle(); } done++; console.error(`step ${done}: cpu ${((t1 - t0) / 1000).toFixed(1)} s`); }
  lines.push(`\n## after ${done} step${done > 1 ? 's' : ''}`);
  if (CPU_ONLY) { mappedReport(Object.keys(PH_MAP), () => null, done); continue; }
  await gpu.sync();
  lines.push(`\n### state\n\n${HEAD}`);
  ['pi', 'theta', 'u', 'Ts', 'q', 'qc', 'ice'].forEach((name, a) => lines.push(row(name, stats(gpu.state[a], cpu.state[a], before[a]))));
  lines.push(row('seaIce.concentration', stats(gpu.seaIce.concentration, cpu.seaIce.concentration)));
  if (cpu.ocean && gpu.ocean) {
    const [co, go] = [cpu.ocean.serialize(), await gpu.ocean.serialize()];
    lines.push(`\n### ocean\n\n${HEAD}`);
    for (const key of ['h', 'u', 'T', 'S', 'eta']) lines.push(row(key, stats(go[key], co[key], oceanBefore[key])));
    const content = (o, key) => Float64Array.from(o[key], (v, n) => v * o.h[n]);
    for (const key of ['T', 'S']) lines.push(row(`h·${key}`, stats(content(go, key), content(co, key), content(oceanBefore, key))));
  }
  if (cpu.land && gpu.land) {
    const [cl, gl] = [cpu.land.serialize(), await gpu.land.serialize()];
    lines.push(`\n### land\n\n${HEAD}`);
    for (const key of Object.keys(cl)) if (gl[key] && cl[key] && cl[key].length === gl[key].length && cl[key].length) lines.push(row(key, stats(gl[key], cl[key])));
  }
  const layout = gpu.gpu.layout.PH, ph = await readBuffer(gpu.gpu.device, gpu.gpu.buffers.PH, gpu.gpu.buffers.PH.size, Float32Array);
  const entries = Object.entries(layout).filter(([k, v]) => k !== 'total' && typeof v === 'number').sort((x, y) => x[1] - y[1]);
  const spans = new Map(entries.map(([name, start], e) => [name, ph.subarray(start, e + 1 < entries.length ? entries[e + 1][1] : layout.total ?? ph.length)]));
  const named = [...spans.keys()].filter((n) => spans.get(n).length >= C);
  mappedReport(named, (n) => {
    let g = spans.get(n);
    if ((n === 'DEPTH' || n === 'MIXTOP') && cpu.surfaceGeopotential) { g = Float32Array.from(g); for (let i = 0; i < C; i++) g[i] += cpu.surfaceGeopotential[i] / 9.80665; }
    return g;
  }, done);
}
const text = lines.join('\n');
if (OUT) writeFileSync(OUT, text + '\n');
console.log(text);
console.log(failures ? `\n${failures} failure${failures > 1 ? 's' : ''}` : '\nno failures');
process.exit(failures ? 1 : 0);
