#!/usr/bin/env node
/*
 * The default-physics digest: a fresh start at N = 16 on createModel's own defaults (so the land's coverGoal is 'water',
 * the model default; 'observed' is a spin-up option), with the observed geography, the ocean atlas's March and the
 * observed sea ice with its concentration, as a spin-up's fresh start takes them, and the clear-sky pass for the cloud
 * effects; stepped STEPS steps of 1350 s on the CPU from 00 UTC on day 0, then reduced to a fixed table of physical rows:
 * global and 30S–30N means of the step-mean fields. A change to a default path that this run executes moves the state
 * hash and, past its first steps, rows; --compare lists the rows that moved against an earlier digest, so a re-pin can
 * cite them.
 *
 * What it cannot see: the held radiation and drags of the page and spin-ups (radiationEvery and dragEvery are 1 at
 * N = 16, against 8 or 16 and 2 or 4 at N = 64 and 128); the ocean's cadence at N >= 64; WGSL-only code (the GPU
 * engine's parity tests cover it); restore, regrid and re-seed paths, since it starts fresh; and any dated input outside
 * 00–06 UTC of day 0. Those need the N = 64 scorecard (scripts/scorecard64.mjs) or their own tests.
 *
 *   node scripts/physicsDigest.mjs [--steps=16] [--out=digest.json] [--compare=earlier.json]
 *
 * Rows are area-weighted means over all cells unless named otherwise: land and sea split by the geography's land mask,
 * sea ice over sea cells (not lakes) by hemisphere, the soil fill over open land (off the ice sheets). Radiation is the
 * model's own mean over its radiation steps (readMeans), precipitation the accumulated rain over the elapsed time, and
 * every other row the mean of the field after each step. "sea-cell surface temperature" is surfaceT over sea cells, the
 * ice's skin on iced cells, not the SST. The state hash fingerprints, at the end, the atmosphere's state arrays, the sea
 * ice, the land, the ocean, the lakes, the boundary layer, the sea skin and the decks (the groups it hashed are listed).
 * Six hours from rest and a dry initial atmosphere are a transient, not a climate: the rows detect a change to the
 * physics and size it; the N = 64 scorecard is what judges the model against observation.
 * test/helpers/physicsDigest.json is the digest of record, with the Node version, the commit and the data/ tree (all of
 * it, not only what the digest reads) it was made on and a hash per group;
 * test/physicsDigest.test.mjs runs the digest and fails on any moved row or state, and a change that moves it re-pins
 * it (--out) and cites the rows --compare lists.
 */
import { readFileSync, writeFileSync, realpathSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { createHash } from 'node:crypto';
import { execSync } from 'node:child_process';

const DAY = 86400, DT = 1350, MOVED = 1e-9;

export function areaMean(mesh, values, mask = null) {
  let a = 0, s = 0;
  for (let i = 0; i < mesh.nCells; i++) if (!mask || mask[i]) { a += mesh.areaCell[i]; s += mesh.areaCell[i] * values[i]; }
  return a > 0 ? s / a : NaN;
}

export function bandMask(mesh, south, north, within = null) {
  return Uint8Array.from({ length: mesh.nCells }, (_, i) => {
    const lat = mesh.latCell[i] * 180 / Math.PI;
    return lat >= south && lat <= north && (!within || within[i]) ? 1 : 0;
  });
}

export function hemisphereArea(mesh, cover, sea, north) {
  let s = 0;
  for (let i = 0; i < mesh.nCells; i++) if (sea[i] && (mesh.latCell[i] > 0) === north) s += mesh.areaCell[i] * cover[i];
  return s / 1e12;
}

export function compareDigests(earlier, later) {
  const before = new Map(earlier.rows.map((r) => [r.name, r])), moved = [], added = [];
  for (const row of later.rows) {
    const old = before.get(row.name);
    if (!old) { added.push(row.name); continue; }
    before.delete(row.name);
    const delta = row.value - old.value, scale = Math.max(Math.abs(row.value), Math.abs(old.value), 1e-30);
    if (!(Math.abs(delta) <= MOVED * scale)) moved.push({ name: row.name, units: row.units, before: old.value, after: row.value, delta, relative: delta / scale });
  }
  moved.sort((x, y) => Math.abs(y.relative) - Math.abs(x.relative));
  return { moved, added, removed: [...before.keys()], sameState: earlier.stateHash === later.stateHash };
}

const typedArrays = (object) => Object.keys(object ?? {}).sort().flatMap((key) => {
  const value = object[key];
  if (ArrayBuffer.isView(value)) return [value];
  if (Array.isArray(value) && value.length && value.every((x) => typeof x === 'number')) return [Float64Array.from(value)];
  return [];
});

export function stateHash(model) {
  const parts = [
    ['atmosphere', model.state], ['sea ice', typedArrays(model.seaIce)], ['land', typedArrays(model.land?.serialize?.())],
    ['ocean', typedArrays(model.ocean?.serialize?.())], ['lakes', typedArrays(model.lake?.serialize?.())], ['boundary layer', typedArrays(model.boundaryLayer)],
    ['sea skin', typedArrays(model.skin)], ['decks', typedArrays({ stratus: model.radiation.stratus, stratusFraction: model.radiation.stratusFraction, mlmTop: model.radiation.mlmTop })],
  ].filter(([, arrays]) => arrays.length);
  const hash = createHash('sha256'), groupHashes = {};
  for (const [name, arrays] of parts) {
    const group = createHash('sha256');
    for (const a of arrays) { const bytes = new Uint8Array(a.buffer, a.byteOffset, a.byteLength); hash.update(bytes); group.update(bytes); }
    groupHashes[name] = group.digest('hex').slice(0, 16);
  }
  return { hash: hash.digest('hex').slice(0, 32), groups: parts.map(([name, arrays]) => `${name} ${arrays.length}`), groupHashes };
}

function gitOf(command, paths) {
  try {
    const run = (cmd) => execSync(cmd, { cwd: fileURLToPath(new URL('..', import.meta.url)), encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] }).trim();
    return `${run(command)}${run(`git status --porcelain --untracked-files=no -- ${paths}`) ? ' with changes' : ''}`;
  } catch { return null; }
}
const commitOf = () => gitOf('git rev-parse --short=12 HEAD', 'js scripts data');
const dataOf = () => gitOf('git rev-parse --short=12 HEAD:data', 'data');

export function formatValue(x) {
  return !Number.isFinite(x) ? String(x) : Math.abs(x) >= 100 ? x.toFixed(2) : Math.abs(x) >= 1 ? x.toFixed(3) : x.toPrecision(4);
}

export async function runDigest({ steps = 16 } = {}) {
  const { Grid } = await import('../js/grid.module.js');
  const { createModel } = await import('../js/model.module.js');
  const { initializeState } = await import('../js/physics/init.module.js');
  const { stateGeography, FRESH_GEOGRAPHY } = await import('../js/geography.module.js');
  const { gridHeader } = await import('../js/stateGrid.module.js');
  const { START_CLIMATOLOGY_FILE, loadClimatology } = await import('../js/ocean/climatology.module.js');
  const { SEA_ICE_START_FILE, loadSeaIceStart } = await import('../js/physics/seaIceStart.module.js');
  const repo = (path) => fileURLToPath(new URL(`../${path}`, import.meta.url));
  const grid = new Grid(16), { topography, geography, land } = stateGeography({ N: 16, ...gridHeader(grid), ...FRESH_GEOGRAPHY });
  const model = createModel(grid, { topography, geography, land, ocean: { climatology: await loadClimatology(repo(START_CLIMATOLOGY_FILE)) }, radiation: { clearSkyPass: true } });
  const { mesh, state } = model, C = mesh.nCells;
  initializeState(model, { geostrophic: !model.surfaceGeopotential, seaIce: await loadSeaIceStart(repo(SEA_ICE_START_FILE)) }).forEach((values, a) => state[a].set(values));
  const landMask = model.geography.land, lake = model.geography.lake;
  for (let i = 0; i < C; i++) if (landMask[i]) state[6][i] = 0;
  model.ocean.initialize(state[3], state[6]);
  model.land.initialize(model.time);
  const fields = ['surfaceT', 'evaporation', 'sensible', 'vapor', 'cloud', 'wind', 'ice', 'lowCover'], sum = Object.fromEntries(fields.map((f) => [f, new Float64Array(C)]));
  const moist = model.moist, radiation = model.radiation;
  for (let n = 0; n < steps; n++) {
    model.step(DT);
    const [pi, , , surfaceT, q] = state;
    for (let i = 0; i < C; i++) {
      sum.surfaceT[i] += surfaceT[i]; sum.evaporation[i] += radiation.evaporation[i]; sum.sensible[i] += radiation.sensibleHeat[i];
      sum.vapor[i] += moist.columnWater(pi, q, i); sum.cloud[i] += model.cloudWater(i);
      sum.wind[i] += model.surface.windSpeed[i]; sum.ice[i] += model.seaIce.cover(i, state[6][i]); sum.lowCover[i] += radiation.lowCover[i];
    }
  }
  radiation.readMeans(model.radiationSteps);
  const mean = Object.fromEntries(fields.map((f) => [f, Float64Array.from(sum[f], (s) => s / steps)]));
  const sea = Uint8Array.from({ length: C }, (_, i) => (!landMask[i] && !(lake && lake[i]) ? 1 : 0)), landOnly = Uint8Array.from(landMask, (x) => (x ? 1 : 0));
  const openLand = Uint8Array.from({ length: C }, (_, i) => (landMask[i] && !(model.geography.iceSheet && model.geography.iceSheet[i]) ? 1 : 0));
  const tropics = bandMask(mesh, -30, 30), mmDay = (x) => x * DAY, C2K = 273.15;
  const precip = Float64Array.from(moist.precipitation, (x) => x / model.time), convective = Float64Array.from(moist.convectivePrecipitation, (x) => x / model.time);
  const summed = radiation.summed, total = (name) => areaMean(mesh, summed[name]);
  const fill = Float64Array.from({ length: C }, (_, i) => (openLand[i] ? model.land.soil[i] / model.land.capacity(i) : 0));
  const rows = [];
  const row = (name, units, value) => rows.push({ name, units, value });
  for (const [label, mask] of [['global', null], ['30S-30N', tropics]]) {
    row(`${label} surface temperature`, '°C', areaMean(mesh, mean.surfaceT, mask) - C2K);
    row(`${label} absorbed solar`, 'W/m²', areaMean(mesh, radiation.meanAbsorbedSolar, mask));
    row(`${label} outgoing longwave`, 'W/m²', areaMean(mesh, radiation.meanOutgoingLongwave, mask));
    row(`${label} shortwave cloud effect`, 'W/m²', areaMean(mesh, radiation.meanShortwaveCloudEffect, mask));
    row(`${label} longwave cloud effect`, 'W/m²', areaMean(mesh, radiation.meanLongwaveCloudEffect, mask));
    row(`${label} precipitation`, 'mm/day', mmDay(areaMean(mesh, precip, mask)));
    row(`${label} convective precipitation`, 'mm/day', mmDay(areaMean(mesh, convective, mask)));
    row(`${label} evaporation`, 'mm/day', mmDay(areaMean(mesh, mean.evaporation, mask)));
    row(`${label} sensible heat`, 'W/m²', areaMean(mesh, mean.sensible, mask));
    row(`${label} water vapor path`, 'kg/m²', areaMean(mesh, mean.vapor, mask));
    row(`${label} cloud water path`, 'g/m²', 1000 * areaMean(mesh, mean.cloud, mask));
    row(`${label} low cloud cover`, 'fraction', areaMean(mesh, mean.lowCover, mask));
    row(`${label} surface wind speed`, 'm/s', areaMean(mesh, mean.wind, mask));
  }
  row('land surface temperature', '°C', areaMean(mesh, mean.surfaceT, landOnly) - C2K);
  row('sea-cell surface temperature', '°C', areaMean(mesh, mean.surfaceT, sea) - C2K);
  row('30S-30N land precipitation', 'mm/day', mmDay(areaMean(mesh, precip, bandMask(mesh, -30, 30, landOnly))));
  row('30S-30N sea precipitation', 'mm/day', mmDay(areaMean(mesh, precip, bandMask(mesh, -30, 30, sea))));
  row('Arctic sea-ice area', 'Mkm²', hemisphereArea(mesh, mean.ice, sea, true));
  row('Antarctic sea-ice area', 'Mkm²', hemisphereArea(mesh, mean.ice, sea, false));
  row('open-land soil fill', 'fraction', areaMean(mesh, fill, openLand));
  row('global net TOA', 'W/m²', areaMean(mesh, radiation.meanAbsorbedSolar) - areaMean(mesh, radiation.meanOutgoingLongwave));
  row('global planetary albedo', 'fraction', total('reflectedSolar') / total('insolation'));
  row('largest wind', 'm/s', state[2].reduce((m, x) => Math.max(m, Math.abs(x)), 0));
  const { hash, groups, groupHashes } = stateHash(model);
  return { N: 16, steps, dt: DT, hours: steps * DT / 3600, node: process.version, commit: commitOf(), data: dataOf(), stateHash: hash, hashed: groups, groupHashes, rows };
}

if (process.argv[1] && realpathSync(process.argv[1]) === realpathSync(fileURLToPath(import.meta.url))) {
  const args = process.argv.slice(2), given = (name) => { const arg = args.find((a) => a.startsWith(`--${name}=`)); return arg ? arg.slice(name.length + 3) : null; };
  if (args.some((a) => !/^--(steps|out|compare)=/.test(a))) throw new Error('usage: node scripts/physicsDigest.mjs [--steps=16] [--out=digest.json] [--compare=earlier.json]');
  const t0 = performance.now(), digest = await runDigest({ steps: Number(given('steps') ?? 16) });
  console.log(`default-physics digest: N = 16, ${digest.steps} steps of ${digest.dt} s (${digest.hours} h) from a fresh start, ${((performance.now() - t0) / 1000).toFixed(0)} s; state ${digest.stateHash}`);
  const width = Math.max(...digest.rows.map((r) => r.name.length));
  for (const r of digest.rows) console.log(`  ${r.name.padEnd(width)}  ${formatValue(r.value).padStart(10)}  ${r.units}`);
  if (given('out')) writeFileSync(given('out'), JSON.stringify(digest, null, 1) + '\n');
  if (given('compare')) {
    const earlier = JSON.parse(readFileSync(given('compare'), 'utf8'));
    if (earlier.steps !== digest.steps) throw new Error(`${given('compare')} ran ${earlier.steps} steps, this digest ${digest.steps}`);
    const { moved, added, removed, sameState } = compareDigests(earlier, digest);
    console.log(`against ${given('compare')} (state ${earlier.stateHash}): ${sameState ? 'the same state' : 'the state differs'}, ${moved.length} of ${digest.rows.length} rows moved`);
    for (const m of moved) console.log(`  ${m.name.padEnd(width)}  ${formatValue(m.before).padStart(10)} → ${formatValue(m.after).padStart(10)}  ${m.units} (${m.delta >= 0 ? '+' : ''}${formatValue(m.delta)}, ${(100 * m.relative).toFixed(3)} %)`);
    if (added.length) console.log(`  new rows: ${added.join(', ')}`);
    if (removed.length) console.log(`  rows no longer made: ${removed.join(', ')}`);
  }
}
