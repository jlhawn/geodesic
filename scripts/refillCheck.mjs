// One CPU day from a saved state on the grid it was made on, judging the
// mixed layer's refill and the coasts where the 8- and 30-pass grids
// disagree about land (the cells either grid calls land and the other sea,
// and the sea cells beside them):
//   [WORKERS=6] [DAYS=1] [RELAX=passes] [OCEAN='{"refillRate":0.00027778}'] node scripts/refillCheck.mjs [runs/eleven64_day1825.bin]
// RELAX puts the state on another grid than its own, for comparison, and
// OCEAN sets the ocean's options. It passes when the refill takes no more
// than 2 m from a coastal column in one ocean step (so no layer loses more
// to it; each layer's largest loss to every process is reported beside
// it), the cells refilled at 5 m/day or more cover at most 0.2 % of the
// sea, the refill takes −4.9 ± 0.5 W/m² from the mixed layer of the sea
// deeper than 75 m, no edge of a coastal cell is clamped at the speed
// limit, and, under refillRate, no mixed layer stays at refillFloor for
// more than 2 days running. Each day it reports the sea thinned below its
// floor, the cells held at refillFloor and the Gulf of Cádiz. Exits 1 when
// any fails.
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createParallelModel } from '../js/parallel.module.js';
import { withCadence } from '../js/cadence.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { topographyFromInt16, createGeography, rebalanceSurfacePressure } from '../js/geography.module.js';
import { savedDeckField, DECK_FIELDS, savedMoistField, MOIST_FIELDS, savedRadiationField, RADIATION_FIELDS } from '../js/physics/regrid.module.js';
import { SPEED_LIMIT, EPS } from '../js/ocean/layered.module.js';
import { stateGrid, gridMesh, STATE_RELAX } from '../js/stateGrid.module.js';

const FILE = process.argv[2] ?? new URL('../runs/eleven64_day1825.bin', import.meta.url).pathname;
const WORKERS = Number(process.env.WORKERS ?? 6), DAYS = Number(process.env.DAYS ?? 1);
const LOSS = 2, FAST = 5, FAST_SHARE = 0.002, HEAT = -4.9, HEAT_TOLERANCE = 0.5, SHELF = 75, HELD_DAYS = 2, AT_FLOOR = 0.01;
const OCEAN = process.env.OCEAN ? JSON.parse(process.env.OCEAN) : {};
const CADIZ = { south: 35, north: 37.5, west: -10, east: -6 };
const t0 = performance.now(), seconds = () => ((performance.now() - t0) / 1000).toFixed(0);

const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const N = saved.N, STEP = 1350 * 16 / N, found = process.env.RELAX ? { mesh: gridMesh(new Grid(N, { relax: Number(process.env.RELAX) })), relax: Number(process.env.RELAX), source: 'RELAX', mismatches: null } : stateGrid(saved, topography);
const options = { topography, terrain: true, levels: savedLevels(saved), ...(Object.keys(OCEAN).length ? { ocean: OCEAN } : {}) };
Object.assign(options, withCadence(options, STEP, N));
const model = await createParallelModel(found.mesh, options, WORKERS);
const { mesh, ocean, state } = model, C = mesh.nCells, E = mesh.nEdges, L = ocean.layers;

['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => state[a].set(saved[name]));
model.seaIce.load(state[6], saved.concentration ?? null);
for (const field of Object.keys(DECK_FIELDS)) model.radiation[field].set(savedDeckField(saved, field, model));
for (const field of Object.keys(MOIST_FIELDS)) model.moist[field].set(savedMoistField(saved, field, model));
for (const field of Object.keys(RADIATION_FIELDS)) model.radiation[field].set(savedRadiationField(saved, field, model));
const restore = (target, values) => { if (target && values && values.length === target.length) target.set(values); };
const bl = model.boundaryLayer;
restore(bl.depth, saved.boundaryDepth); restore(bl.mixingTop, saved.mixingTop); restore(bl.regime, saved.boundaryRegime); restore(bl.buoyancyFlux, saved.boundaryBuoyancy);
restore(model.surface.windSpeed, saved.windSpeed); restore(model.radiation.evaporation, saved.evaporation);
if (model.exchange && !model.exchange.fixed) restore(model.exchange.heat, saved.exchangeHeat);
if (model.exchange) restore(model.exchange.wind, saved.exchangeWind);
const K0 = model.moist.cumulusK0 ?? 0, K = model.core.K;
for (const name of ['cumulusCover', 'cumulusWater']) {
  const target = model.moist[name], values = saved[name];
  if (!target || !values) continue;
  target.fill(0);
  if (values.length === target.length) target.set(values); else if (values.length === (K - K0) * C) target.set(values, K0 * C);
}
restore(model.moist.subcloudVirtual, saved.subcloudVirtual);
restore(model.seaIce.snow, saved.seaSnow); restore(model.seaIce.snowAlbedo, saved.seaSnowAlbedo);
model.time = saved.time;
ocean.load(saved.ocean, state[3], state[6]);
restore(model.seaIce.oceanFlux, saved.ocean.flux);
if (model.land && saved.land) model.land.load(saved.land, state[6]);
model.dragsDue();
model.restartPrecipitation();
model.core.diagnose(state[0], state[1], null, null);
rebalanceSurfacePressure(model.core, state[0], state[1], model.surfaceGeopotential, model.surfaceGeopotential);

const { cellsOnCell, nEdgesOnCell, edgesOnCell, maxEdges, areaCell, latCell, lonCell } = mesh, deg = 180 / Math.PI;
const masks = STATE_RELAX.map((relax) => createGeography(gridMesh(new Grid(N, { relax })), topography).land);
const disputed = Uint8Array.from({ length: C }, (_, i) => (masks.some((land) => !land[i] !== !masks[0][i]) ? 1 : 0));
const coast = [];
for (let i = 0; i < C; i++) {
  if (!ocean.cellOcean[i]) continue;
  let near = disputed[i];
  for (let m = 0; m < nEdgesOnCell[i]; m++) if (disputed[cellsOnCell[maxEdges * i + m]]) near = 1;
  if (near) coast.push(i);
}
const coastEdges = [...new Set(coast.flatMap((i) => Array.from({ length: nEdgesOnCell[i] }, (_, m) => edgesOnCell[maxEdges * i + m])))].filter((e) => ocean.edgeOcean[e]);
const longitude = (i) => ((((lonCell[i] * deg) + 540) % 360) - 180);
const where = (i) => `${i} (${(latCell[i] * deg).toFixed(1)}, ${longitude(i).toFixed(1)})`;
const minimumThickness = OCEAN.minimumThickness ?? 50, capped = (OCEAN.refillRate ?? 0) > 0, refillFloor = Math.min(minimumThickness, OCEAN.refillFloor ?? 5);
const cadiz = [];
for (let i = 0; i < C; i++) if (ocean.cellOcean[i] && latCell[i] * deg >= CADIZ.south && latCell[i] * deg <= CADIZ.north && longitude(i) >= CADIZ.west && longitude(i) <= CADIZ.east) cadiz.push(i);
const held = new Float64Array(C), heldLongest = new Float64Array(C), thin = new Float64Array(C), thinLongest = new Float64Array(C), shallowest = new Float64Array(C).fill(Infinity);
console.log(`${FILE}: N=${N} on ${found.relax} Lloyd passes (${found.source}${found.mismatches === null ? '' : `, ${found.mismatches} cells off the saved coast`}); ${coast.length} coastal sea cells beside the ${disputed.reduce((s, x) => s + x, 0)} cells the grids disagree on; setup ${seconds()} s`);

const column = (i, out) => { for (let k = 0; k < L; k++) out[k] = ocean.h[k * C + i]; return out; };
let before = coast.map((i) => column(i, new Float64Array(L)));
let oceanSteps = 0, clamps = 0, clampAt = null;
const losses = [{ worst: 0, at: null, over: 0 }, { worst: 0, at: null, over: 0 }], refillStep = { worst: 0, at: null, over: 0 };
const sstStart = coast.map((i) => ocean.T0[i]);
ocean.refilled.water.fill(0); ocean.refilled.heat.fill(0);
let refilledBefore = coast.map(() => 0);
const steps = Math.round(DAYS * 86400 / STEP), stepsPerDay = Math.round(86400 / STEP);
for (let n = 0; n < steps; n++) {
  model.step(STEP);
  for (let i = 0; i < C; i++) {
    if (!ocean.cellOcean[i]) continue;
    const h0 = ocean.h[i];
    let below = 0;
    for (let k = 1; k < L; k++) below += ocean.h[k * C + i] - EPS;
    shallowest[i] = Math.min(shallowest[i], h0);
    held[i] = capped && h0 <= refillFloor + AT_FLOOR ? held[i] + STEP : 0;
    heldLongest[i] = Math.max(heldLongest[i], held[i]);
    thin[i] = h0 < Math.min(minimumThickness, h0 + below) - AT_FLOOR ? thin[i] + STEP : 0;
    thinLongest[i] = Math.max(thinLongest[i], thin[i]);
  }
  if ((n + 1) % stepsPerDay === 0) {
    let atFloor = 0, lowest = -1;
    for (let i = 0; i < C; i++) if (ocean.cellOcean[i]) { if (held[i] > 0) atFloor++; if (lowest < 0 || ocean.h[i] < ocean.h[lowest]) lowest = i; }
    const cadizNow = cadiz.map((i) => ocean.h[i].toFixed(1)).join(', ');
    console.log(`day ${((n + 1) * STEP / 86400).toFixed(0)}: thinned ${(100 * ocean.diagnostics().oceanThinned).toFixed(4)} % of the sea, ${atFloor} cells at refillFloor, thinnest mixed layer ${ocean.h[lowest].toFixed(2)} m at ${where(lowest)}; Cádiz mixed layers ${cadizNow || 'none'} m; ${seconds()} s`);
  }
  const after = coast.map((i) => column(i, new Float64Array(L)));
  if (!after.some((col, c) => col.some((x, k) => x !== before[c][k]))) continue;
  oceanSteps++;
  after.forEach((col, c) => col.forEach((x, k) => {
    const loss = before[c][k] - x, kind = losses[k === 0 ? 0 : 1];
    if (loss > kind.worst) { kind.worst = loss; kind.at = `layer ${k} of ${where(coast[c])}`; }
    if (loss > LOSS) kind.over++;
  }));
  coast.forEach((i, c) => {
    const taken = ocean.refilled.water[i] - refilledBefore[c];
    if (taken > refillStep.worst) { refillStep.worst = taken; refillStep.at = where(i); }
    if (taken > LOSS) refillStep.over++;
  });
  refilledBefore = coast.map((i) => ocean.refilled.water[i]);
  for (const e of coastEdges) for (let k = 0; k < L; k++) if (Math.abs(ocean.u[k * E + e]) >= SPEED_LIMIT) { clamps++; clampAt ??= `edge ${e}, layer ${k}`; }
  before = after;
}
const days = steps * STEP / 86400, totalSeconds = steps * STEP;

let seaArea = 0, heat = 0, water = 0, fastArea = 0, fastCells = 0, fastHeat = 0;
for (let i = 0; i < C; i++) {
  if (!ocean.cellOcean[i] || ocean.D[i] < SHELF) continue;
  const a = areaCell[i], rate = ocean.refilled.water[i] / days;
  seaArea += a; heat += a * ocean.refilled.heat[i]; water += a * ocean.refilled.water[i];
  if (rate >= FAST) { fastArea += a; fastCells++; fastHeat += a * ocean.refilled.heat[i]; }
}
const heatFlux = heat / seaArea / totalSeconds, fastShare = fastArea / seaArea;
const coastRates = coast.map((i) => ocean.refilled.water[i] / days), sstEnd = coast.map((i) => ocean.T0[i]);
const hottest = (values) => { let best = 0; values.forEach((x, c) => { if (x > values[best]) best = c; }); return `${(values[best] - 273.15).toFixed(1)} °C at ${where(coast[best])}`; };
const { oceanThinned, oceanLimited } = ocean.diagnostics();
console.log(`${days.toFixed(2)} days, ${oceanSteps} ocean steps, ${seconds()} s`);
console.log(`refill over the sea deeper than ${SHELF} m: ${(water / seaArea / days).toFixed(3)} m/day, ${heatFlux.toFixed(2)} W/m²; ${fastCells} cells at ${FAST} m/day or more cover ${(100 * fastShare).toFixed(3)} % of it and carry ${(100 * fastHeat / heat).toFixed(1)} % of the heat`);
console.log(`coastal cells: refill up to ${Math.max(...coastRates).toFixed(2)} m/day; SST hottest ${hottest(sstStart)} at the start, ${hottest(sstEnd)} at the end`);
console.log(`thinned at the end: ${(100 * oceanThinned).toFixed(4)} % of the sea; clamped edges in the last ocean step: ${oceanLimited}`);
console.log(`largest loss of a coastal layer in one ocean step, to every process: mixed layer ${losses[0].worst.toFixed(2)} m (${losses[0].at ?? 'none'}; ${losses[0].over} cell-steps over ${LOSS} m), interior ${losses[1].worst.toFixed(2)} m (${losses[1].at ?? 'none'}; ${losses[1].over} layer-steps over ${LOSS} m)`);
const longest = (values) => { let best = -1; for (let i = 0; i < C; i++) if (ocean.cellOcean[i] && (best < 0 || values[i] > values[best])) best = i; return best; };
const heldWorst = longest(heldLongest), thinWorst = longest(thinLongest), heldCells = heldLongest.reduce((s, x) => s + (x > 0 ? 1 : 0), 0);
let thinCells = 0, thinArea = 0, totalSea = 0;
for (let i = 0; i < C; i++) {
  if (ocean.cellOcean[i]) totalSea += areaCell[i];
  if (thinLongest[i] > HELD_DAYS * 86400) { thinCells++; thinArea += areaCell[i]; }
}
console.log(`under ${capped ? `refillRate ${(OCEAN.refillRate * 86400).toFixed(2)} m/day, refillFloor ${refillFloor} m` : 'the instant refill'}: ${heldCells} cells reached refillFloor, the longest held ${(heldLongest[heldWorst] / 86400).toFixed(2)} days at ${where(heldWorst)}; the longest thinned below the ${minimumThickness} m floor ${(thinLongest[thinWorst] / 86400).toFixed(2)} days at ${where(thinWorst)} (shallowest ${shallowest[thinWorst].toFixed(2)} m), ${thinCells} cells (${(100 * thinArea / totalSea).toFixed(4)} % of the sea) thinned for more than ${HELD_DAYS} days`);
for (const i of cadiz) console.log(`Cádiz ${where(i)}: depth ${ocean.D[i].toFixed(0)} m, mixed layer ${ocean.h[i].toFixed(2)} m at the end, ${shallowest[i].toFixed(2)} at its shallowest, refilled ${(ocean.refilled.water[i] / days).toFixed(2)} m/day, thinned at most ${(thinLongest[i] / 86400).toFixed(2)} days running, at refillFloor ${(heldLongest[i] / 86400).toFixed(2)}`);
const checks = [
  [`the refill takes no more than ${LOSS} m from a coastal column in an ocean step (worst ${refillStep.worst.toFixed(2)} m, ${refillStep.at ?? 'none'}; ${refillStep.over} cell-steps over)`, refillStep.worst <= LOSS],
  [`cells refilled at ${FAST} m/day or more cover at most ${100 * FAST_SHARE} % of the sea (${(100 * fastShare).toFixed(3)} %)`, fastShare <= FAST_SHARE],
  [`refill heat ${HEAT} ± ${HEAT_TOLERANCE} W/m² (${heatFlux.toFixed(2)})`, Math.abs(heatFlux - HEAT) <= HEAT_TOLERANCE],
  [`no clamped edge at those coasts (${clamps} layer-edge-steps${clampAt ? `, first ${clampAt}` : ''})`, clamps === 0],
  ...(capped ? [[`no mixed layer held at refillFloor for more than ${HELD_DAYS} days running (longest ${(heldLongest[heldWorst] / 86400).toFixed(2)}, ${where(heldWorst)}; Cádiz ${(Math.max(0, ...cadiz.map((i) => heldLongest[i])) / 86400).toFixed(2)})`, heldLongest[heldWorst] <= HELD_DAYS * 86400]] : []),
];
for (const [text, ok] of checks) console.log(`${ok ? 'PASS' : 'FAIL'} ${text}`);
await model.close();
process.exit(checks.every(([, ok]) => ok) ? 0 : 1);
