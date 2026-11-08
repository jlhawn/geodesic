import { test, after } from 'node:test';
import assert from 'node:assert/strict';
import { writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { stateGeography, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { encodeState } from '../js/stateFile.module.js';
import { gridHeader } from '../js/stateGrid.module.js';
import { oceanCells, inertialStepLimit } from '../js/ocean/layered.module.js';
import { withoutOptionGroups } from '../scripts/optionEnvironment.mjs';

const dir = mkdtempSync(join(tmpdir(), 'diagnostic-scripts-'));
after(() => rmSync(dir, { recursive: true, force: true }));

const mesh = new Grid(16), header = { N: 16, ...gridHeader(mesh), ...FRESH_GEOGRAPHY, terrain: true };
const { topography, geography, land } = stateGeography(header);
const model = createModel(mesh, { topography, geography, land });
const { state, ocean } = model, E = model.mesh.nEdges, DEG = 180 / Math.PI;
initializeState(model, {}).forEach((values, a) => state[a].set(values));
for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) state[6][i] = 0;
model.land.initialize();
ocean.initialize(state[3], state[6]);
const edgeAt = (e) => {
  const [x, y, z] = model.mesh.xEdge.subarray(3 * e, 3 * e + 3);
  return { lat: Math.asin(z / Math.hypot(x, y, z)) * DEG, lon: ((Math.atan2(y, x) * DEG) % 360 + 360) % 360 };
};
for (let k = 0; k < ocean.layers; k++) {
  for (let e = 0; e < E; e++) if (ocean.edgeOcean[e]) { const { lat, lon } = edgeAt(e); ocean.u[k * E + e] = Math.cos(lat / DEG) * Math.sin(3 * lon / DEG + k); }
}
const [pi, theta, u, surfaceT, q, qc, ice] = state;
const file = join(dir, 'state16_day0000.bin');
writeFileSync(file, encodeState({ ...header, K: model.core.K, levels: model.core.levels, day: 0, time: 0, pi, theta, u, surfaceT, q, qc, ice, concentration: model.seaIce.concentration, ocean: await ocean.serialize(), land: await model.land.serialize() }));
const run = (script, env) => spawnSync(process.execPath, [new URL(`../scripts/${script}.mjs`, import.meta.url).pathname, file], { env: { ...withoutOptionGroups(), ...env }, encoding: 'utf8', maxBuffer: 1 << 28 });

test('edgeBudget splits the ocean\'s momentum tendency into the terms the hybrid computes: on a state at N=16 they add up to the tendency on every row of every ocean edge, a pinned row\'s relaxation at the run\'s rate included, with the fixed layers shown and the wind, the bottom\'s and the shelves\' drags at work', () => {
  const wet = [], open = [], shelfEdges = [], { shelf } = oceanCells(model.mesh.nCells, model.geography);
  for (let e = 0; e < E; e++) {
    if (!ocean.edgeOcean[e]) continue;
    const { lat, lon } = edgeAt(e), [a, b] = [model.mesh.cellsOnEdge[2 * e], model.mesh.cellsOnEdge[2 * e + 1]];
    wet.push(e);
    if (open.length < 3 && !model.geography.land[a] && !model.geography.land[b] && Math.abs(lat) < 10 && lon > 200 && lon < 250) open.push(e);
    if (shelfEdges.length < 2 && shelf[a] + shelf[b] > 0) shelfEdges.push(e);
  }
  assert.equal(open.length, 3);
  assert.equal(shelfEdges.length, 2);
  const budget = run('edgeBudget', { EDGES: wet.join(',') });
  assert.equal(budget.status, 0, budget.stderr);
  assert.doesNotMatch(budget.stdout, /NaN|Infinity/);
  const dtOcean = ocean.everySteps * 1350, part = dtOcean / Math.max(1, Math.ceil(dtOcean / inertialStepLimit(model.mesh)));
  assert.match(budget.stdout, new RegExp(`pinned rows relax at 1/${Math.max(3600, part).toFixed(0)} s`));
  const rows = budget.stdout.split('\n').filter((line) => /^edge \d+ layer \d+ /.test(line)), pinned = (line) => / pinned,/.test(line);
  const term = (line, name) => Number(new RegExp(`${name} (\\S+)`).exec(line)[1]), key = (line) => /^edge (\d+) layer (\d+) /.exec(line).slice(1).join(':');
  assert.equal(rows.find((line) => term(line, 'rest') !== 0), undefined);
  const active = new Map(rows.filter((line) => !pinned(line)).map((line) => [key(line), line]));
  for (const e of open) for (let k = 0; k < ocean.table.fixed; k++) assert.ok(active.has(`${e}:${k}`), `edge ${e} layer ${k} is not shown as an active row`);
  for (const e of shelfEdges) assert.notEqual(term(active.get(`${e}:0`), 'shelf drag'), 0, active.get(`${e}:0`));
  assert.ok(rows.some(pinned), 'no row is pinned');
  for (const name of ['Coriolis', 'baroclinic', '-grad K', 'bottom drag', 'nu4']) assert.ok([...active.values()].some((line) => term(line, name) !== 0), `${name} is zero on every row`);
  assert.ok(open.some((e) => term(active.get(`${e}:0`), 'stress') !== 0), 'the wind is zero on the open edges');
});

test('radiationProfiles runs the defaults against OTHER\'s gray longwave and Lacis–Hansen gases on the sequential shortwave, on a state at N=16', () => {
  const profiles = run('radiationProfiles', { TIMES: '2' });
  assert.equal(profiles.status, 0, profiles.stderr);
  assert.match(profiles.stdout, /second \{"longwaveScheme":"gray","solarGases":"lacisHansen","aerosol":"flat","shortwave":"sequential"\}/);
  const rows = profiles.stdout.split('\n').filter((line) => /^\s+\d+(\s+-?\d+\.\d+){5}\s*$/.test(line));
  assert.ok(rows.length > 0, profiles.stdout.slice(0, 400));
});
