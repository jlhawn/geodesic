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
import { oceanCells } from '../js/ocean/layered.module.js';

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
const run = (script, env) => spawnSync(process.execPath, [new URL(`../scripts/${script}.mjs`, import.meta.url).pathname, file], { env: { ...process.env, ...env }, encoding: 'utf8' });

test('edgeBudget splits the ocean\'s momentum tendency into the terms the hybrid computes: on a state at N=16 they add up to the tendency in each of the top 20 layers of open equatorial edges, far above the floor, where no drag acts, and in the top layer of edges under an ice shelf, where the shelf\'s drag does', () => {
  const edges = [], shelfEdges = [], { shelf } = oceanCells(model.mesh.nCells, model.geography);
  for (let e = 0; e < E; e++) {
    const { lat, lon } = edgeAt(e), [a, b] = [model.mesh.cellsOnEdge[2 * e], model.mesh.cellsOnEdge[2 * e + 1]];
    if (edges.length < 3 && ocean.edgeOcean[e] && !model.geography.land[a] && !model.geography.land[b] && Math.abs(lat) < 10 && lon > 200 && lon < 250) edges.push(e);
    if (shelfEdges.length < 2 && ocean.edgeOcean[e] && shelf[a] + shelf[b] > 0) shelfEdges.push(e);
  }
  assert.equal(edges.length, 3);
  assert.equal(shelfEdges.length, 2);
  const budget = run('edgeBudget', { EDGES: [...edges, ...shelfEdges].join(',') });
  assert.equal(budget.status, 0, budget.stderr);
  assert.doesNotMatch(budget.stdout, /NaN|Infinity/);
  const lines = budget.stdout.split('\n').filter((line) => /^edge \d+ layer \d+ /.test(line)), of = (e) => lines.filter((line) => line.startsWith(`edge ${e} `));
  const term = (line, name) => Number(new RegExp(`${name} (\\S+)`).exec(line)[1]);
  for (const e of shelfEdges) {
    const [top] = of(e);
    assert.match(top, / layer 0 .* shelf drag -?\d/);
    assert.notEqual(term(top, 'shelf drag'), 0, top);
    assert.equal(Math.abs(term(top, 'drag\\+rest')), 0, top);
  }
  const upper = edges.flatMap(of).filter((line) => Number(/ layer (\d+) /.exec(line)[1]) < 20);
  for (const e of edges) assert.ok(upper.filter((line) => line.startsWith(`edge ${e} `)).length >= ocean.table.fixed, `edge ${e} shows ${upper.filter((line) => line.startsWith(`edge ${e} `)).length} of the top 20 layers`);
  for (const line of upper) assert.equal(Math.abs(term(line, 'drag\\+rest')), 0, line);
  for (const name of ['Coriolis', 'baroclinic', '-grad K']) assert.ok(upper.some((line) => term(line, name) !== 0), `${name} is zero on every edge`);
  assert.ok(upper.some((line) => / layer 0 /.test(line) && term(line, 'stress') !== 0), 'the wind is zero on every edge');
});

test('radiationProfiles runs the defaults against OTHER\'s gray longwave and Lacis–Hansen gases on the sequential shortwave, on a state at N=16', () => {
  const profiles = run('radiationProfiles', { TIMES: '2' });
  assert.equal(profiles.status, 0, profiles.stderr);
  assert.match(profiles.stdout, /second \{"longwaveScheme":"gray","solarGases":"lacisHansen","aerosol":"flat","shortwave":"sequential"\}/);
  const rows = profiles.stdout.split('\n').filter((line) => /^\s+\d+(\s+-?\d+\.\d+){5}\s*$/.test(line));
  assert.ok(rows.length > 0, profiles.stdout.slice(0, 400));
});
