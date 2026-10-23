import { test } from 'node:test';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { mkdtempSync, readFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { readTopography, geographyOptions, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { withCadence } from '../js/cadence.module.js';
import { CLIMATOLOGY_FILE, loadClimatology } from '../js/ocean/climatology.module.js';
import { decodeState } from '../js/stateFile.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }

const N = 16, DT = 1350 * 16 / N, LEVELS = 'bl60', GPU_DAYS = 3;
const nonFinite = (fields) => Object.entries(fields).filter(([, a]) => ArrayBuffer.isView(a) && a.some((v) => !Number.isFinite(v))).map(([name]) => name);

test(`a fresh start on the default geography at N=${N} on ${LEVELS} stays finite for a day on the CPU`, async () => {
  const topography = readTopography(FRESH_GEOGRAPHY.topography), sheets = geographyOptions(FRESH_GEOGRAPHY);
  const ocean = withCadence({ ocean: {} }, DT, N).ocean, radiation = withCadence({ radiation: { clearSkyPass: true } }, DT, N).radiation;
  const climatology = await loadClimatology(new URL(`../${CLIMATOLOGY_FILE}`, import.meta.url).pathname);
  const model = createModel(new Grid(N), { topography, geography: sheets.geography, ocean: { ...ocean, climatology }, radiation, land: { ...sheets.land }, levels: sigmaInterfaces(LEVELS) });
  const { state } = model, C = model.mesh.nCells;
  initializeState(model, { geostrophic: !model.surfaceGeopotential }).forEach((values, a) => state[a].set(values));
  for (let i = 0; i < C; i++) if (model.geography.land[i]) state[6][i] = 0;
  if (model.load) model.load();
  model.ocean.initialize(state[3], state[6]);
  model.land.initialize();
  for (let n = 0; n < 86400 / DT; n++) model.step(DT);
  const atmosphere = nonFinite(Object.fromEntries(['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].map((name, a) => [name, state[a]])));
  assert.deepEqual([...atmosphere, ...nonFinite({ h: model.ocean.h, u: model.ocean.u, Q: model.ocean.Q, W: model.ocean.W, eta: model.ocean.eta })], []);
});

test(`a fresh start on the default geography at N=${N} on ${LEVELS} stays finite for ${GPU_DAYS} days on the GPU, through scripts/spinup.mjs`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const out = mkdtempSync(join(tmpdir(), 'fresh-'));
  try {
    execFileSync(process.execPath, [new URL('../scripts/spinup.mjs', import.meta.url).pathname], { env: { ...process.env, N: String(N), LEVELS, TAG: 'fresh', OUT: out, MINUTES: '100', DAYS: String(GPU_DAYS) }, stdio: 'pipe' });
    const st = await decodeState(new Uint8Array(readFileSync(join(out, `fresh_day${String(GPU_DAYS).padStart(4, '0')}.bin`))));
    assert.equal(st.day, GPU_DAYS);
    assert.deepEqual(nonFinite({ pi: st.pi, theta: st.theta, u: st.u, surfaceT: st.surfaceT, q: st.q, qc: st.qc }), [], 'atmosphere');
    assert.deepEqual(nonFinite({ h: st.ocean.h, u: st.ocean.u, T: st.ocean.T, S: st.ocean.S, eta: st.ocean.eta }), [], 'ocean');
    const iced = [...st.ice.keys()].filter((i) => st.ice[i] > 0);
    assert.ok(iced.some((i) => st.concentration[i] < 0.95), 'the observed March ice keeps cells partly covered');
    assert.ok(iced.filter((i) => st.land.snow[i] > 0).length > 0.8 * iced.length, 'the observed March ice keeps its snow');
  } finally {
    rmSync(out, { recursive: true, force: true });
  }
});
