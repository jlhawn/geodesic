import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createParallelModel } from '../js/parallel.module.js';
import { initializeState } from '../js/physics/init.module.js';

test('serial and parallel engines stay bit-identical with the boundary layer', async () => {
  const serial = createModel(new Grid(16));
  const parallel = await createParallelModel(new Grid(16), {}, 3);
  try {
    const init = initializeState(serial, {});
    for (let a = 0; a < init.length; a++) { serial.state[a].set(init[a]); parallel.state[a].set(init[a]); }
    serial.ocean.initialize(serial.state[3], serial.state[6]);
    parallel.ocean.initialize(parallel.state[3], parallel.state[6]);
    for (let n = 0; n < 6; n++) { serial.step(900); parallel.step(900); }
    for (let a = 0; a < serial.state.length; a++) for (let x = 0; x < serial.state[a].length; x++) assert.equal(parallel.state[a][x], serial.state[a][x], `state ${a}[${x}]`);
  } finally { await parallel.close(); }
});

test('serial and parallel engines stay bit-identical with the mixed-layer deck carrying its height and gate', async () => {
  const radiation = { stratusSubsidence: 0, minimumInversion: 0 };
  const serial = createModel(new Grid(16), { radiation });
  const parallel = await createParallelModel(new Grid(16), { radiation }, 3);
  try {
    const init = initializeState(serial, {});
    for (let a = 0; a < init.length; a++) { serial.state[a].set(init[a]); parallel.state[a].set(init[a]); }
    serial.ocean.initialize(serial.state[3], serial.state[6]);
    parallel.ocean.initialize(parallel.state[3], parallel.state[6]);
    for (let n = 0; n < 6; n++) { serial.step(900); parallel.step(900); }
    let decked = 0;
    for (let i = 0; i < serial.mesh.nCells; i++) if (serial.radiation.mlmTop[i] > 0) decked++;
    assert.ok(decked > 0, 'the deck runs and hands its height to the boundary layer');
    for (let a = 0; a < serial.state.length; a++) for (let x = 0; x < serial.state[a].length; x++) assert.equal(parallel.state[a][x], serial.state[a][x], `state ${a}[${x}]`);
    for (const name of ['mlmHeight', 'mlmGate', 'mlmTop', 'mlmCover', 'mlmWater']) assert.deepEqual(parallel.radiation[name], serial.radiation[name], name);
  } finally { await parallel.close(); }
});
