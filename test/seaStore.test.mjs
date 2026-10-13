import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { openWaterAlbedo } from '../js/physics/ice.module.js';
import { SOLAR_CONSTANT } from '../js/physics/radiation.module.js';
import { sunTrack } from '../explain/insolation.module.js';
import { createSeasonalColumn, slabCycle } from '../explain/oceanColumn.module.js';
import { startingSea, sunSwing, classLabels, columnOptions, createSeaStore, airWarming, RHO_CP, AIR_COLUMN, SLAB, DAMPING, START_DEPTH, TRANSMISSION } from '../explain/seaStore.module.js';

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const DEG = 180 / Math.PI;

test('the starting sea is initializeState\'s, wherever the start leaves no ice', () => {
  const model = createModel(buildMesh(new Grid(4)), { ocean: false });
  const state = initializeState(model, { profile: new Float64Array(model.core.diagnostics.K).fill(300) });
  const surfaceT = state[3], ice = state[6];
  let checked = 0;
  for (let i = 0; i < model.mesh.nCells; i++) {
    if (ice[i] > 0) continue;
    close(startingSea(model.mesh.latCell[i] * DEG).temperature, surfaceT[i], 1e-9, `sea temperature at ${(model.mesh.latCell[i] * DEG).toFixed(2)}°`);
    checked++;
  }
  assert.ok(checked > 100, `${checked} cells checked`);
});

test('the column starts as createOcean\'s analytic start leaves the sea at each latitude up to 45°: the mixed layer 60 m deep, every class at its depth, temperature and salinity', () => {
  const mesh = buildMesh(new Grid(8)), C = mesh.nCells, ocean = createOcean(mesh);
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => startingSea(lat * DEG).temperature);
  ocean.initialize(surfaceT, new Float64Array(C));
  const L = ocean.layers;
  for (const want of [0, 15, 30, 40, 44.9]) {
    let i = 0;
    for (let j = 0; j < C; j++) if (mesh.latCell[j] * DEG <= 45 && Math.abs(mesh.latCell[j] * DEG - want) < Math.abs(mesh.latCell[i] * DEG - want)) i = j;
    const latitude = mesh.latCell[i] * DEG, column = createSeasonalColumn({ ...columnOptions({ latitude }), days: 0 }), o = column.ocean, c = column.cell, Cc = o.h.length / L;
    close(o.h[c], START_DEPTH, 1e-12, 'the mixed layer');
    close(o.S0[c], ocean.S0[i], 1e-12, `the mixed layer's salinity at ${latitude.toFixed(2)}°`);
    let a = 0, b = 0;
    for (let k = 0; k < L - 1 && a < 3000; k++) {
      a += ocean.h[k * C + i]; b += o.h[k * Cc + c];
      close(b, a, 1e-9, `the base of layer ${k} at ${latitude.toFixed(2)}°`);
      if (ocean.h[k * C + i] > 0.02) {
        close(o.Q[k * Cc + c] / o.h[k * Cc + c], ocean.Q[k * C + i] / ocean.h[k * C + i], 1e-9, `layer ${k}'s temperature`);
        close(o.W[k * Cc + c] / o.h[k * Cc + c], ocean.W[k * C + i] / ocean.h[k * C + i], 1e-9, `layer ${k}'s salinity`);
      }
    }
  }
  assert.throws(() => columnOptions({ latitude: 50 }), /45/);
});

test('the labels are createOcean\'s class temperatures', () => {
  const labels = classLabels();
  assert.equal(labels.length, 46);
  close(labels[19] - 273.15, 12.592, 5e-4, 'the 1025.5 class');
  close(labels[0] - 273.15, 29.394, 5e-4, 'the 1019.5 class');
  close(labels[2] - 273.15, 28.592, 5e-4, 'the 1020.5 class');
  close(labels[42] - 273.15, -1.664, 5e-4, 'the 1027.3 class');
});

test('the heating swings with the yearly harmonic of the sunlight a sea absorbs under a sky that passes 70 percent: ±133 W/m² at 45°N, ±95 at 30°, ±49 at 15°, none on the equator', () => {
  assert.equal(TRANSMISSION, 0.7);
  for (const [latitude, swing] of [[45, 132.8], [30, 94.6], [15, 48.9], [0, 0]]) {
    let sum = 0;
    for (let d = 0; d < 365; d++) {
      const { x, z } = sunTrack(d, 288);
      let lit = 0, reflected = 0;
      for (let k = 0; k < 288; k++) { const mu = Math.cos(latitude / DEG) * x[k] + Math.sin(latitude / DEG) * z[k]; if (mu > 0) { lit += mu; reflected += mu * openWaterAlbedo(mu); } }
      sum += 0.7 * SOLAR_CONSTANT * (lit - reflected) / 288 * Math.sin(2 * Math.PI * (d + 0.5) / 365);
    }
    close(sunSwing(latitude), 2 * sum / 365, 1e-9, `the swing at ${latitude}°, summed directly`);
    close(sunSwing(latitude), swing, 0.05, `the swing at ${latitude}°`);
  }
});

test('a year at 45°N under ±133 W/m² and 0.1 N/m²: the wind holds the mixed layer 57 to 69 m deep through the spring and summer, above its 50 m floor and near 45°N\'s Ekman bound, it deepens past 300 m in winter and retreats at the March equinox, its surface swinging 4.0 K against the slab\'s 20.8', () => {
  const t0 = performance.now(), store = createSeaStore({ latitude: 45, heating: 133, wind: 0.1 });
  const batches = [];
  for (let n = 0; !store.done; n++) { store.run(97); batches.push(store.batch()); }
  const seconds = (performance.now() - t0) / 1000, r = store.column.record, N = store.column.samples - 1;
  assert.equal(N, 365);
  const joined = (name) => batches.flatMap((b) => Array.from(b[name]));
  for (const [name, values] of [['time', r.time], ['sst', r.seas.model.sst], ['slab', r.seas.slab.sst], ['depth', r.mixedDepth], ['stored', r.stored]]) assert.deepEqual(joined(name), Array.from(values), `the batches' ${name}`);
  assert.deepEqual(joined('bases'), Array.from(r.bases), 'the batches\' layer bases');
  assert.equal(store.layout.heating, 133);
  close(store.layout.mean - 273.15, 13, 1e-9, 'the mean temperature');

  let low = Infinity, high = 0;
  for (let n = 21; n <= 151; n++) { low = Math.min(low, r.mixedDepth[n]); high = Math.max(high, r.mixedDepth[n]); }
  close(low, 56.60, 0.05, 'the shallowest from day 21 to 151');
  close(high, 68.63, 0.05, 'the deepest from day 21 to 151');
  close(r.mixedDepth[300], 136.3, 2, 'day 300');
  close(r.mixedDepth[330], 272.5, 2, 'day 330');
  close(r.mixedDepth[363], 313.1, 2, 'day 363');
  assert.ok(r.mixedDepth[365] < 0.6 * r.mixedDepth[363], `the retreat: ${r.mixedDepth[363].toFixed(1)} → ${r.mixedDepth[365].toFixed(1)} m`);

  const model = store.column.readouts(N, 'model'), slab = store.column.readouts(N, 'slab');
  close(model.range, 3.96, 0.02, 'the mixed layer\'s swing');
  close(slab.range, 20.85, 0.02, 'the slab\'s swing');
  assert.equal(model.warmestDay, 158);
  assert.equal(slab.warmestDay, 129);
  close(slab.range, slabCycle({ capacity: SLAB, heating: 133, damping: DAMPING }).range, 0.03 * slab.range, 'the slab against its linear answer');

  let most = 0;
  for (let n = 0; n <= N; n++) if (r.stored[n] > r.stored[most]) most = n;
  close(r.stored[most] / 1e9, 1.063, 0.003, 'the most heat in the top 100 m');
  close(airWarming(r.stored[most]), 103.9, 0.3, 'as kelvins of the air column');
  close(AIR_COLUMN / RHO_CP, 2.505, 0.001, 'the air column in meters of sea');
  console.log(`a year in ${seconds.toFixed(1)} s; floor ${low.toFixed(2)}–${high.toFixed(2)} m on days 21–151; ${r.mixedDepth[300].toFixed(1)}, ${r.mixedDepth[330].toFixed(1)}, ${r.mixedDepth[363].toFixed(1)} → ${r.mixedDepth[365].toFixed(1)} m on days 300, 330, 363 and 365; swings ${model.range.toFixed(2)} and ${slab.range.toFixed(2)} K, warmest on days ${model.warmestDay} and ${slab.warmestDay}; the top 100 m stored ${(r.stored[most] / 1e9).toFixed(3)} GJ/m² on day ${most}, ${airWarming(r.stored[most]).toFixed(1)} K of the air column`);
});
