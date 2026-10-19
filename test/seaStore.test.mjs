import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { openWaterAlbedo } from '../js/physics/ice.module.js';
import { SOLAR_CONSTANT } from '../js/physics/radiation.module.js';
import { sunTrack } from '../explain/insolation.module.js';
import { createSeasonalColumn, columnMesh, YEAR_DAYS } from '../explain/oceanColumn.module.js';
import { startingSea, sunSwing, columnOptions, createSeaStore, airWarming, TRANSMISSION } from '../explain/seaStore.module.js';

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const DEG = 180 / Math.PI;
const C = (k) => k - 273.15;

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


test('the start is the model\'s analytic start at the latitude: the twelve surface layers hold T and S stretched by the tilt, from 2 °C and 34.7 psu below', () => {
  const tDeep = 275.15, sDeep = 34.7;
  for (const latitude of [0, 15, 30, 45]) {
    const column = createSeasonalColumn(columnOptions({ latitude, days: 0 })), o = column.ocean, n = columnMesh().nCells, sea = startingSea(latitude), scale = 400 * sea.stretch;
    const decay = (a, b) => scale * (Math.exp(-a / scale) - Math.exp(-b / scale)) / (b - a);
    let z = 0, worst = 0;
    for (let k = 0; k < 12; k++) {
      const h = o.h[k * n + column.cell], f = decay(z, z + h);
      worst = Math.max(worst, Math.abs(o.Q[k * n + column.cell] / h - (tDeep + (sea.temperature - tDeep) * f)), Math.abs(o.W[k * n + column.cell] / h - (sDeep + (sea.salinity - sDeep) * f)));
      z += h;
    }
    console.log(`${latitude}°: surface ${C(sea.temperature).toFixed(2)} °C, ${sea.salinity.toFixed(3)} psu, e-folding ${scale.toFixed(0)} m; the largest difference from the model's layers ${worst.toExponential(1)}`);
    assert.ok(worst < 1e-6, `${latitude}°: ${worst}`);
  }
  close(400 * startingSea(45).stretch, 400, 1e-9, 'e-folding at 45°');
  close(400 * startingSea(0).stretch, 520, 1e-9, 'e-folding at the equator');
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


const day = (run, field, n) => run.record[field][n];
const extremum = (values, from, to, sign) => { let best = from; for (let j = from; j <= to; j++) if (sign * values[j] > sign * values[best]) best = j; return best; };

test('a year at 45°N under 0.1 N/m²: the wind stirs the start to 53 m, spring heating shoals the mixed layer only to 41 m, the surface peaks at 16.7 °C 56 days after the heating where the slab peaks at 23.6 °C after 38, and the top 100 m stores 1.05 GJ/m²', () => {
  const run = createSeaStore({}).column.run(), r = run.record, sst = r.seas.model.sst, slab = r.seas.slab.sst;
  assert.equal(run.samples, YEAR_DAYS + 1);
  close(r.mixedDepth[0], 7.3, 0.5, 'the start\'s mixed layer');
  close(C(sst[0]), 13.0, 0.05, 'the start\'s surface');
  const early = extremum(sst, 0, 40, -1);
  close(C(sst[early]), 12.33, 0.03, 'the surface after the wind has stirred the start'); close(early, 11, 2, 'the day');
  close(r.mixedDepth[10], 53.2, 1, 'the mixed layer on day 10');
  const shallow = extremum(r.mixedDepth, 20, 200, -1);
  close(r.mixedDepth[shallow], 40.9, 0.5, 'the shallowest mixed layer'); close(shallow, 110, 25, 'the day');
  const within = []; for (let d = 20; d <= 200; d++) if (r.mixedDepth[d] < r.mixedDepth[shallow] + 1) within.push(d);
  close(within[0], 72, 2, 'the first day within a meter of the shallowest'); close(within[within.length - 1], 137, 2, 'the last'); assert.equal(within.length, within[within.length - 1] - within[0] + 1);
  const warm = extremum(sst, 0, 365, 1), slabWarm = extremum(slab, 0, 365, 1);
  close(C(sst[warm]), 16.66, 0.05, 'the warmest surface'); close(warm, 147, 2, 'the warmest day');
  close(C(slab[slabWarm]), 23.60, 0.05, 'the slab\'s warmest'); close(slabWarm, 129, 1, 'the slab\'s warmest day');
  const peak = YEAR_DAYS / 4;
  close(warm - peak, 56, 2, 'the lag of the surface'); close(slabWarm - peak, 38, 1, 'the slab\'s lag');
  const cold = extremum(sst, 150, 365, -1);
  close(C(sst[cold]), 10.68, 0.05, 'the coldest surface'); close(cold, 355, 4, 'the coldest day');
  close(Math.max(...sst) - Math.min(...sst), 5.97, 0.1, 'the surface\'s swing'); close(Math.max(...slab) - Math.min(...slab), 20.82, 0.1, 'the slab\'s swing');
  const stored = extremum(r.stored, 0, 365, 1);
  close(r.stored[stored] / 1e9, 1.05, 0.01, 'GJ/m² stored in the top 100 m'); close(stored, 167, 3, 'the day');
  close(airWarming(r.stored[stored]), 103, 1, 'K of the air column'); close(Math.max(...r.seas.slab.heat) / 1e9, 0.447, 0.003, 'the slab\'s most');
  const deepest = extremum(r.boundaryDepth, 150, 365, 1);
  close(r.boundaryDepth[deepest], 154.7, 2, 'the deepest boundary layer'); close(deepest, 347, 5, 'the day');
  close(Math.max(...r.mixedDepth), 150.7, 2, 'the deepest mixed layer');
  close(r.boundaryDepth[365], 85.8, 5, 'the boundary layer at the end'); close(r.mixedDepth[365], 150, 3, 'the mixed layer at the end');
  close(r.mixedDepth[20], 55, 1.5, 'the mixed layer on day 20'); close(r.mixedDepth[350], 150.6, 2.5, 'the mixed layer on day 350'); close(C(r.seas.model.sst[355]), 10.68, 0.05, 'the surface on day 355');
  for (const [d, depth] of [[200, 61], [250, 90], [300, 122]]) close(r.mixedDepth[d], depth, 2.5, `the mixed layer on day ${d}`);
  console.log(`45°N, 0.1 N/m²: mixed layer ${r.mixedDepth[10].toFixed(1)} m on day 10, ${r.mixedDepth[shallow].toFixed(1)} m on day ${shallow}; surface ${C(sst[early]).toFixed(2)} °C on day ${early}, peak ${C(sst[warm]).toFixed(2)} °C on day ${warm}, low ${C(sst[cold]).toFixed(2)} °C on day ${cold}; slab ${C(slab[slabWarm]).toFixed(2)} °C on day ${slabWarm}; stored ${(r.stored[stored] / 1e9).toFixed(3)} GJ/m² on day ${stored}; boundary layer ${r.boundaryDepth[deepest].toFixed(1)} m on day ${deepest}`);
  let turn = 0;
  while (r.seas.model.flux[turn] > 0 || turn < 100) turn++;
  console.log(`the surface flux turns negative on day ${turn}`);
  close(turn, 168, 3, 'the day the column starts to lose heat');
  close(C(run.temperatureAt(120, 40)), 16.3, 0.05, 'on day 120, the water at 40 m'); close(C(run.temperatureAt(120, 60)), 12.33, 0.05, 'and at 60 m'); close(C(run.temperatureAt(0, 60)), 11.4, 0.1, 'the start at 60 m');
  for (const d of [0, 120, 350]) {
    const layers = run.layersAt(d);
    assert.ok(layers.slice(0, 40).every((l) => l.mode < 0.01), `day ${d}: the top 40 layers are held at their minimum thicknesses`);
    assert.ok(layers.slice(40, 46).every((l) => l.mode > 0.99), `day ${d}: the layers from 40 follow their isopycnals`);
    close(layers[40].top, 334, 2, `day ${d}: top of layer 40`); close(layers[40].thickness, 62, 2, `day ${d}: thickness of layer 40`);
    const thick = layers.slice(40, 46).map((l) => l.thickness);
    close(Math.min(...thick), 62, 2, `day ${d}: thinnest isopycnal layer`); close(Math.max(...thick), 148, 3, `day ${d}: thickest`);
  }
});

test('with no wind the same year leaves the mixed layer in the top layers all spring and summer: the surface peaks at 20.3 °C 23 days after the heating, the top 100 m stores 0.60 GJ/m², and convection alone mixes it to 161 m by winter', () => {
  const run = createSeaStore({ wind: 0 }).column.run(), r = run.record, sst = r.seas.model.sst;
  for (let d = 20; d <= 133; d++) assert.ok(r.mixedDepth[d] < 6, `day ${d}: ${r.mixedDepth[d]}`);
  close(Math.min(...r.mixedDepth.slice(20, 200)), 1.8, 0.5, 'the shallowest mixed layer');
  const warm = extremum(sst, 0, 365, 1);
  close(C(sst[warm]), 20.32, 0.05, 'the warmest surface'); close(warm, 114, 2, 'the warmest day'); close(warm - YEAR_DAYS / 4, 23, 2, 'the lag');
  close(Math.max(...sst) - Math.min(...sst), 10.30, 0.1, 'the swing');
  const stored = extremum(r.stored, 0, 365, 1);
  close(r.stored[stored] / 1e9, 0.596, 0.01, 'GJ/m² stored'); close(stored, 155, 3, 'the day'); close(airWarming(r.stored[stored]), 58, 1, 'K of the air column');
  close(Math.max(...r.mixedDepth), 161.5, 2, 'the deepest mixed layer'); close(Math.max(...r.boundaryDepth), 163.3, 2, 'the deepest boundary layer');
  assert.ok(r.mixedDepth[360] < 6 && r.mixedDepth[365] < 6, 'the mixed layer is gone by day 360');
  console.log(`45°N, no wind: peak ${C(sst[warm]).toFixed(2)} °C on day ${warm}, stored ${(r.stored[stored] / 1e9).toFixed(3)} GJ/m², mixed layer to ${Math.max(...r.mixedDepth).toFixed(1)} m, ${r.mixedDepth[365].toFixed(1)} m on day 365`);
});

test('the worker\'s batches carry the record: every field the figure draws, the layers\' minimum thicknesses of the ocean\'s own table', () => {
  const store = createSeaStore({ days: 30 });
  store.run(1000);
  const batch = store.batch(), L = store.layout.layers;
  assert.equal(L, 58);
  assert.equal(batch.time.length, 31);
  for (const name of ['forcing', 'sst', 'slab', 'depth', 'boundary', 'salinity', 'stored', 'heat', 'slabHeat']) assert.equal(batch[name].length, 31, name);
  assert.equal(batch.bases.length, 31 * L); assert.equal(batch.temperatures.length, 31 * L);
  assert.deepEqual(Array.from(store.layout.minimum), Array.from(store.column.ocean.table.minimum, Math.fround));
  close(batch.bases[L - 1], 4000, 1e-3, 'the floor'); close(batch.salinity[0], 34.7, 0.1, 'the top layer\'s salinity');
  assert.equal(store.batch().time.length, 0);
});
