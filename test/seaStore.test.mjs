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

