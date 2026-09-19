import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { saturationHumidity, LATENT_HEAT } from '../js/physics/moist.module.js';

const model = createModel(new Grid(3));
const { core, moist } = model;
const { K, C, dSigma, sigmaMid, cp, g, exnerLayer } = core.diagnostics;

function column(surfaceT, humidity, lapse = 6.5e-3) {
  const pi = new Float64Array(C).fill(101325);
  const theta = new Float64Array(K * C);
  const q = new Float64Array(K * C);
  const qc = new Float64Array(K * C);
  core.diagnose(pi, theta);
  for (let i = 0; i < C; i++) {
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      const p = pi[i] * sigmaMid[k];
      const height = 7000 * Math.log(pi[i] / p);
      const temperature = Math.max(200, surfaceT - lapse * height);
      theta[idx] = temperature / exnerLayer[idx];
      q[idx] = humidity * saturationHumidity(temperature, p);
    }
  }
  return [pi, theta, q, qc];
}

const moistEnthalpy = (pi, theta, q, i) => {
  let sum = 0;
  for (let k = 0; k < K; k++) { const idx = k * C + i; sum += (cp * theta[idx] * exnerLayer[idx] + LATENT_HEAT * q[idx]) * pi[i] * dSigma[k] / g; }
  return sum;
};

test('saturation humidity is near 22 g/kg at 300 K and 1000 hPa and falls with cold and pressure', () => {
  assert.ok(Math.abs(saturationHumidity(300, 1e5) - 0.0223) < 5e-4);
  assert.ok(saturationHumidity(273.15, 1e5) < 0.004 && saturationHumidity(273.15, 1e5) > 0.0035);
  assert.ok(saturationHumidity(280, 5e4) > saturationHumidity(280, 1e5));
});

test('saturation adjustment moves exactly the supersaturation into cloud water and conserves moist enthalpy and water', () => {
  const [pi, theta, q, qc] = column(300, 1.3);
  const before = moistEnthalpy(pi, theta, q, 0), waterBefore = moist.columnWater(pi, q, 0);
  core.diagnoseColumn(0, pi, theta, q, qc);
  const formed = moist.condenseColumn(0, pi, theta, q, qc);
  assert.ok(formed > 0);
  for (let k = 0; k < K; k++) {
    const idx = k * C;
    const qs = saturationHumidity(theta[idx] * exnerLayer[idx], pi[0] * sigmaMid[k]);
    assert.ok(q[idx] <= qs * (1 + 1e-3), `layer ${k} still supersaturated: ${q[idx]} > ${qs}`);
    assert.ok(qc[idx] >= 0);
  }
  assert.ok(Math.abs(moist.columnWater(pi, qc, 0) - formed) < 1e-12 * waterBefore);
  assert.ok(Math.abs(moist.columnWater(pi, q, 0) + moist.columnWater(pi, qc, 0) - waterBefore) < 1e-12 * waterBefore);
  assert.ok(Math.abs(moistEnthalpy(pi, theta, q, 0) - before) < 1e-9 * before);
  // dry the air again: the cloud evaporates back, cooling the layers
  for (let k = 0; k < K; k++) q[k * C] *= 0.5;
  const thetaBefore = Float64Array.from(theta), cloudBefore = moist.columnWater(pi, qc, 0);
  const evaporated = -moist.condenseColumn(0, pi, theta, q, qc);
  assert.ok(evaporated > 0 && evaporated <= cloudBefore * (1 + 1e-12));
  for (let k = 0; k < K; k++) assert.ok(theta[k * C] <= thetaBefore[k * C]);
});

test('autoconversion rains out cloud water above the threshold and conserves water', () => {
  const [pi, theta, q, qc] = column(290, 0.5);
  qc.fill(0);
  qc[(K - 3) * C] = 1e-3;
  qc[(K - 5) * C] = 1e-4;
  const before = moist.columnWater(pi, qc, 0);
  const rain = moist.autoconvertColumn(0, pi, qc, 900);
  assert.ok(rain > 0);
  assert.ok(Math.abs(before - moist.columnWater(pi, qc, 0) - rain) < 1e-15);
  assert.ok(qc[(K - 3) * C] < 1e-3 && qc[(K - 3) * C] > 1.5e-4, 'the thick cloud converts toward the threshold');
  assert.ok(qc[(K - 5) * C] < 1e-4 && qc[(K - 5) * C] > 0.9e-4, 'the thin cloud only decays slowly');
});

test('Betts–Miller convection warms and dries an unstable moist column, conserving enthalpy against its rain', () => {
  const [pi, theta, q, qc] = column(302, 0.9, 7e-3);
  const dt = 450;
  const before = moistEnthalpy(pi, theta, q, 0), waterBefore = moist.columnWater(pi, q, 0);
  core.diagnoseColumn(0, pi, theta, q, qc);
  const top = moist.referenceProfile(0, pi, theta, q);
  assert.ok(top >= 0 && top < K - 1, `expected a convecting column, got top ${top}`);
  const rain = moist.convectColumn(0, pi, theta, q, dt);
  assert.ok(rain > 0, 'deep convection should rain');
  const water = moist.columnWater(pi, q, 0);
  assert.ok(Math.abs(water - (waterBefore - rain)) < 1e-9 * waterBefore, `water ${water} vs ${waterBefore - rain}`);
  assert.ok(Math.abs(moistEnthalpy(pi, theta, q, 0) - before) < 1e-9 * before);
});

test('Betts–Miller leaves a stable dry column alone', () => {
  const [pi, theta, q, qc] = column(280, 0.2, 5e-3);
  const copy = Float64Array.from(theta);
  core.diagnoseColumn(0, pi, theta, q, qc);
  const rain = moist.convectColumn(0, pi, theta, q, 450);
  assert.equal(rain, 0);
  assert.deepEqual(theta, copy);
});

test('transport alone conserves total water', () => {
  const dry = createModel(new Grid(4), { physics: false });
  const init = initializeState(dry, {});
  for (let a = 0; a < init.length; a++) dry.state[a].set(init[a]);
  dry.state[5].set(dry.state[4].map((x) => 0.05 * x));
  const total = () => { const d = dry.diagnostics(); return d.columnWater + d.columnCloud; };
  const water0 = total();
  for (let n = 0; n < 20; n++) dry.step(900);
  const water1 = total();
  assert.ok(Math.abs(water1 - water0) < 1e-10 * water0, `water drifted ${water0} → ${water1}`);
  assert.equal(dry.moist.budget.lost, 0);
});

test('with sources on, the water column change over a day matches evaporation minus precipitation', () => {
  const wet = createModel(new Grid(4));
  const init = initializeState(wet, {});
  for (let a = 0; a < init.length; a++) wet.state[a].set(init[a]);
  const dt = 900, steps = 96;
  const total = (d) => d.columnWater + d.columnCloud;
  const water0 = total(wet.diagnostics());
  let evaporated = 0, precipitated = 0, diag = null;
  for (let n = 0; n < steps; n++) { wet.step(dt); diag = wet.diagnostics(); evaporated += diag.evaporation * dt; precipitated += diag.precipitation * dt; }
  const water1 = total(diag);
  let area = 0; for (let i = 0; i < wet.mesh.nCells; i++) area += wet.mesh.areaCell[i];
  const rained = (wet.moist.budget.condensation + wet.moist.budget.convection) / area;
  const residual = (water1 - water0) - (evaporated - precipitated);
  console.log(`one day at N=4: evaporation ${evaporated.toFixed(3)} mm, precipitation ${precipitated.toFixed(3)} mm (budget ${rained.toFixed(3)}), water ${water0.toFixed(2)} → ${water1.toFixed(2)} mm of which cloud ${(1000 * diag.columnCloud).toFixed(0)} g/m², residual ${residual.toExponential(2)} mm, lost ${(wet.moist.budget.lost / area).toExponential(1)}`);
  assert.ok(evaporated > 0.5 && evaporated < 15, `evaporation ${evaporated} mm/day out of range`);
  assert.ok(Math.abs(rained - precipitated) < 1e-9 * (1 + rained));
  assert.ok(Math.abs(residual) < 0.05 * evaporated, `residual ${residual} vs evaporation ${evaporated}`);
  assert.ok(Number.isFinite(diag.maxWind) && diag.maxWind < 80);
});
