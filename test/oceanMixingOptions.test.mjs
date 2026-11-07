import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean, backgroundDiffusivities, HENYEY_RATIO, EQUATORIAL_BACKGROUND } from '../js/ocean/layered.module.js';

/*
 * The ocean's optional interior closures on the CPU at N = 16, each against
 * its published form and against the default it leaves alone: the
 * latitude dependence of the background diffusivity (backgroundLatitude).
 */
const mesh = buildMesh(new Grid(16)), { nCells: C, latCell } = mesh;
const DT = 1350, DEGREE = Math.PI / 180;
const latitudes = (degrees) => ({ nCells: degrees.length, latCell: Float64Array.from(degrees, (d) => d * DEGREE) });

test('backgroundLatitude henyey is MOM6\'s Henyey scaling: κ_b at 30°, the floor at the equator, 1.68 κ_b at the poles, symmetric and rising with |latitude|', () => {
  const degrees = [-90, -60, -30, -15, -5, -1, 0, 1, 5, 15, 30, 60, 90], n = degrees.length;
  const kappa = backgroundDiffusivities(latitudes(degrees), 'henyey', 1e-5, 1e-6);
  const G = (d) => { const s = Math.max(1e-10, Math.abs(Math.sin(d * DEGREE))), x = HENYEY_RATIO / s; return s * Math.log(x + Math.sqrt(x * x - 1)); };
  for (let m = 0; m < n; m++) {
    const expected = Math.max(1e-6, 1e-5 * G(degrees[m]) / G(30));
    assert.ok(Math.abs(kappa[m] - expected) <= 1e-12 * expected, `${degrees[m]}°: ${kappa[m]} against ${expected}`);
    assert.equal(kappa[m], kappa[n - 1 - m]);
  }
  assert.ok(Math.abs(kappa[10] - 1e-5) < 1e-17, `30°: ${kappa[10]}`);
  assert.equal(kappa[6], 1e-6);
  assert.ok(Math.abs(kappa[12] / 1e-5 - 1.684) < 1e-3, `pole: ${kappa[12]}`);
  for (let m = 7; m < n; m++) assert.ok(kappa[m] >= kappa[m - 1], `${degrees[m]}°`);
  const om4 = backgroundDiffusivities(latitudes([0, 2, 5]), 'henyey', 1.5e-5, 2e-6);
  assert.deepEqual(Array.from(om4.subarray(0, 2)), [2e-6, 2e-6]);
  assert.ok(om4[2] > 3.6e-6 && om4[2] < 3.7e-6, `OM4's 5°: ${om4[2]}`);
});

test('backgroundLatitude equatorial is MOM025\'s: κ_eq within 5°, κ_b poleward of 15°, and a raised cosine between with the mean at 10°', () => {
  const degrees = [0, 3, 5, 7.5, 10, 12.5, 15.5, 20, 45, 90, -10, -90];
  const kappa = backgroundDiffusivities(latitudes(degrees), 'equatorial', 1e-5, 1e-6);
  assert.equal(kappa[0], 1e-6); assert.equal(kappa[1], 1e-6); assert.equal(kappa[2], 1e-6);
  for (const m of [6, 7, 8, 9, 11]) assert.equal(kappa[m], 1e-5, `${degrees[m]}°`);
  assert.ok(Math.abs(kappa[4] - 5.5e-6) < 1e-18, `10°: ${kappa[4]}`);
  assert.equal(kappa[10], kappa[4]);
  const quarter = 1e-6 + 9e-6 * 0.5 * (1 - Math.cos(Math.PI / 4));
  assert.ok(Math.abs(kappa[3] - quarter) < 1e-18, `7.5°: ${kappa[3]}`);
  for (let m = 1; m < 10; m++) assert.ok(kappa[m] >= kappa[m - 1], `${degrees[m]}°`);
  assert.deepEqual(EQUATORIAL_BACKGROUND, { inner: 5, outer: 15 });
});

test('backgroundLatitude uniform holds κ_b in every cell, and an unknown shape or a negative equatorial diffusivity is refused', () => {
  const kappa = backgroundDiffusivities(mesh, 'uniform', 1e-5, 1e-6);
  assert.equal(kappa.length, C);
  for (const v of kappa) assert.equal(v, 1e-5);
  assert.throws(() => backgroundDiffusivities(mesh, 'tanh', 1e-5, 1e-6), /backgroundLatitude is one of/);
  assert.throws(() => createOcean(mesh, { backgroundLatitude: 'step' }), /backgroundLatitude is one of/);
  assert.throws(() => createOcean(mesh, { backgroundLatitude: 'equatorial', equatorialDiffusivity: -1e-6 }), /equatorialDiffusivity/);
});

function stratified(options) {
  const ocean = createOcean(mesh, { everySteps: 1, iceDrift: false, ...options });
  ocean.initialize(Float64Array.from(latCell, (lat) => 271.5 + 30 * Math.cos(lat) ** 2), new Float64Array(C));
  return ocean;
}

test('one column phase under backgroundLatitude equatorial leaves every column poleward of 15° as the uniform background does, bit for bit, and moves a tenth of the heat within 5°', () => {
  const uniform = stratified({}), equatorial = stratified({ backgroundLatitude: 'equatorial' }), L = uniform.layers;
  assert.deepEqual(Array.from(equatorial.Q), Array.from(uniform.Q));
  const before = Float64Array.from(uniform.Q);
  uniform.phases.column(DT); equatorial.phases.column(DT);
  let near = 0, moved = 0, movedLess = 0;
  for (let i = 0; i < C; i++) {
    if (!uniform.cellOcean[i]) continue;
    const degrees = Math.abs(latCell[i]) / DEGREE;
    let a = 0, b = 0;
    for (let k = 0; k < L; k++) {
      const n = k * C + i;
      if (degrees > 15.01) { assert.equal(equatorial.Q[n], uniform.Q[n]); assert.equal(equatorial.h[n], uniform.h[n]); assert.equal(equatorial.W[n], uniform.W[n]); }
      a += Math.abs(uniform.Q[n] - before[n]); b += Math.abs(equatorial.Q[n] - before[n]);
    }
    if (degrees < 5 && a > 0) { near++; moved += a; movedLess += b; }
  }
  assert.ok(near > 0, 'no stratified column within 5°');
  const ratio = movedLess / moved;
  assert.ok(ratio > 0.08 && ratio < 0.12, `within 5° the column moves ${ratio} of the uniform background's heat`);
});
