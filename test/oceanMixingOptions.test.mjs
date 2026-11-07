import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean, backgroundDiffusivities, HENYEY_RATIO, EQUATORIAL_BACKGROUND } from '../js/ocean/layered.module.js';
import { interiorMixing, shearTaper } from '../js/ocean/column.module.js';

/*
 * The ocean's optional interior closures on the CPU at N = 16, each against
 * its published form and against the default it leaves alone: the
 * latitude dependence of the background diffusivity (backgroundLatitude)
 * and the Richardson-number cutoff of the interior shear mixing
 * (shearCutoff).
 */
const mesh = buildMesh(new Grid(16)), { nCells: C, nEdges: E, latCell, latEdge } = mesh;
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

test('shearTaper is Large, McWilliams and Doney\'s interior shape: 1 at Ri 0, (3/4)³ at half the cutoff, 0 at and above it', () => {
  assert.equal(shearTaper(0, 0.7), 1);
  assert.ok(Math.abs(shearTaper(0.35, 0.7) - 0.421875) < 1e-15);
  assert.equal(shearTaper(0.7, 0.7), 0);
  assert.equal(shearTaper(5, 0.7), 0);
  for (let r = 0.01; r < 0.7; r += 0.01) assert.ok(shearTaper(r, 0.7) < shearTaper(r - 0.01, 0.7));
});

test('under shearCutoff the interior keeps PP81\'s background and scales only its shear term by the taper, zero above the cutoff; without it the shear term is PP81\'s', () => {
  const K = 8, h = new Float64Array(K).fill(10), T = Float64Array.from({ length: K }, (_, k) => 293.15 - 0.5 * k), S = new Float64Array(K).fill(35);
  const u = Float64Array.of(1, 0.75, 0.55, 0.4, 0.3, 0.25, 0.24, 0.235), v = new Float64Array(K);
  const scratch = { rho: new Float64Array(K) }, nu0 = 5e-3, nub = 1e-4, kb = 1e-5;
  const run = (extra) => { const kappa = new Float64Array(K), nu = new Float64Array(K); interiorMixing({ h, T, S, u, v, K, shearViscosity: nu0, backgroundViscosity: nub, backgroundDiffusivity: kb, ...extra }, kappa, nu, scratch); return { kappa, nu }; };
  const off = run({}), infinite = run({ shearCutoff: Infinity }), cut = run({ shearCutoff: 0.7 });
  assert.deepEqual(Array.from(infinite.nu), Array.from(off.nu));
  assert.deepEqual(Array.from(infinite.kappa), Array.from(off.kappa));
  let below = 0, above = 0;
  for (let k = 0; k < K - 1; k++) {
    const ri = (Math.sqrt(nu0 / (off.nu[k] - nub)) - 1) / 5;
    if (ri >= 0.7 * (1 + 1e-9)) {
      above++;
      assert.equal(cut.nu[k], nub, `interface ${k} at Ri ${ri}`);
      assert.ok(Math.abs(cut.kappa[k] - (nub / (1 + 5 * ri) + kb)) < 1e-12 * cut.kappa[k], `interface ${k} at Ri ${ri}`);
    } else if (ri < 0.7 * (1 - 1e-6)) {
      below++;
      const ratio = (cut.nu[k] - nub) / (off.nu[k] - nub);
      assert.ok(Math.abs(ratio - (1 - (ri / 0.7) ** 2) ** 3) < 1e-6, `interface ${k} at Ri ${ri}: ${ratio}`);
      assert.ok(Math.abs(cut.kappa[k] - (cut.nu[k] / (1 + 5 * ri) + kb)) < 1e-9 * cut.kappa[k], `interface ${k}`);
    }
  }
  assert.ok(below >= 2 && above >= 2, `${below} interfaces below the cutoff and ${above} above`);
  assert.throws(() => createOcean(mesh, { shearCutoff: 0 }), /shearCutoff/);
});

test('shearCutoff reaches the ocean\'s column phase: under a sheared surface current every interface\'s viscosity is at most the default\'s and some are less', () => {
  const sheared = (options) => {
    const ocean = stratified(options), L = ocean.layers;
    for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) ocean.u[k * E + e] = 0.4 * Math.exp(-k / 3) * Math.cos(latEdge[e]) * Math.sin(3 * e);
    ocean.phases.column(DT);
    return ocean;
  };
  const base = sheared({}), cut = sheared({ shearCutoff: 0.7 });
  let less = 0;
  for (let n = 0; n < base.viscosity.length; n++) {
    assert.ok(cut.viscosity[n] <= base.viscosity[n] * (1 + 1e-12), `entry ${n}: ${cut.viscosity[n]} over ${base.viscosity[n]}`);
    if (cut.viscosity[n] < base.viscosity[n] * (1 - 1e-6)) less++;
  }
  assert.ok(less > 0, 'the cutoff changed no viscosity');
});
