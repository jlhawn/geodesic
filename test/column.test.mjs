import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  HYBRID_LAYERS, labelBoundaries, columnScratch, interfacesOf, layerDensities, isopycnalInterfaces, targetInterfaces, remapColumn, remapMeans, solveColumn,
  sunlightShares, mixedLayerDepth, convectiveAdjustment, boundaryLayer, interiorMixing, mixColumn, H_MIN, GRAVITY,
} from '../js/ocean/column.module.js';
import { seawaterDensity, labelTemperature, SEAWATER } from '../js/ocean/seawater.module.js';
import { transmitted } from '../js/ocean/layered.module.js';
import { ekmanSpiral } from '../explain/ekmanTheory.module.js';
import { pwpColumn } from './helpers/pwp.mjs';

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const within = (a, b, ratio, what) => assert.ok(a < b * ratio && a > b / ratio, `${what}: ${a} against ${b}`);
const RHO_CP = 1025 * 3985, RHO0 = SEAWATER.rho0;
function lcg(seed) {
  let s = seed >>> 0;
  return () => { s = (1664525 * s + 1013904223) >>> 0; return s / 4294967296; };
}
const sum = (a, n = a.length) => { let t = 0; for (let k = 0; k < n; k++) t += a[k]; return t; };
const potentialEnergy = (h, T, S, K) => {
  let z = 0, pe = 0;
  for (let k = 0; k < K; k++) { pe -= GRAVITY * (z + 0.5 * h[k]) * h[k] * (seawaterDensity(T[k], S[k]) - RHO0); z += h[k]; }
  return pe;
};

test('the hybrid table has 58 layers whose labels increase, 3 m at the top, the surface twelve lighter than any class, 10 m minima through the thermocline and a z grid to about 3 km', () => {
  const { densities, minimum, depths, K } = HYBRID_LAYERS;
  assert.equal(K, 58);
  assert.equal(minimum[0], 3);
  assert.equal(depths[0], 0);
  for (let k = 1; k < K; k++) {
    assert.ok(densities[k] > densities[k - 1], `label ${k}`);
    assert.ok(minimum[k] >= minimum[k - 1], `minimum ${k}`);
    close(depths[k + 1], depths[k] + minimum[k], 1e-9, `depth ${k + 1}`);
  }
  assert.ok(densities[11] < 1019.5 && densities[12] === 1019.5);
  assert.ok(depths[12] > 50 && depths[12] < 65, `surface layers to ${depths[12]} m`);
  for (let k = 21; k < 44; k++) assert.equal(minimum[k], 10);
  assert.ok(depths[44] > 350 && depths[44] < 400, `thermocline labels to ${depths[44]} m`);
  assert.ok(depths[K] > 2800 && depths[K] < 3300, `table to ${depths[K]} m`);
});

test('the isopycnal position of a boundary on a linear density ramp is where the ramp crosses it, the surface where the surface is denser, the floor where the bottom is lighter', () => {
  const K = 20, h = new Float64Array(K).fill(10), rho = new Float64Array(K), scratch = columnScratch(K), out = new Float64Array(K + 1);
  for (let k = 0; k < K; k++) rho[k] = 1024 + 0.01 * (10 * k + 5);
  const boundaries = new Float64Array(K + 1).fill(Infinity);
  const cases = [[1, 1023, 0], [2, 1024.25, 25], [3, 1024.5, 50], [4, 1025.0, 100], [5, 1025.55, 155], [6, 1026.0, 200]];
  for (const [k, b] of cases) boundaries[k] = b;
  for (let k = 7; k < K; k++) boundaries[k] = 1030 + k;
  isopycnalInterfaces(h, rho, K, boundaries, out, scratch);
  for (const [k, b, z] of cases) close(out[k], z, 1e-9, `boundary ${b}`);
});

test('the isopycnal generator is a fixed point on uniform layers at their labels, and relax moves an isopycnal interface a fraction of the way', () => {
  const { densities, K } = HYBRID_LAYERS, boundaries = labelBoundaries(densities), random = lcg(7), scratch = columnScratch(K);
  const h = new Float64Array(K), T = new Float64Array(K), S = new Float64Array(K).fill(35), rho = new Float64Array(K);
  for (let k = 0; k < K; k++) { h[k] = 5 + 95 * random(); T[k] = labelTemperature(densities[k], 35); }
  layerDensities(T, S, K, rho);
  for (let k = 1; k < K; k++) assert.ok(rho[k] > rho[k - 1]);
  const current = new Float64Array(K + 1), target = new Float64Array(K + 1), weight = new Float64Array(K + 1).fill(1);
  interfacesOf(h, K, current);
  targetInterfaces('isopycnal', { h, rho, K, depth: current[K], table: HYBRID_LAYERS, boundaries, previousWeight: weight }, target, weight, scratch);
  for (let k = 0; k <= K; k++) close(target[k], current[k], 1e-9, `interface ${k}`);
  for (let k = 1; k < K; k++) assert.equal(weight[k], 1);
  T[20] = labelTemperature(0.5 * (boundaries[21] + densities[21]), 35);
  layerDensities(T, S, K, rho);
  const relaxed = new Float64Array(K + 1), full = new Float64Array(K + 1), w2 = new Float64Array(K + 1);
  targetInterfaces('isopycnal', { h, rho, K, depth: current[K], table: HYBRID_LAYERS, boundaries, previousWeight: weight, relax: 1 }, full, w2, scratch);
  targetInterfaces('isopycnal', { h, rho, K, depth: current[K], table: HYBRID_LAYERS, boundaries, previousWeight: weight, relax: 0.1 }, relaxed, w2, scratch);
  close(full[21], current[20] + H_MIN, 1e-9, 'layer 20, now denser than the boundary, goes to layer 21 but for its minimum');
  close(relaxed[21] - current[21], 0.1 * (current[20] - current[21]), 1e-9, 'a tenth of the way');
  for (const k of [19, 20, 22, 23]) close(relaxed[k], current[k], 1e-9, `interface ${k} stays`);
});

test('the z* generator stretches the depth table to the column and collapses the layers below the floor to the minimum thickness', () => {
  const { K, depths, minimum } = HYBRID_LAYERS, scratch = columnScratch(K), out = new Float64Array(K + 1), weight = new Float64Array(K + 1);
  for (const [D, eta] of [[4000, 0.5], [2000, -0.3], [50, 0]]) {
    const H = D + eta, h = new Float64Array(K).fill(H / K), rho = new Float64Array(K).fill(1026);
    targetInterfaces('zstar', { h, rho, K, depth: D, table: HYBRID_LAYERS, boundaries: labelBoundaries(HYBRID_LAYERS.densities) }, out, weight, scratch);
    close(out[K], H, 1e-9, `${D} m: the column`);
    const spacing = Math.min(H_MIN, H / K);
    for (let k = 0; k < K; k++) {
      const thickness = out[k + 1] - out[k];
      assert.ok(thickness >= spacing - 1e-9, `${D} m, layer ${k}: ${thickness} m`);
      assert.equal(weight[k], 0);
      if (k < K - 1 && depths[k + 1] < D - (K - k - 1) * spacing * D / H) close(thickness, minimum[k] * H / D, 1e-9, `${D} m, layer ${k} stretched`);
    }
  }
});

function warmPool(K, table, D) {
  const h = new Float64Array(K), T = new Float64Array(K), S = new Float64Array(K).fill(35), rho = new Float64Array(K);
  const profile = (z) => 273.15 + (z < 100 ? 28 : 5 + 23 * Math.exp(-(z - 100) / 250));
  let z = 0;
  for (let k = 0; k < K; k++) { h[k] = k < K - 1 ? Math.min(table.minimum[k], D - z - (K - k - 1) * H_MIN) : D - z; T[k] = profile(z + 0.5 * h[k]); z += h[k]; }
  layerDensities(T, S, K, rho);
  return { h, T, S, rho, profile };
}

test('the hybrid generator holds the surface layers at the table, follows the thermocline isopycnals, and collapses the labels denser than the bottom water', () => {
  const table = HYBRID_LAYERS, { K, densities, depths } = table, boundaries = labelBoundaries(densities), scratch = columnScratch(K), D = 4000;
  const { h, rho, profile } = warmPool(K, table, D), out = new Float64Array(K + 1), weight = new Float64Array(K + 1), current = new Float64Array(K + 1);
  interfacesOf(h, K, current);
  targetInterfaces('hybrid', { h, rho, K, depth: D, table, boundaries }, out, weight, scratch);
  close(out[K], D, 1e-9, 'the column');
  for (let k = 1; k <= 12; k++) { close(out[k], depths[k], 1e-9, `surface interface ${k}`); assert.equal(weight[k], 0); }
  const densityAt = (z) => seawaterDensity(profile(z), 35), bottom = densityAt(D);
  let isopycnal = 0;
  for (let k = 13; k < K; k++) {
    assert.ok(out[k] >= out[k - 1] + Math.min(H_MIN, table.minimum[k - 1]) - 1e-9, `layer ${k - 1} keeps its water`);
    if (weight[k] < 1 || boundaries[k] > bottom) continue;
    isopycnal++;
    let low = 0, high = D;
    for (let i = 0; i < 60; i++) { const mid = 0.5 * (low + high); if (densityAt(mid) < boundaries[k]) low = mid; else high = mid; }
    const holding = Math.max(0, current.findIndex((z) => z > out[k]) - 1);
    close(out[k], low, Math.max(4, 0.1 * h[holding]), `isopycnal interface ${k} (label ${densities[k]}) at ${out[k].toFixed(1)} m against the profile's ${low.toFixed(1)} m, in a ${h[holding].toFixed(0)} m layer`);
  }
  assert.ok(isopycnal >= 12, `${isopycnal} isopycnal interfaces`);
  for (let k = 0; k < K; k++) if (densities[k] > bottom + 0.05) close(out[k + 1] - out[k], H_MIN, 1e-9, `label ${densities[k]} is below the floor`);
});

test('in a uniform Southern Ocean column every interface is held at the table and on a 50 m shelf every layer keeps its minimum', () => {
  const table = HYBRID_LAYERS, { K, densities, depths } = table, boundaries = labelBoundaries(densities), scratch = columnScratch(K);
  const out = new Float64Array(K + 1), weight = new Float64Array(K + 1);
  for (const D of [4000, 50]) {
    const water = seawaterDensity(274.15, 34.7), h = new Float64Array(K).fill(D / K), rho = new Float64Array(K).fill(water);
    targetInterfaces('hybrid', { h, rho, K, depth: D, table, boundaries }, out, weight, scratch);
    close(out[K], D, 1e-9, `${D} m: the column`);
    const spacing = Math.min(H_MIN, D / K);
    let held = 0, collapsed = 0;
    for (let k = 1; k < K; k++) {
      assert.ok(out[k] - out[k - 1] >= spacing - 1e-9, `${D} m: layer ${k - 1}`);
      if (boundaries[k] > water) { collapsed++; close(out[k], D - (K - k) * spacing, 1e-9, `${D} m: interface ${k} collapsed to the floor`); continue; }
      held++;
      assert.equal(weight[k], 0, `${D} m: interface ${k} is held`);
      if (D === 4000 && depths[k] < 2000) close(out[k], depths[k], 1e-9, `interface ${k} at the table`);
    }
    assert.ok(held > 40 && collapsed > 5, `${held} held, ${collapsed} collapsed`);
    if (D === 50) assert.ok(out[12] < 30, `the surface layers are squeezed to ${out[12]} m`);
  }
});

function randomGrid(K, random, H) {
  const h = new Float64Array(K);
  let total = 0;
  for (let k = 0; k < K; k++) { h[k] = 0.5 + 4 * random(); total += h[k]; }
  for (let k = 0; k < K; k++) h[k] *= H / total;
  return h;
}

test('the remap conserves the amount, makes no new extremum, and leaves a column on the same grid alone, in both orders', () => {
  const K = 30, random = lcg(11), scratch = columnScratch(K), amounts = new Float64Array(K), means = new Float64Array(K);
  for (const order of [3, 2]) for (let trial = 0; trial < 20; trial++) {
    const hOld = randomGrid(K, random, 100), hNew = randomGrid(K, random, 100), q = Float64Array.from({ length: K }, () => 10 * random());
    let low = Infinity, high = -Infinity, before = 0;
    for (let k = 0; k < K; k++) { low = Math.min(low, q[k]); high = Math.max(high, q[k]); before += hOld[k] * q[k]; }
    remapColumn(hOld, q, K, hNew, amounts, scratch, { order });
    close(sum(amounts), before, 1e-12 * Math.abs(before), `order ${order}, trial ${trial}: the amount`);
    for (let k = 0; k < K; k++) {
      const mean = amounts[k] / hNew[k];
      assert.ok(mean >= low - 1e-12 && mean <= high + 1e-12, `order ${order}, trial ${trial}, layer ${k}: ${mean} outside [${low}, ${high}]`);
    }
    remapMeans(hOld, q, K, hOld, means, scratch, { order });
    for (let k = 0; k < K; k++) close(means[k], q[k], 1e-12 * Math.max(1, Math.abs(q[k])), `order ${order}, trial ${trial}, layer ${k} on the same grid`);
  }
});

test('the parabolic remap reproduces a linear profile exactly away from the boundary layers', () => {
  const K = 30, random = lcg(5), scratch = columnScratch(K), means = new Float64Array(K);
  const hOld = randomGrid(K, random, 200), hNew = randomGrid(K, random, 200), zOld = new Float64Array(K + 1), zNew = new Float64Array(K + 1);
  interfacesOf(hOld, K, zOld);
  interfacesOf(hNew, K, zNew);
  const q = Float64Array.from({ length: K }, (_, k) => 2 + 0.3 * (zOld[k] + 0.5 * hOld[k]));
  remapMeans(hOld, q, K, hNew, means, scratch);
  let checked = 0;
  for (let i = 0; i < K; i++) {
    if (zNew[i] < zOld[1] || zNew[i + 1] > zOld[K - 1]) continue;
    checked++;
    close(means[i], 2 + 0.3 * (zNew[i] + 0.5 * hNew[i]), 1e-9, `layer ${i}`);
  }
  assert.ok(checked >= 20, `${checked} layers checked`);
});

test('on a smooth profile the parabolic remap converges at third order and beats the linear one', () => {
  const H = 400, center = 200, width = 40, exact = (a, b) => width * (Math.log(Math.cosh((b - center) / width)) - Math.log(Math.cosh((a - center) / width))) / (b - a);
  const errors = {};
  for (const order of [3, 2]) for (const N of [40, 80]) {
    const scratch = columnScratch(N), dz = H / N, hOld = new Float64Array(N).fill(dz), hNew = new Float64Array(N), q = new Float64Array(N), means = new Float64Array(N);
    for (let k = 0; k < N; k++) q[k] = exact(k * dz, (k + 1) * dz);
    for (let k = 0; k < N; k++) hNew[k] = k === 0 ? dz * 1.37 : k === N - 1 ? dz * 0.63 : dz;
    remapMeans(hOld, q, N, hNew, means, scratch, { order });
    let error = 0, z = 0;
    for (let k = 0; k < N; k++) { error += Math.abs(means[k] - exact(z, z + hNew[k])) * hNew[k]; z += hNew[k]; }
    errors[`${order}:${N}`] = error / H;
  }
  console.log(`remap L1 error, PPM 40/80 cells ${errors['3:40'].toExponential(2)} / ${errors['3:80'].toExponential(2)}, PLM ${errors['2:40'].toExponential(2)} / ${errors['2:80'].toExponential(2)}`);
  assert.ok(errors['3:40'] / errors['3:80'] > 5, `PPM ratio ${errors['3:40'] / errors['3:80']}`);
  assert.ok(errors['2:40'] / errors['2:80'] > 3, `PLM ratio ${errors['2:40'] / errors['2:80']}`);
  assert.ok(errors['3:80'] < errors['2:80'] / 2, 'PPM beats PLM');
});

test('a step keeps its sharpness through the remap: the layer holding the jump takes the exact fraction and the rest stay at 0 or 1', () => {
  const K = 10, scratch = columnScratch(K), hOld = new Float64Array(K).fill(20), q = Float64Array.from({ length: K }, (_, k) => (k < 5 ? 0 : 1)), hNew = new Float64Array(K).fill(20), means = new Float64Array(K);
  hNew[4] = 27;
  hNew[5] = 13;
  remapMeans(hOld, q, K, hNew, means, scratch);
  for (let k = 0; k < K; k++) {
    if (k === 4) close(means[k], 7 / 27, 1e-12, 'the layer holding the jump');
    else close(means[k], k < 4 ? 0 : 1, 1e-12, `layer ${k}`);
  }
});

test('the implicit solve conserves the column, mixes it uniform under a large diffusivity, and agrees with the two-layer solution', () => {
  const K = 6, random = lcg(3), scratch = columnScratch(K), h = randomGrid(K, random, 60), q = Float64Array.from({ length: K }, () => 5 + 10 * random()), kappa = Float64Array.from({ length: K }, () => 1e-3 * random());
  kappa[K - 1] = 0;
  const before = Float64Array.from(q), total = (a) => { let t = 0; for (let k = 0; k < K; k++) t += h[k] * a[k]; return t; };
  solveColumn(h, q, K, kappa, 3600, scratch);
  close(total(q), total(before), 1e-12 * total(before), 'the column');
  const huge = Float64Array.from(kappa, (x, k) => (k < K - 1 ? 1e6 : 0)), uniform = Float64Array.from(before);
  solveColumn(h, uniform, K, huge, 3600, scratch);
  for (let k = 0; k < K; k++) close(uniform[k], total(before) / sum(h), 1e-6, `layer ${k} uniform`);
  const none = Float64Array.from(before);
  assert.equal(solveColumn(h, none, K, new Float64Array(K), 3600, scratch), false);
  assert.deepEqual(none, before);
  const two = columnScratch(2), h2 = Float64Array.from([10, 30]), q2 = Float64Array.from([20, 10]), k2 = Float64Array.from([2e-3, 0]), dt = 1800;
  const r = dt * k2[0] / (0.5 * (h2[0] + h2[1]));
  const q0 = (h2[0] * q2[0] * (h2[1] + r) + r * h2[1] * q2[1]) / (h2[0] * h2[1] + r * (h2[0] + h2[1]));
  solveColumn(h2, q2, 2, k2, dt, two);
  close(q2[0], q0, 1e-12, 'the top layer');
  close(h2[0] * q2[0] + h2[1] * q2[1], 10 * 20 + 30 * 10, 1e-9, 'the pair');
});

test('the sunlight shares are the Jerlov transmission across each layer and sum to one, the floor keeping what passes it', () => {
  for (const [K, dz] of [[40, 10], [4, 5]]) {
    const z = Float64Array.from({ length: K + 1 }, (_, k) => k * dz), shares = new Float64Array(K);
    sunlightShares(z, K, 'IB', shares);
    close(sum(shares), 1, 1e-12, `${K} layers sum`);
    for (let k = 0; k < K - 1; k++) close(shares[k], transmitted(k * dz) - transmitted((k + 1) * dz), 1e-12, `${K} layers, share ${k}`);
    assert.ok(shares[K - 1] >= transmitted((K - 1) * dz) - transmitted(K * dz), 'the floor keeps the rest');
  }
});

test('the mixed-layer depth interpolates the first excess over the surface density between layer centers', () => {
  const h = Float64Array.from([10, 10, 10, 10, 10]), rho = Float64Array.from([1024, 1024.01, 1024.02, 1024.1, 1024.5]);
  close(mixedLayerDepth(h, rho, 5), 26.25, 1e-9, 'the depth');
  close(mixedLayerDepth(h, Float64Array.from([1024, 1024, 1024, 1024, 1024]), 5), 50, 1e-9, 'a uniform column');
});

test('convective adjustment merges unstable neighbors into uniform slabs, conserves heat and salt, and returns the potential energy the surface slab released', () => {
  const K = 5, scratch = columnScratch(K), h = Float64Array.from([5, 5, 10, 20, 40]), T = Float64Array.from([10, 15, 15, 14, 12], (t) => t + 273.15), S = new Float64Array(K).fill(35);
  const heat = sum(Float64Array.from(h, (x, k) => x * T[k])), before = potentialEnergy(h, T, S, K);
  const result = convectiveAdjustment(h, T, S, K, scratch);
  close(sum(Float64Array.from(h, (x, k) => x * T[k])), heat, 1e-9, 'heat');
  for (let k = 0; k < 4; k++) close(T[k], 13.875 + 273.15, 1e-12, `layer ${k}`);
  close(T[4], 12 + 273.15, 1e-12, 'the bottom layer');
  assert.equal(result.depth, 40);
  assert.ok(result.released > 0, `released ${result.released} J/m²`);
  close(result.released, before - potentialEnergy(h, T, S, K), 1e-9, 'the energy released');
  assert.equal(result.interior, 0);
  const stable = Float64Array.from([20, 18, 16, 14, 12], (t) => t + 273.15), copy = Float64Array.from(stable);
  const nothing = convectiveAdjustment(h, stable, S, K, scratch);
  assert.deepEqual(stable, copy);
  assert.equal(nothing.released, 0);
  assert.equal(nothing.depth, h[0]);
});

function ramp(K, dz, top = 20, gradient = 0.1) {
  const h = new Float64Array(K).fill(dz), T = new Float64Array(K), S = new Float64Array(K).fill(35);
  for (let k = 0; k < K; k++) T[k] = 273.15 + top - gradient * (k + 0.5) * dz;
  return { h, T, S };
}

test('the boundary layer uses the energy it is given: no wind and no convection mix nothing, and under wind the work of mixing matches the energy to within the search', () => {
  const K = 50, { h, T, S } = ramp(K, 2), scratch = columnScratch(K), kappa = new Float64Array(K);
  const still = boundaryLayer({ h, T, S, K, ustar: 0, f: 1e-4, dt: 1350 }, kappa, scratch);
  assert.equal(still.depth, 0);
  assert.ok(kappa.every((x) => x === 0));
  const windy = boundaryLayer({ h, T, S, K, ustar: 0.01, f: 1e-4, dt: 1350 }, kappa, scratch);
  console.log(`ramp under u* 1 cm/s at f 1e-4: boundary layer ${windy.depth.toFixed(2)} m, work ${windy.work.toExponential(3)} of ${windy.energy.toExponential(3)} J/m²`);
  assert.ok(windy.depth > 5 && windy.depth < 60, `depth ${windy.depth}`);
  assert.ok(windy.work <= windy.energy * (1 + 1e-9) && windy.work > 0.99 * windy.energy, `work ${windy.work} of ${windy.energy}`);
  let active = 0;
  for (let k = 0; k < K - 1; k++) if (kappa[k] > 0) { active++; assert.ok((k + 1) * 2 < windy.depth); }
  assert.ok(active > 0);
});

test('the boundary-layer depth is insensitive to the vertical resolution', () => {
  const depths = [1, 4].map((dz) => {
    const K = Math.round(100 / dz), { h, T, S } = ramp(K, dz), scratch = columnScratch(K), kappa = new Float64Array(K);
    return boundaryLayer({ h, T, S, K, ustar: 0.01, f: 1e-4, dt: 1350 }, kappa, scratch).depth;
  });
  console.log(`boundary layer at 1 m / 4 m layers: ${depths.map((d) => d.toFixed(2)).join(' / ')} m`);
  within(depths[0], depths[1], 1.1, 'the depth at 1 m against 4 m');
});

test('wind over a mixed layer entrains at the Kraus–Turner rate: the potential energy raised is the energy supplied and the base deepens as 2 E t / (g Δρ h₀)', () => {
  const K = 150, h = new Float64Array(K).fill(1), T = new Float64Array(K), S = new Float64Array(K).fill(35), scratch = columnScratch(K), kappa = new Float64Array(K);
  for (let k = 0; k < K; k++) T[k] = 273.15 + (k < 50 ? 20 : 15);
  const dt = 1350, steps = 256, ustar = 0.02, start = potentialEnergy(h, T, S, K);
  let supplied = 0, worked = 0;
  for (let n = 0; n < steps; n++) {
    const layer = boundaryLayer({ h, T, S, K, ustar, f: 0, dt, langmuir: 0 }, kappa, scratch);
    supplied += layer.energy;
    worked += layer.work;
    solveColumn(h, T, K, kappa, dt, scratch);
    solveColumn(h, S, K, kappa, dt, scratch);
  }
  const raised = potentialEnergy(h, T, S, K) - start, dRho = seawaterDensity(273.15 + 15, 35) - seawaterDensity(273.15 + 20, 35);
  const expected = 50 + 2 * supplied / (GRAVITY * dRho * 50);
  let base = 0;
  while (base < K && T[base] > 273.15 + 17.5) base++;
  console.log(`Kraus–Turner: supplied ${supplied.toFixed(1)} J/m², worked ${worked.toFixed(1)}, raised ${raised.toFixed(1)}; base at ${base} m against ${expected.toFixed(1)} m`);
  close(worked, supplied, 0.02 * supplied, 'the work against the energy');
  assert.ok(raised < supplied && raised > 0.75 * supplied, `the potential energy raised, ${raised}, is the work less what cabbeling released, against ${supplied}`);
  close(base, expected, 3, 'the base');
});

test('surface cooling without wind convects to a little beyond the non-penetrative depth √(2Q/Γ)', () => {
  const K = 100, { h, T, S } = ramp(K, 1, 20, 0.05), scratch = columnScratch(K), out = { kappa: new Float64Array(K), nu: new Float64Array(K) };
  const heat = sum(Float64Array.from(h, (x, k) => x * T[k]));
  let deepest = 0;
  for (let n = 0; n < 40; n++) {
    const result = mixColumn({ h, T, S, K, heat: -0.5, ustar: 0, f: 1e-4, dt: 1350 }, out, scratch);
    deepest = Math.max(deepest, result.depth);
  }
  close(sum(Float64Array.from(h, (x, k) => x * T[k])), heat - 20, 1e-9, 'heat');
  const rho = layerDensities(T, S, K, new Float64Array(K)), mixed = mixedLayerDepth(h, rho, K, 0.003);
  console.log(`convection after 20 K m of cooling on a 0.05 K/m ramp: mixed to ${mixed.toFixed(1)} m, boundary layer reached ${deepest.toFixed(1)} m, √(2Q/Γ) ${Math.sqrt(2 * 20 / 0.05).toFixed(1)} m`);
  assert.ok(mixed >= Math.sqrt(2 * 20 / 0.05) - 1 && mixed < 1.3 * Math.sqrt(2 * 20 / 0.05), `mixed to ${mixed} m`);
});

test('the interior mixing is the Pacanowski–Philander form of the Richardson number', () => {
  const K = 4, h = new Float64Array(K).fill(10), T = Float64Array.from([20, 18, 16, 14], (t) => t + 273.15), S = new Float64Array(K).fill(35), scratch = columnScratch(K);
  const kappa = new Float64Array(K), nu = new Float64Array(K), u = new Float64Array(K), v = new Float64Array(K);
  interiorMixing({ h, T, S, u, v, K }, kappa, nu, scratch);
  for (let k = 0; k < K - 1; k++) { close(nu[k], 1e-4, 1e-6, `calm viscosity ${k}`); close(kappa[k], 1e-5, 1e-6, `calm diffusivity ${k}`); }
  u[0] = 1;
  interiorMixing({ h, T, S, u, v, K }, kappa, nu, scratch);
  const rho = layerDensities(T, S, K, new Float64Array(K)), ri = GRAVITY * (rho[1] - rho[0]) / RHO0 * 10 / 1;
  close(nu[0], 1e-2 / (1 + 5 * ri) ** 2 + 1e-4, 1e-12, 'the sheared viscosity');
  close(kappa[0], nu[0] / (1 + 5 * ri) + 1e-5, 1e-12, 'the sheared diffusivity');
  T[0] = 273.15 + 10;
  interiorMixing({ h, T, S, u, v, K }, kappa, nu, scratch);
  assert.equal(kappa[0], 0.1);
});

function diurnalColumn({ wind, sun, loss, latitude = 45, days = 30, dt = 900 }) {
  const layers = [];
  for (let k = 0; k < 40; k++) layers.push(1);
  for (let k = 0; k < 20; k++) layers.push(2);
  for (let k = 0; k < 30; k++) layers.push(4);
  const K = layers.length, h = Float64Array.from(layers), T = new Float64Array(K), S = new Float64Array(K).fill(35), u = new Float64Array(K), v = new Float64Array(K);
  const scratch = columnScratch(K), out = { kappa: new Float64Array(K), nu: new Float64Array(K) }, interfaces = new Float64Array(K + 1);
  interfacesOf(h, K, interfaces);
  for (let k = 0; k < K; k++) { const z = interfaces[k] + 0.5 * h[k]; T[k] = 273.15 + (z < 40 ? 20 : 20 - 0.05 * (z - 40)); }
  const f = 2 * 7.292e-5 * Math.sin(latitude * Math.PI / 180), tau = 1.2 * 1.2e-3 * wind * wind, ustar = Math.sqrt(tau / RHO0), damp = 1 / (2 * 86400);
  const c = Math.cos(f * dt), s = Math.sin(f * dt), steps = Math.round(days * 86400 / dt), dawn = [], tops = new Map();
  for (let step = 0; step < steps; step++) {
    const middle = (step + 0.5) * dt, hour = (middle / 3600) % 24;
    mixColumn({ h, T, S, u, v, K, heat: -loss * dt / RHO_CP, light: sun(middle) * dt / RHO_CP, ustar, f, dt, options: { langmuir: 0 } }, out, scratch);
    for (let k = 0; k < K; k++) { const uu = u[k], vv = v[k]; u[k] = (uu * c + vv * s) * (1 - damp * dt); v[k] = (vv * c - uu * s) * (1 - damp * dt); }
    u[0] += tau / (RHO0 * h[0]) * dt;
    solveColumn(h, u, K, out.nu, dt, scratch);
    solveColumn(h, v, K, out.nu, dt, scratch);
    const day = Math.floor(middle / 86400);
    if (day < days - 10) continue;
    if (Math.abs(hour - 5.875) < 0.2) {
      let k = 0;
      while (k + 1 < K && T[0] - T[k + 1] < 0.02) k++;
      dawn.push(interfaces[k + 1]);
    }
    const [low, high] = tops.get(day) ?? [Infinity, -Infinity];
    tops.set(day, [Math.min(low, T[0]), Math.max(high, T[0])]);
  }
  const median = (a) => [...a].sort((x, y) => x - y)[Math.floor(a.length / 2)];
  return { dawn: median(dawn), range: median([...tops.values()].map(([low, high]) => high - low)) };
}

test('under a daily cycle the column physics puts the dawn depth within a factor 1.5 or 2.5 m of a resolved Price–Weller–Pinkel column at 2 m/s, and deeper under an 8 m/s wind, which PWP\'s shear criterion does not feel', () => {
  const MEAN_SUN = 250, sun = (t) => { const hour = (t / 3600) % 24; return hour > 6 && hour < 18 ? MEAN_SUN * Math.PI * Math.sin(Math.PI * (hour - 6) / 12) : 0; };
  const cases = [[2, 0.6], [2, 0.8], [8, 0.7]].map(([wind, ratio]) => {
    const ours = diurnalColumn({ wind, sun, loss: ratio * MEAN_SUN }), pwp = pwpColumn({ wind, sun, loss: ratio * MEAN_SUN, latitude: 45 });
    console.log(`U ${wind} m/s, N/S ${ratio}: dawn depth ${ours.dawn.toFixed(1)} m against PWP's ${pwp.dawn.toFixed(1)} m; daily range of the top layer ${ours.range.toFixed(2)} K against PWP's ${pwp.range.toFixed(2)} K`);
    return { wind, ratio, ours, pwp };
  });
  for (const { wind, ratio, ours, pwp } of cases) {
    if (wind > 4) assert.ok(ours.dawn >= pwp.dawn && ours.dawn < 2.5 * pwp.dawn, `U ${wind}, N/S ${ratio}: the dawn depth ${ours.dawn} against PWP's ${pwp.dawn}`);
    else assert.ok(Math.abs(ours.dawn - pwp.dawn) <= 2.5 || (ours.dawn < 1.5 * pwp.dawn && ours.dawn > pwp.dawn / 1.5), `U ${wind}, N/S ${ratio}: the dawn depth ${ours.dawn} against PWP's ${pwp.dawn}`);
    within(ours.range, pwp.range, 2.5, `U ${wind}, N/S ${ratio}: the daily range`);
  }
});

test('a steady stress on a column of constant viscosity spins up the Ekman spiral: the transport τ/(ρ₀ f) at right angles and the surface current 45° to the wind', () => {
  const K = 64, dz = 1.5, h = new Float64Array(K).fill(dz), u = new Float64Array(K), v = new Float64Array(K), nu = new Float64Array(K).fill(0.02), scratch = columnScratch(K);
  nu[K - 1] = 0;
  const tau = 0.1, f = 1e-4, dt = 60, epsilon = f / 200, a = 0.5 * f * dt, steps = Math.round(120 * 86400 / dt);
  const re = (1 - a * a) / (1 + a * a), im = -2 * a / (1 + a * a);
  for (let n = 0; n < steps; n++) {
    for (let k = 0; k < K; k++) { const uu = u[k], vv = v[k]; u[k] = (uu * re - vv * im) * (1 - epsilon * dt); v[k] = (uu * im + vv * re) * (1 - epsilon * dt); }
    u[0] += tau / (RHO0 * h[0]) * dt;
    solveColumn(h, u, K, nu, dt, scratch);
    solveColumn(h, v, K, nu, dt, scratch);
  }
  const spiral = ekmanSpiral({ tau, f, viscosity: 0.02, rho: RHO0 });
  let east = 0, north = 0;
  for (let k = 0; k < K; k++) { east += h[k] * u[k]; north += h[k] * v[k]; }
  console.log(`Ekman: transport ${east.toExponential(3)} / ${north.toExponential(3)} m²/s against 0 / ${spiral.transport[1].toExponential(3)}; surface ${(Math.atan2(v[0], u[0]) * 180 / Math.PI).toFixed(2)}° at ${Math.hypot(u[0], v[0]).toFixed(4)} m/s`);
  close(east, 0, 0.015 * Math.abs(spiral.transport[1]), 'eastward transport');
  close(north, spiral.transport[1], 0.015 * Math.abs(spiral.transport[1]), 'northward transport');
  for (const k of [0, 6, 13, 26]) {
    const [ue, ve] = spiral.at(-(k + 0.5) * dz);
    close(u[k], ue, 0.02 * spiral.surface, `east at ${(k + 0.5) * dz} m`);
    close(v[k], ve, 0.02 * spiral.surface, `north at ${(k + 0.5) * dz} m`);
  }
  const [u0, v0] = spiral.at(-0.5 * dz);
  close(Math.atan2(v[0], u[0]) * 180 / Math.PI, Math.atan2(v0, u0) * 180 / Math.PI, 1, 'the surface angle');
});
