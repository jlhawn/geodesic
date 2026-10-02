import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { spongeGeometry, spongeRates, dampEddies, spongeSigmaFor, SPONGE, lidFrictionRates, lidFrictionFor, LID_FRICTION } from '../js/dynamics/sponge.module.js';
import { P0, sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};

function zonalField(mesh, wind) {
  return Float64Array.from({ length: mesh.nEdges }, (_, e) => {
    const x = mesh.xEdge[3 * e], y = mesh.xEdge[3 * e + 1], z = mesh.xEdge[3 * e + 2], lat = Math.asin(z / Math.hypot(x, y, z)), lon = Math.atan2(y, x);
    const [east, north] = wind(lat, lon), n = mesh.nEdge.subarray(3 * e, 3 * e + 3);
    return east * (-Math.sin(lon) * n[0] + Math.cos(lon) * n[1]) + north * (-Math.sin(lat) * Math.cos(lon) * n[0] - Math.sin(lat) * Math.sin(lon) * n[1] + Math.cos(lat) * n[2]);
  });
}
const rms = (a) => Math.sqrt(a.reduce((s, x) => s + x * x, 0) / a.length);

// The axial angular momentum per unit mass of a layer's wind, summed over the cells: Σ A a cos φ u_east.
function angularMomentum(mesh, u) {
  const v = cellVector(mesh, u);
  let sum = 0;
  for (let i = 0; i < mesh.nCells; i++) sum += mesh.areaCell[i] * mesh.radius * Math.cos(mesh.latCell[i]) * (-Math.sin(mesh.lonCell[i]) * v[3 * i] + Math.cos(mesh.lonCell[i]) * v[3 * i + 1]);
  return sum;
}

test('the sponge keeps the zonal mean and damps the zonally asymmetric wind at its rate', () => {
  const lines = [];
  for (const N of [16, 32]) {
    const mesh = buildMesh(new Grid(N)), geometry = spongeGeometry(mesh), means = new Float64Array(2 * geometry.bands);
    const mean = zonalField(mesh, (lat) => [30 * Math.cos(lat) + 80 * Math.exp(-(((lat + 1.05) / 0.25) ** 2)) - 40 * Math.exp(-(((lat - 0.5) / 0.3) ** 2)), 3 * Math.sin(2 * lat)]);
    const wave = zonalField(mesh, (lat, lon) => [20 * Math.cos(lat) ** 2 * Math.cos(2 * lon), 15 * Math.cos(lat) ** 2 * Math.sin(2 * lon + 0.4)]);
    const rate = 1 / 86400, dt = 3600, keep = 1 / (1 + rate * dt);
    const alone = Float64Array.from(mean);
    dampEddies(mesh, geometry, alone, rate, dt, means);
    const kept = rms(alone.map((x, e) => x - mean[e]));
    const both = Float64Array.from(mean, (x, e) => x + wave[e]);
    dampEddies(mesh, geometry, both, rate, dt, means);
    const waveLeft = both.map((x, e) => x - mean[e]), expected = rms(wave) * keep;
    const torque = (angularMomentum(mesh, both) - angularMomentum(mesh, Float64Array.from(mean, (x, e) => x + wave[e]))) / (rate * dt * angularMomentum(mesh, mean));
    const offset = kept / (1 - keep);
    lines.push(`N=${N} (${geometry.bands} bands): a zonal flow of rms ${rms(mean).toFixed(1)} m/s is ${offset.toFixed(3)} m/s rms from its reconstructed zonal mean, ${(100 * offset / rms(mean)).toFixed(2)} %, and moves ${kept.toExponential(1)} m/s rms in an hour at a one-day rate (the wave ${(rms(wave) - rms(waveLeft)).toFixed(3)} m/s, expected ${(rms(wave) - expected).toFixed(3)}); the axial angular momentum changes by ${torque.toExponential(1)} of what a Rayleigh drag at the same rate would remove`);
    assert.ok(offset < (N === 16 ? 0.02 : 0.006) * rms(mean), `the zonal flow is ${offset} m/s rms from its reconstruction`);
    assert.ok(Math.abs(rms(waveLeft) - expected) < 0.02 * (rms(wave) - expected), `the wave left ${rms(waveLeft)} m/s rms, expected ${expected}`);
    assert.ok(Math.abs(torque) < 0.01, `the sponge's torque is ${torque} of a Rayleigh drag's`);
  }
  console.log(lines.join('\n'));
});

test('the lid friction damps the zonal mean at its rate, removing the axial angular momentum a Rayleigh drag on the zonal flow would, and leaves the eddies to the sponge', () => {
  const mesh = buildMesh(new Grid(32)), geometry = spongeGeometry(mesh), means = new Float64Array(2 * geometry.bands);
  const mean = zonalField(mesh, (lat) => [30 * Math.cos(lat) + 80 * Math.exp(-(((lat + 1.05) / 0.25) ** 2)), 0]);
  const wave = zonalField(mesh, (lat, lon) => [20 * Math.cos(lat) ** 2 * Math.cos(2 * lon), 15 * Math.cos(lat) ** 2 * Math.sin(2 * lon + 0.4)]);
  const rate = 1 / (3 * 86400), dt = 3600, keep = 1 / (1 + rate * dt);
  const both = Float64Array.from(mean, (x, e) => x + wave[e]);
  dampEddies(mesh, geometry, both, 0, dt, means, rate);
  const torque = (angularMomentum(mesh, both) - angularMomentum(mesh, Float64Array.from(mean, (x, e) => x + wave[e]))) / ((keep - 1) * angularMomentum(mesh, mean));
  const waveLeft = rms(both.map((x, e) => x - keep * mean[e] - wave[e]));
  console.log(`N=32, a three-day lid friction over an hour: the axial angular momentum changes by ${torque.toFixed(4)} of a Rayleigh drag's on the zonal flow, the wave by ${waveLeft.toExponential(1)} m/s rms against the mean's ${((1 - keep) * rms(mean)).toFixed(3)}`);
  assert.ok(Math.abs(torque - 1) < 0.01, `the friction's torque is ${torque} of a Rayleigh drag's on the zonal flow`);
  assert.ok(waveLeft < 0.05 * (1 - keep) * rms(mean), `the wave moved by ${waveLeft} m/s rms`);
});

test('the lid friction profiles averaged over each layer\'s mass reach only the layers above 65 km', () => {
  for (const name of ['bl36', 'bl34']) {
    const levels = sigmaInterfaces(name);
    for (const [profile, points] of Object.entries(LID_FRICTION)) {
      const rates = lidFrictionRates(levels, profile);
      assert.ok(rates[0] > 0 && rates.slice(1).every((r) => r === 0), `${name} ${profile}: ${Array.from(rates.slice(0, 3))}`);
      assert.ok(rates[0] < 1 / (points[points.length - 1][1] * 86400), `${name} ${profile}: the layer's rate under the profile's largest`);
    }
    assert.ok(lidFrictionRates(levels, null).every((r) => r === 0));
  }
  const bl36 = sigmaInterfaces('bl36'), at = (p) => 7e3 * Math.log(101325 / p);
  const top = bl36[1] * 101325, below = (p) => (at(p) < 65e3 ? 1 : 0);
  let share = 0;
  for (let n = 0; n < 1000; n++) share += 1 - below((n + 0.5) * top / 1000);
  share /= 1000;
  const flat = lidFrictionRates(bl36, [[65e3, 4]]);
  assert.ok(Math.abs(flat[0] * 4 * 86400 - share) < 1e-3, `a uniform four-day rate above 65 km over ${(100 * share).toFixed(1)} % of the top layer's mass: ${flat[0] * 4 * 86400}`);
  assert.equal(lidFrictionFor('bl34'), null);
  assert.equal(lidFrictionFor('cam26'), null);
  assert.throws(() => lidFrictionRates(bl36, 'none of these'));
});

test('the sponge falls linearly in σ from its rate at the top to zero at its depth', () => {
  const rates = spongeRates(Float64Array.from([0.001, 0.005, 0.009, 0.011]), 0.01, 2);
  assert.deepEqual(Array.from(rates, (r) => +(r * 2 * 86400).toFixed(6)), [0.9, 0.5, 0.1, 0]);
  assert.ok(spongeRates(Float64Array.from([0.001]), 0.01, 0).every((r) => r === 0));
});

test('the sponge begins at 78 Pa on bl36 and at σ 0.005 on the other grids, and the model takes its grid\'s onset', () => {
  assert.equal(spongeSigmaFor('bl36'), 78 / 101325);
  for (const name of ['bl34', 'cam26', null]) assert.equal(spongeSigmaFor(name), SPONGE.sigma);
  for (const name of ['bl34', 'bl36']) {
    const levels = sigmaInterfaces(name), model = createModel(new Grid(4), { ocean: false, levels });
    const expected = spongeRates(model.core.sigmaMid, spongeSigmaFor(name), SPONGE.days);
    assert.deepEqual(Array.from(model.core.spongeRates), Array.from(expected));
    for (let k = 0; k < model.core.K; k++) if (model.core.sigmaMid[k] * P0 > (name === 'bl36' ? 78 : 0.005 * P0)) assert.equal(model.core.spongeRates[k], 0);
  }
});

test('the model returns the kinetic energy the sponge and the lid friction remove as heat', () => {
  const model = createModel(new Grid(4), { ocean: false, nu4Hours: Infinity, divergenceDamping: 0, levels: sigmaInterfaces('bl36'), surface: { spongeSigma: 0.05, spongeDays: 1, lidFriction: 'rind' } });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const { state, mesh, core } = model;
  const { K, C, E, dSigma, g, cp } = core.diagnostics;
  assert.ok(core.spongeRates[0] > 0 && core.spongeRates[K - 1] === 0 && core.spongeMeanRates[0] > 0);
  let seed = 3;
  const rnd = () => { seed = (seed * 16807) % 2147483647; return seed / 2147483647; };
  for (let x = 0; x < state[2].length; x++) state[2][x] += 15 * (rnd() - 0.5);
  state[0].fill(P0);
  core.diagnose(state[0], state[1], state[4], state[5]);
  const kinetic = () => {
    let sum = 0;
    for (let k = 0; k < K; k++) for (let e = 0; e < E; e++) sum += 0.25 * mesh.dcEdge[e] * mesh.dvEdge[e] * (state[0][mesh.cellsOnEdge[2 * e]] + state[0][mesh.cellsOnEdge[2 * e + 1]]) * dSigma[k] / g * state[2][k * E + e] ** 2;
    return sum;
  };
  const kineticBefore = kinetic(), thetaBefore = Float64Array.from(state[1]);
  model.phases.closure(0, K, 3600);
  const lost = kineticBefore - kinetic();
  model.phases.dissipate(0, C);
  let heat = 0;
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) heat += mesh.areaCell[i] * state[0][i] * dSigma[k] / g * cp * (state[1][k * C + i] - thetaBefore[k * C + i]) * core.diagnostics.exnerLayer[k * C + i];
  assert.ok(lost > 0);
  assert.ok(Math.abs(heat - lost) < 1e-9 * lost, `${heat} J of heat against ${lost} J of kinetic energy removed`);
});

test('twenty GPU steps with the sponge, and with the lid friction, track twenty CPU steps, and the top treatment moves the flow by far more than the engines differ', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const run = (sponge, friction = null) => {
    const model = createModel(new Grid(8), { physics: false, nu4Hours: Infinity, divergenceDamping: 0, core: sponge ? { spongeRates: sponge, spongeMeanRates: friction } : {} });
    const init = initializeState(model, {});
    for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
    const E = model.mesh.nEdges;
    for (let k = 0; k < 4; k++) model.state[2].set(zonalField(model.mesh, (lat, lon) => [40 * Math.cos(lat) + 15 * Math.cos(lat) ** 2 * Math.cos(3 * lon), 10 * Math.cos(lat) ** 2 * Math.sin(3 * lon)]), k * E);
    return model;
  };
  const probe = run(null), sponge = spongeRates(probe.core.sigmaMid, 0.02, 0.5);
  for (const friction of [null, Float64Array.from(sponge, (r, k) => (k < 2 ? 1 / (2 * 86400) : 0))]) {
    const model = run(sponge, friction), free = run(null);
    const meanTheta = Float64Array.from({ length: model.core.K }, (_, k) => { let s = 0, a = 0; for (let i = 0; i < model.mesh.nCells; i++) { s += model.mesh.areaCell[i] * model.state[1][k * model.mesh.nCells + i]; a += model.mesh.areaCell[i]; } return s / a; });
    const gpu = await createGpuCore(model.mesh, { dragCoefficient: 0, topDragDays: 0, spongeRates: sponge, spongeMeanRates: friction, gravityWaves: false, referenceTheta: meanTheta });
    gpu.upload(model.state);
    const dt = 900;
    for (let n = 0; n < 20; n++) { model.step(dt); free.step(dt); await gpu.step(dt); }
    const got = await gpu.download();
    let uDiff = 0, uMoved = 0;
    for (let x = 0; x < got[2].length; x++) { uDiff = Math.max(uDiff, Math.abs(model.state[2][x] - got[2][x])); uMoved = Math.max(uMoved, Math.abs(model.state[2][x] - free.state[2][x])); }
    console.log(`20 steps at N=8 with a half-day sponge over σ < 0.02 (${Array.from(sponge).filter((r) => r > 0).length} layers)${friction ? ' and a two-day friction on the top two layers\' zonal mean' : ''}: the top treatment moves the wind by up to ${uMoved.toFixed(2)} m/s, the engines differ by ${uDiff.toExponential(1)}`);
    assert.ok(uDiff < 0.01 * uMoved, `engines differ by ${uDiff} m/s, the top treatment moved the wind by ${uMoved}`);
  }
});
