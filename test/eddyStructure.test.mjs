import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { createOcean, eddyOptions, eddyStratification, eddyStructureFactors, EDDY_STRUCTURE_FLOOR } from '../js/ocean/layered.module.js';

/*
 * The eddy diffusivity's local first-mode speed and its stratification
 * structure in the CPU ocean (ocean/layered.module.js): the column
 * arithmetic, and at N=16 the eddy phase on a sea with a continent.
 */
const g = 9.81, rho0 = 1025, DT = 1350;
const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));

function ocean(options = {}) {
  const model = createModel(new Grid(16), { topography });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  const cpu = createOcean(model.mesh, { geography: model.geography, vertical: 'hybrid', everySteps: 1, iceDrift: false, ...options });
  cpu.initialize(Float64Array.from(model.state[3]), Float64Array.from(model.state[6]));
  return { cpu, mesh: model.mesh };
}
const totals = ({ cpu, mesh }) => {
  const C = mesh.nCells, sums = [0, 0, 0];
  for (let n = 0; n < cpu.h.length; n++) { const a = mesh.areaCell[n % C]; sums[0] += a * cpu.h[n]; sums[1] += a * cpu.Q[n]; sums[2] += a * cpu.W[n]; }
  return sums;
};

test('a column of uniform N has the first-mode speed N times the depth between its first and last layer centers over π', () => {
  const n = 40, dz = 100, N2 = 1e-5, h = new Float64Array(n).fill(dz), nsq = new Float64Array(n);
  const rho = Float64Array.from({ length: n }, (_, k) => 1025 + rho0 * N2 * dz * k / g);
  const speed = eddyStratification(h, rho, n, g, rho0, nsq);
  assert.ok(Math.abs(speed - Math.sqrt(N2) * (n - 1) * dz / Math.PI) < 1e-9 * speed, `${speed}`);
  for (let k = 0; k < n - 1; k++) assert.ok(Math.abs(nsq[k] / N2 - 1) < 1e-9, `interface ${k}: ${nsq[k]}`);
});

test('classes held at their minimum add nothing to the speed and take the buoyancy frequency of the water around them', () => {
  const thick = 20, dz = 100, N2 = 2e-5, step = rho0 * N2 * dz / g, nsq = new Float64Array(thick + 3);
  const h = [], rho = [];
  for (let k = 0; k < thick; k++) {
    h.push(dz); rho.push(1025 + step * k);
    if (k === 9) for (let m = 1; m <= 3; m++) { h.push(0.5); rho.push(1025 + step * (k + m / 4)); }
  }
  const speed = eddyStratification(h, rho, h.length, g, rho0, nsq), across = dz + 1.5;
  const expected = (Math.sqrt(N2) * dz * (thick - 2) + Math.sqrt(g * step * across / rho0)) / Math.PI;
  assert.ok(Math.abs(speed / expected - 1) < 1e-9, `${speed} against ${expected}`);
  for (const k of [9, 10, 11, 12]) assert.ok(Math.abs(nsq[k] / (g * step / (rho0 * across)) - 1) < 1e-9, `interface ${k}: ${nsq[k]}`);
  assert.ok(Math.abs(nsq[8] / N2 - 1) < 1e-9 && Math.abs(nsq[13] / N2 - 1) < 1e-9);
});

test('the structure is one down to the strongest N² below the taper depth and N²/N²_peak beneath it, never under the floor', () => {
  const n = 12, h = new Float64Array(n).fill(50), out = new Float64Array(n);
  const nsq = Float64Array.from([0, 3e-4, 3e-5, 8e-5, 1e-4, 6e-5, 3e-5, 1e-5, 5e-6, 1e-6, 1e-7, 0]);
  eddyStructureFactors(h, nsq, n, 150, 0.1, out);
  const expected = [1, 1, 1, 1, 1, 0.6, 0.3, 0.1, 0.1, 0.1, 0.1, 1];
  for (let k = 0; k < n; k++) assert.ok(Math.abs(out[k] - expected[k]) < 1e-12, `interface ${k}: ${out[k]}`);
  eddyStructureFactors(h, new Float64Array(n), n, 150, 0.1, out);
  assert.deepEqual(Array.from(out), new Array(n).fill(1));
});

test('eddyOptions takes a speed or local, a positive power, uniform or stratification and a floor from 0 to 1', () => {
  assert.deepEqual(eddyOptions({}), { local: false, shape: false, columns: false, waveSpeed: 2, power: 2, floor: EDDY_STRUCTURE_FLOOR });
  assert.deepEqual(eddyOptions({ eddyWaveSpeed: 'local', eddyResolutionPower: 4 }), { local: true, shape: false, columns: true, waveSpeed: 2, power: 4, floor: EDDY_STRUCTURE_FLOOR });
  for (const bad of [{ eddyWaveSpeed: 'global' }, { eddyWaveSpeed: 0 }, { eddyResolutionPower: 0 }, { eddyStructure: 'n2' }, { eddyStructureFloor: 1.5 }]) assert.throws(() => eddyOptions(bad), JSON.stringify(bad));
});

test('the default eddy options step the same bits as their explicit values', () => {
  const a = ocean(), b = ocean({ eddyWaveSpeed: 2, eddyResolutionPower: 2, eddyStructure: 'uniform' });
  a.cpu.phases.eddy(DT); b.cpu.phases.eddy(DT);
  for (const key of ['h', 'Q', 'W']) assert.ok(Buffer.from(a.cpu[key].buffer).equals(Buffer.from(b.cpu[key].buffer)), key);
});

test('with the local speed and the stratification structure the eddy phase keeps the water, heat and salt', () => {
  const run = ocean({ eddyWaveSpeed: 'local', eddyResolutionPower: 4, eddyStructure: 'stratification' }), before = totals(run);
  for (let n = 0; n < 4; n++) run.cpu.phases.eddy(DT);
  const after = totals(run);
  for (let q = 0; q < 3; q++) assert.ok(Math.abs(after[q] - before[q]) <= 1e-12 * Math.abs(before[q]), `${['water', 'heat', 'salt'][q]}: ${before[q]} to ${after[q]}`);
});

test('the local speed gives equatorial edges less of κ₀ than midlatitude ones, and every edge some', () => {
  const { cpu, mesh } = ocean({ eddyWaveSpeed: 'local', eddyResolutionPower: 4 }), kappa = cpu.eddyDiffusivities();
  let low = 0, nLow = 0, mid = 0, nMid = 0;
  for (let e = 0; e < mesh.nEdges; e++) {
    if (!cpu.edgeOcean[e]) continue;
    assert.ok(kappa[e] > 0 && kappa[e] <= 1000, `edge ${e}: ${kappa[e]}`);
    const lat = Math.abs(mesh.latEdge[e]) * 180 / Math.PI;
    if (lat < 5) { low += kappa[e]; nLow++; } else if (lat > 35 && lat < 55) { mid += kappa[e]; nMid++; }
  }
  assert.ok(nLow > 0 && nMid > 0 && low / nLow < mid / nMid, `equator ${low / nLow}, midlatitudes ${mid / nMid}`);
});

test('the stratification structure cuts the eddy fluxes in the weakly stratified deep water', () => {
  const uniform = ocean(), shaped = ocean({ eddyStructure: 'stratification' });
  uniform.cpu.phases.eddy(DT); shaped.cpu.phases.eddy(DT);
  const deep = ({ cpu, mesh }) => {
    const E = mesh.nEdges, C = mesh.nCells, flux = new Float64Array(cpu.shared.eddyFlux), L = cpu.layers;
    let sum = 0;
    for (let e = 0; e < E; e++) {
      const a = mesh.cellsOnEdge[2 * e];
      let z = 0;
      for (let k = 0; k < L; k++) { z += cpu.h[k * C + a]; if (z > 1500) sum += Math.abs(flux[k * E + e]); }
    }
    return sum;
  };
  const u = deep(uniform), s = deep(shaped);
  assert.ok(u > 0 && s < 0.5 * u, `uniform ${u}, stratification ${s}`);
});
