import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { boundaryLayer, mixColumn, columnScratch, surfaceBuoyancy, langmuirInverse, RH18_ENERGY, MSTAR_SCHEMES, LANGMUIR_SCHEMES, LANGMUIR_RESCALE } from '../js/ocean/column.module.js';
import { SEAWATER } from '../js/ocean/seawater.module.js';

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const K = 60, ustar = 0.0098, dt = 2700, f = 1.3e-4;
const column = () => ({ h: new Float64Array(K).fill(5), T: new Float64Array(K).fill(272), S: Float64Array.from({ length: K }, (_, k) => 34 + (k < 20 ? 0 : 0.02 * (k - 19))), kappa: new Float64Array(K), scratch: columnScratch(K) });
const energy = (depth, buoyancy) => {
  const { cN1, cN2, cN3, cS1, cS2, cC, cap } = RH18_ENERGY, u3 = ustar ** 3, decay = cN2 * Math.exp(cN3 * depth * f / ustar), neutral = cN1 * decay / (1 + decay);
  const loss = Math.max(0, -buoyancy), gain = Math.max(0, buoyancy);
  const mstar = Math.min(cap, neutral * (1 - cC * loss / (loss + 2 * neutral * u3 / depth)) + cS1 * (gain * gain * depth / (u3 * ustar * ustar * f)) ** cS2);
  return SEAWATER.rho0 * dt * mstar * u3;
};

test('with mstarScheme rh18 the wind\'s energy is Reichl and Hallberg\'s (2018) m*N C + m*S: the neutral fit, reduced under a destabilizing surface buoyancy flux and added to under a stabilizing one', () => {
  const depths = {};
  for (const buoyancy of [0, -4e-8, -2e-7, 4e-8, 1e-7]) {
    const { h, T, S, kappa, scratch } = column();
    const layer = boundaryLayer({ h, T, S, K, ustar, f, dt, langmuir: 0, mstarScheme: 'rh18', buoyancy }, kappa, scratch);
    console.log(`B ${buoyancy}: ${layer.depth.toFixed(2)} m, ${layer.energy.toFixed(4)} J/m²`);
    close(layer.energy, energy(layer.depth, buoyancy), 1e-9 * layer.energy, `energy at B ${buoyancy}`);
    close(layer.work, layer.energy, 1e-4 * layer.energy, `the work of mixing equals the energy at B ${buoyancy}`);
    depths[buoyancy] = layer.depth;
  }
  assert.ok(depths[-2e-7] < depths[0], `strong cooling takes the wind's share down: ${depths[-2e-7]} m against ${depths[0]} m neutral`);
  assert.ok(depths[1e-7] > depths[0], `warming adds the stable term: ${depths[1e-7]} m against ${depths[0]} m neutral`);
  assert.deepEqual(MSTAR_SCHEMES, ['fixed', 'rh18']);
});

test('the fixed scheme takes no account of the surface buoyancy flux, and mixColumn forms that flux from the heat, light and salt it is given only under rh18, taking the flux it is handed otherwise', () => {
  const heat = 0.004, light = 0.001, salt = -0.01;
  const run = (extra, options) => { const { h, T, S, scratch } = column(); return mixColumn({ h, T, S, K, ustar, f, dt, ...extra, options: { langmuir: 0, ...options } }, { kappa: new Float64Array(K), nu: new Float64Array(K) }, scratch).depth; };
  const flux = surfaceBuoyancy(272, 34, heat + light, salt, dt);
  assert.ok(flux > 0, `warming and freshening are stabilizing: ${flux}`);
  const fluxes = { heat, light, salt };
  assert.equal(run(fluxes, { mstarScheme: 'rh18' }), run({ ...fluxes, buoyancy: flux }, { mstarScheme: 'rh18' }));
  assert.notEqual(run(fluxes, { mstarScheme: 'rh18' }), run({ ...fluxes, buoyancy: 0 }, { mstarScheme: 'rh18' }));
  assert.equal(run({ ...fluxes, buoyancy: flux }, {}), run({ ...fluxes, buoyancy: -flux }, {}));
  assert.equal(run(fluxes, {}), run({ ...fluxes, buoyancy: 0 }, {}));
});

test('the ocean keeps the heat and salt its surface gave it over the step, the light that passed the top layer and the brine of ice grown included, and leaves them out of the fixed m*', () => {
  const mesh = buildMesh(new Grid(16)), C = mesh.nCells, DT = 1350;
  const stress = new Float64Array(mesh.nEdges).fill(0.02);
  const give = (options) => {
    const ocean = createOcean(mesh, { everySteps: 1, ...options }), start = new Float64Array(C).fill(285), ice = new Float64Array(C), flux = new Float64Array(C);
    ocean.initialize(start, ice);
    const warm = Float64Array.from(start, (t, i) => t + (mesh.latCell[i] > 0 ? 0.5 : -0.5)), grown = Float64Array.from(ice, (v, i) => (mesh.latCell[i] < -1.2 ? 0.1 : 0));
    ocean.accumulate(new Float64Array(C).fill(1e-5), new Float64Array(C), DT);
    ocean.advance(warm, grown, flux, stress, DT);
    return ocean;
  };
  const ocean = give({ mstarScheme: 'rh18' }), rest = give({ mstarScheme: 'fixed', langmuirScheme: 'factor' });
  let warmed = 0, cooled = 0, brined = 0;
  for (let i = 0; i < C; i++) {
    if (!ocean.cellOcean[i] || ocean.shelf?.[i]) continue;
    const south = mesh.latCell[i] < -1.2;
    assert.equal(rest.heatIn[i], 0);
    if (south) { assert.ok(ocean.saltIn[i] > 2 && ocean.saltIn[i] < 3.5, `brine ${ocean.saltIn[i]} psu m at cell ${i}`); brined++; continue; }
    const side = mesh.latCell[i] > 0 ? 1 : -1;
    const thick = ocean.heatIn[i] / (side * 0.5);
    assert.ok(thick > 2.5 && thick < 3.6, `the heat given to cell ${i} spreads over its ${thick} m top layer`);
    close(ocean.saltIn[i], 35 * 0.0135 / 1000, 4e-5, `evaporated salt at cell ${i}`);
    if (side > 0) warmed++; else cooled++;
  }
  assert.ok(warmed > 20 && cooled > 20 && brined > 5, `${warmed} warmed, ${cooled} cooled, ${brined} brined cells`);
  assert.throws(() => createOcean(mesh, { everySteps: 1, mstarScheme: 'bulk' }), /mstarScheme is one of fixed, rh18/);
});

test('with langmuirScheme rescale the energy is the rh18 m* times min(5, 1 + waves 0.447 La^-1.33) of the wind\'s surface-layer Langmuir number at the depth reached, and langmuir 0 or waves 0 leave the enhancement out', () => {
  const run = (extra) => { const { h, T, S, kappa, scratch } = column(); return boundaryLayer({ h, T, S, K, ustar, f, dt, mstarScheme: 'rh18', ...extra }, kappa, scratch); };
  const enhance = (depth, waves) => Math.min(LANGMUIR_RESCALE.cap, 1 + waves * LANGMUIR_RESCALE.coefficient * (1 / Math.sqrt(langmuirInverse(ustar, Math.max(depth, 1)))) ** LANGMUIR_RESCALE.exponent);
  const plain = run({ langmuir: 0 }), off = run({ langmuir: 0, langmuirScheme: 'rescale' }), bare = run({ langmuirScheme: 'rescale', waves: 0 });
  assert.equal(off.depth, plain.depth);
  assert.equal(bare.depth, plain.depth);
  for (const waves of [1, 0.5]) {
    const layer = run({ langmuirScheme: 'rescale', waves }), factor = enhance(layer.depth, waves);
    console.log(`waves ${waves}: ${layer.depth.toFixed(2)} m against ${plain.depth.toFixed(2)} m without, enhancement ${factor.toFixed(3)}`);
    assert.ok(factor > 1 + 0.2 * waves && factor < LANGMUIR_RESCALE.cap, `the enhancement at ${layer.depth} m is ${factor}`);
    close(layer.energy, energy(layer.depth, 0) * factor, 1e-9 * layer.energy, `energy at waves ${waves}`);
    close(layer.work, layer.energy, 1e-4 * layer.energy, `the work of mixing equals the energy at waves ${waves}`);
    assert.ok(layer.depth > plain.depth, `the enhancement deepens the layer: ${layer.depth} m against ${plain.depth} m`);
  }
  assert.ok(run({ langmuirScheme: 'rescale', waves: 0.5 }).depth < run({ langmuirScheme: 'rescale' }).depth);
  assert.deepEqual(LANGMUIR_SCHEMES, ['factor', 'rescale']);
  const mesh = buildMesh(new Grid(16));
  assert.throws(() => createOcean(mesh, { everySteps: 1, mstarScheme: 'fixed', langmuirScheme: 'rescale' }), /langmuirScheme 'rescale' multiplies the m\* of mstarScheme 'rh18'/);
  assert.throws(() => createOcean(mesh, { everySteps: 1, mstarScheme: 'rh18', langmuirScheme: 'additive' }), /langmuirScheme is one of factor, rescale/);
});
