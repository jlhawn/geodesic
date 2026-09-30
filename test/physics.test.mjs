import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore, P0, CP_DRY, R_DRY } from '../js/dynamics/sigmaCore.module.js';
import { curl, divergence } from '../js/dynamics/operators.module.js';
import { createRadiation, sunDirection, AXIAL_TILT, DAY, YEAR, waterVaporAbsorptivity, adiabaticWaterLapse, inversionStrength, entrainmentIndex, ringMean } from '../js/physics/radiation.module.js';
import { smoothCells } from '../js/levels.module.js';
import { LATENT_HEAT, saturationHumidity, liftingCondensationLevel } from '../js/physics/moist.module.js';
import { createSurface } from '../js/physics/surface.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { createMixedLayer, dycomsLongwave } from '../js/physics/mixedLayer.module.js';

const N = +(process.env.PHYSICS_TEST_N ?? 6);
const grid = new Grid(N);
const mesh = buildMesh(grid);
const core = createSigmaCore(mesh);
const { K, C, E } = core.diagnostics;
const EPS = 1e-12;
const OVERCAST = { cloudCover: 'overcast' };

function random(seed) {
  let s = seed >>> 0;
  return () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
}

function sampleState(seed) {
  const rnd = random(seed);
  const pi = Float64Array.from({ length: C }, () => P0 * (0.95 + 0.1 * rnd()));
  const theta = new Float64Array(K * C);
  for (let i = 0; i < C; i++) for (let k = 0; k < K; k++) theta[k * C + i] = 280 + 250 * (1 - core.sigmaMid[k]) ** 2 + 5 * (rnd() - 0.5);
  const surfaceT = Float64Array.from({ length: C }, () => 270 + 30 * rnd());
  const u = Float64Array.from({ length: K * E }, () => 20 * (rnd() - 0.5));
  return [pi, theta, u, surfaceT];
}

test('the sun follows the equinox, the solstice and the daily rotation', () => {
  const s0 = sunDirection(0);
  assert.ok(Math.abs(s0[0] - 1) < EPS && Math.abs(s0[1]) < EPS && Math.abs(s0[2]) < EPS);
  const solstice = sunDirection(YEAR / 4);
  assert.ok(Math.abs(solstice[2] - Math.sin(AXIAL_TILT)) < 1e-9);
  const noonLater = sunDirection(DAY / 2);
  assert.ok(noonLater[0] < -0.9999 && Math.abs(noonLater[1]) < 1e-9);
});

test('radiation column: layer and surface fluxes sum to absorbed solar minus outgoing longwave', () => {
  const radiation = createRadiation(mesh, core);
  const [pi, theta, u, surfaceT] = sampleState(3);
  core.diagnose(pi, theta);
  radiation.setTime(0.3 * DAY);
  let worst = 0;
  for (let i = 0; i < C; i++) {
    const surfaceFlux = radiation.column(i, pi[i], theta, surfaceT[i], 5);
    let layers = 0, scale = 0;
    for (let k = 0; k < K; k++) { layers += radiation.layerFlux[k]; scale += Math.abs(radiation.layerFlux[k]); }
    const residual = layers + surfaceFlux - (radiation.budget.absorbedSolar - radiation.budget.outgoingLongwave);
    worst = Math.max(worst, Math.abs(residual) / (scale + Math.abs(surfaceFlux)));
    assert.ok(radiation.budget.outgoingLongwave > 0 && radiation.budget.outgoingLongwave < 600);
  }
  assert.ok(worst < EPS);
  const pi0 = new Float64Array(C).fill(P0);
  core.diagnose(pi0, theta);
  radiation.column(0, P0, theta, 288, 5);
  let transmitted = 1;
  for (let k = 0; k < K; k++) transmitted *= 1 - radiation.emissivity[k];
  assert.ok(Math.abs(transmitted - Math.exp(-radiation.opticalDepth(mesh.latCell[0]))) < 1e-12);
  assert.ok(radiation.opticalDepth(0) > radiation.opticalDepth(Math.PI / 2));
});

test('radiation warms the column where it absorbs and cools the slab where it emits', () => {
  const radiation = createRadiation(mesh, core);
  const [pi, theta, u, surfaceT] = sampleState(5);
  surfaceT.fill(320);
  core.diagnose(pi, theta);
  radiation.setTime(DAY / 4);
  const out = [new Float64Array(C), new Float64Array(K * C), new Float64Array(K * E), new Float64Array(C)];
  const wind = new Float64Array(C).fill(5);
  radiation.apply([pi, theta, u, surfaceT], out, wind);
  let night = -1;
  for (let i = 0; i < C; i++) if (radiation.insolation(i) === 0) { night = i; break; }
  assert.ok(night >= 0);
  assert.ok(radiation.surfaceFlux[night] < 0);
  assert.ok(out[1].every(Number.isFinite) && radiation.surfaceFlux.every(Number.isFinite));
});

test('convective adjustment leaves a stable column and conserves enthalpy', () => {
  const surface = createSurface(mesh, core);
  const [pi, theta] = sampleState(7);
  for (let i = 0; i < C; i++) for (let k = 0; k < K; k++) theta[k * C + i] = 300 + 40 * core.sigmaMid[k] + 20 * Math.sin(3 * k + i);
  core.diagnose(pi, theta);
  const { exnerLayer, dSigma } = core.diagnostics;
  const enthalpy = (i) => { let h = 0; for (let k = 0; k < K; k++) h += CP_DRY * theta[k * C + i] * exnerLayer[k * C + i] * dSigma[k]; return h; };
  const before = Float64Array.from({ length: C }, (_, i) => enthalpy(i));
  const mixes = surface.convectiveAdjustment(pi, theta);
  assert.ok(mixes > 0);
  for (let i = 0; i < C; i++) {
    assert.ok(Math.abs(enthalpy(i) - before[i]) / before[i] < EPS);
    for (let k = 0; k < K - 1; k++) assert.ok(theta[(k + 1) * C + i] <= theta[k * C + i] * (1 + 1e-9));
  }
});

test('surface drag only removes kinetic energy', () => {
  const surface = createSurface(mesh, core);
  const state = sampleState(11);
  const [pi, theta, u] = state;
  core.diagnose(pi, theta);
  surface.lowestWindSpeed(u);
  const out = [new Float64Array(C), new Float64Array(K * C), new Float64Array(K * E), new Float64Array(C)];
  surface.apply(state, out);
  for (let k = 0; k < K; k++) {
    let work = 0;
    for (let e = 0; e < E; e++) work += mesh.dcEdge[e] * mesh.dvEdge[e] * u[k * E + e] * out[2][k * E + e];
    if (k === K - 1) assert.ok(work < 0); else assert.equal(work, 0);
  }
});

function cellKineticWeights(m, pi, u, k, dSigma, g) {
  let sum = 0;
  for (let e = 0; e < m.nEdges; e++) {
    const a = m.cellsOnEdge[2 * e], b = m.cellsOnEdge[2 * e + 1];
    sum += 0.25 * m.dcEdge[e] * m.dvEdge[e] * (pi[a] + pi[b]) * dSigma[k] / g * u[e] * u[e];
  }
  return sum;
}

test('the heat the drags return equals the kinetic energy they remove', () => {
  const surface = createSurface(mesh, core, { topDragDays: 5, topSigma: 0.05 });
  const state = sampleState(13);
  const [pi, theta, u] = state;
  core.diagnose(pi, theta);
  const { exnerLayer, dSigma, g } = core.diagnostics;
  surface.lowestWindSpeed(u);
  const out = [new Float64Array(C), new Float64Array(K * C), new Float64Array(K * E), new Float64Array(C)];
  surface.applyLayers(state, out);
  let kineticRate = 0;
  for (let k = 0; k < K; k++) for (let e = 0; e < E; e++) {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
    kineticRate += 0.5 * mesh.dcEdge[e] * mesh.dvEdge[e] * (pi[a] + pi[b]) * dSigma[k] / g * u[k * E + e] * out[2][k * E + e];
  }
  surface.heatLayers(state, out);
  let heatRate = 0;
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) heatRate += mesh.areaCell[i] * pi[i] * dSigma[k] / g * CP_DRY * out[1][k * C + i] * exnerLayer[k * C + i];
  assert.ok(kineticRate < 0);
  assert.ok(Math.abs(heatRate + kineticRate) < 1e-10 * -kineticRate, `heat ${heatRate} W against kinetic ${kineticRate} W`);
});

test('the closure and the boundary-layer mixing return the kinetic energy they remove as heat', () => {
  const model = createModel(new Grid(4), { ocean: false });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const { state, mesh: m, core: c } = model;
  const { K: layers, C: cells, E: edges, dSigma, g, cp } = c.diagnostics;
  const rnd = random(5);
  for (let x = 0; x < state[2].length; x++) state[2][x] += 15 * (rnd() - 0.5);
  model.surface.lowestWindSpeed(state[2]);
  model.phases.physics(0, cells, 900, model.totals);
  state[0].fill(P0);
  c.diagnose(state[0], state[1], state[4], state[5]);
  const energy = () => {
    let sum = 0;
    for (let k = 0; k < layers; k++) {
      sum += cellKineticWeights(m, state[0], state[2].subarray(k * edges, (k + 1) * edges), k, dSigma, g);
      for (let i = 0; i < cells; i++) sum += m.areaCell[i] * state[0][i] * dSigma[k] / g * cp * state[1][k * cells + i] * c.diagnostics.exnerLayer[k * cells + i];
    }
    return sum;
  };
  const kinetic = () => { let sum = 0; for (let k = 0; k < layers; k++) sum += cellKineticWeights(m, state[0], state[2].subarray(k * edges, (k + 1) * edges), k, dSigma, g); return sum; };
  const before = energy(), kineticBefore = kinetic();
  model.phases.closure(0, layers, 900);
  model.phases.mixMomentum(0, edges, 900);
  const lost = kineticBefore - kinetic();
  model.phases.dissipate(0, cells);
  assert.ok(lost > 0);
  assert.ok(Math.abs(energy() - before) < 1e-9 * lost, `energy changed by ${energy() - before} J against ${lost} J of kinetic energy removed`);
});

test('divergence damping takes kinetic energy from the divergent flow alone and the model returns it as heat: mass, vorticity and total energy are unchanged', () => {
  const model = createModel(new Grid(4), { ocean: false, nu4Hours: Infinity, divergenceDamping: 0.05 });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const { state, mesh: m, core: c } = model;
  const { K: layers, C: cells, E: edges, V: vertices, dSigma, g, cp } = c.diagnostics;
  assert.equal(c.nu4, 0);
  assert.throws(() => createSigmaCore(m, { divergenceDamping: 0.05 }));
  const rnd = random(7);
  for (let x = 0; x < state[2].length; x++) state[2][x] += 15 * (rnd() - 0.5);
  model.surface.lowestWindSpeed(state[2]);
  model.phases.physics(0, cells, 900, model.totals);
  state[0].fill(P0);
  c.diagnose(state[0], state[1], state[4], state[5]);
  const layerKinetic = (k) => cellKineticWeights(m, state[0], state[2].subarray(k * edges, (k + 1) * edges), k, dSigma, g);
  const energy = () => {
    let sum = 0;
    for (let k = 0; k < layers; k++) {
      sum += layerKinetic(k);
      for (let i = 0; i < cells; i++) sum += m.areaCell[i] * state[0][i] * dSigma[k] / g * cp * state[1][k * cells + i] * c.diagnostics.exnerLayer[k * cells + i];
    }
    return sum;
  };
  const kinetic = () => { let sum = 0; for (let k = 0; k < layers; k++) sum += layerKinetic(k); return sum; };
  const fields = () => {
    const vorticity = new Float64Array(layers * vertices), div = new Float64Array(layers * cells);
    for (let k = 0; k < layers; k++) {
      curl(m, state[2].subarray(k * edges, (k + 1) * edges), vorticity.subarray(k * vertices, (k + 1) * vertices));
      divergence(m, state[2].subarray(k * edges, (k + 1) * edges), div.subarray(k * cells, (k + 1) * cells));
    }
    return { vorticity, div };
  };
  const rms = (x) => Math.sqrt(x.reduce((a, b) => a + b * b, 0) / x.length);
  const before = energy(), kineticBefore = kinetic(), pi = Float64Array.from(state[0]), start = fields();
  model.phases.closure(0, layers, 900);
  const lost = kineticBefore - kinetic(), end = fields();
  model.phases.dissipate(0, cells);
  let turned = 0;
  for (let x = 0; x < start.vorticity.length; x++) turned = Math.max(turned, Math.abs(end.vorticity[x] - start.vorticity[x]));
  const scale = start.vorticity.reduce((a, b) => Math.max(a, Math.abs(b)), 0);
  console.log(`one step at c = 0.05: the divergence falls from rms ${rms(start.div).toExponential(2)} to ${rms(end.div).toExponential(2)} /s, ${(100 * lost / kineticBefore).toFixed(2)} % of the kinetic energy goes to heat, the vorticity moves by at most ${(turned / scale).toExponential(1)} of its largest value`);
  assert.ok(lost > 0 && rms(end.div) < 0.9 * rms(start.div), `kinetic energy lost ${lost}, divergence ${rms(start.div)} → ${rms(end.div)}`);
  assert.ok(turned < 1e-12 * scale, `vorticity moved by ${turned} against ${scale}`);
  assert.deepEqual(state[0], pi);
  assert.ok(Math.abs(energy() - before) < 1e-9 * lost, `energy changed by ${energy() - before} J against ${lost} J of kinetic energy removed`);
});

test('the assembled model steps a uniform atmosphere without blowing up', () => {
  const model = createModel(new Grid(4));
  const init = initializeState(model, { seedAmplitude: 0, geostrophic: false });
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  model.state[0].fill(P0);
  for (let n = 0; n < 20; n++) model.step(600);
  const d = model.diagnostics();
  assert.ok(Number.isFinite(d.maxWind) && d.maxWind < 30);
  assert.ok(Math.abs(d.mass - P0) / P0 < 1e-12);
  assert.ok(d.absorbedSolar > 150 && d.absorbedSolar < 350, `absorbed solar ${d.absorbedSolar}`);
  assert.ok(d.outgoingLongwave > 100 && d.outgoingLongwave < 400);
});

test('cloud water reflects sunlight and closes the window: lower OLR, higher reflection, exact closure', () => {
  const radiation = createRadiation(mesh, core, OVERCAST);
  const [pi, theta, u, surfaceT] = sampleState(3);
  const q = new Float64Array(K * C), qc = new Float64Array(K * C);
  core.diagnose(pi, theta, q, qc);
  radiation.setTime(0);
  let day = -1;
  for (let i = 0; i < C; i++) if (radiation.insolation(i) > 300) { day = i; break; }
  assert.ok(day >= 0);
  const clear = { ...evaluate(radiation, day, pi, theta, surfaceT, q, qc) };
  for (let k = 12; k < 16; k++) qc[k * C + day] = 5e-4;
  const cloudy = { ...evaluate(radiation, day, pi, theta, surfaceT, q, qc) };
  assert.ok(cloudy.reflected > clear.reflected + 50, `reflected ${clear.reflected} → ${cloudy.reflected}`);
  assert.ok(cloudy.olr < clear.olr - 20, `OLR ${clear.olr} → ${cloudy.olr}`);
  assert.ok(cloudy.closure < EPS && clear.closure < EPS);
});

test('with cloudCover: pdf a cloudy layer covers the part of a uniform total-water distribution of half-width (1 − RHc) qs above saturation, RHc 0.85 inside the boundary layer and 0.8 above: its shortwave is the cover-weighted blend of the clear column and the column with W / f, its emissivity f (1 − exp(−a W / f)), and the column closes', () => {
  const pdf = createRadiation(mesh, core), overcast = createRadiation(mesh, core, OVERCAST);
  pdf.setTime(0); overcast.setTime(0);
  const noon = brightest(pdf), [pi, theta, , surfaceT] = sampleState(3), q = new Float64Array(K * C), qc = new Float64Array(K * C);
  pi[noon] = P0;
  core.diagnose(pi, theta, q, qc);
  const { exnerLayer, sigmaMid, geopotential, g } = core.diagnostics, k = K - 4, idx = k * C + noon;
  const qs = saturationHumidity(theta[idx] * exnerLayer[idx], P0 * sigmaMid[k]);
  q[idx] = qs; qc[idx] = 0.02 * qs;
  const height = (geopotential[idx] - geopotential[(K - 1) * C + noon]) / g;
  const run = (r, depth, water = qc) => {
    const flux = r.column(noon, P0, theta, surfaceT[noon], 5, r.opticalDepth(mesh.latCell[noon]), r.insolation(noon), q[(K - 1) * C + noon], q, water, 0.07, 0.07, 1, 1.5e-3, 0, depth);
    let layers = 0;
    for (let n = 0; n < K; n++) layers += r.layerFlux[n];
    return { ...r.budget, closure: Math.abs(layers + flux + LATENT_HEAT * r.budget.evaporation - (r.budget.absorbedSolar - r.budget.outgoingLongwave)) };
  };
  for (const [depth, rhc] of [[height + 100, 0.85], [0.5 * height, 0.8]]) {
    const f = 0.5 + 0.02 / (2 * (1 - rhc)), thick = Float64Array.from(qc, (x) => x / f);
    const partial = run(pdf, depth), inCloud = run(overcast, depth, thick), full = run(overcast, depth);
    console.log(`cloud of ${(1000 * qc[idx]).toFixed(2)} g/kg in a saturated layer ${height.toFixed(0)} m up, RHc ${rhc}: cover ${f.toFixed(3)}, reflectance ${partial.cloudReflectance.toFixed(4)} against ${full.cloudReflectance.toFixed(4)} overcast`);
    assert.ok(Math.abs(partial.cloudReflectance - f * inCloud.cloudReflectance) < 1e-12, `reflectance ${partial.cloudReflectance} against ${f} × ${inCloud.cloudReflectance}`);
    assert.ok(partial.cloudReflectance < 0.8 * full.cloudReflectance, 'thin cloud in a partly humid layer reflects less than overcast');
    assert.ok(partial.outgoingLongwave > full.outgoingLongwave, 'and closes less of the window');
    assert.ok(partial.closure < 1e-9 * partial.absorbedSolar, `closure ${partial.closure}`);
  }
  q[idx] = qs * 1.25;
  const wet = run(pdf, 0), wetOvercast = run(overcast, 0);
  assert.equal(wet.cloudReflectance, wetOvercast.cloudReflectance, 'total water at qs + w and beyond is overcast');
});

test('the surface sees the direct beam in clear sky and diffuse light under thick cloud', () => {
  const radiation = createRadiation(mesh, core, OVERCAST);
  const [pi, theta, u, surfaceT] = sampleState(3);
  const q = new Float64Array(K * C), qc = new Float64Array(K * C);
  core.diagnose(pi, theta, q, qc);
  radiation.setTime(0);
  let day = -1;
  for (let i = 0; i < C; i++) if (radiation.insolation(i) > 300) { day = i; break; }
  const absorbed = (direct, diffuse) => {
    radiation.column(day, pi[day], theta, surfaceT[day], 5, radiation.opticalDepth(mesh.latCell[day]), radiation.insolation(day), q[(K - 1) * C + day], q, qc, direct, diffuse);
    return radiation.budget.absorbedSolar;
  };
  const clearBright = absorbed(0.3, 0.06), clearDark = absorbed(0.06, 0.06);
  assert.ok(clearBright / clearDark < 0.85, `clear sky: bright direct albedo absorbs ${clearBright / clearDark} of the dark`);
  for (let k = 12; k < 16; k++) qc[k * C + day] = 5e-4;
  const cloudyBright = absorbed(0.3, 0.06), cloudyDark = absorbed(0.06, 0.06);
  assert.ok(cloudyBright / cloudyDark > 0.97, `thick cloud: bright direct albedo absorbs ${cloudyBright / cloudyDark} of the dark`);
  assert.ok(cloudyDark < 0.5 * clearDark, 'the cloud reflects most of the beam');
});

test('water vapour absorbs sunlight by the Lacis–Hansen curve: a humid column takes a tenth or more of the beam, most of it low down, and dry air none', () => {
  assert.equal(waterVaporAbsorptivity(0), 0);
  assert.ok(Math.abs(waterVaporAbsorptivity(1) - 0.099) < 0.002 && Math.abs(waterVaporAbsorptivity(5) - 0.154) < 0.002);
  for (let y = 0.01; y < 20; y *= 1.5) assert.ok(waterVaporAbsorptivity(1.5 * y) > waterVaporAbsorptivity(y));
  const radiation = createRadiation(mesh, core), none = createRadiation(mesh, core, { vaporAbsorption: 0 });
  const [pi, theta, u, surfaceT] = sampleState(3);
  const q = new Float64Array(K * C), qc = new Float64Array(K * C);
  core.diagnose(pi, theta, q, qc);
  radiation.setTime(0); none.setTime(0);
  let day = -1;
  for (let i = 0; i < C; i++) if (radiation.insolation(i) > 1200) { day = i; break; }
  assert.ok(day >= 0);
  const dry = evaluate(radiation, day, pi, theta, surfaceT, q, qc), dryNone = evaluate(none, day, pi, theta, surfaceT, q, qc);
  assert.equal(dry.absorbed, dryNone.absorbed, 'dry air absorbs nothing');
  let water = 0;
  for (let k = 0; k < K; k++) {
    const idx = k * C + day, p = pi[day] * core.sigmaMid[k];
    q[idx] = 0.8 * saturationHumidity(theta[idx] * core.diagnostics.exnerLayer[idx], p);
    water += q[idx] * pi[day] * core.diagnostics.dSigma[k] / 9.80665;
  }
  const humid = evaluate(radiation, day, pi, theta, surfaceT, q, qc), humidNone = evaluate(none, day, pi, theta, surfaceT, q, qc);
  const extra = humid.layers.map((f, k) => f - humidNone.layers[k]);
  const taken = extra.reduce((a, b) => a + b, 0), beam = radiation.insolation(day) * 0.97;
  assert.ok(water > 20, `column water ${water} kg/m²`);
  assert.ok(taken > 0.1 * beam && taken < 0.25 * beam, `vapour takes ${taken} of ${beam} W/m²`);
  assert.ok(humid.surfaceShortwave < humidNone.surfaceShortwave - 0.9 * taken, 'what the vapour takes does not reach the surface');
  assert.ok(humid.absorbed >= humidNone.absorbed, 'the planet absorbs no less: the vapour takes light the surface would have absorbed or reflected');
  assert.ok(humid.closure < EPS);
  assert.ok(extra.every((f) => f >= -1e-9));
  const lower = extra.slice(K / 2).reduce((a, b) => a + b, 0);
  assert.ok(lower > 0.5 * taken, `the lower half of the column takes ${lower} of ${taken}`);
});

test('cloud water absorbs 1 − exp(−0.4 m²/kg × W) of the sunlight that meets it — 3.9 % at 100 g/m², 15 % at 400 g/m² — heating the cloudy layers in proportion to their water, and incident = reflected + atmosphere + surface', () => {
  const radiation = createRadiation(mesh, core, OVERCAST), scattering = createRadiation(mesh, core, { cloudSolarAbsorption: 0, ...OVERCAST });
  radiation.setTime(0); scattering.setTime(0);
  const noon = brightest(radiation), beam = radiation.insolation(noon), incident = beam * 0.97;
  const pi = new Float64Array(C).fill(P0), theta = sampleState(3)[1], q = new Float64Array(K * C), qc = new Float64Array(K * C);
  core.diagnose(pi, theta, q, qc);
  const { dSigma, g } = core.diagnostics, shares = new Map([[20, 0.2], [21, 0.3], [22, 0.5]]);
  const overcast = (path) => {
    qc.fill(0);
    for (const [k, share] of shares) qc[k * C + noon] = share * path / (P0 * dSigma[k] / g);
    const run = (r, sun) => {
      const surface = r.column(noon, P0, theta, 290, 5, r.opticalDepth(mesh.latCell[noon]), sun, 0, q, qc, 0.07);
      return { surface, layers: Array.from(r.layerFlux), ...r.budget };
    };
    const day = run(radiation, beam), night = run(radiation, 0), scattered = run(scattering, beam);
    const heating = day.layers.map((f, k) => f - scattered.layers[k]);
    const atmosphere = day.layers.reduce((a, f, k) => a + f - night.layers[k], 0), surface = day.surface - night.surface;
    return { ...day, heating, absorbed: heating.reduce((a, b) => a + b, 0), closure: Math.abs(beam - day.reflectedSolar - atmosphere - surface) / beam, scattered };
  };
  const thin = overcast(0.1), thick = overcast(0.4);
  console.log(`at noon (${beam.toFixed(0)} W/m²) a dry column's cloud absorbs ${(100 * thin.absorbed / incident).toFixed(2)} % of the ${incident.toFixed(0)} W/m² reaching it at 100 g/m² and ${(100 * thick.absorbed / incident).toFixed(2)} % at 400 g/m²; reflected ${thin.reflectedSolar.toFixed(1)} against ${thin.scattered.reflectedSolar.toFixed(1)} W/m² when it only scatters; closure ${thin.closure.toExponential(1)}, ${thick.closure.toExponential(1)}`);
  assert.ok(Math.abs(thin.absorbed / incident - 0.039) < 0.002, `100 g/m² absorbs ${thin.absorbed / incident}`);
  assert.ok(Math.abs(thick.absorbed / incident - 0.15) < 0.005, `400 g/m² absorbs ${thick.absorbed / incident}`);
  for (const c of [thin, thick]) {
    assert.ok(Math.abs(c.cloudSolar - c.absorbed) < 1e-9 * c.absorbed && Math.abs(c.atmosphereSolar - c.scattered.atmosphereSolar - c.absorbed) < 1e-9 * c.absorbed);
    c.heating.forEach((h, k) => {
      if (shares.has(k)) assert.ok(Math.abs(h - shares.get(k) * c.absorbed) < 1e-9 * c.absorbed, `layer ${k} takes ${h / c.absorbed} of the cloud's absorption`);
      else assert.equal(h, 0, `layer ${k} holds no cloud and takes none`);
    });
    assert.ok(c.reflectedSolar < c.scattered.reflectedSolar && c.surfaceShortwave < c.scattered.surfaceShortwave);
    assert.ok(c.closure < 1e-9, `closure ${c.closure}`);
  }
});

function deckColumns(options = {}) {
  const radiation = createRadiation(mesh, core, { stratus: true, mixedLayerDeck: false, ...OVERCAST, ...options }), off = createRadiation(mesh, core, { stratus: false, ...OVERCAST });
  radiation.setTime(0); off.setTime(0);
  let noon = 0;
  for (let i = 0; i < C; i++) if (radiation.insolation(i) > radiation.insolation(noon)) noon = i;
  const pi = new Float64Array(C).fill(P0), q = new Float64Array(K * C), qc = new Float64Array(K * C), unstable = new Float64Array(K * C);
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) unstable[k * C + i] = 297 + 16 * (1 - core.sigmaMid[k]) + 600 * Math.max(0, 0.6 - core.sigmaMid[k]) ** 2;
  core.diagnose(pi, unstable, q, qc);
  const { exnerLayer } = core.diagnostics;
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) q[k * C + i] = 0.7 * saturationHumidity(unstable[k * C + i] * exnerLayer[k * C + i], P0 * core.sigmaMid[k]);
  const bottom = (K - 1) * C + noon, top = radiation.stabilityLayer * C + noon;
  const withStability = (lts) => { const theta = Float64Array.from(unstable); theta[top] = theta[bottom] + lts; return theta; };
  const run = (r, theta, surfaceT, openSea, depth = 1200, cloud = qc) => {
    core.diagnose(pi, theta, q, cloud);
    const flux = r.column(noon, P0, theta, surfaceT, 5, r.opticalDepth(mesh.latCell[noon]), r.insolation(noon), q[bottom], q, cloud, 0.07, 0.06, 1, 1.5e-3, openSea, depth);
    let layers = 0, scale = 0;
    for (let k = 0; k < K; k++) { layers += r.layerFlux[k]; scale += Math.abs(r.layerFlux[k]); }
    const latent = LATENT_HEAT * r.budget.evaporation;
    const closure = Math.abs(layers + flux + latent - (r.budget.absorbedSolar - r.budget.outgoingLongwave)) / (scale + Math.abs(flux) + latent);
    const { stabilityIndex, ...budget } = r.budget;
    return { ...budget, flux, closure, layers: Array.from(r.layerFlux) };
  };
  const inversion = (theta) => {
    core.diagnose(pi, theta, q, qc);
    const { g, kappa, geopotential } = core.diagnostics;
    const lowerT = theta[bottom] * exnerLayer[bottom], upperT = theta[top] * exnerLayer[top];
    const lcl = liftingCondensationLevel(lowerT, q[bottom], P0 * core.sigmaMid[K - 1], kappa);
    const T = 0.5 * (lowerT + upperT), qs = saturationHumidity(T, 85000), vaporR = R_DRY / 0.622;
    const moist = g / CP_DRY * (1 + LATENT_HEAT * qs / (R_DRY * T)) / (1 + LATENT_HEAT * LATENT_HEAT * qs / (CP_DRY * vaporR * T * T));
    return theta[top] - theta[bottom] - (g / CP_DRY - moist) * ((geopotential[top] - geopotential[bottom]) / g - Math.max(0, CP_DRY * (lowerT - lcl.temperature) / g));
  };
  const cover = (theta) => Math.min(1, Math.max(0, 0.19 + 0.08 * (inversion(theta) - 1)));
  return { radiation, off, noon, q, qc, unstable, bottom, top, withStability, run, inversion, cover, airTemperature: unstable[bottom] * exnerLayer[bottom] };
}

test('a stable column over warm open sea carries a stratocumulus deck that reflects the noon sun and lowers the OLR; land, cold sea and stratus: false carry none', () => {
  const { radiation, off, unstable, withStability, run, inversion, cover } = deckColumns();
  const stable = withStability(20);
  const deck = run(radiation, stable, 298, 1), clear = run(radiation, unstable, 298, 1);
  console.log(`LTS 20 K (EIS ${inversion(stable).toFixed(1)} K) against ${(unstable[radiation.stabilityLayer * C] - unstable[(K - 1) * C]).toFixed(1)} K (EIS ${inversion(unstable).toFixed(1)} K) over a 298 K sea under a 1200 m boundary layer: deck ${deck.stratus.toFixed(4)} kg/m² on ${deck.stratusFraction.toFixed(3)} of the cell; noon absorbed solar ${deck.absorbedSolar.toFixed(0)} against ${clear.absorbedSolar.toFixed(0)} W/m², OLR ${deck.outgoingLongwave.toFixed(1)} against ${clear.outgoingLongwave.toFixed(1)}`);
  assert.ok(deck.stratus > 0.02 && clear.stratus === 0 && clear.stratusFraction === 0, `deck ${deck.stratus} against ${clear.stratus} kg/m²`);
  assert.ok(deck.stratusFraction > 0 && Math.abs(deck.stratusFraction - cover(stable)) < 1e-9, `cover ${deck.stratusFraction}`);
  assert.ok(deck.absorbedSolar < clear.absorbedSolar - 60, `absorbed solar ${deck.absorbedSolar} against ${clear.absorbedSolar}`);
  assert.ok(deck.outgoingLongwave < clear.outgoingLongwave, `OLR ${deck.outgoingLongwave} against ${clear.outgoingLongwave}`);
  assert.ok(deck.closure < EPS);
  const half = run(radiation, stable, 298, 0.5);
  assert.ok(Math.abs(half.stratusFraction - 0.5 * deck.stratusFraction) < 1e-15 && half.stratus === deck.stratus, 'half the cell under ice halves the cover and keeps the deck');
  assert.equal(run(radiation, stable, 298, 0).stratus, 0, 'land and full ice carry none');
  assert.equal(run(radiation, stable, 275, 1).stratus, 0, 'a sea colder than 5 °C carries none');
  const before = run(radiation, stable, 298, 0);
  assert.deepEqual(run(off, stable, 298, 1), before, 'stratus: false is the column without a deck');
  assert.deepEqual(run(off, unstable, 298, 1), run(radiation, unstable, 298, 0));
  assert.deepEqual(run(radiation, stable, 275, 1), run(off, stable, 275, 1), 'a deck with no cover leaves the column bit-identical');
});

test('the deck is two independent columns: its shortwave is the cover-weighted mean of the overcast and clear columns, the overcast one is the column with the deck water condensed in its layer, and the column closes', () => {
  const { radiation, off, noon, qc, withStability, run } = deckColumns();
  const part = run(radiation, withStability(18), 298, 1), full = run(radiation, withStability(30), 298, 1), none = run(radiation, withStability(30), 298, 0);
  const f = part.stratusFraction;
  assert.ok(f > 0.2 && f < 0.8 && full.stratusFraction === 1 && none.stratusFraction === 0, `covers ${f}, ${full.stratusFraction}, ${none.stratusFraction}`);
  assert.ok(part.stratus > 0 && part.stratus === full.stratus);
  for (const name of ['reflectedSolar', 'absorbedSolar', 'surfaceShortwave', 'surfaceDirect', 'cloudReflectance']) {
    const mean = f * full[name] + (1 - f) * none[name];
    assert.ok(Math.abs(part[name] - mean) < 1e-9, `${name} ${part[name]} against the mean ${mean}`);
  }
  console.log(`${f.toFixed(3)} cover of a ${part.stratus.toFixed(4)} kg/m² deck at noon reflects ${part.reflectedSolar.toFixed(2)} W/m², the weighted mean of ${full.reflectedSolar.toFixed(2)} overcast and ${none.reflectedSolar.toFixed(2)} clear; closure ${part.closure.toExponential(1)}`);
  assert.ok(part.closure < EPS && full.closure < EPS, `closure ${part.closure}, ${full.closure}`);
  const condensed = Float64Array.from(qc), layer = radiation.stratusLayer;
  condensed[layer * C + noon] = full.stratus / (P0 * core.diagnostics.dSigma[layer] / core.diagnostics.g);
  const overcast = run(off, withStability(30), 298, 1, 1200, condensed);
  for (const name of ['reflectedSolar', 'absorbedSolar', 'surfaceShortwave', 'surfaceDirect', 'outgoingLongwave', 'surfaceFlux']) {
    assert.ok(Math.abs(full[name] - overcast[name]) <= 1e-9 * Math.max(1, Math.abs(overcast[name])), `${name} ${full[name]} against the condensed deck's ${overcast[name]}`);
  }
  full.layers.forEach((f, k) => assert.ok(Math.abs(f - overcast.layers[k]) <= 1e-9 * Math.max(1, Math.abs(f)), `layer ${k}: ${f} against ${overcast.layers[k]}`));
});

test('the deck\'s own water absorbs in the deck layer in the same two-column blend: half the cover absorbs half as much', () => {
  const { radiation, withStability, run } = deckColumns();
  const scattering = createRadiation(mesh, core, { stratus: true, mixedLayerDeck: false, cloudSolarAbsorption: 0 });
  scattering.setTime(0);
  const layer = radiation.stratusLayer, stable = withStability(30);
  const absorbed = (openSea) => {
    const lit = run(radiation, stable, 298, openSea), dark = run(scattering, stable, 298, openSea);
    const heating = lit.layers.map((f, k) => f - dark.layers[k]);
    heating.forEach((h, k) => { if (k !== layer) assert.equal(h, 0, `layer ${k}`); });
    return { ...lit, deckLayer: heating[layer] };
  };
  const full = absorbed(1), half = absorbed(0.5);
  console.log(`a ${(1000 * full.stratus).toFixed(1)} g/m² deck at noon absorbs ${full.deckLayer.toFixed(2)} W/m² in its layer under full cover, ${half.deckLayer.toFixed(2)} under half`);
  assert.ok(full.stratusFraction === 1 && half.stratusFraction === 0.5 && half.stratus === full.stratus);
  assert.ok(full.deckLayer > 0 && Math.abs(half.deckLayer - 0.5 * full.deckLayer) < 1e-12 * full.deckLayer, `half ${half.deckLayer} against full ${full.deckLayer}`);
  assert.ok(Math.abs(full.cloudSolar - full.deckLayer) < 1e-12 * full.deckLayer && full.closure < EPS && half.closure < EPS);
});

test('the deck is as thick as the boundary layer above its condensation level, with the adiabatic water of that thickness scaled down', () => {
  const { g, kappa } = core.diagnostics, lapse = adiabaticWaterLapse(290, 95000, CP_DRY, R_DRY, g);
  assert.ok(Math.abs(lapse - 2.44e-6) < 0.02e-6, `Γ_l ${lapse} kg/m³ per m`);
  const { radiation, off, q, bottom, withStability, run, airTemperature } = deckColumns();
  q[bottom] = 0.8 * saturationHumidity(airTemperature, P0 * core.sigmaMid[K - 1]);
  const lcl = liftingCondensationLevel(airTemperature, q[bottom], P0 * core.sigmaMid[K - 1], kappa);
  const base = CP_DRY * (airTemperature - lcl.temperature) / g;
  const stable = withStability(20);
  const shallow = run(radiation, stable, 298, 1, 500), deep = run(radiation, stable, 298, 1, 1200), under = run(radiation, stable, 298, 1, base - 20);
  const adiabatic = 0.15 * 0.5 * adiabaticWaterLapse(lcl.temperature, lcl.pressure, CP_DRY, R_DRY, g) * (1200 - base) ** 2;
  console.log(`80 % humid air condenses ${base.toFixed(0)} m up: a 500 m boundary layer holds ${(1000 * shallow.stratus).toFixed(2)} g/m², a 1200 m one ${(1000 * deep.stratus).toFixed(1)} g/m² (a sub-adiabatic 0.15 of ${(1000 * adiabatic / 0.15).toFixed(0)}), one topped below the base none`);
  assert.ok(base > 200 && base < 500, `condensation level ${base} m`);
  assert.ok(shallow.stratus > 0 && deep.stratus >= 3 * shallow.stratus, `deck ${shallow.stratus} at 500 m, ${deep.stratus} at 1200 m`);
  assert.ok(Math.abs(deep.stratus - Math.min(0.15, adiabatic)) < 1e-12, `deck ${deep.stratus} against ${adiabatic}`);
  assert.equal(under.stratus, 0);
  assert.equal(under.stratusFraction, 0);
  assert.deepEqual(under, run(off, stable, 298, 1, base - 20), 'a boundary layer below its condensation level leaves the column bit-identical');
  assert.equal(run(createRadiation(mesh, core, { stratusScale: 10, mixedLayerDeck: false }), stable, 298, 1, 1200).stratus, 0.15, 'the deck water stops at stratusWaterMax');
});

test('the deck follows the estimated inversion strength, which tells a capping inversion from a free troposphere warmed with the sea', () => {
  const { radiation, q, unstable, bottom, top, withStability, run, inversion } = deckColumns();
  const { g, kappa, exnerLayer, geopotential } = core.diagnostics;
  const humid = (theta) => { q[bottom] = 0.8 * saturationHumidity(theta[bottom] * exnerLayer[bottom], P0 * core.sigmaMid[K - 1]); return theta; };
  const stable = humid(withStability(18));
  const expected = inversion(stable);
  const lowerT = stable[bottom] * exnerLayer[bottom], lcl = liftingCondensationLevel(lowerT, q[bottom], P0 * core.sigmaMid[K - 1], kappa);
  const depth = (geopotential[top] - geopotential[bottom]) / g - CP_DRY * (lowerT - lcl.temperature) / g;
  const eis = inversionStrength(stable[top] - stable[bottom], lowerT, stable[top] * exnerLayer[top], depth, CP_DRY, R_DRY, g);
  const deck = run(radiation, stable, 298, 1);
  assert.ok(Math.abs(eis - expected) < 1e-9 && Math.abs(radiation.budget.stabilityIndex - expected) < 1e-9, `EIS ${eis} (column ${radiation.budget.stabilityIndex}) against ${expected}`);
  assert.ok(eis > 3 && eis < 7, `EIS ${eis} K at LTS 18 K`);
  assert.ok(Math.abs(deck.stratusFraction - (0.19 + 0.08 * (expected - 1))) < 1e-9, `cover ${deck.stratusFraction}`);
  const column = (shift, surfaceT) => {
    const theta = Float64Array.from(unstable, (t) => t + shift);
    theta[top] = theta[bottom] + 16;
    humid(theta);
    return { theta, inversion: inversion(theta), ...run(radiation, theta, surfaceT, 1) };
  };
  const warm = column(4, 302), cool = column(-6, 292);
  console.log(`LTS 18 K over a 298 K sea at 80 % humidity: EIS ${eis.toFixed(2)} K, cover ${deck.stratusFraction.toFixed(3)}; LTS 16 K over a 302 K sea: EIS ${warm.inversion.toFixed(2)} K, cover ${warm.stratusFraction.toFixed(3)}; over a 292 K sea: EIS ${cool.inversion.toFixed(2)} K, cover ${cool.stratusFraction.toFixed(3)}`);
  assert.ok(Math.abs(warm.theta[top] - warm.theta[bottom] - (cool.theta[top] - cool.theta[bottom])) < 1e-9);
  assert.ok(cool.inversion > warm.inversion + 2, `EIS ${cool.inversion} over the cool sea against ${warm.inversion} over the warm one`);
  assert.ok(cool.stratusFraction > warm.stratusFraction && warm.stratusFraction > 0, `cover ${cool.stratusFraction} against ${warm.stratusFraction}`);
  for (const c of [warm, cool]) assert.ok(Math.abs(c.stratusFraction - (0.19 + 0.08 * (c.inversion - 1))) < 1e-9);
});

test('with stratusIndex: \'ectei\' the deck follows the entrainment index: a dry layer above thins it, a layer as humid as the lowest leaves it as EIS has it, and the default is EIS', () => {
  const { radiation, q, bottom, top, withStability, run, inversion } = deckColumns();
  const entraining = createRadiation(mesh, core, { stratusIndex: 'ectei', mixedLayerDeck: false }), explicit = createRadiation(mesh, core, { stratusIndex: 'eis', mixedLayerDeck: false });
  entraining.setTime(0); explicit.setTime(0);
  const { exnerLayer } = core.diagnostics;
  const theta = withStability(24);
  q[bottom] = 0.8 * saturationHumidity(theta[bottom] * exnerLayer[bottom], P0 * core.sigmaMid[K - 1]);
  const covers = (upperQ) => {
    q[top] = upperQ;
    const eis = run(radiation, theta, 298, 1), eisIndex = radiation.budget.stabilityIndex;
    const ectei = run(entraining, theta, 298, 1), ecteiIndex = entraining.budget.stabilityIndex;
    const same = run(explicit, theta, 298, 1), sameIndex = explicit.budget.stabilityIndex;
    return { eis, ectei, same, eisIndex, ecteiIndex, sameIndex, expected: inversion(theta) };
  };
  const dry = covers(0.2 * q[bottom]), moist = covers(q[bottom]);
  const gap = 0.23 * LATENT_HEAT / CP_DRY * 0.8 * q[bottom];
  console.log(`LTS 24 K over a 298 K sea, 80 % humid below: EIS ${dry.eisIndex.toFixed(2)} K covers ${dry.eis.stratusFraction.toFixed(3)}; under a 700 hPa layer at 0.2 of the lowest humidity ECTEI ${dry.ecteiIndex.toFixed(2)} K covers ${dry.ectei.stratusFraction.toFixed(3)}; as humid as the lowest ECTEI ${moist.ecteiIndex.toFixed(2)} K against EIS ${moist.eisIndex.toFixed(2)} K`);
  assert.ok(Math.abs(dry.eisIndex - dry.expected) < 1e-9, `EIS ${dry.eisIndex} against ${dry.expected}`);
  assert.ok(Math.abs(dry.ecteiIndex - (dry.expected - gap)) < 1e-9 && Math.abs(entrainmentIndex(dry.expected, q[bottom], 0.2 * q[bottom], CP_DRY) - (dry.expected - gap)) < 1e-9, `ECTEI ${dry.ecteiIndex} against ${dry.expected - gap}`);
  assert.ok(Math.abs(dry.ectei.stratusFraction - (0.19 + 0.08 * (dry.ecteiIndex - 1))) < 1e-9, `cover ${dry.ectei.stratusFraction}`);
  assert.ok(dry.ectei.stratusFraction > 0 && dry.ectei.stratusFraction < dry.eis.stratusFraction - 0.2, `cover ${dry.ectei.stratusFraction} under ECTEI against ${dry.eis.stratusFraction} under EIS`);
  assert.ok(Math.abs(moist.ecteiIndex - moist.eisIndex) < 1e-9 && Math.abs(moist.ectei.stratusFraction - moist.eis.stratusFraction) < 1e-9, `ECTEI ${moist.ecteiIndex} against EIS ${moist.eisIndex}`);
  for (const c of [dry, moist]) {
    assert.deepEqual(c.same, c.eis, 'stratusIndex: \'eis\' is the default');
    assert.equal(c.sameIndex, c.eisIndex);
  }
  assert.throws(() => createRadiation(mesh, core, { stratusIndex: 'lts' }));
});

const REDIAGNOSED = { prognosticHeight: false, gateMemory: 0 };

function mixedLayerColumn(sinking = 0.4) {
  const pi = new Float64Array(C).fill(P0), theta = new Float64Array(K * C), q = new Float64Array(K * C), qc = new Float64Array(K * C);
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) {
    const sigma = core.sigmaMid[k], mixed = sigma > 0.9;
    theta[k * C + i] = mixed ? 289 : 298 + 60 * (0.9 - sigma);
    q[k * C + i] = mixed ? 9e-3 : 1.5e-3 * Math.min(1, sigma / 0.8) ** 3;
  }
  core.diagnose(pi, theta, q, qc);
  const { g, geopotential, piSigmaDot } = core.diagnostics;
  for (let k = 1; k < K; k++) piSigmaDot.fill(sinking * (1 - core.levels[k]), k * C, (k + 1) * C);
  const depth = Float64Array.from({ length: C }, (_, i) => geopotential[(K - 3) * C + i] / g + 100);
  return { pi, theta, q, qc, depth, mixedDepth: (i) => depth[i] - geopotential[(K - 1) * C + i] / g };
}

function mixedLayerRun(r, noon, column, openSea = 1, dt = 900, air = column.theta) {
  const { q, qc, mixedDepth } = column, bottom = (K - 1) * C + noon;
  const flux = r.column(noon, P0, air, 292, 5, r.opticalDepth(mesh.latCell[noon]), r.insolation(noon), q[bottom], q, qc, 0.07, 0.06, 1, 1.5e-3, openSea, mixedDepth(noon), dt);
  let layers = 0, scale = 0;
  for (let k = 0; k < K; k++) { layers += r.layerFlux[k]; scale += Math.abs(r.layerFlux[k]); }
  const latent = LATENT_HEAT * r.budget.evaporation;
  return { ...r.budget, closure: Math.abs(layers + flux + latent - (r.budget.absorbedSolar - r.budget.outgoingLongwave)) / (scale + Math.abs(flux) + latent) };
}

function brightest(r) {
  let noon = 0;
  for (let i = 0; i < C; i++) if (r.insolation(i) > r.insolation(noon)) noon = i;
  return noon;
}

test('the mixed-layer deck on a stable column over a warm sea carries the water path and cover of the mixed-layer model, advanced one step, in the same two-column blend', () => {
  const shadow = createRadiation(mesh, core, { mixedLayerDeck: true, subsidenceMemory: 1e-9, ...REDIAGNOSED }), off = createRadiation(mesh, core, { stratus: false });
  shadow.setTime(0); off.setTime(0);
  const noon = brightest(shadow), column = mixedLayerColumn(), { pi, theta, q, qc, depth } = column;
  const run = (r, openSea, dt = 900, air = theta) => mixedLayerRun(r, noon, column, openSea, dt, air);
  const deck = run(shadow, 1), clear = run(off, 1), start = run(shadow, 1, 0), half = run(shadow, 0.5);
  console.log(`289 K mixed layer under a 298 K free troposphere, 100 m above its third layer, over a 292 K sea: the mixed-layer model holds ${(1000 * deck.mlmWater).toFixed(1)} g/m² (${(1000 * start.mlmWater).toFixed(1)} before its step) on ${deck.mlmCover} of the cell, entraining ${(1000 * deck.mlmEntrainment).toFixed(2)} mm/s; noon absorbed solar ${deck.absorbedSolar.toFixed(0)} against ${clear.absorbedSolar.toFixed(0)} W/m²`);
  assert.ok(deck.mlmWater > 0.03 && deck.mlmWater < 0.15 && deck.mlmCover === 1, `water ${deck.mlmWater}, cover ${deck.mlmCover}`);
  assert.ok(deck.mlmEntrainment > 1e-3 && deck.mlmEntrainment < 1e-2, `entrainment ${deck.mlmEntrainment}`);
  assert.ok(deck.stratus === deck.mlmWater && deck.stratusFraction === deck.mlmCover, 'the deck is the mixed-layer model\'s');
  assert.ok(deck.mlmWater !== start.mlmWater && Math.abs(deck.mlmWater - start.mlmWater) < 0.1 * start.mlmWater, 'one step adjusts the water path');
  assert.ok(half.stratusFraction === 0.5 * deck.stratusFraction && half.stratus === deck.stratus, 'half the cell under ice halves the cover');
  assert.ok(deck.absorbedSolar < clear.absorbedSolar - 300, `absorbed solar ${deck.absorbedSolar} against ${clear.absorbedSolar}`);
  assert.ok(deck.closure < EPS && Number.isFinite(deck.stabilityIndex));
  const flat = Float64Array.from(theta);
  for (let k = 0; k < K; k++) if (core.sigmaMid[k] > 0.8) flat[k * C + noon] = 289;
  core.diagnose(pi, flat, q, qc);
  const uncapped = run(shadow, 1, 900, flat);
  assert.ok(uncapped.stratus === 0 && uncapped.stratusFraction === 0 && uncapped.mlmCover === 0, 'a layer without an inversion carries no deck');
  core.diagnose(pi, theta, q, qc);
  const out = [new Float64Array(C), new Float64Array(K * C), new Float64Array(K * E), new Float64Array(C), new Float64Array(K * C)];
  shadow.apply([pi, theta, new Float64Array(K * E), new Float64Array(C).fill(292), q, qc], out, new Float64Array(C).fill(5), null, 0, C, null, null, null, new Float64Array(C).fill(1), depth, 900);
  for (let i = 0; i < C; i++) {
    assert.ok(shadow.mlmCover[i] === 1 && shadow.mlmWater[i] > 0.03 && shadow.mlmEntrainment[i] > 1e-3, `cell ${i}`);
    assert.ok(shadow.stratus[i] === Math.min(0.15, shadow.mlmWater[i]) && shadow.stratusFraction[i] === shadow.mlmCover[i]);
  }
  core.diagnostics.piSigmaDot.fill(0);
});

test('the mixed-layer deck needs subsidence and a capping inversion: a column under ascent or under a 1 K jump has none, and falls back to no deck', () => {
  const shadow = createRadiation(mesh, core, { mixedLayerDeck: true, subsidenceMemory: 1e-9, ...REDIAGNOSED });
  shadow.setTime(0);
  const noon = brightest(shadow), empty = { mlmCover: 0, mlmWater: 0, mlmEntrainment: 0, stratus: 0, stratusFraction: 0 };
  const pick = (b) => ({ mlmCover: b.mlmCover, mlmWater: b.mlmWater, mlmEntrainment: b.mlmEntrainment, stratus: b.stratus, stratusFraction: b.stratusFraction });
  const sinking = mixedLayerColumn(0.4);
  assert.ok(mixedLayerRun(shadow, noon, sinking).mlmCover === 1, 'a sinking, capped column carries the deck');
  const rising = mixedLayerColumn(-0.4);
  assert.deepEqual(pick(mixedLayerRun(shadow, noon, rising)), empty, 'a rising column carries none');
  const column = mixedLayerColumn(0.4), above = (K - 4) * C + noon;
  const mlm = createMixedLayer({ cp: CP_DRY, R: R_DRY, g: core.diagnostics.g, latentHeat: LATENT_HEAT, referencePressure: P0, cloudLevels: 8 });
  const jump = (thetaAbove) => mlm.diagnose({ h: column.depth[noon], thetaL: 289, qt: 9e-3 }, { surfacePressure: P0, sensibleHeat: 0, evaporation: 0, thetaLAbove: thetaAbove, qtAbove: column.q[above], subsidence: () => 0, radiation: dycomsLongwave() }).virtualJump;
  const slope = jump(301) - jump(300), withJump = (target) => { const theta = Float64Array.from(column.theta); theta[above] = 300 + (target - jump(300)) / slope; return theta; };
  const weak = withJump(1), strong = withJump(3);
  console.log(`first layer above h at ${weak[above].toFixed(2)} K gives Δθ_v ${jump(weak[above]).toFixed(2)} K, at ${strong[above].toFixed(2)} K ${jump(strong[above]).toFixed(2)} K`);
  core.diagnose(column.pi, weak, column.q, column.qc);
  assert.deepEqual(pick(mixedLayerRun(shadow, noon, column, 1, 900, weak)), empty, 'a 1 K jump carries none');
  core.diagnose(column.pi, strong, column.q, column.qc);
  assert.ok(mixedLayerRun(shadow, noon, column, 1, 900, strong).mlmCover > 0, 'a 3 K jump carries the deck');
  core.diagnostics.piSigmaDot.fill(0);
});

test('the regime test reads the subsidence averaged over subsidenceMemory: a column that turns from rising to sinking gains its deck once the running mean passes the threshold', () => {
  const instant = createRadiation(mesh, core, { mixedLayerDeck: true, subsidenceMemory: 1e-9, ...REDIAGNOSED }), memory = createRadiation(mesh, core, { mixedLayerDeck: true, ...REDIAGNOSED });
  instant.setTime(0); memory.setTime(0);
  const { subsidenceMemory, stratusSubsidence } = memory.deckGates;
  assert.equal(subsidenceMemory, 2 * DAY);
  const noon = brightest(memory), column = mixedLayerColumn(0.4), dt = 3600, keep = Math.exp(-dt / subsidenceMemory);
  mixedLayerRun(instant, noon, column, 1, dt);
  const sinking = instant.mlmSubsidence[noon];
  memory.mlmSubsidence[noon] = -sinking;
  let first = -1, expected = -sinking;
  for (let n = 1; n <= 48; n++) {
    const deck = mixedLayerRun(memory, noon, column, 1, dt);
    expected = expected * keep + sinking * (1 - keep);
    assert.ok(Math.abs(memory.mlmSubsidence[noon] - expected) < 1e-12 * Math.abs(sinking), `step ${n}: ${memory.mlmSubsidence[noon]} against ${expected}`);
    assert.equal(deck.mlmCover > 0, expected <= -stratusSubsidence, `step ${n}`);
    if (first < 0 && deck.mlmCover > 0) first = n;
  }
  console.log(`a column that rose at ${(-1000 * sinking).toFixed(2)} mm/s and now sinks as fast: its ${subsidenceMemory / DAY}-day mean passes ${(-1000 * stratusSubsidence).toFixed(2)} mm/s and the deck appears after ${first} hourly steps`);
  assert.ok(first > 1 && first < 48);
  core.diagnostics.piSigmaDot.fill(0);
});

test('the deck reads πσ̇ averaged over the cell and its neighbours, twice over, at the two interfaces bracketing h: its subsidence is that of the field smoothCells smooths twice, and grid-scale noise in πσ̇ barely reaches it', () => {
  const rnd = random(11), field = Float64Array.from({ length: C }, () => rnd() - 0.5);
  for (let passes = 0; passes <= 2; passes++) {
    let smoothed = Float64Array.from(field);
    for (let p = 0; p < passes; p++) smoothed = smoothCells(mesh, smoothed, new Float64Array(C));
    for (let i = 0; i < C; i++) assert.equal(ringMean(mesh, field, 0, i, passes), smoothed[i], `cell ${i}, ${passes} passes`);
  }
  assert.throws(() => createRadiation(mesh, core, { subsidenceSmoothing: 3 }));
  const read = createRadiation(mesh, core, { subsidenceMemory: 1e-9, ...REDIAGNOSED }), unsmoothed = createRadiation(mesh, core, { subsidenceMemory: 1e-9, subsidenceSmoothing: 0, ...REDIAGNOSED });
  read.setTime(0); unsmoothed.setTime(0);
  const column = mixedLayerColumn(0.4), { piSigmaDot } = core.diagnostics;
  for (let k = 1; k < K; k++) for (let i = 0; i < C; i++) piSigmaDot[k * C + i] += 0.4 * (rnd() - 0.5);
  const raw = Float64Array.from(piSigmaDot), presmoothed = Float64Array.from(raw);
  for (let k = 1; k < K; k++) presmoothed.set(smoothCells(mesh, smoothCells(mesh, raw.subarray(k * C, (k + 1) * C), new Float64Array(C)), new Float64Array(C)), k * C);
  const w = (r) => Float64Array.from({ length: C }, (_, i) => { mixedLayerRun(r, i, column, 1, 900); return r.mlmSubsidence[i]; });
  const smoothedRead = w(read), noisy = w(unsmoothed);
  piSigmaDot.set(presmoothed);
  const readOfSmoothed = w(unsmoothed);
  const spread = (x) => { const m = x.reduce((a, b) => a + b, 0) / C; return Math.sqrt(x.reduce((a, b) => a + (b - m) ** 2, 0) / C); };
  let worst = 0;
  for (let i = 0; i < C; i++) worst = Math.max(worst, Math.abs(smoothedRead[i] - readOfSmoothed[i]) / Math.abs(readOfSmoothed[i]));
  console.log(`under ±0.2 Pa/s of cell-to-cell noise on a 0.04 Pa/s sink the deck's subsidence spreads over cells by ${(1000 * spread(noisy)).toFixed(2)} mm/s unsmoothed and ${(1000 * spread(smoothedRead)).toFixed(2)} mm/s as it reads it; against the presmoothed field it agrees to ${worst.toExponential(1)}`);
  assert.ok(worst < 1e-12, `the deck's reading against the presmoothed field: ${worst}`);
  assert.ok(spread(smoothedRead) < 0.3 * spread(noisy), `spread ${spread(smoothedRead)} against ${spread(noisy)}`);
  piSigmaDot.fill(0);
});

function boundaryLayerTop(column, i) {
  const { g, cp, geopotential, exnerLower, exnerLayer } = core.diagnostics, bottom = (K - 1) * C + i, thetaV = core.arrays.thetaV;
  const surface = geopotential[bottom] - cp * thetaV[bottom] * (exnerLower[bottom] - exnerLayer[bottom]);
  return { height: column.mixedDepth(i) + (geopotential[bottom] - surface) / g, offset: surface / g };
}

test('the deck carries its inversion height: from the boundary-layer top it deepens step after step by its own dh/dt, while the re-diagnosed deck starts from that top each step; mlmTop hands the height to the boundary layer', () => {
  const carried = createRadiation(mesh, core, { subsidenceMemory: 1e-9 }), rediagnosed = createRadiation(mesh, core, { subsidenceMemory: 1e-9, ...REDIAGNOSED });
  carried.setTime(0); rediagnosed.setTime(0);
  const noon = brightest(carried), column = mixedLayerColumn(0.4), top = boundaryLayerTop(column, noon), dt = 900;
  let previous = top.height, fixed = null;
  for (let n = 1; n <= 16; n++) {
    const deck = mixedLayerRun(carried, noon, column, 1, dt), again = mixedLayerRun(rediagnosed, noon, column, 1, dt);
    const h = carried.mlmHeight[noon];
    assert.ok(deck.mlmCover === 1 && again.mlmCover === 1, `step ${n}`);
    assert.ok(h > previous && h - previous < dt * 0.02, `step ${n}: ${previous} → ${h}`);
    assert.ok(Math.abs(deck.mlmTop - (h + top.offset)) < 1e-9 * h && again.mlmTop === 0, `step ${n}: mlmTop ${deck.mlmTop}`);
    fixed ??= rediagnosed.mlmHeight[noon];
    assert.equal(rediagnosed.mlmHeight[noon], fixed, `step ${n}: the re-diagnosed deck starts from the boundary-layer top`);
    previous = h;
  }
  console.log(`over 4 h of a steady sinking column the carried inversion climbs from the boundary-layer top at ${top.height.toFixed(1)} m to ${previous.toFixed(1)} m; the re-diagnosed deck ends each step at ${fixed.toFixed(1)} m`);
  assert.ok(previous > fixed + 10, `carried ${previous}, re-diagnosed ${fixed}`);
  carried.mlmHeight[noon] = 0.5 * top.height;
  mixedLayerRun(carried, noon, column, 1, dt);
  assert.ok(carried.mlmHeight[noon] >= top.height, `a carried height below the boundary-layer top starts from the top: ${carried.mlmHeight[noon]}`);
  core.diagnostics.piSigmaDot.fill(0);
});

test('the gates switch the deck through their running mean: a standing deck outlives failing gates by ln 2 × gateMemory and a new one waits about as long, deepening meanwhile no further than the inversion ceiling, and without its deck the carried height relaxes toward the boundary-layer top over heightMemory', () => {
  const r = createRadiation(mesh, core, { subsidenceMemory: 1e-9 });
  r.setTime(0);
  const noon = brightest(r), dt = 3600, fresh = 1 - Math.exp(-dt / DAY), relaxed = Math.exp(-dt / DAY);
  let gate = 0.5, column = mixedLayerColumn(0.4);
  const top = boundaryLayerTop(column, noon);
  for (let n = 1; n <= 48; n++) {
    const deck = mixedLayerRun(r, noon, column, 1, dt);
    gate += (1 - gate) * fresh;
    assert.ok(Math.abs(r.mlmGate[noon] - gate) < 1e-12 && deck.mlmCover > 0 && deck.mlmTop > 0, `sinking step ${n}: gate ${r.mlmGate[noon]} against ${gate}`);
  }
  const standing = r.mlmHeight[noon];
  column = mixedLayerColumn(-0.4);
  let lastOn = 0, switchedOff = 0;
  for (let n = 1; n <= 48; n++) {
    const before = r.mlmHeight[noon], deck = mixedLayerRun(r, noon, column, 1, dt);
    gate -= gate * fresh;
    assert.ok(Math.abs(r.mlmGate[noon] - gate) < 1e-12, `rising step ${n}`);
    assert.equal(deck.mlmTop > 0, gate > 0.5, `rising step ${n}: gate ${gate}`);
    if (gate > 0.5) { lastOn = n; switchedOff = r.mlmHeight[noon]; continue; }
    const expected = top.height + (before - top.height) * relaxed;
    assert.ok(Math.abs(r.mlmHeight[noon] - expected) < 1e-9 * expected && deck.stratusFraction === 0, `rising step ${n}: height ${r.mlmHeight[noon]} against ${expected}`);
  }
  const fallen = r.mlmHeight[noon];
  column = mixedLayerColumn(0.4);
  let formed = 0;
  for (let n = 1; n <= 48 && !formed; n++) {
    const deck = mixedLayerRun(r, noon, column, 1, dt);
    gate += (1 - gate) * fresh;
    assert.equal(deck.mlmTop > 0, gate > 0.5, `sinking again, step ${n}`);
    if (deck.mlmTop > 0) formed = n;
  }
  console.log(`hourly steps: after 48 h of sinking the carried inversion stands at ${standing.toFixed(0)} m over a ${top.height.toFixed(0)} m boundary layer; under ascent the deck stays ${lastOn} h, ending at ${switchedOff.toFixed(0)} m, and the height relaxes to ${fallen.toFixed(0)} m by 48 h; sinking again, the deck returns after ${formed} h`);
  assert.ok(lastOn >= 12 && lastOn <= 24 && formed >= 12 && formed <= 36, `stayed ${lastOn} h, returned after ${formed} h`);
  const settled = top.height + (switchedOff - top.height) * Math.exp(-(48 - lastOn) * dt / DAY);
  const freeTroposphere = core.diagnostics.geopotential[(K - 4) * C + noon] / core.diagnostics.g - top.offset;
  assert.ok(switchedOff > standing && switchedOff <= freeTroposphere - 1 + 1e-9, `under ascent the deck deepens to ${switchedOff} m, no further than 1 m under the midpoint of the first free-tropospheric layer at ${freeTroposphere} m`);
  assert.ok(standing > top.height + 20 && Math.abs(fallen - settled) < 1e-9 * settled, `standing ${standing}, fallen ${fallen} against ${settled}`);
  core.diagnostics.piSigmaDot.fill(0);
});

function modelDigest(radiation) {
  const model = createModel(new Grid(4), { ocean: { eddyDiffusivity: 0 }, divergenceDamping: 0, ...(radiation ? { radiation } : {}) });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  for (let n = 0; n < 12; n++) model.step(900);
  const hash = createHash('sha256');
  for (const a of [...model.state, model.radiation.stratus, model.radiation.stratusFraction, model.radiation.surfaceFlux, model.radiation.outgoing]) hash.update(new Uint8Array(a.buffer, a.byteOffset, a.byteLength));
  return { digest: hash.digest('hex').slice(0, 32), radiation: model.radiation };
}

test('with mixedLayerDeck: false and the purely scattering clouds of cloudSolarAbsorption: 0, cloudScattering: 55 the model is bit-identical to the engine before the mixed-layer deck; by default the deck follows the mixed-layer model', () => {
  const before = '2c1fe52a35a17a762501cb6b3e31ab46';
  assert.equal(modelDigest({ mixedLayerDeck: false, cloudSolarAbsorption: 0, cloudScattering: 55, ...OVERCAST }).digest, before);
  assert.notEqual(modelDigest({ mixedLayerDeck: false }).digest, before, 'by default cloud water absorbs sunlight');
  const fresh = modelDigest();
  assert.equal(fresh.digest, modelDigest({ mixedLayerDeck: true }).digest);
  assert.notEqual(fresh.digest, before);
  assert.ok(fresh.radiation.stratusFraction.every((f) => f === 0) && fresh.radiation.mlmCover.every((f) => f === 0), 'twelve steps from rest build no 2 K inversion, and no empirical deck stands in');
  const shadow = modelDigest({ mixedLayerDeck: true, stratusSubsidence: 0, minimumInversion: 0 });
  const { mlmCover, mlmWater, mlmEntrainment, stratusFraction } = shadow.radiation;
  let covered = 0;
  for (let i = 0; i < mlmCover.length; i++) {
    assert.ok(mlmCover[i] >= 0 && mlmCover[i] <= 1 && mlmWater[i] >= 0 && mlmEntrainment[i] >= 0 && Number.isFinite(mlmWater[i] + mlmEntrainment[i]));
    assert.ok(stratusFraction[i] <= mlmCover[i]);
    if (stratusFraction[i] > 0) covered++;
  }
  assert.ok(covered > 0, 'some cells carry a mixed-layer deck');
});

test('the mixed layer feels the sunlight the column absorbs in the deck\'s layer: with the purely scattering clouds of cloudSolarAbsorption: 0, cloudScattering: 55 it feels none and the engine is bit-identical to the deck before it absorbed sunlight, with stratusSolar: false it feels none while the column absorbs', () => {
  const forced = { stratusSubsidence: 0, minimumInversion: 0, subsidenceSmoothing: 0, subsidenceMemory: 10 * DAY }, scatteringOnly = { cloudSolarAbsorption: 0, cloudScattering: 55, ...OVERCAST };
  assert.equal(modelDigest({ stratusSolar: false, ...scatteringOnly }).digest, '359b2a50159e3dfa0236098bac19547d');
  assert.equal(modelDigest({ ...forced, ...REDIAGNOSED, stratusSolar: false, ...scatteringOnly }).digest, '9587f5400c7a4b789fa02c0c4e41015c');
  assert.equal(modelDigest({ ...forced, ...REDIAGNOSED, ...scatteringOnly }).digest, '9587f5400c7a4b789fa02c0c4e41015c');
  assert.notEqual(modelDigest({ ...forced, ...REDIAGNOSED }).digest, '9587f5400c7a4b789fa02c0c4e41015c');
  assert.notEqual(modelDigest({ ...forced, ...scatteringOnly }).digest, '9587f5400c7a4b789fa02c0c4e41015c', 'the carried height and the gate\'s memory change the deck');
  const shadow = createRadiation(mesh, core, { subsidenceMemory: 1e-9, ...REDIAGNOSED }), dark = createRadiation(mesh, core, { subsidenceMemory: 1e-9, stratusSolar: false, ...REDIAGNOSED });
  const scattering = createRadiation(mesh, core, { subsidenceMemory: 1e-9, cloudSolarAbsorption: 0, ...REDIAGNOSED });
  shadow.setTime(0); dark.setTime(0); scattering.setTime(0);
  const noon = brightest(shadow), column = mixedLayerColumn(), lit = mixedLayerRun(shadow, noon, column), unlit = mixedLayerRun(dark, noon, column), half = mixedLayerRun(shadow, noon, column, 0.5);
  console.log(`at noon the mixed layer of the stable column holds ${(1000 * lit.mlmWater).toFixed(2)} g/m² after its step with the ${lit.mlmSolar.toFixed(2)} W/m² its deck's layer absorbs, ${(1000 * unlit.mlmWater).toFixed(2)} without, entraining ${(1000 * lit.mlmEntrainment).toFixed(3)} and ${(1000 * unlit.mlmEntrainment).toFixed(3)} mm/s`);
  assert.ok(lit.mlmCover === 1 && lit.mlmSolar > 0 && lit.mlmSolar === lit.cloudSolar, `the mixed layer feels ${lit.mlmSolar} W/m², the column absorbs ${lit.cloudSolar}`);
  assert.ok(half.mlmSolar === lit.mlmSolar && Math.abs(half.cloudSolar - 0.5 * lit.cloudSolar) < 1e-12 * lit.cloudSolar, 'the mixed layer feels the power per unit deck area');
  assert.ok(unlit.mlmSolar === 0 && unlit.cloudSolar > 0, 'with stratusSolar: false the column still absorbs');
  assert.deepEqual([scattering.budget.mlmSolar, mixedLayerRun(scattering, noon, column).mlmWater], [0, unlit.mlmWater]);
  assert.ok(lit.mlmWater < unlit.mlmWater && lit.mlmWater > 0.9 * unlit.mlmWater, `water ${lit.mlmWater} against ${unlit.mlmWater}`);
  assert.ok(lit.closure < EPS && unlit.closure < EPS);
  core.diagnostics.piSigmaDot.fill(0);
});

function evaluate(radiation, i, pi, theta, surfaceT, q, qc) {
  const flux = radiation.column(i, pi[i], theta, surfaceT[i], 5, radiation.opticalDepth(mesh.latCell[i]), radiation.insolation(i), q[(K - 1) * C + i], q, qc, 0.07);
  let layers = 0, scale = 0;
  for (let k = 0; k < K; k++) { layers += radiation.layerFlux[k]; scale += Math.abs(radiation.layerFlux[k]); }
  const latent = LATENT_HEAT * radiation.budget.evaporation;
  const residual = layers + flux + latent - (radiation.budget.absorbedSolar - radiation.budget.outgoingLongwave);
  return { reflected: radiation.budget.reflectedSolar, olr: radiation.budget.outgoingLongwave, absorbed: radiation.budget.absorbedSolar, surfaceShortwave: radiation.budget.surfaceShortwave, layers: Array.from(radiation.layerFlux), closure: Math.abs(residual) / (scale + Math.abs(flux) + latent) };
}
