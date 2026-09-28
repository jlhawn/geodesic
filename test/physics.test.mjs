import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore, P0, CP_DRY, R_DRY } from '../js/dynamics/sigmaCore.module.js';
import { createRadiation, sunDirection, AXIAL_TILT, DAY, YEAR, waterVaporAbsorptivity, adiabaticWaterLapse, inversionStrength, entrainmentIndex } from '../js/physics/radiation.module.js';
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
  const radiation = createRadiation(mesh, core);
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

test('the surface sees the direct beam in clear sky and diffuse light under thick cloud', () => {
  const radiation = createRadiation(mesh, core);
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

function deckColumns(options = {}) {
  const radiation = createRadiation(mesh, core, { stratus: true, ...options }), off = createRadiation(mesh, core, { stratus: false });
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
  assert.equal(run(createRadiation(mesh, core, { stratusScale: 10 }), stable, 298, 1, 1200).stratus, 0.15, 'the deck water stops at stratusWaterMax');
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
  const entraining = createRadiation(mesh, core, { stratusIndex: 'ectei' }), explicit = createRadiation(mesh, core, { stratusIndex: 'eis' });
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
  const shadow = createRadiation(mesh, core, { mixedLayerDeck: true, subsidenceMemory: 1e-9 }), off = createRadiation(mesh, core, { stratus: false });
  shadow.setTime(0); off.setTime(0);
  const noon = brightest(shadow), column = mixedLayerColumn(), { pi, theta, q, qc, depth } = column;
  const run = (r, openSea, dt = 900, air = theta) => mixedLayerRun(r, noon, column, openSea, dt, air);
  const deck = run(shadow, 1), clear = run(off, 1), start = run(shadow, 1, 0), half = run(shadow, 0.5);
  console.log(`289 K mixed layer under a 298 K free troposphere, 100 m above its third layer, over a 292 K sea: the mixed-layer model holds ${(1000 * deck.mlmWater).toFixed(1)} g/m² (${(1000 * start.mlmWater).toFixed(1)} before its step) on ${deck.mlmCover} of the cell, entraining ${(1000 * deck.mlmEntrainment).toFixed(2)} mm/s; noon absorbed solar ${deck.absorbedSolar.toFixed(0)} against ${clear.absorbedSolar.toFixed(0)} W/m²`);
  assert.ok(deck.mlmWater > 0.03 && deck.mlmWater < 0.15 && deck.mlmCover === 1, `water ${deck.mlmWater}, cover ${deck.mlmCover}`);
  assert.ok(deck.mlmEntrainment > 1e-3 && deck.mlmEntrainment < 1e-2, `entrainment ${deck.mlmEntrainment}`);
  assert.ok(deck.stratus === deck.mlmWater && deck.stratusFraction === deck.mlmCover, 'the deck is the mixed-layer model\'s');
  assert.ok(deck.mlmWater !== start.mlmWater && Math.abs(deck.mlmWater - start.mlmWater) < 0.05 * start.mlmWater, 'one step adjusts the water path');
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
  const shadow = createRadiation(mesh, core, { mixedLayerDeck: true, subsidenceMemory: 1e-9 });
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

test('the regime test reads the subsidence averaged over subsidenceMemory: a column that starts sinking gains its deck only once the running mean passes the floor', () => {
  const instant = createRadiation(mesh, core, { mixedLayerDeck: true, subsidenceMemory: 1e-9 }), memory = createRadiation(mesh, core, { mixedLayerDeck: true });
  instant.setTime(0); memory.setTime(0);
  const noon = brightest(memory), column = mixedLayerColumn(0.4), dt = 3600, keep = Math.exp(-dt / (3 * DAY));
  mixedLayerRun(instant, noon, column, 1, dt);
  const sinking = instant.mlmSubsidence[noon];
  let first = -1;
  for (let n = 1; n <= 24; n++) {
    const deck = mixedLayerRun(memory, noon, column, 1, dt);
    const expected = sinking * (1 - keep ** n);
    assert.ok(Math.abs(memory.mlmSubsidence[noon] - expected) < 1e-12 * Math.abs(sinking), `step ${n}: ${memory.mlmSubsidence[noon]} against ${expected}`);
    assert.equal(deck.mlmCover > 0, expected <= -3e-4, `step ${n}`);
    if (first < 0 && deck.mlmCover > 0) first = n;
  }
  console.log(`under a steady ${(1000 * sinking).toFixed(2)} mm/s the 3-day mean passes −0.3 mm/s and the deck appears after ${first} hourly steps`);
  assert.ok(first > 1 && first < 24);
  core.diagnostics.piSigmaDot.fill(0);
});

function modelDigest(radiation) {
  const model = createModel(new Grid(4), radiation ? { radiation } : {});
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  for (let n = 0; n < 12; n++) model.step(900);
  const hash = createHash('sha256');
  for (const a of [...model.state, model.radiation.stratus, model.radiation.stratusFraction, model.radiation.surfaceFlux, model.radiation.outgoing]) hash.update(new Uint8Array(a.buffer, a.byteOffset, a.byteLength));
  return { digest: hash.digest('hex').slice(0, 32), radiation: model.radiation };
}

test('with mixedLayerDeck: false the model is bit-identical to the engine before the mixed-layer deck; with it the deck follows the mixed-layer model', () => {
  const before = '1a7a2bd66c1edabedde1875fd6868c4a';
  assert.equal(modelDigest().digest, before);
  assert.equal(modelDigest({ mixedLayerDeck: false }).digest, before);
  const fresh = modelDigest({ mixedLayerDeck: true });
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

function evaluate(radiation, i, pi, theta, surfaceT, q, qc) {
  const flux = radiation.column(i, pi[i], theta, surfaceT[i], 5, radiation.opticalDepth(mesh.latCell[i]), radiation.insolation(i), q[(K - 1) * C + i], q, qc, 0.07);
  let layers = 0, scale = 0;
  for (let k = 0; k < K; k++) { layers += radiation.layerFlux[k]; scale += Math.abs(radiation.layerFlux[k]); }
  const latent = LATENT_HEAT * radiation.budget.evaporation;
  const residual = layers + flux + latent - (radiation.budget.absorbedSolar - radiation.budget.outgoingLongwave);
  return { reflected: radiation.budget.reflectedSolar, olr: radiation.budget.outgoingLongwave, absorbed: radiation.budget.absorbedSolar, surfaceShortwave: radiation.budget.surfaceShortwave, layers: Array.from(radiation.layerFlux), closure: Math.abs(residual) / (scale + Math.abs(flux) + latent) };
}
