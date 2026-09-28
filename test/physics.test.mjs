import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore, P0, CP_DRY, R_DRY } from '../js/dynamics/sigmaCore.module.js';
import { createRadiation, sunDirection, AXIAL_TILT, DAY, YEAR, waterVaporAbsorptivity, adiabaticWaterLapse } from '../js/physics/radiation.module.js';
import { LATENT_HEAT, saturationHumidity, liftingCondensationLevel } from '../js/physics/moist.module.js';
import { createSurface } from '../js/physics/surface.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';

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
    const flux = r.column(noon, P0, theta, surfaceT, 5, r.opticalDepth(mesh.latCell[noon]), r.insolation(noon), q[bottom], q, cloud, 0.07, 0.06, 1, 1.5e-3, openSea, depth);
    let layers = 0, scale = 0;
    for (let k = 0; k < K; k++) { layers += r.layerFlux[k]; scale += Math.abs(r.layerFlux[k]); }
    const latent = LATENT_HEAT * r.budget.evaporation;
    const closure = Math.abs(layers + flux + latent - (r.budget.absorbedSolar - r.budget.outgoingLongwave)) / (scale + Math.abs(flux) + latent);
    return { ...r.budget, flux, closure, layers: Array.from(r.layerFlux) };
  };
  return { radiation, off, noon, q, qc, unstable, bottom, withStability, run, airTemperature: unstable[bottom] * exnerLayer[bottom] };
}

test('a stable column over warm open sea carries a stratocumulus deck that reflects the noon sun and lowers the OLR; land, cold sea and stratus: false carry none', () => {
  const { radiation, off, unstable, withStability, run } = deckColumns();
  const stable = withStability(20);
  const deck = run(radiation, stable, 298, 1), clear = run(radiation, unstable, 298, 1);
  console.log(`LTS 20 K against ${(unstable[radiation.stabilityLayer * C] - unstable[(K - 1) * C]).toFixed(1)} K over a 298 K sea under a 1200 m boundary layer: deck ${deck.stratus.toFixed(4)} kg/m² on ${deck.stratusFraction.toFixed(3)} of the cell; noon absorbed solar ${deck.absorbedSolar.toFixed(0)} against ${clear.absorbedSolar.toFixed(0)} W/m², OLR ${deck.outgoingLongwave.toFixed(1)} against ${clear.outgoingLongwave.toFixed(1)}`);
  assert.ok(deck.stratus > 0.02 && clear.stratus === 0 && clear.stratusFraction === 0, `deck ${deck.stratus} against ${clear.stratus} kg/m²`);
  assert.ok(Math.abs(deck.stratusFraction - (0.057 * 20 - 0.556)) < 1e-9, `cover ${deck.stratusFraction}`);
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
  const half = run(radiation, withStability(1.056 / 0.057), 298, 1), full = run(radiation, withStability(30), 298, 1), none = run(radiation, withStability(30), 298, 0);
  assert.ok(Math.abs(half.stratusFraction - 0.5) < 1e-12 && full.stratusFraction === 1 && none.stratusFraction === 0, `covers ${half.stratusFraction}, ${full.stratusFraction}, ${none.stratusFraction}`);
  assert.ok(half.stratus > 0 && half.stratus === full.stratus);
  for (const name of ['reflectedSolar', 'absorbedSolar', 'surfaceShortwave', 'surfaceDirect', 'cloudReflectance']) {
    const mean = 0.5 * (full[name] + none[name]);
    assert.ok(Math.abs(half[name] - mean) < 1e-9, `${name} ${half[name]} against the mean ${mean}`);
  }
  console.log(`half cover of a ${half.stratus.toFixed(4)} kg/m² deck at noon reflects ${half.reflectedSolar.toFixed(2)} W/m², the mean of ${full.reflectedSolar.toFixed(2)} overcast and ${none.reflectedSolar.toFixed(2)} clear; closure ${half.closure.toExponential(1)}`);
  assert.ok(half.closure < EPS && full.closure < EPS, `closure ${half.closure}, ${full.closure}`);
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

function evaluate(radiation, i, pi, theta, surfaceT, q, qc) {
  const flux = radiation.column(i, pi[i], theta, surfaceT[i], 5, radiation.opticalDepth(mesh.latCell[i]), radiation.insolation(i), q[(K - 1) * C + i], q, qc, 0.07);
  let layers = 0, scale = 0;
  for (let k = 0; k < K; k++) { layers += radiation.layerFlux[k]; scale += Math.abs(radiation.layerFlux[k]); }
  const latent = LATENT_HEAT * radiation.budget.evaporation;
  const residual = layers + flux + latent - (radiation.budget.absorbedSolar - radiation.budget.outgoingLongwave);
  return { reflected: radiation.budget.reflectedSolar, olr: radiation.budget.outgoingLongwave, absorbed: radiation.budget.absorbedSolar, surfaceShortwave: radiation.budget.surfaceShortwave, layers: Array.from(radiation.layerFlux), closure: Math.abs(residual) / (scale + Math.abs(flux) + latent) };
}
