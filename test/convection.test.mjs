import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { ADJUST_KERNELS } from '../js/gpu/physics.gpu.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { MID_LEVEL, saturationHumidity, cloudSaturation, liquidFraction, LATENT_HEAT, FUSION_HEAT, MOIST_DEFAULTS, COUPLED_REGIME, DECK_CLOSED, BECHTOLD, SUBCLOUD_LAYERS, DEEP_CLOUD_DEPTH, IFS_ENTRAINMENT, TEST_PARCEL, DEPARTURE, IFS_PRECIPITATION, fogConstants } from '../js/physics/moist.module.js';
import { syntheticTopography, surfaceGeopotential } from '../js/geography.module.js';
import { VIRTUAL_FACTOR } from '../js/dynamics/sigmaCore.module.js';
import { REGIME } from '../js/physics/boundaryLayer.module.js';
import { PREVIOUS_MOIST, PREVIOUS_BERGERON, PREVIOUS_PCAPE_MINIMUM } from './helpers/previousDefaults.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }

const levels = sigmaInterfaces('bl34');
const build = (moist = {}, N = 16) => createModel(new Grid(N), { ocean: false, levels, moist });
const FOG_TOPOGRAPHY = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 : -4000));

/*
 * Jordan's (1958) mean West Indies sounding for the hurricane season:
 * hPa, °C and the mixing ratio in g/kg, with a stratosphere above.
 */
const JORDAN = [[1015, 26.3, 18.9], [1000, 25.4, 18.4], [950, 22.2, 16.6], [900, 19.6, 14.4], [850, 17.3, 12.4], [800, 14.8, 10.5], [750, 12.0, 8.8],
  [700, 8.9, 7.4], [650, 5.6, 6.0], [600, 1.8, 4.7], [550, -2.3, 3.6], [500, -6.9, 2.6], [450, -11.9, 1.8], [400, -17.3, 1.2], [350, -24.4, 0.66],
  [300, -32.4, 0.35], [250, -42.1, 0.14], [200, -54.0, 0.04], [175, -60.2, 0.02], [150, -67.2, 0.01], [125, -73.0, 0.005], [100, -73.3, 0.003],
  [70, -69, 0.003], [50, -62, 0.003], [30, -55, 0.003], [10, -45, 0.003], [1, -20, 0.003]];
function jordan(p) {
  if (p >= 100 * JORDAN[0][0]) return { T: 273.15 + JORDAN[0][1], q: 1e-3 * JORDAN[0][2] };
  const x = Math.log(p);
  for (let n = 0; n < JORDAN.length - 1; n++) {
    const [p0, t0, q0] = JORDAN[n], [p1, t1, q1] = JORDAN[n + 1];
    if ((p <= 100 * p0 && p >= 100 * p1) || n === JORDAN.length - 2) {
      const f = Math.max(0, (x - Math.log(100 * p0)) / (Math.log(100 * p1) - Math.log(100 * p0)));
      return { T: 273.15 + t0 + f * (t1 - t0), q: 1e-3 * (q0 + f * (q1 - q0)) };
    }
  }
  return null;
}

function saturatedAt(energy, p, guess) {
  let t = guess;
  for (let n = 0; n < 20; n++) {
    const qs = saturationHumidity(t, p);
    t -= (1004.64 * t + LATENT_HEAT * qs - energy) / (1004.64 + LATENT_HEAT * LATENT_HEAT * qs / (461.5 * t * t));
  }
  return t;
}

/*
 * A stratocumulus-topped column over a 26 °C sea: a mixed layer at
 * θ 297 K and 16 g/kg to its condensation level, saturated along its
 * moist adiabat above that up to the inversion at 1.3 km, whose θv
 * jump is 1.5 K, and dry air above (3 g/kg, at most 30 % humid) whose
 * θ rises 5 K/km to a 210 K stratosphere.
 */
function stratocumulus(z, p, below) {
  const kappa = 0.2857, exner = (pressure) => Math.pow(pressure / 1e5, kappa);
  const theta = 297, qt = 0.016;
  if (z < 1300) {
    const T = theta * exner(p);
    if (qt < saturationHumidity(T, p)) return { T, q: qt, qc: 0 };
    const t = saturatedAt(1004.64 * T + LATENT_HEAT * qt, p, T - 2);
    const qs = saturationHumidity(t, p);
    return { T: t, q: qs, qc: qt - qs };
  }
  const thetaV = below.thetaTop * (1 + 0.608 * below.qTop - below.qcTop) + 1.5;
  const thetaAbove = thetaV / (1 + 0.608 * 0.003) + 5e-3 * (z - 1300);
  const T = Math.max(210, thetaAbove * exner(p));
  return { T, q: Math.min(0.003, 0.3 * saturationHumidity(T, p)), qc: 0 };
}

/*
 * Column i of `model` set from profile(z, p, memo) (heights above the
 * surface), iterated so the heights are the core's.
 */
function place(model, i, surfacePressure, profile) {
  const { K, C, sigmaMid, exnerLayer, geopotential, g } = model.core.diagnostics;
  const [pi, theta, , , q, qc] = model.state;
  pi[i] = surfacePressure;
  for (let k = 0; k < K; k++) { theta[k * C + i] = 300; q[k * C + i] = 0; qc[k * C + i] = 0; }
  for (let pass = 0; pass < 4; pass++) {
    model.core.diagnoseColumn(i, pi, theta, q, qc);
    const memo = {};
    for (let k = K - 1; k >= 0; k--) {
      const idx = k * C + i, p = pi[i] * sigmaMid[k], z = geopotential[idx] / g;
      const air = profile(z, p, memo);
      theta[idx] = air.T / exnerLayer[idx];
      q[idx] = air.q;
      qc[idx] = air.qc ?? 0;
      if (z < 1300) { memo.thetaTop = theta[idx]; memo.qTop = q[idx]; memo.qcTop = qc[idx]; }
    }
  }
  model.core.diagnoseColumn(i, pi, theta, q, qc);
}
const jordanColumn = (model, i) => place(model, i, 101500, (z, p) => { const air = jordan(p); return { T: air.T, q: Math.min(air.q, saturationHumidity(air.T, p)) }; });

function snapshot(model, i) {
  const { K, C } = model.core.diagnostics, [, theta, , , q, qc] = model.state;
  return { theta: Float64Array.from({ length: K }, (_, k) => theta[k * C + i]), q: Float64Array.from({ length: K }, (_, k) => q[k * C + i]), qc: Float64Array.from({ length: K }, (_, k) => qc[k * C + i]) };
}
function budget(model, i) {
  const { K, C, dSigma, g, cp, exnerLayer } = model.core.diagnostics, [pi, theta, , , q, qc] = model.state;
  let enthalpy = 0, water = 0, heat = 0, vapour = 0;
  for (let k = 0; k < K; k++) {
    const idx = k * C + i, mass = pi[i] * dSigma[k] / g;
    heat += cp * theta[idx] * exnerLayer[idx] * mass;
    vapour += q[idx] * mass;
    enthalpy += (cp * theta[idx] * exnerLayer[idx] + LATENT_HEAT * q[idx]) * mass;
    water += (q[idx] + qc[idx]) * mass;
  }
  return { enthalpy, water, heat, vapour };
}
const setDepth = (model, i, height) => { const { K, C, geopotential, g } = model.core.diagnostics; model.boundaryLayer.depth[i] = geopotential[(K - 1) * C + i] / g + height; };

test('a stratocumulus-topped column under a 1.5 K inversion with dry air above lifts no deep plume and never rains', () => {
  const model = build(), { moist } = model, [pi, theta, , , q, qc] = model.state;
  place(model, 0, 101500, stratocumulus);
  setDepth(model, 0, 1300);
  model.boundaryLayer.buoyancyFlux[0] = 1e-4;
  model.boundaryLayer.friction[0] = 0.25;
  model.radiation.mlmGate[0] = 0.3;
  const { K, C } = model.core.diagnostics;
  let cloudy = 0;
  for (let k = 0; k < K; k++) if (qc[k * C] > 0) cloudy++;
  assert.ok(cloudy >= 2, `a deck of ${cloudy} layers`);
  let shallow = 0;
  for (let n = 0; n < 40; n++) {
    assert.equal(moist.plumeColumn(0, pi, theta, q, qc, 600), 0, `step ${n} rains`);
    assert.ok(!moist.deep.deep, `step ${n} lifts a deep plume`);
    if (moist.cumulusBaseFlux[0] > 0) shallow++;
  }
  console.log(`stratocumulus: cloud in ${cloudy} layers; a shallow plume on ${shallow} of 40 steps`);
});

test('the convection trace alone takes none of the condensation\'s heating', () => {
  const model = build(), { moist } = model, [pi, theta, , , q, qc] = model.state;
  const { K, C } = model.core.diagnostics;
  place(model, 0, 101500, stratocumulus);
  setDepth(model, 0, 1300);
  let cloudy = -1;
  for (let k = 0; k < K; k++) if (qc[k * C] > 0) cloudy = k;
  q[cloudy * C] += 1e-3;
  const before = snapshot(model, 0);
  moist.trace.convection = new Float64Array(K * C);
  moist.adjust(model.state, 0, 1, 600);
  assert.ok(theta[cloudy * C] > before.theta[cloudy], 'the supersaturated layer condenses and warms');
  assert.ok(moist.trace.convection.every((x) => x === 0), 'nothing is charged to convection');
});

/*
 * A trade-cumulus column: a mixed layer at θ 299 K and 17 g/kg to
 * 600 m, a conditionally unstable cloud layer (the Jordan sounding
 * warmed toward it) to 1.8 km, a 6 K inversion and dry air above.
 */
function tradeCumulus(z, p) {
  const kappa = 0.2857, exner = (pressure) => Math.pow(pressure / 1e5, kappa);
  if (z < 600) return { T: 299 * exner(p), q: 0.017 };
  const air = jordan(p);
  if (z < 1800) return { T: air.T - 1.5, q: Math.min(0.9 * saturationHumidity(air.T - 1.5, p), air.q + 0.002) };
  const T = air.T + 6;
  return { T, q: Math.min(0.002, 0.2 * saturationHumidity(T, p)) };
}

/*
 * A trade-wind column over a 26 °C sea: a mixed layer at θ 298.5 K and
 * 16 g/kg to 600 m, whose air is buoyant from its condensation level with
 * little inhibition below it, a conditionally unstable cloud layer whose
 * θ rises 3 K/km at 90 % humidity to the inversion at 1.3 km, a 1.5 K θv
 * jump there and dry air above (3 g/kg, at most 30 % humid) whose θ rises
 * 5 K/km to a 210 K stratosphere.
 */
function tradeWind(z, p) {
  const kappa = 0.2857, exner = (pressure) => Math.pow(pressure / 1e5, kappa);
  if (z < 600) return { T: 298.5 * exner(p), q: 0.016 };
  if (z < 1300) {
    const T = (298.5 + 3e-3 * (z - 600)) * exner(p);
    return { T, q: Math.min(0.9 * saturationHumidity(T, p), 0.016 - 3e-6 * (z - 600)) };
  }
  const thetaV = (298.5 + 3e-3 * 700) * (1 + 0.608 * (0.016 - 3e-6 * 700)) + 1.5;
  const T = Math.max(210, (thetaV / (1 + 0.608 * 0.003) + 5e-3 * (z - 1300)) * exner(p));
  return { T, q: Math.min(0.003, 0.3 * saturationHumidity(T, p)) };
}
function tradeWindColumn(options = {}, { buoyancy = 4e-4, gate = 0.3 } = {}) {
  const model = build(options);
  place(model, 0, 101500, tradeWind);
  setDepth(model, 0, 600);
  model.boundaryLayer.buoyancyFlux[0] = buoyancy;
  model.boundaryLayer.friction[0] = 0.25;
  model.radiation.mlmGate[0] = gate;
  return model;
}

test('a trade-wind column lifts a cumulus plume that tops in the inversion and detrains there, rains nothing, dries the boundary layer and moistens the inversion layer, with column enthalpy and water exact', () => {
  const model = tradeWindColumn(), { moist } = model, [pi, theta, , , q, qc] = model.state;
  const { K, C, sigmaMid, exnerLayer, geopotential, g, dSigma } = model.core.diagnostics, dt = 600;
  const before = budget(model, 0), snap = snapshot(model, 0);
  const rain = moist.cumulusColumn(0, pi, theta, q, qc, dt);
  const after = budget(model, 0), { top, source, inhibition } = moist.cumulus, base = moist.cumulusBaseFlux[0];
  const height = (k) => geopotential[k * C] / g - geopotential[(K - 1) * C] / g;
  const upper = (k) => height(k) + 0.5 * (height(k - 1) - height(k)), lower = (k) => height(k) - 0.5 * (height(k) - height(k + 1));
  const dq = (k) => (q[k * C] - snap.q[k]) / dt * 86400 * 1000;
  let boundary = 0, boundaryBefore = 0;
  for (let k = source; k < K; k++) { boundary += q[k * C] * pi[0] * dSigma[k] / g; boundaryBefore += snap.q[k] * pi[0] * dSigma[k] / g; }
  console.log(`trade wind: base mass flux ${base.toFixed(4)} kg/m²/s, inhibition ${inhibition.toFixed(2)} J/kg, source layers ${source}–${K - 1} (to ${upper(source).toFixed(0)} m), top layer ${top} at ${(pi[0] * sigmaMid[top] / 100).toFixed(0)} hPa (${lower(top).toFixed(0)}–${upper(top).toFixed(0)} m); moistening ${Array.from({ length: K - top }, (_, n) => `${(pi[0] * sigmaMid[top + n] / 100).toFixed(0)}: ${dq(top + n).toFixed(1)}`).join(', ')} g/kg/d; boundary layer water ${((boundary - boundaryBefore) / dt * 86400).toFixed(2)} kg/m²/d; enthalpy off by ${((after.enthalpy - before.enthalpy) / before.enthalpy).toExponential(1)}, water by ${((after.water - before.water) / before.water).toExponential(1)}`);
  assert.ok(base > 0.005 && base < 0.1, `base mass flux ${base}`);
  assert.ok(lower(top) < 1300 && upper(top) > 1300, `the plume tops at ${lower(top)}–${upper(top)} m`);
  assert.equal(rain, 0);
  assert.ok(dq(top) > 0, `the inversion layer moistens by ${dq(top)} g/kg/d`);
  assert.ok(boundary < boundaryBefore, 'the boundary layer dries');
  for (let k = 0; k < top; k++) assert.ok(theta[k * C] === snap.theta[k] && q[k * C] === snap.q[k], `layer ${k} above the top changed`);
  assert.ok(Math.abs(after.enthalpy - before.enthalpy) < 1e-12 * before.enthalpy, `enthalpy ${before.enthalpy} → ${after.enthalpy}`);
  assert.ok(Math.abs(after.water - before.water) < 1e-12 * before.water, `water ${before.water} → ${after.water}`);
  let covered = 0;
  for (let k = top; k < source; k++) if (moist.cumulusCover[k * C] > 0) covered++;
  assert.ok(covered > 0 && moist.cumulusCover.every((f) => f >= 0 && f <= 1), 'the plume layers have a cumulus fraction');
  const full = tradeWindColumn(), [, , , , fq, fqc] = full.state;
  const whole = budget(full, 0);
  full.moist.adjust(full.state, 0, 1, dt);
  const adjusted = budget(full, 0);
  assert.equal(full.moist.rain[0], 0, 'nothing rains through the adjustment');
  assert.ok(fqc[top * C] + fq[top * C] > snap.q[top] + snap.qc[top], 'the inversion layer gains water through the adjustment');
  assert.ok(Math.abs(adjusted.enthalpy - whole.enthalpy) < 1e-12 * whole.enthalpy && Math.abs(adjusted.water - whole.water) < 1e-12 * whole.water, 'enthalpy and water exact through the adjustment');
});

test('with cumulusMemory the cumulus cloud of a plume that fires every other step settles between on and off: its cover and cover × water relax toward each step\'s with e^(−Δt/τ), to half the plume\'s on average, and with 0 it is the step\'s own', () => {
  const dt = 600, tau = MOIST_DEFAULTS.cumulusMemory, keep = Math.exp(-dt / tau), steps = 24;
  const run = (options) => {
    const model = tradeWindColumn(options), { moist } = model, start = model.state.map((a) => Float64Array.from(a)), subcloud = Float64Array.from(moist.subcloudVirtual);
    const { K, C } = model.core.diagnostics, covers = [], paths = [], waters = [];
    for (let n = 0; n < steps; n++) {
      model.state.forEach((a, j) => a.set(start[j]));
      moist.subcloudVirtual.set(subcloud);
      model.boundaryLayer.buoyancyFlux[0] = n % 2 === 0 ? 4e-4 : 0;
      moist.adjust(model.state, 0, 1, dt);
      covers.push(Float64Array.from({ length: K }, (_, k) => moist.cumulusCover[k * C]));
      paths.push(Float64Array.from({ length: K }, (_, k) => moist.cumulusCover[k * C] * moist.cumulusWater[k * C]));
      waters.push(Float64Array.from({ length: K }, (_, k) => moist.cumulusWater[k * C]));
    }
    return { covers, paths, waters, K };
  };
  assert.equal(MOIST_DEFAULTS.cumulusMemory, 1800);
  const fresh = run({ cumulusMemory: 0 }), kept = run({}), { K } = fresh;
  const layers = Array.from({ length: K }, (_, k) => k).filter((k) => fresh.covers[0][k] > 0);
  assert.ok(layers.length > 0, 'the plume makes a cloud');
  for (let n = 0; n < steps; n++) for (let k = 0; k < K; k++) assert.equal(fresh.covers[n][k], n % 2 === 0 ? fresh.covers[0][k] : 0, `cumulusMemory 0: step ${n} layer ${k} is the plume's own`);
  let worst = 0;
  for (const k of layers) {
    let cover = 0, path = 0;
    for (let n = 0; n < steps; n++) {
      cover = fresh.covers[n][k] + (cover - fresh.covers[n][k]) * keep;
      path = fresh.paths[n][k] + (path - fresh.paths[n][k]) * keep;
      worst = Math.max(worst, Math.abs(kept.covers[n][k] - cover) / cover, Math.abs(kept.paths[n][k] - path) / path);
    }
    const x = fresh.covers[0][k], on = kept.covers[steps - 2][k], off = kept.covers[steps - 1][k];
    console.log(`layer ${k}: the plume's cover ${x.toFixed(4)} on alternate steps; with τ ${tau} s and Δt ${dt} s the cloud settles at ${on.toFixed(4)} after an on step and ${off.toFixed(4)} after an off step (x/(1 + e^(−Δt/τ)) ${(x / (1 + keep)).toFixed(4)}, mean ${(0.5 * (on + off)).toFixed(4)}); its water ${(1000 * kept.waters[steps - 1][k]).toFixed(4)} against the plume's ${(1000 * fresh.waters[0][k]).toFixed(4)} g/kg`);
    assert.ok(Math.abs(on - x / (1 + keep)) < 1e-3 * x && Math.abs(off - keep * x / (1 + keep)) < 1e-3 * x, `layer ${k}: ${on}, ${off} against ${x / (1 + keep)}, ${keep * x / (1 + keep)}`);
    assert.ok(Math.abs(0.5 * (on + off) - 0.5 * x) < 1e-3 * x, `layer ${k}: mean ${0.5 * (on + off)} against ${0.5 * x}`);
    assert.ok(Math.abs(kept.waters[steps - 1][k] - fresh.waters[0][k]) < 1e-12 * fresh.waters[0][k], `layer ${k}: the cloud keeps the plume's in-cloud water`);
    assert.ok(Math.abs(kept.covers[0][k] - (1 - keep) * x) < 1e-15, `layer ${k}: the first step from no cloud`);
  }
  assert.ok(worst < 1e-12, `the cover and path follow X ← X' + (X − X') e^(−Δt/τ) to ${worst}`);
});

test('with cumulusRain the plume rains its condensate above the threshold, keeping column enthalpy exact and water with the rain', () => {
  const model = tradeWindColumn({ cumulusRain: 2e-4 }), { moist } = model, [pi, theta, , , q, qc] = model.state;
  const before = budget(model, 0);
  const rain = moist.cumulusColumn(0, pi, theta, q, qc, 600);
  const after = budget(model, 0);
  assert.ok(rain > 0, `rain ${rain}`);
  assert.ok(Math.abs(after.enthalpy - before.enthalpy) < 1e-12 * before.enthalpy, `enthalpy ${before.enthalpy} → ${after.enthalpy}`);
  assert.ok(Math.abs(after.water + rain - before.water) < 1e-12 * before.water, `water ${before.water} → ${after.water} + ${rain}`);
});

test('a stable surface, an active deck or a condensation level above the shallow top gives no cumulus base mass flux and leaves the column untouched', () => {
  for (const [options, forcing, why] of [[{}, { buoyancy: -1e-4 }, 'a stable surface'], [{}, { buoyancy: 0 }, 'no surface buoyancy flux'], [{}, { gate: 0.8 }, 'an active deck'], [{ shallowTop: 980e2 }, {}, 'a condensation level above the shallow top']]) {
    const model = tradeWindColumn(options, forcing), [pi, theta, , , q, qc] = model.state, snap = snapshot(model, 0);
    assert.equal(model.moist.cumulusColumn(0, pi, theta, q, qc, 600), 0, why);
    assert.equal(model.moist.cumulusBaseFlux[0], 0, why);
    assert.deepEqual(snapshot(model, 0), snap, `${why}: the column is untouched`);
    assert.ok(model.moist.cumulusCover.every((f) => f === 0), `${why}: no cumulus cover`);
  }
  const halfOpen = tradeWindColumn({}, { gate: 0.55 }), open = tradeWindColumn();
  for (const m of [halfOpen, open]) { const [pi, theta, , , q, qc] = m.state; m.moist.cumulusColumn(0, pi, theta, q, qc, 600); }
  assert.ok(Math.abs(halfOpen.moist.cumulusBaseFlux[0] - 0.5 * open.moist.cumulusBaseFlux[0]) < 1e-12, 'half the base mass flux on the deck gate\'s ramp');
});

test('no column convects under an active deck, and an undecided gate lets it', () => {
  const model = build(), { moist } = model, [pi, theta, , , q, qc] = model.state;
  jordanColumn(model, 0);
  setDepth(model, 0, 500);
  model.boundaryLayer.buoyancyFlux[0] = 4e-4;
  model.boundaryLayer.friction[0] = 0.25;
  model.radiation.mlmGate[0] = 0.8;
  const before = snapshot(model, 0);
  assert.equal(moist.plumeColumn(0, pi, theta, q, qc, 600), 0);
  assert.equal(moist.cumulusBaseFlux[0], 0);
  assert.deepEqual(snapshot(model, 0), before);
  model.radiation.mlmGate[0] = 0.5;
  assert.ok(moist.plumeColumn(0, pi, theta, q, qc, 600) > 0 && moist.deep.deep, 'an undecided gate lets it rain');
});

test('without the fog scheme autoconversion stays out of the lowest two layers, or with autoconversionFloor: boundaryLayer out of the boundary layer, and rain evaporates only into cloud-free layers; an unknown floor is refused', () => {
  const off = { fogDroplets: null, fogDeposition: 0 };
  const lowest = build(off), [pl, tl, , , ql, qcl] = lowest.state, KL = lowest.core.K, CL = lowest.core.diagnostics.C;
  jordanColumn(lowest, 0);
  for (let k = 0; k < KL; k++) qcl[k * CL] = 0;
  qcl[(KL - 1) * CL] = 2e-3; qcl[(KL - 2) * CL] = 2e-3;
  assert.equal(lowest.moist.autoconvertColumn(0, pl, tl, ql, qcl, 600), 0, 'cloud in the lowest two layers does not rain');
  assert.ok(qcl[(KL - 1) * CL] === 2e-3 && qcl[(KL - 2) * CL] === 2e-3);
  qcl[(KL - 3) * CL] = 2e-3;
  assert.ok(lowest.moist.autoconvertColumn(0, pl, tl, ql, qcl, 600) > 0 || qcl[(KL - 3) * CL] < 2e-3, 'the third layer converts');
  assert.throws(() => build({ autoconversionFloor: 'surface' }), /'lowest', 'boundaryLayer', 'mixingTop' or 'none'/);
  const model = build({ autoconversionFloor: 'boundaryLayer', ...off }), { moist } = model, [pi, theta, , , q, qc] = model.state;
  const { K, C, geopotential, g } = model.core.diagnostics;
  jordanColumn(model, 0);
  setDepth(model, 0, 1000);
  const layerAt = (height) => { let best = K - 1; for (let k = 0; k < K; k++) if (Math.abs(geopotential[k * C] / g - height) < Math.abs(geopotential[best * C] / g - height)) best = k; return best; };
  const low = layerAt(500), high = layerAt(4000), cloudBelow = layerAt(2500), dry = layerAt(3000);
  for (let k = 0; k < K; k++) q[k * C] = 0.5 * saturationHumidity(theta[k * C] * model.core.diagnostics.exnerLayer[k * C], pi[0] * model.core.diagnostics.sigmaMid[k]);
  qc[low * C] = 2e-3;
  const qBefore = Float64Array.from({ length: K }, (_, k) => q[k * C]);
  assert.equal(moist.autoconvertColumn(0, pi, theta, q, qc, 600), 0, 'cloud inside the boundary layer does not rain');
  assert.equal(qc[low * C], 2e-3);
  qc[high * C] = 2e-3;
  qc[cloudBelow * C] = 1e-5;
  moist.autoconvertColumn(0, pi, theta, q, qc, 600);
  assert.ok(q[dry * C] > qBefore[dry], 'the rain evaporates into a cloud-free layer below its cloud');
  assert.equal(q[cloudBelow * C], qBefore[cloudBelow], 'and not into a cloudy one');
});

/*
 * The fog scheme of a floored layer k of column i holding `liquid`, by the
 * formulas of the moist module's header: what settles or deposits out of
 * it, what then drizzles, and what it keeps.
 */
function fogStep(model, i, k, liquid, dt, { droplets = MOIST_DEFAULTS.fogDroplets[1], drag = 0, efficiency = 0, limit = MOIST_DEFAULTS.fogDepositionLimit, settles = true } = {}) {
  const { R, g, sigmaMid, dSigma, exnerLayer, C } = model.core.diagnostics, [pi, theta] = model.state, idx = k * C + i;
  const temperature = theta[idx] * exnerLayer[idx], pressure = pi[i] * sigmaMid[k], mass = pi[i] * dSigma[k] / g;
  const settle = 7.32e5 * Math.pow(droplets * 1e6, -2 / 3), drizzle = 1350 * 1.5 * Math.pow(droplets, -1.79);
  let courant = 0;
  if (settles) courant = settle * Math.pow(pressure / (R * temperature) * liquid, 2 / 3) * dt * sigmaMid[k] * g / (R * temperature * dSigma[k]);
  const deposition = Math.min(efficiency * drag, limit * pressure / (R * temperature)) * dt * g / (pi[i] * dSigma[k]);
  courant += deposition;
  const leaving = liquid * courant / (1 + courant), kept = liquid - leaving;
  const converted = settles && kept > 0 ? -kept * Math.expm1(-drizzle * Math.pow(kept, 1.47) * dt) : 0;
  return { courant, deposition, leaving, kept, converted, mass, settle, drizzle };
}
const STOKES_GAMMA = 1000 * 9.81 / (18 * 1.8e-5) * 13 * 12 * Math.pow(6 * 40320 / (Math.PI * 1000 * 39916800), 2 / 3);

test('in the lowest two layers fog settles into the layer below and out of the lowest to the ground and drizzles, and the third layer converts by Kessler as before', () => {
  const dt = 600, model = build(), { moist } = model, [pi, theta, , , q, qc] = model.state, { K, C, g, dSigma } = model.core.diagnostics;
  jordanColumn(model, 0);
  for (let k = 0; k < K; k++) qc[k * C] = 0;
  model.boundaryLayer.surfaceDrag.fill(0);
  qc[(K - 1) * C] = 2e-3; qc[(K - 2) * C] = 2e-3;
  const before = budget(model, 0);
  const upper = fogStep(model, 0, K - 2, 2e-3, dt), settling = upper.leaving * upper.mass, arriving = 2e-3 + settling * g / (pi[0] * dSigma[K - 1]);
  const lower = fogStep(model, 0, K - 1, arriving, dt);
  const rain = moist.autoconvertColumn(0, pi, theta, q, qc, dt);
  const expected = upper.mass * upper.converted + lower.leaving * lower.mass + lower.mass * lower.converted;
  const after = budget(model, 0);
  console.log(`2 g/kg in the lowest two layers over ${dt} s at N ${MOIST_DEFAULTS.fogDroplets[1]}/cm³: the upper settles ${(1e3 * upper.leaving).toFixed(4)} and drizzles ${(1e3 * upper.converted).toFixed(5)} g/kg (Courant ${upper.courant.toFixed(3)}), the lowest ${(1e3 * lower.leaving).toFixed(4)} and ${(1e3 * lower.converted).toFixed(5)} g/kg; ${(1e3 * rain).toFixed(3)} g/m² reach the ground, ${(1e3 * moist.falling.fog).toFixed(3)} of them as fog`);
  assert.ok(Math.abs(rain - expected) <= 1e-15 * expected, `rain ${rain} against ${expected}`);
  assert.ok(Math.abs(qc[(K - 2) * C] - (2e-3 - upper.leaving - upper.converted)) <= 1e-15 * 2e-3 && Math.abs(qc[(K - 1) * C] - (arriving - lower.leaving - lower.converted)) <= 1e-15 * 2e-3);
  assert.ok(Math.abs(after.water + rain - before.water) <= 1e-15 * before.water, `water ${before.water} → ${after.water} + ${rain}`);
  assert.equal(moist.falling.moved, true);
  assert.equal(moist.falling.fog, lower.leaving * lower.mass);
  const third = (fog) => {
    const m = build(fog), [p3, t3, , , q3, c3] = m.state, K3 = m.core.K, C3 = m.core.diagnostics.C;
    jordanColumn(m, 0);
    for (let k = 0; k < K3; k++) c3[k * C3] = 0;
    c3[(K3 - 3) * C3] = 2e-3;
    m.moist.autoconvertColumn(0, p3, t3, q3, c3, dt);
    return c3[(K3 - 3) * C3];
  };
  const kept = third({}), unfogged = third({ fogDroplets: null, fogDeposition: 0 });
  assert.ok(kept < 2e-3 && kept === unfogged, `the third layer converts by Kessler alone: ${kept} against ${unfogged}`);
});

test('with autoconversionFloor boundaryLayer the fog scheme settles cloud down through the boundary layer and out of it, and leaves the layers above alone', () => {
  const dt = 600, model = build({ autoconversionFloor: 'boundaryLayer' }), { moist } = model, [pi, theta, , , q, qc] = model.state;
  const { K, C, geopotential, g, dSigma, exnerLayer, sigmaMid } = model.core.diagnostics;
  jordanColumn(model, 0);
  setDepth(model, 0, 1000);
  model.boundaryLayer.surfaceDrag.fill(0);
  for (let k = 0; k < K; k++) { qc[k * C] = 0; q[k * C] = 0.5 * saturationHumidity(theta[k * C] * exnerLayer[k * C], pi[0] * sigmaMid[k]); }
  const height = (k) => geopotential[k * C] / g - geopotential[(K - 1) * C] / g;
  let top = K;
  for (let k = K - 1; k >= 0; k--) if (height(k) < 700) { qc[k * C] = 2e-3 * (1 + 0.1 * (K - k)); top = k; }
  const was = snapshot(model, 0), before = budget(model, 0), mass = (k) => pi[0] * dSigma[k] / g;
  const rain = moist.autoconvertColumn(0, pi, theta, q, qc, dt), after = budget(model, 0), now = snapshot(model, 0);
  let above = 0, below = 0;
  for (let k = 0; k < K; k++) {
    above += (was.qc[k] + was.q[k]) * mass(k); below += (now.qc[k] + now.q[k]) * mass(k);
    assert.ok(below <= above * (1 + 1e-15), `layer ${k}: the water above its base grew from ${above} to ${below}`);
    if (height(k) > 1100) assert.ok(now.theta[k] === was.theta[k] && now.q[k] === was.q[k] && now.qc[k] === was.qc[k], `layer ${k} above the boundary layer changed`);
  }
  for (let k = top; k < K; k++) assert.ok(now.qc[k] > 0 && now.qc[k] < was.qc[k] + 2e-3, `layer ${k}: ${now.qc[k]}`);
  assert.ok(rain > 0 && Math.abs(after.water + rain - before.water) <= 1e-15 * before.water);
  assert.ok(moist.falling.fog > 0 && moist.falling.fog < rain);
});

test('with autoconversionFloor mixingTop the cloud between the boundary-layer depth and the mixing top is floored too, settling and drizzling instead of converting over the lifetime, and the layers above the mixing top are left alone', () => {
  const dt = 600, run = (floor) => {
    const model = build({ autoconversionFloor: floor }), { moist } = model, [pi, theta, , , q, qc] = model.state;
    const { K, C, geopotential, g, exnerLayer, sigmaMid } = model.core.diagnostics;
    jordanColumn(model, 0);
    setDepth(model, 0, 700);
    model.boundaryLayer.mixingTop[0] = model.boundaryLayer.depth[0] + 600;
    model.boundaryLayer.surfaceDrag.fill(0);
    const height = (k) => geopotential[k * C] / g - geopotential[(K - 1) * C] / g;
    for (let k = 0; k < K; k++) { q[k * C] = 0.5 * saturationHumidity(theta[k * C] * exnerLayer[k * C], pi[0] * sigmaMid[k]); qc[k * C] = height(k) < 1500 ? 2e-4 : 0; }
    const was = snapshot(model, 0);
    moist.autoconvertColumn(0, pi, theta, q, qc, dt);
    return { was, now: snapshot(model, 0), height, K };
  };
  const mixed = run('mixingTop'), layer = run('boundaryLayer'), lifetime = 1 - Math.exp(-dt / 3600);
  let checked = 0;
  for (let k = 0; k < mixed.K; k++) {
    const z = mixed.height(k);
    if (z > 1500) { assert.equal(mixed.now.qc[k], mixed.was.qc[k], `layer ${k} above the mixing top changed`); continue; }
    if (z > 750 && z < 1200) {
      checked++;
      assert.ok(layer.was.qc[k] - layer.now.qc[k] > 0.9 * lifetime * layer.was.qc[k], `boundaryLayer: layer ${k} converts over the lifetime`);
      assert.ok(mixed.now.qc[k] > 0 && mixed.was.qc[k] - mixed.now.qc[k] < 0.5 * lifetime * mixed.was.qc[k], `mixingTop: layer ${k} keeps its cloud but for settling and drizzle (${mixed.was.qc[k]} to ${mixed.now.qc[k]})`);
    }
  }
  assert.ok(checked >= 2, `${checked} layers between the depth and the mixing top`);
});

test('fog settles implicitly and conservatively: a floored layer keeps 1/(1 + c) of its liquid, at the mass-weighted Stokes speed of a gamma distribution, never emptying a layer it feeds, with moist enthalpy exact, and only from the liquid share while ice falls', () => {
  for (const [N, cm] of [[150, 1.16], [60, 2.14]]) {
    const { settleLand, settleSea } = fogConstants([N, N]);
    assert.equal(settleLand, settleSea);
    assert.ok(Math.abs(STOKES_GAMMA - 7.32e5) < 1e-3 * 7.32e5, `the derivation gives ${STOKES_GAMMA}`);
    const speed = settleLand * Math.pow(0.3e-3, 2 / 3);
    assert.ok(Math.abs(100 * speed - cm) < 5e-3 && Math.abs(speed - STOKES_GAMMA * Math.pow(N * 1e6, -2 / 3) * Math.pow(0.3e-3, 2 / 3)) < 1e-3 * speed, `${N}/cm³: ${100 * speed} cm/s at 0.3 g/m³`);
  }
  const dt = 600, model = build(), { moist } = model, [pi, theta, , , q, qc] = model.state, { K, C, sigmaMid, exnerLayer } = model.core.diagnostics;
  jordanColumn(model, 0);
  model.boundaryLayer.surfaceDrag.fill(0);
  for (let k = 0; k < K; k++) qc[k * C] = 0;
  qc[(K - 1) * C] = 1e-3;
  const one = fogStep(model, 0, K - 1, 1e-3, dt);
  moist.autoconvertColumn(0, pi, theta, q, qc, dt);
  assert.ok(Math.abs(qc[(K - 1) * C] + one.converted - 1e-3 / (1 + one.courant)) <= 1e-15 * 1e-3, 'the lowest layer keeps 1/(1 + c)');

  const stack = build({ autoconversionFloor: 'boundaryLayer' }), [ps, ts, , , qs, cs] = stack.state;
  jordanColumn(stack, 0);
  setDepth(stack, 0, 1500);
  stack.boundaryLayer.surfaceDrag.fill(0);
  const { geopotential, g } = stack.core.diagnostics, fed = [];
  for (let k = 0; k < K; k++) { cs[k * C] = 0; if (geopotential[k * C] / g - geopotential[(K - 1) * C] / g < 1400) { cs[k * C] = 1e-3; fed.push(k); } }
  const start = budget(stack, 0);
  let ground = 0;
  for (const step of [60, 600, 6000, 60000]) {
    ground += stack.moist.autoconvertColumn(0, ps, ts, qs, cs, step);
    for (let k = 0; k < K; k++) assert.ok(cs[k * C] >= 0 && Number.isFinite(cs[k * C]), `layer ${k}: ${cs[k * C]} after ${step} s`);
    for (const k of fed) assert.ok(cs[k * C] > 0, `layer ${k} emptied after ${step} s`);
  }
  assert.ok(fed.length > 8 && Math.abs(budget(stack, 0).water + ground - start.water) <= 1e-14 * start.water);

  const dryBelow = build(), [pd, td, , , qd, cd] = dryBelow.state;
  jordanColumn(dryBelow, 0);
  dryBelow.boundaryLayer.surfaceDrag.fill(0);
  for (let k = 0; k < K; k++) cd[k * C] = 0;
  qd[(K - 2) * C] = saturationHumidity(td[(K - 2) * C] * exnerLayer[(K - 2) * C], pd[0] * sigmaMid[K - 2]);
  qd[(K - 1) * C] *= 0.7;
  cd[(K - 2) * C] = 3e-3;
  dryBelow.core.diagnoseColumn(0, pd, td, qd, cd);
  const enthalpy = budget(dryBelow, 0), vapourBelow = qd[(K - 1) * C];
  const fell = dryBelow.moist.autoconvertColumn(0, pd, td, qd, cd, dt);
  assert.ok(dryBelow.moist.falling.moved && cd[(K - 1) * C] > 0, 'fog settles into the subsaturated lowest layer');
  dryBelow.moist.condenseColumn(0, pd, td, qd, cd);
  const settled = budget(dryBelow, 0);
  assert.ok(qd[(K - 1) * C] > vapourBelow && cd[(K - 1) * C] === 0, 'and evaporates there in the second adjustment');
  assert.ok(Math.abs(settled.enthalpy - enthalpy.enthalpy) <= 1e-14 * enthalpy.enthalpy, `moist enthalpy ${enthalpy.enthalpy} → ${settled.enthalpy}`);
  assert.ok(Math.abs(settled.water + fell - enthalpy.water) <= 1e-15 * enthalpy.water);

  for (const iceFall of [MOIST_DEFAULTS.iceFall, null]) {
    const cold = (fog) => {
      const m = build({ iceFall, ...fog }), [pc, tc, , , qq, cc] = m.state;
      place(m, 0, 101500, (z, p) => { const T = 255 + 0.004 * z; return { T, q: cloudSaturation(T, p, true, 273.15, 235.15, {}).qs, qc: 0 }; });
      m.boundaryLayer.surfaceDrag.fill(0);
      cc[(K - 1) * C] = 2e-3;
      m.core.diagnoseColumn(0, pc, tc, qq, cc);
      const share = liquidFraction(tc[(K - 1) * C] * m.core.diagnostics.exnerLayer[(K - 1) * C], 273.15, 235.15);
      const liquid = iceFall === null ? 2e-3 : share * 2e-3, expected = fogStep(m, 0, K - 1, liquid, dt);
      m.moist.autoconvertColumn(0, pc, tc, qq, cc, dt);
      return { qc: cc[(K - 1) * C], share, expected };
    };
    const on = cold({}), off = cold({ fogDroplets: null, fogDeposition: 0 });
    console.log(`supercooled fog at liquid share ${on.share.toFixed(3)} with iceFall ${iceFall}: the fog takes ${(1e3 * (off.qc - on.qc)).toFixed(4)} g/kg of the ${(2 * (iceFall === null ? 1 : on.share)).toFixed(3)} g/kg it acts on`);
    assert.ok(on.share > 0.3 && on.share < 0.9);
    assert.ok(Math.abs(off.qc - on.qc - on.expected.leaving - on.expected.converted) <= 1e-15 * 2e-3, `the fog takes ${off.qc - on.qc} against ${on.expected.leaving + on.expected.converted}`);
  }
});

test('the floored layers drizzle by Khairoutdinov and Kogan, negligibly at fog\'s water and not at all without the fog scheme', () => {
  const dt = 600, model = build(), { K, C } = model.core.diagnostics;
  jordanColumn(model, 0);
  model.boundaryLayer.surfaceDrag.fill(0);
  const run = (amount, m = model) => {
    const [p, t, , , qq, cc] = m.state;
    for (let k = 0; k < K; k++) cc[k * C] = 0;
    cc[(K - 1) * C] = amount;
    const expected = fogStep(m, 0, K - 1, amount, dt);
    const rain = m.moist.autoconvertColumn(0, p, t, qq, cc, dt);
    return { rain, expected, qc: cc[(K - 1) * C] };
  };
  const thick = run(2e-3), { drizzleSea } = fogConstants(MOIST_DEFAULTS.fogDroplets);
  assert.equal(thick.expected.drizzle, drizzleSea);
  const kk = thick.expected.kept * (1 - Math.exp(-drizzleSea * Math.pow(thick.expected.kept, 1.47) * dt));
  assert.ok(Math.abs(thick.expected.kept - thick.qc - kk) <= 1e-12 * kk && kk > 0, `drizzle ${thick.expected.kept - thick.qc} against ${kk}`);
  assert.ok(Math.abs(thick.rain - thick.expected.mass * (thick.expected.leaving + kk)) <= 1e-12 * thick.rain);
  const thin = run(3e-4), share = thin.expected.converted / thin.expected.kept;
  console.log(`over ${dt} s at N 60/cm³ drizzle takes ${(100 * kk / thick.expected.kept).toFixed(3)} % of 2 g/kg and ${(100 * share).toFixed(3)} % of 0.3 g/kg`);
  assert.ok(share > 0 && share < 1e-2, `0.3 g/kg drizzles ${share} of itself`);
  const none = build({ fogDroplets: null, fogDeposition: 0 });
  jordanColumn(none, 0);
  const plain = run(2e-3, none);
  assert.ok(plain.rain === 0 && plain.qc === 2e-3);
});

test('fog stays bounded: a lowest layer fed 200 mm/d holds under 10 g/kg, one fed 50 g/m²/h sits at the scalar fixed point near observed fog water, and 75 g/kg falls below 1 g/kg within six steps over the sea and seven over land, with water exact', () => {
  const model = createModel(new Grid(16), { topography: FOG_TOPOGRAPHY, ocean: false, levels: sigmaInterfaces('bl36') }), { moist } = model, { K, C, g, dSigma } = model.core.diagnostics;
  const [pi, theta, , , q, qc] = model.state, land = model.geography.land, iceSheet = model.geography.iceSheet;
  const landCell = [...Array(C).keys()].find((i) => land[i] && !iceSheet[i] && Math.abs(model.mesh.latCell[i]) < 0.8), seaCell = [...Array(C).keys()].find((i) => !land[i] && Math.abs(model.mesh.latCell[i]) < 0.8);
  const rows = [];
  for (const [name, i, N] of [['land', landCell, MOIST_DEFAULTS.fogDroplets[0]], ['sea', seaCell, MOIST_DEFAULTS.fogDroplets[1]]]) {
    jordanColumn(model, i);
    const rho = pi[i] * model.core.diagnostics.sigmaMid[K - 1] / (model.core.diagnostics.R * theta[(K - 1) * C + i] * model.core.diagnostics.exnerLayer[(K - 1) * C + i]);
    for (const dt of [337.5, 168.75]) {
      const mass = pi[i] * dSigma[K - 1] / g, drag = rho * 0.0012 * 3;
      const feed = (rate, start, steps, depositing) => {
        model.boundaryLayer.surfaceDrag[i] = depositing ? drag : 0;
        for (let k = 0; k < K; k++) qc[k * C + i] = 0;
        qc[(K - 1) * C + i] = start;
        const before = budget(model, i);
        let most = start, ground = 0, given = 0;
        const trace = [];
        for (let n = 0; n < steps; n++) {
          qc[(K - 1) * C + i] += rate * dt / mass; given += rate * dt;
          ground += moist.autoconvertColumn(i, pi, theta, q, qc, dt);
          most = Math.max(most, qc[(K - 1) * C + i]);
          trace.push(qc[(K - 1) * C + i]);
          assert.ok(Number.isFinite(qc[(K - 1) * C + i]));
        }
        assert.ok(Math.abs(budget(model, i).water + ground - given - before.water) <= 1e-14 * before.water, 'water exact');
        return { most, last: qc[(K - 1) * C + i], trace };
      };
      const days = Math.round(2 * 86400 / dt);
      const heavy = feed(200 / 86400, 0, days, true);
      assert.ok(heavy.most < 10e-3, `${name} at ${dt} s: 200 mm/d holds up to ${heavy.most}`);
      const moderate = feed(50e-3 / 3600, 0, days, false);
      let fixed = 0;
      for (let n = 0; n < 4 * days; n++) { const l = fixed + 50e-3 / 3600 * dt / mass, step = fogStep(model, i, K - 1, l, dt, { droplets: N }); fixed = l - step.leaving - step.converted; }
      assert.ok(Math.abs(moderate.last - fixed) < 0.05 * fixed, `${name} at ${dt} s: ${moderate.last} against the fixed point ${fixed}`);
      const deposited = feed(50e-3 / 3600, 0, days, true);
      const decay = feed(0, 75e-3, 7, true);
      if (dt === 337.5) assert.ok(decay.trace[name === 'sea' ? 5 : 6] < 1e-3, `${name}: 75 g/kg decays to ${decay.trace.join(', ')}`);
      rows.push(`${name} dt ${dt}: 200 mm/d max ${(1e3 * heavy.most).toFixed(2)}, 50 g/m²/h ${(1e3 * moderate.last).toFixed(3)} (fixed point ${(1e3 * fixed).toFixed(3)}; with deposition at C_D 0.0012 W 3 ${(1e3 * deposited.last).toFixed(3)}), 75 g/kg after six and seven steps ${(1e3 * decay.trace[5]).toFixed(2)} and ${(1e3 * decay.trace[6]).toFixed(2)} g/kg`);
    }
  }
  console.log(rows.join('\n'));
});

test('fog deposits onto the surface at E C_D W up to fogDepositionLimit, with the continental droplet number on land and the sea\'s on ice sheets, and not without the implicit drag', () => {
  const dt = 600, make = (moist = {}, surface = {}) => createModel(new Grid(16), { topography: FOG_TOPOGRAPHY, ocean: false, levels, moist, surface });
  const model = make(), { K, C } = model.core.diagnostics, land = model.geography.land, iceSheet = model.geography.iceSheet;
  const landCell = [...Array(C).keys()].find((i) => land[i] && !iceSheet[i] && Math.abs(model.mesh.latCell[i]) < 0.8), sheetCell = [...Array(C).keys()].find((i) => iceSheet[i]);
  assert.ok(landCell !== undefined && sheetCell !== undefined);
  const lowest = (m, i, drag, amount = 1e-3) => {
    const [pi, theta, , , q, qc] = m.state;
    jordanColumn(m, i);
    for (let k = 0; k < K; k++) qc[k * C + i] = 0;
    qc[(K - 1) * C + i] = amount;
    m.boundaryLayer.surfaceDrag[i] = drag;
    const rain = m.moist.autoconvertColumn(i, pi, theta, q, qc, dt);
    return { rain, qc: qc[(K - 1) * C + i], fog: m.moist.falling.fog };
  };
  const drag = 0.06, [N_LAND, N_SEA] = MOIST_DEFAULTS.fogDroplets;
  const deposited = lowest(model, landCell, drag), expected = fogStep(model, landCell, K - 1, 1e-3, dt, { droplets: N_LAND, drag, efficiency: 1 });
  const { R, g, dSigma } = model.core.diagnostics, pi = model.state[0];
  assert.ok(Math.abs(expected.deposition - drag * dt * g / (pi[landCell] * dSigma[K - 1])) <= 1e-15 * expected.deposition, 'c_d = E surfaceDrag Δt g/(π Δσ)');
  assert.ok(Math.abs(deposited.qc - (expected.kept - expected.converted)) <= 1e-15 * 1e-3 && deposited.fog === expected.leaving * expected.mass, `kept ${deposited.qc} against ${expected.kept - expected.converted}`);
  const settledOnly = fogStep(model, landCell, K - 1, 1e-3, dt, { droplets: N_LAND });
  console.log(`1 g/kg over a continental cell for ${dt} s: settling alone takes ${(1e3 * settledOnly.leaving).toFixed(4)} g/kg (c ${settledOnly.courant.toFixed(3)}), with ρ C_D W ${drag} kg/m²/s of deposition ${(1e3 * expected.leaving).toFixed(4)} (c_d ${expected.deposition.toFixed(3)})`);
  const off = make({ fogDeposition: 0 }), noDeposit = lowest(off, landCell, drag);
  assert.ok(Math.abs(noDeposit.qc - (settledOnly.kept - settledOnly.converted)) <= 1e-15 * 1e-3, 'fogDeposition 0 deposits nothing');
  const explicit = make({}, { implicitDrag: false }), noDrag = lowest(explicit, landCell, drag);
  assert.ok(Math.abs(noDrag.qc - (settledOnly.kept - settledOnly.converted)) <= 1e-15 * 1e-3, 'without the implicit drag nothing deposits');
  const capped = lowest(model, landCell, 10), limit = fogStep(model, landCell, K - 1, 1e-3, dt, { droplets: N_LAND, drag: 10, efficiency: 1 });
  const temperature = model.state[1][(K - 1) * C + landCell] * model.core.diagnostics.exnerLayer[(K - 1) * C + landCell], rho = pi[landCell] * model.core.diagnostics.sigmaMid[K - 1] / (R * temperature);
  assert.ok(10 > 0.25 * rho && Math.abs(limit.deposition - 0.25 * rho * dt * g / (pi[landCell] * dSigma[K - 1])) <= 1e-15 * limit.deposition && Math.abs(capped.qc - (limit.kept - limit.converted)) <= 1e-15 * 1e-3, 'the limit binds');
  const sheet = lowest(model, sheetCell, 0), sheetExpected = fogStep(model, sheetCell, K - 1, 1e-3, dt, { droplets: N_SEA }), sheetAsLand = fogStep(model, sheetCell, K - 1, 1e-3, dt, { droplets: N_LAND });
  assert.ok(Math.abs(sheet.qc - (sheetExpected.kept - sheetExpected.converted)) <= 1e-15 * 1e-3 && Math.abs(sheet.qc - (sheetAsLand.kept - sheetAsLand.converted)) > 1e-8, 'an ice-sheet cell settles at the sea\'s droplet number');
});

test('fog reaching the ground follows the surface\'s phase rule: over cold land it is snow and its fusion heat warms the lowest layer, over warm sea it is rain', () => {
  const dt = 600, model = createModel(new Grid(16), { topography: FOG_TOPOGRAPHY, ocean: false, levels, boundaryLayer: false, moist: { cumulusMemory: 0 } });
  const { moist, seaIce } = model, { K, C, g, cp, dSigma, exnerLayer } = model.core.diagnostics, [pi, theta] = model.state, land = model.geography.land, iceSheet = model.geography.iceSheet;
  const landCell = [...Array(C).keys()].find((i) => land[i] && !iceSheet[i] && Math.abs(model.mesh.latCell[i]) < 0.8), seaCell = [...Array(C).keys()].find((i) => !land[i] && Math.abs(model.mesh.latCell[i]) < 0.8);
  const fog = (i, surfaceT) => place(model, i, 100000, (height, p) => { const z = height - (model.surfaceGeopotential ? model.surfaceGeopotential[i] / g : 0), T = surfaceT + 0.03 * Math.min(z, 300) - 0.0065 * Math.max(0, z - 300); return { T, q: (z < 150 ? 1 : 0.6) * cloudSaturation(T, p, true, 273.15, 235.15, {}).qs, qc: z < 100 ? 1.5e-3 : 0 }; });
  for (const [i, surfaceT, frozen] of [[landCell, 265, true], [seaCell, 285, false]]) {
    fog(i, surfaceT);
    model.state[3][i] = surfaceT;
    const kept = model.state.map((a) => Float64Array.from(a)), snowBefore = land[i] ? model.land.snow[i] : seaIce.snow[i];
    moist.adjust(model.state, i, i + 1, dt);
    const amount = moist.rain[i], fogged = moist.falling.fog, thetaMoist = theta[(K - 1) * C + i], plumeSnow = moist.convectiveSnow[i];
    model.state.forEach((a, n) => a.set(kept[n]));
    model.phases.adjust(i, i + 1, dt);
    const snowAfter = land[i] ? model.land.snow[i] : seaIce.snow[i];
    const heating = frozen ? seaIce.latentHeatFusion * (amount - plumeSnow) * g / (cp * pi[i] * dSigma[K - 1] * exnerLayer[(K - 1) * C + i]) : 0;
    console.log(`${frozen ? 'supercooled fog at 265 K over land' : 'fog at 285 K over the sea'}: ${(1e3 * amount).toFixed(3)} g/m² reach the ground, ${(1e3 * fogged).toFixed(3)} of them as fog; snow +${(1e3 * (snowAfter - snowBefore)).toFixed(3)} g/m², fusion heating ${heating.toExponential(2)} K of θ`);
    assert.ok(fogged > 0 && amount >= fogged);
    if (frozen) assert.ok(snowAfter - snowBefore === amount || Math.abs(snowAfter - snowBefore - amount) <= 1e-15 * (snowBefore + amount), `snow rose by ${snowAfter - snowBefore} against ${amount}`);
    else assert.equal(snowAfter, snowBefore);
    assert.ok(Math.abs(theta[(K - 1) * C + i] - (thetaMoist + heating)) <= 1e-12, `θ ${theta[(K - 1) * C + i]} against ${thetaMoist} + ${heating}`);
  }
});

test('the fog options are checked on both engines', async () => {
  const refused = [{ fogDroplets: 'x' }, { fogDroplets: [0, 60] }, { fogDroplets: [150] }, { fogDroplets: [NaN, 60] }, { fogDeposition: -1 }, { fogDeposition: NaN }, { fogDepositionLimit: 0 }];
  for (const options of refused) assert.throws(() => build(options), /fogD/, `CPU: ${JSON.stringify(options)}`);
  assert.throws(() => build({ autoconversionFloor: 'surface' }), /autoconversionFloor/);
  assert.doesNotThrow(() => build({ fogDroplets: null, fogDeposition: 0 }));
  if (!gpuAvailable) return;
  const { physicsConstants } = await import('../js/gpu/physics.gpu.js');
  const { PHYSICS_DEFAULTS, createGpuCore } = await import('../js/gpu/core.gpu.js');
  for (const options of refused) assert.throws(() => physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, ...options }), /fogD/, `GPU: ${JSON.stringify(options)}`);
  await assert.rejects(() => createGpuCore(build().mesh, { levels, physics: { autoconversionFloor: 'surface' } }), /autoconversionFloor/);
});

test('with upperCloudLifetime cloud water in the layers above the shallow top converts over that lifetime and the cloud below over cloudLifetime (with no ice falling, all of it as liquid)', () => {
  const dt = 600, model = build({ upperCloudLifetime: 1800, iceFall: null, ...PREVIOUS_BERGERON }), plain = build({ iceFall: null, ...PREVIOUS_BERGERON }), { K, C, sigmaMid } = model.core.diagnostics;
  for (const m of [model, plain]) jordanColumn(m, 0);
  const pressure = (k) => model.state[0][0] * sigmaMid[k];
  let upper = 0, lower = K - 3;
  while (pressure(upper) < 400e2) upper++;
  while (pressure(lower) < 800e2) lower++;
  assert.ok(pressure(upper) < MOIST_DEFAULTS.shallowTop && pressure(lower) > MOIST_DEFAULTS.shallowTop && lower < K - 2);
  const after = (m) => {
    const [pi, theta, , , q, qc] = m.state;
    for (let k = 0; k < K; k++) qc[k * C] = 0;
    qc[upper * C] = 1e-4; qc[lower * C] = 1e-4;
    m.moist.autoconvertColumn(0, pi, theta, q, qc, dt);
    return [qc[upper * C], qc[lower * C]];
  };
  const [up, low] = after(model), [plainUp, plainLow] = after(plain);
  assert.ok(Math.abs(up - 1e-4 * Math.exp(-dt / 1800)) < 1e-18, `upper ${up}`);
  assert.ok(Math.abs(low - 1e-4 * Math.exp(-dt / MOIST_DEFAULTS.cloudLifetime)) < 1e-18, `lower ${low}`);
  assert.equal(low, plainLow);
  assert.ok(Math.abs(plainUp - 1e-4 * Math.exp(-dt / MOIST_DEFAULTS.cloudLifetime)) < 1e-18, `upper without the option ${plainUp}`);
});

test('stratiform cloud converts over stratiformLifetime: below the mixing top of a coupled column, above it by the EIS share, over sea ice by its cover; cloud the plume left and cloud below the top of a surface-driven column over cloudLifetime', () => {
  const dt = 600, model = build(), plain = build({ stratiformLifetime: null }), { K, C, sigmaMid, geopotential, g } = model.core.diagnostics;
  for (const m of [model, plain]) jordanColumn(m, 0);
  let k = K - 3;
  while (model.state[0][0] * sigmaMid[k] < 900e2) k++;
  assert.ok(k < K - 2);
  const z = geopotential[k * C] / g, short = MOIST_DEFAULTS.cloudLifetime, long = MOIST_DEFAULTS.stratiformLifetime;
  const after = (m, { regime, inside, share = 0, iced = 0, plumeTop = 0 }) => {
    const [pi, theta, , , q, qc] = m.state;
    for (let j = 0; j < K; j++) qc[j * C] = 0;
    qc[k * C] = 1e-4;
    m.boundaryLayer.regime[0] = regime;
    m.boundaryLayer.mixingTop[0] = inside ? z + 100 : z - 100;
    m.radiation.stratiform[0] = share;
    m.moist.cumulusBaseFlux[0] = plumeTop > 0 ? 0.01 : 0;
    m.moist.cumulusTop[0] = plumeTop;
    m.moist.autoconvertColumn(0, pi, theta, q, qc, dt, null, iced);
    return qc[k * C];
  };
  const converts = (lifetime) => 1e-4 * Math.exp(-dt / lifetime);
  const cases = [
    ['coupled, below the mixing top', { regime: REGIME.COUPLED, inside: true }, long],
    ['surface-driven, below the mixing top', { regime: REGIME.SURFACE, inside: true, share: 1 }, short],
    ['decoupled, below the mixing top', { regime: REGIME.DECOUPLED, inside: true }, short],
    ['above the mixing top, EIS share 0.5', { regime: REGIME.SURFACE, inside: false, share: 0.5 }, short + 0.5 * (long - short)],
    ['above the mixing top, EIS share 0', { regime: REGIME.COUPLED, inside: false }, short],
    ['surface-driven over ice of cover 0.6', { regime: REGIME.SURFACE, inside: true, iced: 0.6 }, short + 0.6 * (long - short)],
    ['coupled under a plume that topped above the layer', { regime: REGIME.COUPLED, inside: true, iced: 1, plumeTop: 0.5 * model.state[0][0] }, short],
    ['coupled with a plume that topped below the layer', { regime: REGIME.COUPLED, inside: true, plumeTop: model.state[0][0] }, long],
  ];
  for (const [label, column, lifetime] of cases) {
    const left = after(model, column);
    assert.ok(Math.abs(left - converts(lifetime)) < 1e-18, `${label}: ${left} against ${converts(lifetime)} over ${lifetime} s`);
    assert.ok(Math.abs(after(plain, column) - converts(short)) < 1e-18, `${label} without the option`);
  }
  console.log(`a layer at ${(model.state[0][0] * sigmaMid[k] / 100).toFixed(0)} hPa keeps ${(100 * converts(long) / 1e-4).toFixed(2)} % of its cloud water over ${dt} s as stratiform cloud, ${(100 * converts(short) / 1e-4).toFixed(2)} % otherwise`);
});

function plumeColumn(options = {}, profile = null, N = 16) {
  const model = build(options, N);
  if (profile) place(model, 0, 101500, profile); else jordanColumn(model, 0);
  setDepth(model, 0, 500);
  model.boundaryLayer.buoyancyFlux[0] = 4e-4;
  model.boundaryLayer.friction[0] = 0.25;
  model.radiation.mlmGate[0] = 0.3;
  return model;
}

test('the Jordan sounding lifts a deep plume that rains, heats most between 400 and 500 hPa above its cloud-base layer, cools the subcloud layer through its downdraft under the threshold closure, and keeps column enthalpy and water exact', () => {
  const model = plumeColumn({ capeClosure: 'threshold', plumeSourceDepth: 'boundaryLayer', plumeEntrainmentLaw: 'gregory' }), { moist } = model, [pi, theta, , , q, qc] = model.state;
  const { K, C, sigmaMid, geopotential, g } = model.core.diagnostics, dt = 600;
  moist.trace.convection = new Float64Array(K * C);
  const before = budget(model, 0);
  moist.adjust(model.state, 0, 1, dt);
  const after = budget(model, 0), { deep, falling } = moist, rain = moist.rain[0];
  const hPa = (k) => pi[0] * sigmaMid[k] / 100, rate = (k) => moist.trace.convection[k * C] / dt * 86400;
  let peak = deep.top;
  for (let k = deep.top; k < deep.base; k++) if (rate(k) > rate(peak)) peak = k;
  let lowest = 0, mass = 0;
  for (let k = 0; k < K; k++) if (geopotential[k * C] / g - geopotential[(K - 1) * C] / g < 100) { lowest += rate(k) * pi[0] * sigmaMid[k]; mass += pi[0] * sigmaMid[k]; }
  console.log(`Jordan (1958), plume: top ${(moist.cumulusTop[0] / 100).toFixed(0)} hPa, CAPE ${deep.cape.toFixed(0)} J/kg, base flux ${deep.baseFlux.toFixed(4)} kg/m²/s with a downdraft of ${deep.downdraft.toFixed(4)} from ${hPa(deep.start).toFixed(0)} hPa; rain ${(rain * 86400 / dt).toFixed(1)} mm/d after ${(falling.evaporated * 86400 / dt).toFixed(1)} evaporating below cloud base; heating above the cloud-base layer peaks at ${hPa(peak).toFixed(0)} hPa (${rate(peak).toFixed(1)} K/d), ${rate(deep.base).toFixed(1)} K/d in that layer at ${hPa(deep.base).toFixed(0)} hPa, ${(lowest / mass).toFixed(1)} K/d over the lowest 100 m; enthalpy off by ${((after.enthalpy - before.enthalpy) / before.enthalpy).toExponential(1)}, water by ${((after.water + rain - before.water) / before.water).toExponential(1)}`);
  assert.ok(deep.deep && moist.cumulusTop[0] < 300e2, `top ${moist.cumulusTop[0]}`);
  assert.ok(rain > 0 && moist.convectivePrecipitation[0] === rain && moist.largeScalePrecipitation[0] === 0, `rain ${rain}`);
  assert.ok(hPa(peak) > 400 && hPa(peak) < 500, `heating peaks at ${hPa(peak)} hPa`);
  assert.ok(deep.downdraft > 0 && deep.start > deep.top && deep.start < deep.base, 'a downdraft from between the top and cloud base');
  for (let k = deep.base + 1; k < K; k++) assert.ok(rate(k) < 0, `subcloud layer ${k} cools by ${rate(k)} K/d`);
  assert.ok(Math.abs(after.enthalpy - before.enthalpy) < 1e-12 * before.enthalpy, `enthalpy ${before.enthalpy} → ${after.enthalpy}`);
  assert.ok(Math.abs(after.water + rain - before.water) < 1e-12 * before.water, `water ${before.water} → ${after.water} + ${rain}`);
});

test('the trade-wind column\'s cumulus base flux at Grant\'s 0.03 is half that at 0.06, with column enthalpy and water exact', () => {
  const dt = 600, flux = (options) => {
    const model = tradeWindColumn(options), [pi, theta, , , q, qc] = model.state, before = budget(model, 0);
    model.moist.cumulusColumn(0, pi, theta, q, qc, dt);
    const after = budget(model, 0);
    return { base: model.moist.cumulusBaseFlux[0], enthalpy: (after.enthalpy - before.enthalpy) / before.enthalpy, water: (after.water - before.water) / before.water };
  };
  const grant = flux({}), double = flux({ cumulusClosure: 0.06 });
  console.log(`trade wind: base mass flux ${grant.base.toFixed(5)} kg/m²/s at c ${MOIST_DEFAULTS.cumulusClosure}, ${double.base.toFixed(5)} at 0.06; enthalpy ${grant.enthalpy.toExponential(1)}, water ${grant.water.toExponential(1)}`);
  assert.equal(MOIST_DEFAULTS.cumulusClosure, 0.03);
  assert.ok(Math.abs(grant.base - 0.5 * double.base) < 1e-15 * double.base, `${grant.base} against half of ${double.base}`);
  assert.ok(Math.abs(grant.enthalpy) < 1e-15 && Math.abs(grant.water) < 1e-15, 'enthalpy and water');
});

test('a trade-wind column lifts exactly the shallow cumulus plume and rains nothing', () => {
  const plume = tradeWindColumn(), shallow = tradeWindColumn(), dt = 600;
  const [pi, theta, , , q, qc] = plume.state, [sp, st, , , sq, sqc] = shallow.state;
  const rain = plume.moist.plumeColumn(0, pi, theta, q, qc, dt), shallowRain = shallow.moist.cumulusColumn(0, sp, st, sq, sqc, dt);
  assert.ok(!plume.moist.deep.deep && plume.moist.cumulusBaseFlux[0] > 0, 'a shallow plume');
  assert.equal(rain, 0);
  assert.equal(shallowRain, 0);
  assert.deepEqual(snapshot(plume, 0), snapshot(shallow, 0));
  assert.deepEqual(plume.moist.cumulusCover, shallow.moist.cumulusCover);
  assert.equal(plume.moist.cumulusBaseFlux[0], shallow.moist.cumulusBaseFlux[0]);
  assert.equal(plume.moist.cumulusTop[0], shallow.moist.cumulusTop[0]);
  const full = tradeWindColumn();
  full.moist.adjust(full.state, 0, 1, dt);
  assert.equal(full.moist.rain[0], 0, 'nothing rains through the adjustment');
});

test('a drier free troposphere entrains the plume to a lower top, and a dry enough one keeps it shallow', () => {
  const run = (factor) => {
    const model = plumeColumn({ capeClosure: 'threshold', plumeSourceDepth: 'boundaryLayer', plumeEntrainmentLaw: 'gregory', plumeConversion: 'zhangMcFarlane', plumeCape: 70 }, (z, p) => { const air = jordan(p); return { T: air.T, q: Math.min(air.q * (p < 850e2 ? factor : 1), saturationHumidity(air.T, p)) }; });
    model.moist.adjust(model.state, 0, 1, 600);
    return { top: model.moist.cumulusTop[0], deep: model.moist.deep.deep, cape: model.moist.deep.cape, rain: model.moist.rain[0] };
  };
  const moist = run(1), drier = run(0.8), dry = run(0.6), driest = run(0.4);
  console.log(`free-tropospheric humidity × 1, 0.8, 0.6, 0.4: plume tops ${[moist, drier, dry, driest].map((r) => (r.top / 100).toFixed(0)).join(', ')} hPa, CAPE ${[moist, drier, dry, driest].map((r) => r.cape.toFixed(0)).join(', ')} J/kg`);
  assert.ok(moist.deep && drier.deep && dry.deep && !driest.deep);
  assert.ok(moist.top < drier.top && drier.top < dry.top && dry.top < driest.top, 'the top sinks as the air dries');
  assert.ok(moist.rain > drier.rain && drier.rain > dry.rain && driest.rain === 0);
});

test('under the threshold closure the deep base flux relaxes the CAPE toward plumeCape over plumeRelaxation, and plumeClosure: maximum gives it at least the shallow closure', () => {
  const flux = (options) => { const model = plumeColumn({ capeClosure: 'threshold', plumeSourceDepth: 'boundaryLayer', plumeEntrainmentLaw: 'gregory', ...options }), [pi, theta, , , q, qc] = model.state; model.moist.plumeColumn(0, pi, theta, q, qc, 600); return { ...model.moist.deep }; };
  const hour = flux({}), twoHours = flux({ plumeRelaxation: 7200 }), lower = flux({ plumeCape: 0 });
  const cape0 = MOIST_DEFAULTS.plumeCape;
  assert.ok(Math.abs(hour.baseFlux - (hour.cape - cape0) / (3600 * hour.consumption)) < 1e-12 * hour.baseFlux, `base flux ${hour.baseFlux}`);
  assert.ok(Math.abs(twoHours.baseFlux - 0.5 * hour.baseFlux) < 1e-12 * hour.baseFlux, 'twice the relaxation time, half the flux');
  assert.ok(Math.abs(lower.baseFlux / hour.baseFlux - hour.cape / (hour.cape - cape0)) < 1e-9, 'the flux scales with the CAPE above plumeCape');
  const shallowOnly = flux({ plumeCape: 1e6 }), maximum = flux({ plumeCape: 1e6, plumeClosure: 'maximum' });
  assert.ok(!shallowOnly.deep && maximum.deep && maximum.baseFlux > 0, 'with CAPE below plumeCape the deep plume runs only under plumeClosure: maximum');
});

function subcloudVirtual(model, i, heating, dt) {
  const { K, C, exnerLayer } = model.core.diagnostics, [, theta, , , q, qc] = model.state, KL = model.moist.subcloudLayers;
  for (let k = K - KL; k < K; k++) {
    const idx = k * C + i;
    model.moist.subcloudVirtual[(k - K + KL) * C + i] = theta[idx] * exnerLayer[idx] * (1 + VIRTUAL_FACTOR * q[idx] - qc[idx]) - heating / 86400 * dt;
  }
}
function closureByHand(model, i) {
  const { moist } = model, { deep } = moist, { K, C, cp, g, R, sigmaMid, dSigma, geopotential, exnerLayer, exnerLower } = model.core.diagnostics, thetaV = model.core.arrays.thetaV, [pi, theta, , , q, qc] = model.state;
  const upper = (k) => (geopotential[k * C + i] + cp * thetaV[k * C + i] * (exnerLayer[k * C + i] - exnerLower[(k - 1) * C + i])) / g;
  let pcape = 0, weighted = 0, thickness = 0;
  for (let k = 0; k < K; k++) if (moist.plumeCounted[k]) pcape += moist.plumeBuoyancy[k] / g * pi[i] * dSigma[k];
  for (let k = deep.top; k < deep.base; k++) { const depth = upper(k) - upper(k + 1); weighted += depth * Math.sqrt(0.5 * (moist.plumeSpeed[k] + moist.plumeSpeed[k + 1])); thickness += depth; }
  const H = upper(deep.top) - upper(deep.base), speed = weighted / thickness, scale = 1 + BECHTOLD.resolution * Math.sqrt(model.mesh.areaCell[i]) / BECHTOLD.reference;
  const b = (K - 1) * C + i, ground = (geopotential[b] - cp * thetaV[b] * (exnerLower[b] - exnerLayer[b])) / g;
  return { pcape, H, speed, scale, tau: Math.min(BECHTOLD.longest, Math.max(BECHTOLD.shortest, scale * H / speed)), baseHeight: upper(deep.base) - ground };
}

test('under the Bechtold closure the Jordan plume removes its PCAPE over α_x H / w̄, the base flux scales as 1/α_x with the cell size, and the column keeps its enthalpy and water', () => {
  const dt = 600, runs = {};
  for (const N of [32, 64]) {
    const model = plumeColumn({}, null, N), { moist } = model, [pi, theta, , , q, qc] = model.state;
    const before = budget(model, 0);
    moist.adjust(model.state, 0, 1, dt);
    const after = budget(model, 0), rain = moist.rain[0], deep = { ...moist.deep }, hand = closureByHand(model, 0);
    const flux = Math.max(0, deep.pcape - deep.pcapeBoundary) / (deep.tau * deep.consumptionP);
    runs[N] = { deep, hand, flux };
    console.log(`N=${N} (dx ${(Math.sqrt(model.mesh.areaCell[0]) / 1e3).toFixed(0)} km, α_x ${hand.scale.toFixed(3)}): Jordan's plume tops at ${(moist.cumulusTop[0] / 100).toFixed(0)} hPa, CAPE ${deep.cape.toFixed(0)} J/kg, PCAPE ${deep.pcape.toFixed(2)} Pa (by hand ${hand.pcape.toFixed(2)}), F_P ${deep.consumptionP.toExponential(4)}, H ${deep.depth.toFixed(0)} m, w̄ ${deep.speed.toFixed(3)} m/s, τ ${deep.tau.toFixed(0)} s, base flux ${deep.baseFlux.toFixed(5)} kg/m²/s, rain ${(rain * 86400 / dt).toFixed(1)} mm/d; enthalpy off by ${((after.enthalpy - before.enthalpy) / before.enthalpy).toExponential(1)}, water by ${((after.water + rain - before.water) / before.water).toExponential(1)}`);
    assert.ok(deep.deep && deep.pcapeBoundary === 0, 'a deep plume with no saved state, so no boundary-layer forcing');
    assert.ok(Math.abs(deep.pcape - hand.pcape) < 1e-12 * hand.pcape, `PCAPE ${deep.pcape} against ${hand.pcape}`);
    assert.ok(Math.abs(deep.depth - hand.H) < 1e-9 * hand.H && Math.abs(deep.speed - hand.speed) < 1e-12 * hand.speed, 'H and w̄');
    assert.ok(Math.abs(deep.tau - hand.tau) < 1e-12 * hand.tau && hand.tau > BECHTOLD.shortest && hand.tau < BECHTOLD.longest, `τ ${deep.tau} against ${hand.tau}, within the bounds`);
    assert.ok(Math.abs(deep.baseFlux - flux) < 1e-12 * flux, `base flux ${deep.baseFlux} against ${flux}`);
    assert.ok(Math.abs(after.enthalpy - before.enthalpy) < 1e-15 * before.enthalpy && Math.abs(after.water + rain - before.water) < 1e-15 * before.water, 'enthalpy and water');
  }
  const ratio = runs[64].deep.baseFlux / runs[32].deep.baseFlux, expected = runs[32].hand.scale / runs[64].hand.scale;
  assert.ok(Math.abs(ratio / expected - 1) < 1e-9, `halving dx scales the flux by ${ratio}, α_x(dx)/α_x(dx/2) ${expected}`);
});

test('under the Bechtold closure a warming subcloud layer over the sea takes its PCAPE_bl over the advective time z_base / ū_bl, and over land over the turnover time, which holds a plume under strong heating back unless the PCAPE minimum keeps it firing', () => {
  const dt = 600, model = plumeColumn({}, null, 16), { moist, mesh } = model, [pi, theta, u, , q, qc] = model.state;
  const { K, C, dSigma } = model.core.diagnostics, E = mesh.nEdges, lon = mesh.lonCell[0], east = [-Math.sin(lon), Math.cos(lon), 0];
  for (let k = 0; k < K; k++) for (let e = 0; e < E; e++) u[k * E + e] = 5 * (east[0] * mesh.nEdge[3 * e] + east[1] * mesh.nEdge[3 * e + 1] + east[2] * mesh.nEdge[3 * e + 2]);
  subcloudVirtual(model, 0, 2, dt);
  moist.plumeColumn(0, pi, theta, q, qc, dt, u);
  const deep = { ...moist.deep }, hand = closureByHand(model, 0);
  let below = 0;
  for (let k = Math.max(deep.base, K - SUBCLOUD_LAYERS); k < K; k++) below += pi[0] * dSigma[k];
  const expected = hand.baseHeight / 5 * 2 / 86400 * below;
  console.log(`sea: subcloud +2 K/d over ${(below / 100).toFixed(0)} hPa, ū_bl ${deep.boundaryWind.toFixed(3)} m/s (5 by construction), z_base ${hand.baseHeight.toFixed(0)} m: PCAPE_bl ${deep.pcapeBoundary.toFixed(3)} Pa against ${expected.toFixed(3)} by hand, PCAPE ${deep.pcape.toFixed(1)} Pa`);
  assert.ok(Math.abs(deep.pcapeBoundary / expected - 1) < 0.01, `PCAPE_bl ${deep.pcapeBoundary} against ${expected}`);
  const fire = (heating, options = PREVIOUS_PCAPE_MINIMUM) => {
    const land = plumeColumn({ ...options, land: Uint8Array.from({ length: 16 * 16 * 10 + 2 }, (_, i) => (i === 0 ? 1 : 0)) }, null, 16), [lp, lt, , , lq, lqc] = land.state;
    subcloudVirtual(land, 0, heating, dt);
    land.moist.plumeColumn(0, lp, lt, lq, lqc, dt);
    return { ...land.moist.deep };
  };
  const heated = fire(10), calm = fire(0), floored = fire(10, {});
  console.log(`land: under +10 K/d below cloud base τ_bl ${heated.boundaryTime.toFixed(0)} s (H / w̄) gives PCAPE_bl ${heated.pcapeBoundary.toFixed(0)} Pa against PCAPE ${heated.pcape.toFixed(0)}: deep ${heated.deep}; with no heating base flux ${calm.baseFlux.toFixed(4)} kg/m²/s`);
  assert.ok(!heated.deep && heated.pcapeBoundary > heated.pcape && Math.abs(heated.boundaryTime - heated.depth / heated.speed) < 1e-9 * heated.boundaryTime, 'strong heating holds the land plume back over the turnover time');
  assert.ok(calm.deep && calm.baseFlux > 0 && calm.pcapeBoundary === 0, 'the calm land column fires');
  assert.ok(floored.deep && floored.baseFlux > 0 && floored.pcapeBoundary > floored.pcape, 'under the default PCAPE minimum the heated land plume fires');
});

test('a cooling subcloud layer gives no boundary-layer part under pcapeBoundary positive, so the land plume at night removes its PCAPE over τ as with no tendency; signed lets the cooling add to the PCAPE', () => {
  const dt = 600, fire = (heating, options = {}) => {
    const land = plumeColumn({ land: Uint8Array.from({ length: 16 * 16 * 10 + 2 }, (_, i) => (i === 0 ? 1 : 0)), ...options }, null, 16), [lp, lt, , , lq, lqc] = land.state;
    subcloudVirtual(land, 0, heating, dt);
    land.moist.plumeColumn(0, lp, lt, lq, lqc, dt);
    return { ...land.moist.deep };
  };
  const calm = fire(0), night = fire(-10), signed = fire(-10, { pcapeBoundary: 'signed' }), day = fire(10), signedDay = fire(10, { pcapeBoundary: 'signed' });
  console.log(`land under −10 K/d below cloud base: PCAPE_bl ${night.pcapeBoundary} Pa and base flux ${night.baseFlux.toFixed(5)} kg/m²/s (${calm.baseFlux.toFixed(5)} with no tendency); signed PCAPE_bl ${signed.pcapeBoundary.toFixed(0)} Pa against PCAPE ${signed.pcape.toFixed(0)}, base flux ${signed.baseFlux.toFixed(5)}`);
  assert.equal(MOIST_DEFAULTS.pcapeBoundary, 'positive');
  assert.ok(night.deep && night.pcapeBoundary === 0 && night.baseFlux === calm.baseFlux, 'a cooling subcloud layer leaves the closure as with no tendency');
  assert.ok(signed.pcapeBoundary < -signed.pcape && signed.baseFlux > 2 * calm.baseFlux, 'signed, the cooling adds to the PCAPE');
  assert.deepEqual(day, signedDay, 'a warming subcloud layer is the same under both');
});

test('the deep plume leaves the lowest 50 hPa with the IFS surface-flux excess, the plain 50 hPa mean when the fluxes vanish, and the shallow plume keeps the mixed layer', () => {
  const dt = 600, run = (fluxes, options = {}) => {
    const model = plumeColumn(options), { moist } = model, [pi, theta, , , q, qc] = model.state;
    if (fluxes) { model.radiation.sensibleHeat[0] = fluxes[0]; model.radiation.evaporation[0] = fluxes[1] / LATENT_HEAT; }
    const before = budget(model, 0);
    moist.adjust(model.state, 0, 1, dt);
    const after = budget(model, 0), rain = moist.rain[0];
    return { model, deep: { ...moist.deep }, state: snapshot(model, 0), enthalpy: (after.enthalpy - before.enthalpy) / before.enthalpy, water: (after.water + rain - before.water) / before.water };
  };
  const plain = run(null), still = run([0, 0]), forced = run([10, 130]), boundary = run([10, 130], { plumeSourceDepth: 'boundaryLayer' }), convective = run([10, 130], { excessVelocity: 'convective' });
  const night = run([-30, 1e-6 * LATENT_HEAT]), dew = run([5, -5e-6 * LATENT_HEAT]);
  const { K, C, sigmaMid, dSigma, cp, g, R, geopotential, exnerLayer, exnerLower } = forced.model.core.diagnostics;
  const model = plumeColumn(), [pi, theta, , , q, qc] = model.state;
  let mass = 0, energy = 0, water = 0;
  for (let k = K - 1; k >= 0 && (k === K - 1 || pi[0] * sigmaMid[k] >= pi[0] - MOIST_DEFAULTS.cumulusSourceDepth); k--) {
    const idx = k * C, dp = pi[0] * dSigma[k], T = theta[idx] * exnerLayer[idx];
    mass += dp; energy += dp * (cp * T + geopotential[idx] - LATENT_HEAT * qc[idx]); water += dp * (q[idx] + qc[idx]);
  }
  const b = (K - 1) * C, T1 = theta[b] * exnerLayer[b], density = pi[0] * sigmaMid[K - 1] / (R * T1), z1 = cp * model.core.arrays.thetaV[b] * (exnerLower[b] - exnerLayer[b]) / g;
  const velocity = 1.2 * Math.cbrt(0.1 ** 3 + 1.5 * g * z1 * 0.4 / T1 * (10 / (density * cp) + 0.61 * T1 * 130 / (density * LATENT_HEAT))), dT = Math.min(3, 1.5 * 10 / (density * cp * velocity)), dq = Math.min(2e-3, 1.5 * 130 / (density * LATENT_HEAT * velocity));
  const mixed = Math.max(Math.cbrt(4e-4 * (model.boundaryLayer.depth[0] - geopotential[b] / g)), 0.25), mixedT = Math.min(3, 1.5 * 10 / (density * cp * mixed)), mixedQ = Math.min(2e-3, 1.5 * 130 / (density * LATENT_HEAT * mixed));
  console.log(`Jordan, the lowest 50 hPa (${(mass / 100).toFixed(1)} hPa of layers) under 10 W/m² sensible and 130 latent, w* of eq. 6.20 at the lowest layer's ${z1.toFixed(1)} m ${velocity.toFixed(3)} m/s: excess ${forced.deep.excessT.toFixed(4)} K and ${(1e3 * forced.deep.excessQ).toFixed(4)} g/kg (by hand ${dT.toFixed(4)}, ${(1e3 * dq).toFixed(4)}; with the shallow closure's w* ${mixed.toFixed(3)} m/s ${convective.deep.excessT.toFixed(4)} K and ${(1e3 * convective.deep.excessQ).toFixed(4)} g/kg); CAPE ${plain.deep.cape.toFixed(1)} J/kg plain, ${forced.deep.cape.toFixed(1)} with the excess (${convective.deep.cape.toFixed(1)} with the shallow closure's w*), ${boundary.deep.cape.toFixed(1)} from the boundary layer with it; enthalpy ${forced.enthalpy.toExponential(1)}, water ${forced.water.toExponential(1)}`);
  assert.deepEqual(still.state, plain.state, 'with no surface fluxes the plume leaves with the plain 50 hPa mean');
  assert.equal(still.deep.cape, plain.deep.cape);
  assert.ok(Math.abs(plain.deep.sourceS - energy / mass) < 1e-12 * energy / mass && Math.abs(plain.deep.sourceQ - water / mass) < 1e-15 && plain.deep.sourceMass === mass, 'the source is the lowest 50 hPa');
  assert.ok(Math.abs(forced.deep.excessT - dT) < 1e-12 * dT && Math.abs(forced.deep.excessQ - dq) < 1e-12 * dq, 'the excess of IFS eqs. 6.19-6.20');
  assert.ok(Math.abs(convective.deep.excessT - mixedT) < 1e-12 * mixedT && Math.abs(convective.deep.excessQ - mixedQ) < 1e-12 * mixedQ, 'with the shallow closure\'s w*');
  assert.ok(night.deep.excessT === 0 && night.deep.excessQ === 0, 'no excess under a downward buoyancy flux');
  assert.ok(dew.deep.excessT > 0 && dew.deep.excessQ === 0, 'no negative part under dew with an upward buoyancy flux');
  assert.ok(Math.abs(forced.deep.sourceS - (energy / mass + cp * dT)) < 1e-12 * energy / mass && Math.abs(forced.deep.sourceQ - (water / mass + dq)) < 1e-15, 'the source with its excess');
  assert.ok(forced.deep.cape > plain.deep.cape, 'the excess raises the CAPE');
  assert.ok(Math.abs(forced.enthalpy) < 1e-15 && Math.abs(forced.water) < 1e-15, 'enthalpy and water');
});

/*
 * A nocturnal inversion under a moist residual layer that is conditionally
 * unstable: from 283 K at the ground the air warms by `strength` K to the cap
 * at `cap` m (80 % humid), above which it cools at 7 K/km to a 205 K floor
 * (90 % humid), up to `warm` K warmer throughout.
 */
const nocturnalSounding = (cap, strength, warm = 0) => (z, p) => {
  const T = warm + (z < cap ? 283 + strength * z / cap : Math.max(205, 283 + strength - 7e-3 * (z - cap)));
  return { T, q: (z < cap ? 0.8 : 0.9) * saturationHumidity(T, p) };
};

/*
 * Random columns of five kinds (Jordan's sounding warmed or cooled, trade
 * cumulus, stratocumulus, trade wind, a cold dry Jordan) with random
 * boundary layers, deck gates, regimes, mixing tops, sea ice, saved
 * subcloud T_v and surface fluxes, all in single precision.
 */
function randomColumns(model, { nocturnal = false } = {}) {
  const { moist, core, mesh } = model, C = mesh.nCells, { K } = core.diagnostics;
  const [pi, theta, , surfaceT, q, qc] = model.state;
  let seed = 12345;
  const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  const kinds = new Int32Array(C);
  const buoyancy = new Float64Array(C), friction = new Float64Array(C);
  let fogSeed = 99;
  const fogRandom = () => { fogSeed = (fogSeed * 1103515245 + 12345) % 2147483648; return fogSeed / 2147483648; };
  let nightSeed = 2468;
  const nightRandom = () => { nightSeed = (nightSeed * 1103515245 + 12345) % 2147483648; return nightSeed / 2147483648; };
  for (let i = 0; i < C; i++) {
    const kind = Math.floor(random() * 5), warm = 3 * (random() - 0.5), wet = 1 + 0.15 * (random() - 0.3), surface = 100500 + 1500 * random();
    kinds[i] = kind;
    buoyancy[i] = random() < 0.15 ? -2e-4 * random() : 6e-4 * random();
    friction[i] = 0.1 + 0.3 * random();
    const night = nocturnal ? { cooled: nightRandom() < 0.4, cap: 300 + 600 * nightRandom(), strength: 4 + 6 * nightRandom(), flux: -2e-5 * nightRandom() } : null;
    if (kind === 0 && night && night.cooled) {
      place(model, i, surface, nocturnalSounding(night.cap, night.strength, warm));
      buoyancy[i] = night.flux;
    } else if (kind === 0) place(model, i, surface, (z, p) => { const air = jordan(p); const T = air.T + warm * (z < 800 ? 1 : 0.5); return { T, q: Math.min(wet * air.q, 0.95 * saturationHumidity(T, p)) }; });
    else if (kind === 1) place(model, i, surface, (z, p) => { const air = tradeCumulus(z, p); return { T: air.T + warm, q: Math.min(wet * air.q, 0.95 * saturationHumidity(air.T + warm, p)) }; });
    else if (kind === 2) place(model, i, surface, stratocumulus);
    else if (kind === 4) place(model, i, surface, (z, p) => { const air = tradeWind(z, p); return { T: air.T + warm, q: Math.min(wet * air.q, 0.95 * saturationHumidity(air.T + warm, p)) }; });
    else place(model, i, surface, (z, p) => { const air = jordan(p); const T = air.T - 5 + warm; return { T, q: Math.min(0.6 * air.q, 0.95 * saturationHumidity(T, p)) }; });
    surfaceT[i] = 300;
    const { geopotential, g } = core.diagnostics;
    model.boundaryLayer.depth[i] = geopotential[(K - 1) * C + i] / g + 200 + 1300 * random();
    const gate = random();
    model.radiation.mlmGate[i] = gate < 0.2 ? 0.7 : gate < 0.35 ? 0.5 + (gate - 0.2) / 0.15 * (DECK_CLOSED - 0.5) : 0.3;
    for (let k = 0; k < K; k++) if (random() < 0.08) qc[k * C + i] = 1e-3 * random();
    if (fogRandom() < 0.3) {
      for (const k of [K - 2, K - 1]) {
        const x = k * C + i;
        q[x] = saturationHumidity(theta[x] * core.diagnostics.exnerLayer[x], pi[i] * core.diagnostics.sigmaMid[k]);
        qc[x] = 1e-5 * Math.pow(2e3, fogRandom());
      }
    }
  }
  for (let x = 0; x < model.state[2].length; x++) model.state[2][x] = 20 * (random() - 0.5) + 10 * Math.sin(x / mesh.nEdges);
  let layerSeed = 777;
  const layerRandom = () => { layerSeed = (layerSeed * 1103515245 + 12345) % 2147483648; return layerSeed / 2147483648; };
  const { regime, mixingTop } = model.boundaryLayer, { stratiform } = model.radiation, { concentration } = model.seaIce, ice = model.state[6];
  for (let i = 0; i < C; i++) {
    const { geopotential, g } = core.diagnostics;
    regime[i] = Math.floor(4 * layerRandom());
    mixingTop[i] = Math.fround(geopotential[(K - 1) * C + i] / g + 300 + 2200 * layerRandom());
    stratiform[i] = Math.fround(layerRandom() < 0.5 ? 0 : layerRandom());
    if (layerRandom() < 0.15) { ice[i] = 0.5 + layerRandom(); concentration[i] = Math.fround(layerRandom() < 0.3 ? 0 : layerRandom()); }
  }
  for (const a of model.state) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  const KL = moist.subcloudLayers, saved = moist.subcloudVirtual, { exnerLayer } = core.diagnostics;
  for (let i = 0; i < C; i++) {
    core.diagnoseColumn(i, pi, theta, q, qc);
    const unset = layerRandom() < 0.1;
    for (let k = K - KL; k < K; k++) { const x = k * C + i; saved[(k - K + KL) * C + i] = unset ? 0 : Math.fround(theta[x] * exnerLayer[x] * (1 + VIRTUAL_FACTOR * q[x] - qc[x]) + 0.1 * (layerRandom() - 0.5)); }
  }
  for (const a of [model.boundaryLayer.depth, model.radiation.mlmGate, buoyancy, friction]) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  model.boundaryLayer.buoyancyFlux.set(buoyancy);
  model.boundaryLayer.friction.set(friction);
  const sensible = Float64Array.from({ length: C }, () => Math.fround(60 * layerRandom() - 10)), evaporation = Float64Array.from({ length: C }, () => Math.fround(6e-5 * layerRandom()));
  model.radiation.sensibleHeat.set(sensible); model.radiation.evaporation.set(evaporation);
  return { buoyancy, friction, regime, mixingTop, stratiform, concentration, saved, sensible, evaporation };
}

test('the IFS entrainment: on Jordan\'s column the plume entrains 1.75e-3 (1.3 − RH)(q_s/q_s,base)³ where the layer below is buoyant, detrains 0.75e-4 (1.6 − RH), and its mass flux grows by exp((ε − δ)Δz) while buoyant and falls by the organised detrainment above, with column enthalpy and water exact; a drier free troposphere entrains more and tops lower', () => {
  const dt = 600, model = plumeColumn({}), { moist } = model, [pi, theta, , , q, qc] = model.state;
  const { K, C, sigmaMid } = model.core.diagnostics, O = MOIST_DEFAULTS, sat = { qs: 0, slope: 0, liquid: 1 };
  const before = budget(model, 0);
  moist.adjust(model.state, 0, 1, dt);
  const after = budget(model, 0), rain = moist.rain[0], d = moist.deep;
  const fresh = plumeColumn({}), [, ft, , , fq] = fresh.state, { exnerLayer } = fresh.core.diagnostics;
  fresh.core.diagnoseColumn(0, fresh.state[0], ft, fq, fresh.state[5]);
  fresh.moist.condenseColumn(0, fresh.state[0], ft, fq, fresh.state[5]);
  fresh.core.diagnoseColumn(0, fresh.state[0], ft, fq, fresh.state[5]);
  const held = Float64Array.from({ length: K }, (_, k) => ft[k * C] * exnerLayer[k * C]), vapour = Float64Array.from({ length: K }, (_, k) => fq[k * C]);
  const T = (k) => held[k], pk = (k) => fresh.state[0][0] * sigmaMid[k];
  fresh.moist.plumeColumn(0, fresh.state[0], ft, fq, fresh.state[5], dt);
  const f = fresh.moist, fd = f.deep, qsBase = cloudSaturation(T(fd.base - 1), pk(fd.base - 1), O.iceSaturation, O.liquidTemperature, O.iceTemperature, sat).qs;
  let worstE = 0, worstD = 0, worstM = 0, growing = 0, falling = 0;
  for (let k = fd.base - 1; k > fd.top; k--) {
    const qs = cloudSaturation(T(k), pk(k), O.iceSaturation, O.liquidTemperature, O.iceTemperature, sat).qs, rh = Math.min(1, vapour[k] / qs);
    const eps = f.plumeBuoyancy[k + 1] > 0 ? IFS_ENTRAINMENT.entrainment * (IFS_ENTRAINMENT.humidity - rh) * (qs / qsBase) ** 3 : 0;
    const expected = k === fd.base - 1 ? f.plumeEntrained[k] : eps;
    worstE = Math.max(worstE, Math.abs(f.plumeEntrained[k] - expected) / Math.max(1e-12, expected));
    worstD = Math.max(worstD, Math.abs(f.plumeDetrained[k] - IFS_ENTRAINMENT.detrainment * (IFS_ENTRAINMENT.detrainmentHumidity - rh)));
  }
  const { geopotential, g, cp: heat, exnerLower } = fresh.core.diagnostics, thetaV = fresh.core.arrays.thetaV;
  const upper = (k) => (geopotential[k * C] + heat * thetaV[k * C] * (exnerLayer[k * C] - exnerLower[(k - 1) * C])) / g;
  for (let k = fd.base - 1; k > fd.top; k--) {
    const dz = upper(k) - upper(k + 1), M = f.cumulusFlux, ratio = f.plumeBuoyancy[k] > 0 ? Math.exp((f.plumeEntrained[k] - f.plumeDetrained[k]) * dz) : Math.exp(-f.plumeDetrained[k] * dz) * Math.min(1, (IFS_ENTRAINMENT.detrainmentHumidity - f.envHumidity[k]) * Math.sqrt(f.plumeSpeed[k] / f.plumeSpeed[k + 1]));
    worstM = Math.max(worstM, Math.abs(M[k] - M[k + 1] * ratio) / M[k + 1]);
    if (M[k] > M[k + 1]) growing++; else falling++;
  }
  console.log(`Jordan under the IFS entrainment: mass flux off its recurrence by ${worstM.toExponential(1)}, growing through ${growing} layers and falling through ${falling}; top ${(moist.cumulusTop[0] / 100).toFixed(0)} hPa, cloud base ${(pi[0] * model.core.diagnostics.levels[d.base] / 100).toFixed(0)} hPa, CAPE ${d.cape.toFixed(0)} J/kg, base flux ${d.baseFlux.toFixed(4)} kg/m²/s, rain ${(rain * 86400 / dt).toFixed(1)} mm/d; ε at the first cloud layers ${Array.from({ length: 4 }, (_, n) => (1e4 * f.plumeEntrained[fd.base - 1 - n]).toFixed(2)).join(', ')} 1e-4/m; ε off its formula by ${worstE.toExponential(1)}, δ by ${worstD.toExponential(1)}; enthalpy ${((after.enthalpy - before.enthalpy) / before.enthalpy).toExponential(1)}, water ${((after.water + rain - before.water) / before.water).toExponential(1)}`);
  assert.equal(O.plumeEntrainmentLaw, 'ifs');
  assert.ok(d.deep && d.baseFlux > 0, 'a deep plume');
  assert.ok(worstE < 1e-12 && worstD < 1e-18 && worstM < 1e-9, `ε ${worstE}, δ ${worstD}, M ${worstM}`);
  assert.ok(growing > 3 && falling > 0, `${growing} growing, ${falling} falling`);
  assert.ok(Math.abs(after.enthalpy - before.enthalpy) < 1e-15 * before.enthalpy && Math.abs(after.water + rain - before.water) < 1e-15 * before.water, 'enthalpy and water');
  const run = (factor) => {
    const m = plumeColumn({}, (z, p) => { const air = jordan(p); return { T: air.T, q: Math.min(air.q * (p < 850e2 ? factor : 1), saturationHumidity(air.T, p)) }; });
    m.moist.adjust(m.state, 0, 1, 600);
    return { top: m.moist.deep.top >= 0 ? m.moist.cumulusTop[0] : Infinity, entrained: m.moist.plumeEntrained.reduce((a, b) => a + b, 0), deep: m.moist.deep.deep };
  };
  const moistAir = run(1), drier = run(0.8), dry = run(0.6);
  console.log(`free-tropospheric humidity × 1, 0.8, 0.6: tops ${[moistAir, drier, dry].map((r) => (r.top / 100).toFixed(0)).join(', ')} hPa, summed ε ${[moistAir, drier, dry].map((r) => (1e4 * r.entrained).toFixed(1)).join(', ')} 1e-4/m`);
  assert.ok(moistAir.entrained < drier.entrained && drier.entrained < dry.entrained, 'drier air entrains more');
  assert.ok(moistAir.top <= drier.top && drier.top <= dry.top, 'and the top sinks');
});

test('a column convects deep where its plume\'s cloud is deeper than 200 hPa and then as one type only, the deep plume alone, as convectionType top with plumeClosure cape does where both call it deep; elsewhere the shallow plume alone; column enthalpy and water exact', () => {
  const dt = 900, run = (options) => {
    const model = build(options), C = model.mesh.nCells;
    randomColumns(model);
    const [pi, theta, , , q, qc] = model.state, { moist, core } = model, { levels } = core.diagnostics;
    const out = [];
    for (let i = 0; i < C; i++) {
      core.diagnoseColumn(i, pi, theta, q, qc);
      const before = budget(model, i);
      const rain = moist.plumeColumn(i, pi, theta, q, qc, dt);
      const after = budget(model, i), d = moist.deep;
      out.push({ state: snapshot(model, i), deep: d.deep, top: d.top, base: d.base, depth: d.top >= 0 ? pi[i] * (levels[d.base] - levels[d.top]) : 0, shallowTop: d.top >= 0 && levels[d.top] * 1e5 < MOIST_DEFAULTS.shallowTop, flux: moist.cumulusBaseFlux[i], deepFlux: d.baseFlux, enthalpy: Math.abs(after.enthalpy - before.enthalpy) / before.enthalpy, water: Math.abs(after.water + rain - before.water) / before.water });
    }
    return out;
  };
  const depth = run({ convectionType: 'cloudDepth' }), topCape = run({ convectionType: 'top', plumeClosure: 'cape' }), topSeparate = run({ convectionType: 'top' });
  let deep = 0, same = 0, onlyDepth = 0, onlyTop = 0, alone = 0, both = 0, worst = 0;
  depth.forEach((c, i) => {
    worst = Math.max(worst, c.enthalpy, c.water);
    if (c.deep) {
      deep++;
      assert.ok(c.depth > DEEP_CLOUD_DEPTH, `column ${i}: deep with a cloud of ${c.depth} Pa`);
      assert.equal(c.flux, c.deepFlux, `column ${i}: the deep plume runs alone`);
      if (c.deepFlux > 0) alone++;
    }
    const t = topCape[i];
    if (c.top >= 0 && c.depth > DEEP_CLOUD_DEPTH && t.shallowTop) { same++; assert.deepEqual(c.state, t.state, `column ${i}: deep by both rules, as top with cape`); }
    if (c.top >= 0 && c.depth > DEEP_CLOUD_DEPTH && !t.shallowTop) onlyDepth++;
    if (c.top >= 0 && !(c.depth > DEEP_CLOUD_DEPTH) && t.shallowTop) { onlyTop++; assert.ok(!c.deep, `column ${i}: shallow by the cloud's depth`); }
    if (topSeparate[i].deep && topSeparate[i].flux > topSeparate[i].deepFlux) both++;
  });
  console.log(`${depth.length} random columns: ${deep} deep by the cloud's depth (${alone} with a flux, all alone), ${same} deep by both rules and as top with cape, ${onlyDepth} deep by the depth only, ${onlyTop} by the 700 hPa top only; under top with separate ${both} columns run both plumes; enthalpy and water within ${worst.toExponential(1)}`);
  assert.ok(deep > depth.length / 30 && same > depth.length / 30 && alone > depth.length / 30 && both > depth.length / 100, `${deep} deep, ${same} alike, ${alone} alone, ${both} with both under top`);
  assert.ok(worst < 1e-15, `enthalpy and water ${worst}`);
});

/*
 * Plume air of liquid-ice static energy `sl` = c_p T + g z − L q_l −
 * (L + L_f) q_i and total water `qt` by bisection: T, its condensate l =
 * q_t − q_s,mix(T) of cloudSaturation's mix and the ice share (1 − α(T)) l.
 */
function mixedState(sl, qt, z, pressure, O, cp, g) {
  const sat = { qs: 0, slope: 0, liquid: 1 };
  const held = (T) => qt - cloudSaturation(T, pressure, true, O.liquidTemperature, O.iceTemperature, sat).qs;
  const dry = (sl - g * z) / cp;
  if (!(held(dry) > 0)) return { T: dry, l: 0, ice: 0 };
  let lo = dry, hi = dry + 100;
  for (let n = 0; n < 200; n++) {
    const T = 0.5 * (lo + hi), l = Math.max(0, held(T)), a = liquidFraction(T, O.liquidTemperature, O.iceTemperature);
    if (cp * T + g * z - LATENT_HEAT * l - FUSION_HEAT * (1 - a) * l > sl) hi = T; else lo = T;
  }
  const T = 0.5 * (lo + hi), l = Math.max(0, held(T));
  return { T, l, ice: (1 - liquidFraction(T, O.liquidTemperature, O.iceTemperature)) * l };
}

/*
 * A cold column over a 265 K surface: 9 K/km to 205 K, 90 % humid in the
 * lowest 500 m and 85 % above.
 */
function coldColumn(z, p) {
  const T = Math.max(205, 265 - 9e-3 * z);
  return { T, q: (z < 500 ? 0.9 : 0.85) * saturationHumidity(T, p) };
}

/*
 * The IFS's first-guess deep updraught by hand (Cy43r1 §6.4, eqs 6.18–6.21)
 * on column 0 as plumeColumn sees it: from the deep source's s_l and q_t at
 * its top interface at 1 m/s, mixing at 0.4 · 1.75e-3 (q_s/q_s,lowest)³,
 * half the condensate removed at each upper interface, w² by the IFS form,
 * its zero inside the top layer found by bisection and the pressure there
 * interpolated in ln p; with the interface at which its cloud first exceeds
 * 200 hPa.
 */
function handParcel(model, departure = -1) {
  const { K, C, levels, sigmaMid, geopotential, g, cp, R, exnerLayer, exnerLower } = model.core.diagnostics, thetaV = model.core.arrays.thetaV;
  const [pi, theta, , , q, qc] = model.state, O = MOIST_DEFAULTS, sat = { qs: 0, slope: 0, liquid: 1 }, L = LATENT_HEAT;
  const qsOf = (T, p) => cloudSaturation(T, p, O.iceSaturation, O.liquidTemperature, O.iceTemperature, sat).qs;
  const T = (k) => theta[k * C] * exnerLayer[k * C], p = (k) => pi[0] * sigmaMid[k], zMid = (k) => geopotential[k * C] / g;
  const zUp = (k) => (geopotential[k * C] + cp * thetaV[k * C] * (exnerLayer[k * C] - exnerLower[(k - 1) * C])) / g;
  const envS = (k) => cp * T(k) + g * zMid(k) - L * Math.max(0, qc[k * C]), envQ = (k) => Math.max(0, q[k * C]) + Math.max(0, qc[k * C]);
  const state = (sl, qt, z, pressure) => mixedState(sl, qt, z, pressure, O, cp, g);
  const bottom = K - 1;
  let source = bottom;
  while (source > 0 && p(source - 1) >= pi[0] - O.cumulusSourceDepth) source--;
  let sl, qt;
  if (departure >= 0) { source = departure; sl = envS(departure) + cp * 0.2; qt = envQ(departure) + 1e-4; }
  else {
    const twin = build({}); twin.state.forEach((a, n) => a.set(model.state[n])); twin.boundaryLayer.depth.set(model.boundaryLayer.depth); twin.radiation.mlmGate.set(model.radiation.mlmGate);
    twin.core.diagnoseColumn(0, twin.state[0], twin.state[1], twin.state[4], twin.state[5]);
    twin.moist.plumeColumn(0, twin.state[0], twin.state[1], twin.state[4], twin.state[5], 600);
    sl = twin.moist.deep.sourceS; qt = twin.moist.deep.sourceQ;
  }
  let w2 = 1, base = 0, deepAt = NaN;
  const surface = qsOf(T(bottom), p(bottom));
  for (let k = source - 1; k > 0; k--) {
    const lower = zUp(k + 1), depth = zUp(k) - lower, eps = 0.4 * 1.75e-3 * (qsOf(T(k), p(k)) / surface) ** 3, m = (1 + 1.875 * 0.506) * eps;
    const half = Math.exp(-eps * (zMid(k) - lower)), ms = envS(k) + (sl - envS(k)) * half, mq = envQ(k) + (qt - envQ(k)) * half, mid = state(ms, mq, zMid(k), p(k));
    if (!(base > 0) && mid.l > 0) base = pi[0] * levels[k + 1];
    const Tv = T(k) * (1 + 0.608 * Math.max(0, q[k * C]) - Math.max(0, qc[k * C])), B = g * (mid.T * (1 + 0.608 * (mq - mid.l) - mid.l) - Tv) / Tv;
    const w2at = (d) => (m > 0 ? w2 * Math.exp(-2 * m * d) + (B / 3 / m) * (1 - Math.exp(-2 * m * d)) : w2 + 2 * B / 3 * d);
    if (!(w2at(depth) > 0)) {
      if (!(base > 0)) return { base: 0, top: NaN, deep: false, deepAt };
      let a = 0, b = depth;
      for (let n = 0; n < 200; n++) { const c = 0.5 * (a + b); if (w2at(c) > 0) a = c; else b = c; }
      const top = pi[0] * levels[k + 1] * Math.pow(levels[k] / levels[k + 1], 0.5 * (a + b) / depth);
      return { base, top, deep: base - top > DEEP_CLOUD_DEPTH, deepAt, layer: k };
    }
    w2 = w2at(depth);
    if (base > 0 && !Number.isFinite(deepAt) && base - pi[0] * levels[k] > DEEP_CLOUD_DEPTH) deepAt = pi[0] * levels[k];
    const full = Math.exp(-eps * depth);
    sl = envS(k) + (sl - envS(k)) * full; qt = envQ(k) + (qt - envQ(k)) * full;
    const upper = state(sl, qt, zUp(k), pi[0] * levels[k]);
    qt -= 0.5 * upper.l; sl += L * 0.5 * upper.l + FUSION_HEAT * 0.5 * upper.ice;
  }
  return { base, top: pi[0] * levels[1], deep: base > 0 && base - pi[0] * levels[1] > DEEP_CLOUD_DEPTH, deepAt };
}

test('the IFS test parcel types the column: from the deep source at 1 m/s, mixing at 0.4 · 1.75e-3 (q_s/q_s,lowest)³ and keeping half its condensate, its cloud from its first cloudy layer to where its w² vanishes inside a layer; deep beyond 200 hPa, as by hand on Jordan\'s column, a drier one and the trade-wind columns', () => {
  const dt = 600, cases = [
    ['Jordan', plumeColumn({})],
    ['Jordan, free troposphere × 0.3', plumeColumn({}, (z, p) => { const air = jordan(p); return { T: air.T, q: Math.min(air.q * (p < 850e2 ? 0.3 : 1), saturationHumidity(air.T, p)) }; })],
    ['Jordan, free troposphere × 0.1', plumeColumn({}, (z, p) => { const air = jordan(p); return { T: air.T, q: Math.min(air.q * (p < 850e2 ? 0.1 : 1), saturationHumidity(air.T, p)) }; })],
    ['trade wind', tradeWindColumn()],
    ['trade cumulus', plumeColumn({}, tradeCumulus)],
    ['cold column', plumeColumn({}, coldColumn)],
  ];
  assert.equal(MOIST_DEFAULTS.convectionType, 'testParcel');
  assert.deepEqual(TEST_PARCEL, { entrainment: 0.4, removal: 0.5 });
  const lines = [];
  let deepSeen = 0, shallowSeen = 0;
  for (const [label, model] of cases) {
    const [pi, theta, , , q, qc] = model.state;
    model.core.diagnoseColumn(0, pi, theta, q, qc);
    const hand = handParcel(model), before = budget(model, 0);
    const rain = model.moist.plumeColumn(0, pi, theta, q, qc, dt);
    const after = budget(model, 0), d = model.moist.deep;
    lines.push(`${label}: test cloud ${(d.parcelBase / 100).toFixed(1)}–${(d.parcelTop / 100).toFixed(1)} hPa (by hand ${(hand.base / 100).toFixed(1)}–${((hand.deep ? hand.deepAt : hand.top) / 100).toFixed(1)}${hand.deep ? `, w² vanishing at ${(hand.top / 100).toFixed(1)}` : ''}), ${d.parcelDeep ? 'deep' : 'shallow'}; ${d.deep ? `the deep plume tops at ${(model.moist.cumulusTop[0] / 100).toFixed(0)} hPa` : model.moist.cumulusBaseFlux[0] > 0 ? `the shallow plume tops at ${(model.moist.cumulusTop[0] / 100).toFixed(0)} hPa` : 'no plume'}`);
    assert.equal(d.parcelDeep, hand.deep, `${label}: typed alike`);
    assert.ok(Math.abs(d.parcelBase - hand.base) < 1e-6, `${label}: base ${d.parcelBase} against ${hand.base}`);
    const expectedTop = hand.deep ? hand.deepAt : hand.top;
    assert.ok(Math.abs(d.parcelTop - expectedTop) < 1e-6 * expectedTop, `${label}: top ${d.parcelTop} against ${expectedTop}`);
    if (hand.deep) deepSeen++; else shallowSeen++;
    if (!d.parcelDeep) assert.ok(!d.deep, `${label}: shallow by the test parcel, no deep plume`);
    const frozen = d.deep ? model.moist.convectiveFrozen.reduce((x, y) => x + y, 0) : 0;
    assert.ok(Math.abs(after.enthalpy - before.enthalpy - FUSION_HEAT * frozen) < 1e-15 * before.enthalpy && Math.abs(after.water + rain - before.water) < 1e-15 * before.water, `${label}: enthalpy (less the fusion heat of the snow still falling) and water`);
  }
  console.log(lines.join('\n'));
  assert.ok(deepSeen >= 1 && shallowSeen >= 1, `${deepSeen} deep and ${shallowSeen} shallow`);
});

test('under the test parcel a column whose test cloud is deeper than 200 hPa runs the deep plume alone and every other column the shallow plume alone, with column enthalpy and water exact on random columns', () => {
  const dt = 900, model = build({}), C = model.mesh.nCells;
  randomColumns(model);
  const [pi, theta, , , q, qc] = model.state, { moist, core } = model;
  let typed = 0, fired = 0, shallowTyped = 0, shallowRun = 0, worst = 0, deeperThanPlume = 0;
  for (let i = 0; i < C; i++) {
    core.diagnoseColumn(i, pi, theta, q, qc);
    const before = budget(model, i), rain = moist.plumeColumn(i, pi, theta, q, qc, dt), after = budget(model, i), d = moist.deep;
    worst = Math.max(worst, Math.abs(after.enthalpy - before.enthalpy) / before.enthalpy, Math.abs(after.water + rain - before.water) / before.water);
    if (d.parcelDeep) { typed++; assert.ok(d.parcelBase - d.parcelTop > DEEP_CLOUD_DEPTH, `column ${i}: deep with a test cloud of ${d.parcelBase - d.parcelTop} Pa`); }
    else if (d.parcelBase > 0) { shallowTyped++; assert.ok(!(d.parcelBase - d.parcelTop > DEEP_CLOUD_DEPTH), `column ${i}`); }
    if (d.deep) {
      fired++;
      assert.ok(d.parcelDeep, `column ${i}: the deep plume runs only where the test parcel is deep`);
      assert.equal(moist.cumulusBaseFlux[i], d.baseFlux, `column ${i}: the deep plume runs alone`);
      if (!(pi[i] * (model.core.diagnostics.levels[d.base] - model.core.diagnostics.levels[d.top]) > DEEP_CLOUD_DEPTH)) deeperThanPlume++;
    } else if (moist.cumulusBaseFlux[i] > 0) shallowRun++;
  }
  console.log(`${C} random columns: ${typed} typed deep by the test parcel (${fired} fire the deep plume, ${deeperThanPlume} of them with the plume's own cloud no deeper than 200 hPa), ${shallowTyped} typed shallow, ${shallowRun} run the shallow plume; enthalpy and water within ${worst.toExponential(1)}`);
  assert.ok(typed > C / 30 && fired > C / 30 && shallowTyped > C / 30 && shallowRun > C / 30, `${typed}, ${fired}, ${shallowTyped}, ${shallowRun}`);
  assert.ok(worst < 1e-15, `enthalpy and water ${worst}`);
});

function nocturnalColumn(options = {}, { flux = -1e-5, cap = 600, strength = 10 } = {}) {
  const model = plumeColumn({ stratocumulusInversion: null, ...options }, nocturnalSounding(cap, strength));
  model.boundaryLayer.buoyancyFlux[0] = flux;
  return model;
}
function surfaceStack(model) {
  const { K, sigmaMid } = model.core.diagnostics, [pi] = model.state;
  let stack = K - 1;
  while (stack > 0 && pi[0] * sigmaMid[stack - 1] >= pi[0] - MOIST_DEFAULTS.cumulusSourceDepth) stack--;
  return stack;
}
function handDepartures(model) {
  const { sigmaMid } = model.core.diagnostics, [pi] = model.state, tried = [];
  for (let kd = surfaceStack(model) - 1; kd > 0 && pi[0] - pi[0] * sigmaMid[kd] <= 350e2; kd--) tried.push([kd, handParcel(model, kd).deep]);
  return tried;
}
function stepColumn(model, dt) {
  const [pi, theta, , , q, qc] = model.state;
  model.core.diagnoseColumn(0, pi, theta, q, qc);
  const before = budget(model, 0), rain = model.moist.plumeColumn(0, pi, theta, q, qc, dt), after = budget(model, 0), { deep } = model.moist;
  const frozen = deep.deep ? model.moist.convectiveFrozen.reduce((x, y) => x + y, 0) : 0;
  return { deep: { ...deep }, flux: Float64Array.from(model.moist.cumulusFlux), state: snapshot(model, 0), rain, enthalpy: Math.abs(after.enthalpy - before.enthalpy - FUSION_HEAT * frozen) / before.enthalpy, water: Math.abs(after.water + rain - before.water) / before.water };
}

test("plumeDeparture 'elevated' roots the deep plume of a column whose surface test parcel is shallow in the lowest layer above the 50 hPa stack, within 350 hPa of the ground, whose own test parcel is deep: that layer's air plus 0.2 K and 1e-4 kg/kg, no PCAPE_bl, the updraught drawn from that layer alone, with the closure by hand and column enthalpy and water exact", () => {
  const dt = 300;
  assert.equal(MOIST_DEFAULTS.plumeDeparture, 'surface');
  assert.deepEqual(DEPARTURE, { depth: 350e2, temperature: 0.2, humidity: 1e-4, reference: 'lowest' });
  const off = nocturnalColumn(), on = nocturnalColumn({ plumeDeparture: 'elevated' });
  const tried = handDepartures(on), lowest = tried.findIndex(([, deep]) => deep), kd = tried[lowest][0];
  const air = snapshot(on, 0), surface = stepColumn(off, dt), elevated = stepColumn(on, dt), { moist } = on, d = elevated.deep;
  const { K, C, sigmaMid, dSigma, cp, g, geopotential, exnerLayer } = on.core.diagnostics, [pi] = on.state;
  console.log(`nocturnal sounding: surface test ${surface.deep.parcelDeep ? 'deep' : 'shallow'}, hand test parcels from ${tried.map(([k, deep]) => `${k}${deep ? '+' : '-'}`).join(' ')}; elevated departure ${d.departure} (${((pi[0] - pi[0] * sigmaMid[d.departure]) / 100).toFixed(0)} hPa above the ground), PCAPE ${d.pcape.toFixed(1)} Pa, PCAPE_bl ${d.pcapeBoundary}, base flux ${d.baseFlux.toFixed(5)} kg/m²/s, top ${(moist.cumulusTop[0] / 100).toFixed(0)} hPa; enthalpy ${elevated.enthalpy.toExponential(1)}, water ${elevated.water.toExponential(1)}`);
  assert.ok(!surface.deep.parcelDeep && !surface.deep.deep && surface.deep.departure === -1 && surface.deep.rooted, 'with the surface rule the cooled column is shallow');
  assert.ok(lowest >= 0 && tried.slice(lowest + 1).some(([, deep]) => deep), 'a deep departure with deeper ones above it');
  assert.ok(tried.slice(0, lowest).every(([, deep]) => !deep), 'none of the layers below it is deep');
  assert.equal(d.departure, kd, 'the lowest deep departure');
  assert.ok(d.deep && !d.rooted && d.parcelDeep, 'a deep plume rooted above the boundary layer');
  const envS = cp * air.theta[kd] * exnerLayer[kd * C] + geopotential[kd * C] - LATENT_HEAT * air.qc[kd], envQ = Math.max(0, air.q[kd]) + air.qc[kd];
  assert.ok(Math.abs(d.sourceS - (envS + cp * 0.2)) < 1e-13 * envS && d.sourceQ === envQ + 1e-4 && d.excessT === 0.2 && d.excessQ === 1e-4 && d.sourceMass === pi[0] * dSigma[kd], 'the layer\'s air plus the excess');
  assert.equal(d.pcapeBoundary, 0);
  const hand = closureByHand(on, 0), gate = Math.min(1, Math.max(0, 0.5 + (MOIST_DEFAULTS.inhibitionThreshold - d.inhibition) / MOIST_DEFAULTS.inhibitionThreshold));
  const expected = Math.min(gate * hand.pcape / (hand.tau * d.consumptionP), MOIST_DEFAULTS.cumulusBoundaryLoss * d.sourceMass / (g * dt));
  assert.ok(Math.abs(d.pcape - hand.pcape) < 1e-12 * hand.pcape && Math.abs(d.tau - hand.tau) < 1e-12 * hand.tau, 'PCAPE and τ');
  assert.ok(Math.abs(d.baseFlux - expected) < 1e-12 * expected, `base flux ${d.baseFlux} against PCAPE / (τ F_P) or the cap ${expected}`);
  assert.ok(elevated.flux[kd] === 1 && elevated.flux.slice(kd + 1).every((x) => x === 0) && elevated.flux[kd - 1] > 0, 'the updraught is one at the top of the layer and none below it');
  assert.ok(elevated.enthalpy < 1e-15 && elevated.water < 1e-15, 'enthalpy and water');
  const dry = nocturnalColumn({ plumeDeparture: 'elevated', downdraftShare: 0 }), kept = snapshot(dry, 0), after = stepColumn(dry, dt);
  for (let k = kd + 1; k < K; k++) assert.ok(after.state.theta[k] === kept.theta[k] && after.state.q[k] === kept.q[k], `layer ${k} below the source stays put without a downdraft`);
  assert.ok(after.state.theta[kd] !== kept.theta[kd] && after.deep.deep, 'the source layer moves');
  const full = nocturnalColumn({ plumeDeparture: 'elevated' });
  full.moist.adjust(full.state, 0, 1, dt);
  console.log(`rain below the elevated cloud base at ${(pi[0] * sigmaMid[d.base] / 100).toFixed(0)} hPa: ${(full.moist.falling.evaporated * 86400 / dt).toFixed(2)} mm/d evaporates, ${(full.moist.convectivePrecipitation[0] * 86400 / dt).toFixed(2)} mm/d reaches the ground`);
  assert.ok(full.moist.falling.evaporated > 0 && full.moist.convectivePrecipitation[0] < elevated.rain, 'rain evaporates below the cloud base into the cooled air');
});

test("the departure search stops 350 hPa above the ground, at the layer midpoint: a limit 1 Pa short of the lowest deep departure's height loses the plume, a wider one still finds the lowest", () => {
  const dt = 300, model = () => nocturnalColumn({ plumeDeparture: 'elevated' });
  const reference = model(), kd = handDepartures(reference).find(([, deep]) => deep)[0], [pi] = reference.state, { sigmaMid } = reference.core.diagnostics, height = pi[0] - pi[0] * sigmaMid[kd];
  const found = (limit) => { const kept = DEPARTURE.depth; DEPARTURE.depth = limit; try { return stepColumn(model(), dt).deep; } finally { DEPARTURE.depth = kept; } };
  const at = found(height), below = found(height - 1), wide = found(1e6);
  assert.ok(at.departure === kd && below.departure === -1 && !below.deep, `limit ${height}: departure ${at.departure}; 1 Pa less: ${below.departure}`);
  assert.equal(wide.departure, kd, 'a wider limit still finds the lowest');
  assert.equal(DEPARTURE.depth, 350e2);
  assert.ok(height <= 350e2, `the lowest deep departure lies ${height} Pa above the ground`);
});

test("plumeDeparture 'elevatedStable' searches only where the surface buoyancy flux is not positive: a column over a heating surface steps as under 'surface', one over a cooling or neutral surface as under 'elevated'", () => {
  const dt = 300, run = (departure, flux) => stepColumn(nocturnalColumn({ plumeDeparture: departure }, { flux }), dt);
  for (const flux of [4e-4, 1e-12]) {
    const stable = run('elevatedStable', flux), surface = run('surface', flux), elevated = run('elevated', flux);
    assert.deepEqual(stable.state, surface.state, `flux ${flux}: elevatedStable steps as surface`);
    assert.ok(stable.deep.departure === -1 && elevated.deep.departure >= 0 && elevated.deep.deep, `flux ${flux}: stable ${stable.deep.departure}, elevated ${elevated.deep.departure}`);
  }
  for (const flux of [-1e-5, 0, -1e-12]) {
    const stable = run('elevatedStable', flux), elevated = run('elevated', flux);
    assert.ok(stable.deep.departure >= 0 && stable.deep.deep, `flux ${flux}: departure ${stable.deep.departure}`);
    assert.deepEqual(stable.state, elevated.state, `flux ${flux}: elevatedStable steps as elevated`);
    assert.deepEqual(stable.deep, elevated.deep);
  }
});

test("under a +10 K/d subcloud warming over land the elevated plume has no PCAPE_bl and fires, where the same warming holds a surface-rooted Jordan plume back without the PCAPE minimum, and a column whose surface test is deep steps alike under 'surface' and 'elevated'", () => {
  const dt = 300, land = () => Uint8Array.from({ length: 16 * 16 * 10 + 2 }, (_, i) => (i === 0 ? 1 : 0));
  const heated = (model) => { subcloudVirtual(model, 0, 10, dt); return stepColumn(model, dt); };
  const night = heated(nocturnalColumn({ plumeDeparture: 'elevated', land: land() })), calm = stepColumn(nocturnalColumn({ plumeDeparture: 'elevated', land: land() }), dt);
  assert.ok(night.deep.deep && !night.deep.rooted && night.deep.pcapeBoundary === 0 && night.deep.baseFlux === calm.deep.baseFlux, 'the heated elevated plume fires as the calm one does');
  const held = heated(plumeColumn({ ...PREVIOUS_PCAPE_MINIMUM, plumeDeparture: 'elevated', land: land() })), alone = heated(plumeColumn({ ...PREVIOUS_PCAPE_MINIMUM, land: land() }));
  assert.ok(!held.deep.deep && held.deep.pcapeBoundary > held.deep.pcape && held.deep.rooted, 'the surface-rooted plume is held back');
  assert.deepEqual(held.state, alone.state);
  assert.deepEqual({ ...held.deep, typeMargin: 0 }, { ...alone.deep, typeMargin: 0 });
});

test("with plumeDeparture 'elevated' the random columns whose search finds nothing step exactly as under 'surface', those it finds keep column enthalpy and water exact, and no departure lies beyond 350 hPa", () => {
  const dt = 900, surface = build({ stratocumulusInversion: null }), elevated = build({ stratocumulusInversion: null, plumeDeparture: 'elevated' }), C = surface.mesh.nCells, { K, sigmaMid } = surface.core.diagnostics;
  randomColumns(surface, { nocturnal: true }); randomColumns(elevated, { nocturnal: true });
  let found = 0, fired = 0, same = 0, worst = 0, lowest = Infinity, highest = -Infinity;
  for (let i = 0; i < C; i++) {
    const s = surface.state, e = elevated.state;
    surface.core.diagnoseColumn(i, s[0], s[1], s[4], s[5]); elevated.core.diagnoseColumn(i, e[0], e[1], e[4], e[5]);
    const rainS = surface.moist.plumeColumn(i, s[0], s[1], s[4], s[5], dt), before = budget(elevated, i), rainE = elevated.moist.plumeColumn(i, e[0], e[1], e[4], e[5], dt), after = budget(elevated, i), d = elevated.moist.deep;
    worst = Math.max(worst, Math.abs(after.enthalpy - before.enthalpy - FUSION_HEAT * (d.deep ? elevated.moist.convectiveFrozen.reduce((x, y) => x + y, 0) : 0)) / before.enthalpy, Math.abs(after.water + rainE - before.water) / before.water);
    if (d.departure >= 0) {
      found++; if (d.deep) fired++;
      const above = e[0][i] - e[0][i] * sigmaMid[d.departure];
      lowest = Math.min(lowest, above); highest = Math.max(highest, above);
      assert.ok(above <= 350e2 && !d.rooted && d.parcelDeep, `column ${i}: departure ${d.departure} at ${above} Pa`);
    } else {
      same++;
      assert.equal(rainE, rainS, `column ${i}: rain`);
      for (let k = 0; k < K; k++) assert.ok(e[1][k * C + i] === s[1][k * C + i] && e[4][k * C + i] === s[4][k * C + i] && e[5][k * C + i] === s[5][k * C + i], `column ${i}, layer ${k}`);
    }
  }
  console.log(`${C} random columns with a nocturnal pass: ${found} found a departure (${fired} fire the deep plume), between ${(lowest / 100).toFixed(0)} and ${(highest / 100).toFixed(0)} hPa above the ground; ${same} step as under 'surface'; enthalpy and water within ${worst.toExponential(1)}`);
  assert.ok(found > C / 30 && fired > C / 40 && same > C / 2, `${found} found, ${fired} fire, ${same} alike`);
  assert.ok(worst < 1e-15, `enthalpy and water ${worst}`);
});

test("a plumeDeparture other than 'surface', 'elevated' or 'elevatedStable', and a searching one without the test parcel or the Bechtold closure, are refused on both engines", async () => {
  const refused = [{ plumeDeparture: 'deep' }, { plumeDeparture: true }, { plumeDeparture: 'elevated', convectionType: 'cloudDepth' }, { plumeDeparture: 'elevatedStable', convectionType: 'top' }, { plumeDeparture: 'elevated', capeClosure: 'threshold' }];
  for (const options of refused) assert.throws(() => build(options), /plumeDeparture/, `CPU: ${JSON.stringify(options)}`);
  for (const plumeDeparture of ['surface', 'elevated', 'elevatedStable']) build({ plumeDeparture });
  build({ capeClosure: 'threshold' });
  if (!gpuAvailable) return;
  const { physicsConstants } = await import('../js/gpu/physics.gpu.js');
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  for (const options of refused) assert.throws(() => physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, ...options }), /plumeDeparture/, `GPU: ${JSON.stringify(options)}`);
  for (const plumeDeparture of ['surface', 'elevated', 'elevatedStable']) physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, plumeDeparture });
});

test('the mixed-phase plume\'s air holds its liquid-ice static energy: its temperature, condensate and ice as by bisection over the phase ramp, liquid above it and all ice below it', () => {
  const model = build({}), O = MOIST_DEFAULTS, { moist } = model, liquid = build({ plumePhase: 'liquid' }).moist, { g, cp } = model.core.diagnostics, sat = { qs: 0, slope: 0, liquid: 1 };
  assert.equal(O.plumePhase, 'mixed');
  let worstT = 0, worstL = 0, worstI = 0, n = 0, warm = 0;
  for (const p of [900e2, 700e2, 550e2, 450e2, 350e2, 250e2, 180e2]) for (const T of [205, 230, 240, 250, 260, 268, 272, 275, 290]) for (const l of [0, 2e-4, 1e-3, 4e-3]) {
    const z = 7400 * Math.log(101500 / p), qs = cloudSaturation(T, p, true, O.liquidTemperature, O.iceTemperature, sat).qs, a = sat.liquid;
    const sl = cp * T + g * z - LATENT_HEAT * a * l - (LATENT_HEAT + FUSION_HEAT) * (1 - a) * l, qt = qs + l;
    const hand = mixedState(sl, qt, z, p, O, cp, g), air = moist.plumeAir(sl, qt, z, p, T - 1);
    worstT = Math.max(worstT, Math.abs(air.T - hand.T)); worstL = Math.max(worstL, Math.abs(air.liquid - hand.l)); worstI = Math.max(worstI, Math.abs(air.ice - hand.ice)); n++;
    if (T > O.liquidTemperature) { const plain = liquid.plumeAir(sl, qt, z, p, T - 1); assert.ok(Math.abs(plain.T - air.T) < 1e-9 && air.ice === 0, `above the ramp the liquid plume's ${plain.T} against ${air.T}`); warm++; }
    if (T < O.iceTemperature && l > 0) assert.ok(Math.abs(air.ice - air.liquid) < 1e-15, 'below the ramp all of it is ice');
  }
  console.log(`${n} states over 180-900 hPa and 205-290 K: temperature within ${worstT.toExponential(1)} K of the bisection, condensate ${worstL.toExponential(1)}, ice ${worstI.toExponential(1)} kg/kg; ${warm} warm states as the liquid plume's`);
  assert.ok(worstT < 1e-5 && worstL < 1e-8 && worstI < 1e-8, `T ${worstT}, l ${worstL}, ice ${worstI}`);
});

test('the mixed-phase deep plume on Jordan\'s column gains the fusion heat of its frozen condensate, rising higher with more CAPE than the liquid plume, and its frozen rain melts in the first layer at 273.15 K or warmer below where it formed; column enthalpy and water exact', () => {
  const dt = 600, run = (options) => {
    const model = plumeColumn(options), [pi, theta, , , q, qc] = model.state, { K, C, exnerLayer } = model.core.diagnostics;
    model.core.diagnoseColumn(0, pi, theta, q, qc);
    const T = Float64Array.from({ length: K }, (_, k) => theta[k * C] * exnerLayer[k * C]), before = budget(model, 0);
    const rain = model.moist.plumeColumn(0, pi, theta, q, qc, dt), after = budget(model, 0), m = model.moist, d = m.deep;
    let made = 0, melted = 0, highest = -1, meltAt = -1;
    for (let k = d.top + 1; k < d.base; k++) if (m.plumeFrozen[k] > 0) { made += m.cumulusFlux[k] * m.plumeFrozen[k]; if (highest < 0) highest = k; }
    for (let k = 0; k < K; k++) if (m.convectiveMelted[k] > 0) { melted += m.convectiveMelted[k]; if (meltAt < 0) meltAt = k; }
    let firstWarm = -1;
    for (let k = highest + 1; k < K && highest >= 0; k++) if (T[k] >= MOIST_DEFAULTS.liquidTemperature) { firstWarm = k; break; }
    return { top: m.cumulusTop[0], cape: d.cape, flux: d.baseFlux, rain, made, melted, meltAt, firstWarm, meltP: meltAt >= 0 ? pi[0] * model.core.diagnostics.sigmaMid[meltAt] : NaN, meltT: meltAt >= 0 ? T[meltAt] : NaN, enthalpy: (after.enthalpy - before.enthalpy) / before.enthalpy, water: (after.water + rain - before.water) / before.water };
  };
  const mixed = run({}), liquid = run({ plumePhase: 'liquid' });
  console.log(`Jordan, deep plume: mixed phase top ${(mixed.top / 100).toFixed(0)} hPa, CAPE ${mixed.cape.toFixed(0)} J/kg, base flux ${mixed.flux.toFixed(4)}; liquid top ${(liquid.top / 100).toFixed(0)} hPa, CAPE ${liquid.cape.toFixed(0)} J/kg, base flux ${liquid.flux.toFixed(4)} kg/m²/s; frozen rain per unit base flux ${mixed.made.toExponential(3)}, all melting at ${(mixed.meltP / 100).toFixed(0)} hPa (${mixed.meltT.toFixed(2)} K); enthalpy ${mixed.enthalpy.toExponential(1)}, water ${mixed.water.toExponential(1)}`);
  assert.ok(mixed.top < liquid.top && mixed.cape > liquid.cape, 'higher, with more CAPE');
  assert.ok(mixed.made > 0 && Math.abs(mixed.melted - mixed.made) < 1e-15 * mixed.made, `made ${mixed.made}, melted ${mixed.melted}`);
  assert.ok(mixed.meltAt === mixed.firstWarm, `melts in layer ${mixed.meltAt}, the first warm one below ${mixed.firstWarm}`);
  for (const r of [mixed, liquid]) assert.ok(Math.abs(r.enthalpy) < 1e-15 && Math.abs(r.water) < 1e-15, `enthalpy ${r.enthalpy}, water ${r.water}`);
});

test('over a 265 K surface the mixed-phase plume\'s snow reaches the ground frozen and the surface adds no second fusion heat for it: the adjustment changes column enthalpy by L_f times the precipitation under either phase, water exact', () => {
  const dt = 900, run = (options) => {
    const model = plumeColumn(options, coldColumn);
    model.state[3][0] = 268;
    const before = budget(model, 0);
    model.phases.adjust(0, 1, dt);
    const after = budget(model, 0), rain = model.moist.rain[0];
    return { deep: model.moist.deep.deep, rain, convective: model.moist.convectivePrecipitation[0], snow: model.moist.convectiveSnow[0], enthalpy: (after.enthalpy - before.enthalpy - FUSION_HEAT * rain) / before.enthalpy, water: (after.water + rain - before.water) / before.water };
  };
  const mixed = run({}), liquid = run({ plumePhase: 'liquid' });
  console.log(`cold column: mixed phase ${mixed.rain.toExponential(3)} kg/m² of precipitation, ${mixed.snow.toExponential(3)} of it the deep plume's snow; liquid ${liquid.rain.toExponential(3)}; enthalpy less L_f P ${mixed.enthalpy.toExponential(1)} / ${liquid.enthalpy.toExponential(1)}, water ${mixed.water.toExponential(1)} / ${liquid.water.toExponential(1)}`);
  assert.ok(mixed.deep && mixed.snow > 0 && mixed.snow <= mixed.convective, `snow ${mixed.snow} of ${mixed.convective}`);
  assert.equal(liquid.snow, 0);
  for (const r of [mixed, liquid]) assert.ok(Math.abs(r.enthalpy) < 1e-15 && Math.abs(r.water) < 1e-15, `enthalpy ${r.enthalpy}, water ${r.water}`);
});

/*
 * A column over a surface at Ts whose air keeps the mixing ratio of
 * saturation at the cloud base, `base` metres up, below it (lapse 9.7 K/km)
 * and is `humid` saturated above it (lapse 8.5 K/km to 205 K): the plume's
 * base lies above the 0 °C level, with unsaturated air below it.
 */
const coldBase = (Ts, base, humid) => (z, p) => {
  const top = Ts - 9.7e-3 * base;
  if (z < base) return { T: Ts - 9.7e-3 * z, q: 0.98 * saturationHumidity(top, 101500 * Math.pow(top / Ts, 3.5)) };
  const T = Math.max(205, top - 8.5e-3 * (z - base));
  return { T, q: humid * saturationHumidity(T, p) };
};
const COLD_BASES = [276, 278, 280, 283].flatMap((Ts) => [800, 1200, 1600].flatMap((base) => [0.8, 0.95].map((humid) => [Ts, base, humid])));

test('a deep plume whose frozen rain sublimates in the unsaturated air below its cloud base before it melts, or is taken by the downdraft, changes column enthalpy by L_f times the snow reaching the ground only, on the CPU exactly and on the GPU to single precision, alike on both engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const dt = 900, model = build({}), C = model.mesh.nCells, { K, exnerLayer } = model.core.diagnostics, { moist } = model;
  for (let i = 0; i < C; i++) {
    place(model, i, 101500, coldBase(...COLD_BASES[i % COLD_BASES.length]));
    setDepth(model, i, 500);
    model.boundaryLayer.buoyancyFlux[i] = 4e-4; model.boundaryLayer.friction[i] = 0.25; model.radiation.mlmGate[i] = 0.3;
  }
  for (const a of model.state) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  for (const a of [model.boundaryLayer.depth, model.radiation.mlmGate]) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  const [pi, theta, , , q, qc] = model.state;
  for (let i = 0; i < C; i++) model.core.diagnoseColumn(i, pi, theta, q, qc);
  const { createGpuCore } = await import('../js/gpu/core.gpu.js');
  const gpu = await createGpuCore(model.mesh, { levels, physics: {} }), { device, buffers, kernels, layout } = gpu;
  gpu.upload(model.state);
  gpu.uploadPhysics({ mlmGate: model.radiation.mlmGate, concentration: model.seaIce.concentration });
  for (const [field, values] of [['DEPTH', model.boundaryLayer.depth], ['BUOY', model.boundaryLayer.buoyancyFlux], ['USTAR', model.boundaryLayer.friction], ['MIXTOP', model.boundaryLayer.mixingTop]]) device.queue.writeBuffer(buffers.PH, 4 * layout.PH[field], Float32Array.from(values));
  device.queue.writeBuffer(buffers.P, 0, Float32Array.from([dt, 0, 1, 0, 0, 0, 0, 0]));
  const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
  const group = device.createBindGroup({ layout: kernels.adjustMix.getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) });
  for (const name of ADJUST_KERNELS) { pass.setPipeline(kernels[name]); pass.setBindGroup(0, group); pass.dispatchWorkgroups(Math.ceil(C / 64)); }
  pass.end();
  device.queue.submit([encoder.finish()]);
  const gpuState = await gpu.download(), ph = await gpu.downloadPhysics(), gpuBefore = model.state.map((a) => Float64Array.from(a));
  const twin = build({});
  twin.state.forEach((a, n) => a.set(model.state[n]));
  for (const field of ['depth', 'buoyancyFlux', 'friction']) twin.boundaryLayer[field].set(model.boundaryLayer[field]);
  twin.radiation.mlmGate.set(model.radiation.mlmGate);
  twin.phases.adjust(0, C, dt);
  const enthalpyOf = (state, i) => { const [p0, th, , , qv] = state, { dSigma, g, cp } = model.core.diagnostics; let h = 0; for (let k = 0; k < K; k++) { const x = k * C + i; h += (cp * th[x] * exnerLayer[x] + LATENT_HEAT * qv[x]) * p0[i] * dSigma[k] / g; } return h; };
  let sublimating = 0, worstCpu = 0, worstGpu = 0, worstTheta = 0, deep = 0;
  for (let i = 0; i < C; i++) {
    const before = budget(model, i);
    moist.adjust(model.state, i, i + 1, dt);
    const after = budget(model, i), d = moist.deep, rain = moist.rain[i];
    if (!d.deep) continue;
    deep++;
    let melt = -1;
    for (let k = 0; k < K; k++) if (moist.convectiveMelted[k] > 0) melt = k;
    if (melt > d.base) sublimating++;
    worstCpu = Math.max(worstCpu, Math.abs(after.enthalpy - before.enthalpy - FUSION_HEAT * moist.convectiveSnow[i]) / before.enthalpy, Math.abs(after.water + rain - before.water) / before.water);
    const freezing = gpuState[1][(K - 1) * C + i] * exnerLayer[(K - 1) * C + i] < MOIST_DEFAULTS.liquidTemperature;
    worstGpu = Math.max(worstGpu, Math.abs(enthalpyOf(gpuState, i) - enthalpyOf(gpuBefore, i) - (freezing ? FUSION_HEAT * ph.STEPRAIN[i] : 0)) / FUSION_HEAT);
    for (let k = 0; k < K; k++) worstTheta = Math.max(worstTheta, Math.abs(twin.state[1][k * C + i] - gpuState[1][k * C + i]));
  }
  console.log(`${C} cold-based columns: ${deep} deep, ${sublimating} melting below their cloud base; CPU enthalpy less L_f times the snow and water within ${worstCpu.toExponential(1)}; GPU enthalpy off L_f times the precipitation onto freezing air by at most ${worstGpu.toExponential(1)} kg/m² of L_f; θ apart by ${worstTheta.toExponential(1)} K`);
  assert.ok(deep > C / 3 && sublimating > C / 3, `${deep} deep, ${sublimating} melting below the base`);
  assert.ok(worstCpu < 1e-15, `CPU ${worstCpu}`);
  assert.ok(worstGpu < 1e-3 && worstTheta < 1e-3, `GPU ${worstGpu} kg/m² of L_f, θ ${worstTheta} K`);
});

test('the IFS updraught conversion: on Jordan\'s column each layer rains l (1 − exp(−a Δz)) with a = c0/(0.75 w)(1 − exp(−(l/l_crit)²)) and the Bergeron–Findeisen factor below 268.16 K, nothing where l is at most 0.3 g/kg, and the plume detrains more of what it rains than under the previous conversion, above 400 hPa; a drier free troposphere still tops lower', () => {
  const dt = 600, P = IFS_PRECIPITATION, O = MOIST_DEFAULTS, model = plumeColumn({}), [pi, theta, , , q, qc] = model.state, { K, C, levels } = model.core.diagnostics;
  assert.equal(O.plumeConversion, 'sundqvist');
  model.core.diagnoseColumn(0, pi, theta, q, qc);
  const before = budget(model, 0), rain = model.moist.plumeColumn(0, pi, theta, q, qc, dt), after = budget(model, 0), m = model.moist, d = m.deep;
  const { geopotential, g, cp, exnerLayer, exnerLower } = model.core.diagnostics, thetaV = model.core.arrays.thetaV;
  const zUp = (k) => (geopotential[k * C] + cp * thetaV[k * C] * (exnerLayer[k * C] - exnerLower[(k - 1) * C])) / g;
  let worst = 0, dry = 0, raining = 0, made = 0, detrained = 0, weighted = 0, cold = 0;
  for (let k = d.base - 1; k > d.top; k--) {
    const l = m.plumeCarried[k] + m.plumeRain[k], T = m.plumeInterfaceT[k], depth = zUp(k) - zUp(k + 1);
    if (!(l > P.seaThreshold)) { assert.equal(m.plumeRain[k], 0, `layer ${k}: no conversion at ${l}`); dry++; continue; }
    const alpha = liquidFraction(T, O.liquidTemperature, O.iceTemperature), bf = T < P.bergeron ? 1 + 0.5 * Math.sqrt(Math.min(P.bergeron - T, P.bergeron - P.ice)) : 1;
    if (bf > 1) cold++;
    const w = Math.min(P.speed, Math.max(O.plumeVelocity, Math.sqrt(m.plumeSpeed[k]))), a = P.conversion * (P.liquidFactor * alpha + 1 - alpha) * bf / (P.velocityScale * w) * (1 - Math.exp(-((l * bf / P.critical) ** 2)));
    const expected = l * (1 - Math.exp(-a * depth));
    worst = Math.max(worst, Math.abs(m.plumeRain[k] - expected) / expected); raining++;
  }
  const detrainment = (mm, dd) => {
    let rained = 0, out = 0, at = 0;
    for (let k = dd.top; k < dd.base; k++) {
      if (k > dd.top) rained += mm.cumulusFlux[k] * mm.plumeRain[k];
      const leaving = Math.max(0, mm.cumulusFlux[k + 1] - mm.cumulusFlux[k]) * mm.plumeCarried[k + 1];
      out += leaving; at += leaving * pi[0] * levels[k + 1];
    }
    return { share: out / rained, centroid: at / out, rained, out };
  };
  ({ rained: made, out: detrained } = detrainment(m, d));
  const centroid = detrainment(m, d).centroid;
  const old = plumeColumn({ plumeConversion: 'zhangMcFarlane' });
  old.core.diagnoseColumn(0, old.state[0], old.state[1], old.state[4], old.state[5]);
  old.moist.plumeColumn(0, old.state[0], old.state[1], old.state[4], old.state[5], dt);
  const previous = detrainment(old.moist, old.moist.deep);
  console.log(`Jordan under the IFS conversion: rain made in ${raining} layers off the analytic integral by at most ${worst.toExponential(1)} (${cold} below 268.16 K), none in ${dry} at or below 0.3 g/kg; detrained condensate ${(detrained / made).toFixed(3)} of the rain made, centred at ${(centroid / 100).toFixed(0)} hPa (previous conversion ${previous.share.toFixed(3)} at ${(previous.centroid / 100).toFixed(0)} hPa); top ${(m.cumulusTop[0] / 100).toFixed(0)} hPa, CAPE ${d.cape.toFixed(0)} J/kg, rain ${(rain * 86400 / dt).toFixed(1)} mm/d; enthalpy ${((after.enthalpy - before.enthalpy) / before.enthalpy).toExponential(1)}, water ${((after.water + rain - before.water) / before.water).toExponential(1)}`);
  assert.ok(raining > 3 && cold > 0, `${raining} raining, ${cold} cold`);
  assert.ok(worst < 1e-12, `rain off the integral by ${worst}`);
  assert.ok(detrained / made > 2 * previous.share && centroid < 400e2, `detrained ${detrained / made} at ${centroid} against ${previous.share}`);
  assert.ok(Math.abs(after.enthalpy - before.enthalpy) < 1e-15 * before.enthalpy && Math.abs(after.water + rain - before.water) < 1e-15 * before.water, 'enthalpy and water');
  const top = (factor) => {
    const c = plumeColumn({}, (z, p) => { const air = jordan(p); return { T: air.T, q: Math.min(air.q * (p < 850e2 ? factor : 1), saturationHumidity(air.T, p)) }; });
    c.moist.adjust(c.state, 0, 1, dt);
    return c.moist.deep.deep ? c.moist.cumulusTop[0] : Infinity;
  };
  const tops = [1, 0.8, 0.6, 0.4].map(top);
  console.log(`free-tropospheric humidity × 1, 0.8, 0.6, 0.4: tops ${tops.map((t) => (t / 100).toFixed(0)).join(', ')} hPa`);
  for (let n = 1; n < tops.length; n++) assert.ok(tops[n] >= tops[n - 1], `top ${tops[n]} under ${tops[n - 1]}`);
});

async function parity(options, { momentum = false, nocturnal = false, land = false, ascent = false, grid = levels, ground = null, columns = null } = {}) {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js');
  const landMask = land ? Uint8Array.from({ length: 16 * 16 * 10 + 2 }, (_, i) => (i % 3 === 0 ? 1 : 0)) : null;
  const model = createModel(new Grid(16), { ...(ground ? { topography: syntheticTopography(90, 180, ground) } : {}), ocean: false, levels: grid, moist: landMask ? { ...options, land: landMask } : options }), { moist, core, mesh } = model, C = mesh.nCells, { K } = core.diagnostics, dt = 900;
  const [pi, theta, , surfaceT, q, qc] = model.state;
  const fields = randomColumns(model, { nocturnal });
  if (columns) columns(model, fields);
  const { buoyancy, friction, regime, mixingTop, stratiform, concentration, saved, sensible, evaporation } = fields;
  if (ascent) {
    let ascentSeed = 777;
    const ascentRandom = () => { ascentSeed = (ascentSeed * 1103515245 + 12345) % 2147483648; return ascentSeed / 2147483648; };
    for (let k = 1; k < K; k++) for (let i = 0; i < C; i++) core.diagnostics.piSigmaDot[k * C + i] = Math.fround(0.6 * ascentRandom() - 0.4);
    if (core.diagnostics.omega) for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) core.diagnostics.omega[k * C + i] = Math.fround(0.6 * ascentRandom() - 0.4);
  }
  let cloudSeed = 4242;
  const cloudRandom = () => { cloudSeed = (cloudSeed * 1103515245 + 12345) % 2147483648; return cloudSeed / 2147483648; };
  for (let x = moist.cumulusK0 * C; x < K * C; x++) if (cloudRandom() < 0.3) { moist.cumulusCover[x] = Math.fround(0.05 * cloudRandom()); moist.cumulusWater[x] = Math.fround(1e-3 * cloudRandom()); }
  const priorCover = Float64Array.from(moist.cumulusCover), priorWater = Float64Array.from(moist.cumulusWater);
  const phis = ground ? surfaceGeopotential(mesh, model.geography) : null, above = (h, i) => h - (phis ? phis[i] / core.diagnostics.g : 0);
  const gpu = await createGpuCore(mesh, { levels: grid, physics: options, ...(phis ? { surfaceGeopotential: phis } : {}) });
  const { device, buffers, kernels, layout } = gpu;
  gpu.upload(model.state);
  gpu.uploadPhysics({ mlmGate: model.radiation.mlmGate, concentration, cumulusCover: priorCover.subarray(moist.cumulusK0 * C), cumulusWater: priorWater.subarray(moist.cumulusK0 * C) });
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.REGIME, Float32Array.from(regime));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.MIXTOP, Float32Array.from(mixingTop, above));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.STRAT, Float32Array.from(stratiform));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.DEPTH, Float32Array.from(model.boundaryLayer.depth, above));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.BUOY, Float32Array.from(buoyancy));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.USTAR, Float32Array.from(friction));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.SUBTV, Float32Array.from(saved));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.SH, Float32Array.from(sensible));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.EVAP, Float32Array.from(evaporation));
  const landUsed = landMask ?? (ground ? model.geography.land : null);
  if (landUsed) device.queue.writeBuffer(buffers.PH, 4 * layout.PH.LAND, Float32Array.from(landUsed));
  let dragSeed = 31;
  const dragRandom = () => { dragSeed = (dragSeed * 1103515245 + 12345) % 2147483648; return dragSeed / 2147483648; };
  model.boundaryLayer.surfaceDrag.set(Float64Array.from({ length: C }, () => (momentum ? 0 : Math.fround(0.08 * dragRandom()))));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.SDRAG, Float32Array.from(model.boundaryLayer.surfaceDrag));
  if (ascent) device.queue.writeBuffer(buffers.D, 4 * layout.D.PSD, Float32Array.from(core.diagnostics.piSigmaDot));
  if (ascent && core.diagnostics.omega) device.queue.writeBuffer(buffers.D, 4 * layout.D.OMG, Float32Array.from(core.diagnostics.omega));
  device.queue.writeBuffer(buffers.P, 0, Float32Array.from([dt, 0, 1, 0, 0, 0, 0, 0]));
  const group = device.createBindGroup({ layout: kernels.adjustMix.getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) });
  const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
  for (const name of ADJUST_KERNELS) { pass.setPipeline(kernels[name]); pass.setBindGroup(0, group); pass.dispatchWorkgroups(Math.ceil(C / 64)); }
  if (momentum) {
    pass.setPipeline(kernels.mixMomentum);
    pass.setBindGroup(0, device.createBindGroup({ layout: kernels.mixMomentum.getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) }));
    pass.dispatchWorkgroups(Math.ceil(mesh.nEdges / 64));
  }
  pass.end();
  device.queue.submit([encoder.finish()]);
  const after = await gpu.download(), ph = await gpu.downloadPhysics();
  const response = new Float64Array(C), warmed = new Float64Array(C), warmedFlux = new Float64Array(C), thetaResponse = new Float64Array(K * C);
  let windResponse = null;
  {
    const kept = model.state.map((a) => Float64Array.from(a)), keptSaved = Float64Array.from(saved), keptSnow = Float64Array.from(model.seaIce.snow);
    for (let x = 0; x < saved.length; x++) saved[x] *= 1 + 2 ** -23;
    model.phases.adjust(0, C, dt);
    for (let i = 0; i < C; i++) response[i] = moist.convectivePrecipitation[i] + moist.largeScalePrecipitation[i];
    model.state.forEach((a, n) => a.set(kept[n])); saved.set(keptSaved); model.seaIce.snow.set(keptSnow); moist.cumulusCover.set(priorCover); moist.cumulusWater.set(priorWater);
    for (const a of [moist.precipitation, moist.convectivePrecipitation, moist.largeScalePrecipitation, moist.fogPrecipitation]) a.fill(0);
    for (let x = 0; x < theta.length; x++) theta[x] *= 1 + 2 ** -23;
    model.phases.adjust(0, C, dt);
    for (let i = 0; i < C; i++) { warmed[i] = moist.convectivePrecipitation[i] + moist.largeScalePrecipitation[i]; warmedFlux[i] = moist.cumulusBaseFlux[i]; }
    thetaResponse.set(theta);
    if (momentum) { model.phases.mixMomentum(0, mesh.nEdges, dt); windResponse = Float64Array.from(model.state[2]); }
    model.state.forEach((a, n) => a.set(kept[n])); saved.set(keptSaved); model.seaIce.snow.set(keptSnow); moist.cumulusCover.set(priorCover); moist.cumulusWater.set(priorWater);
    for (const a of [moist.precipitation, moist.convectivePrecipitation, moist.largeScalePrecipitation, moist.fogPrecipitation]) a.fill(0);
  }
  moist.trace.convection = new Float64Array(K * C);
  const threshold = new Set();
  const typedAtEdge = new Set();
  const departure = new Int16Array(C).fill(-1);
  let midLevel = 0, deepened = 0, continued = 0, sentBack = 0;
  const continuedColumns = [];
  const downdraftAtEdge = new Set(), { cp, exnerLayer } = core.diagnostics, phiLayer = core.diagnostics.geopotential;
  for (let i = 0; i < C; i++) {
    core.diagnoseColumn(i, pi, theta, q, qc);
    const energy = Float64Array.from({ length: K }, (_, k) => cp * theta[k * C + i] * exnerLayer[k * C + i] + phiLayer[k * C + i] + LATENT_HEAT * q[k * C + i]);
    model.phases.adjust(i, i + 1, dt);
    if (moist.deep.conversionMargin < 1e-3) threshold.add(i);
    if (moist.deep.typeMargin < 1e-3) typedAtEdge.add(i);
    departure[i] = moist.deep.departure;
    if (moist.deep.midLevel) midLevel++;
    if (moist.deep.deepened) deepened++;
    if (moist.deep.deepened && moist.deep.midLevel) sentBack++;
    if (moist.deep.deepened && needsContinuation(moist.deep, pi[i], core.diagnostics.levels, core.diagnostics.sigmaMid)) { continued++; continuedColumns.push(i); }
    const { start, top, base } = moist.deep;
    if (start >= 0) {
      let gap = Infinity;
      for (let k = top + 1; k < base; k++) if (k !== start) gap = Math.min(gap, energy[k] - energy[start]);
      if (gap < 1e-5 * energy[start]) downdraftAtEdge.add(i);
    }
  }
  for (let i = 0; i < C; i++) { const now = moist.convectivePrecipitation[i] + moist.largeScalePrecipitation[i]; response[i] = Math.max(Math.abs(response[i] - now), Math.abs(warmed[i] - now)); }
  for (let x = 0; x < K * C; x++) thetaResponse[x] = Math.abs(thetaResponse[x] - theta[x] * (1 + 2 ** -23));
  if (momentum) {
    const E = mesh.nEdges, { dSigma, g } = core.diagnostics, u = model.state[2], before = Float64Array.from(u);
    model.phases.mixMomentum(0, E, dt);
    let moved = 0, worstU = 0, worstColumn = 0, scale = 0, windOver = -Infinity, windAllowed = 0, besideTied = 0;
    const tied = (i) => downdraftAtEdge.has(i) || threshold.has(i) || typedAtEdge.has(i);
    for (let e = 0; e < E; e++) {
      const columnMass = 0.5 * (pi[mesh.cellsOnEdge[2 * e]] + pi[mesh.cellsOnEdge[2 * e + 1]]), apart = tied(mesh.cellsOnEdge[2 * e]) || tied(mesh.cellsOnEdge[2 * e + 1]);
      if (apart) besideTied++;
      let was = 0, now = 0, changed = false;
      for (let k = 0; k < K; k++) {
        const x = k * E + e, m = columnMass * dSigma[k] / g;
        was += m * before[x]; now += m * u[x]; scale = Math.max(scale, Math.abs(m * before[x]));
        if (u[x] !== before[x]) changed = true;
        if (apart) continue;
        worstU = Math.max(worstU, Math.abs(u[x] - after[2][x]));
        windOver = Math.max(windOver, Math.abs(u[x] - after[2][x]) - (1e-3 + 2 * Math.abs(windResponse[x] - u[x])));
        windAllowed = Math.max(windAllowed, 2 * Math.abs(windResponse[x] - u[x]));
      }
      if (changed) moved++;
      worstColumn = Math.max(worstColumn, Math.abs(now - was));
    }
    console.log(`  momentum: the plume moves the wind on ${moved} of ${E} edges; each edge's column momentum changes by at most ${worstColumn.toExponential(1)} kg/m/s against layer momenta up to ${scale.toFixed(0)}; the engines' winds differ by ${worstU.toExponential(1)} m/s off the ${besideTied} edges beside a column left out below (twice the CPU's response to one f32 ulp of θ at most ${windAllowed.toExponential(1)} m/s)`);
    assert.ok(moved > E / 10, `${moved} edges moved`);
    assert.ok(besideTied <= E / 25, `${besideTied} edges beside columns left out`);
    assert.ok(worstColumn < 1e-12 * scale, `column momentum changed by ${worstColumn}`);
    assert.ok(windOver < 0, `winds differ by ${worstU} m/s, ${windOver} above 1e-3 m/s and twice the response to one ulp of θ`);
  }
  const { sigmaMid } = core.diagnostics;
  let deep = 0, shallow = 0, decked = 0, still = 0, opening = 0, flips = 0, residues = 0, worstTheta = 0, thetaOver = -Infinity, thetaAllowed = 0, worstQ = 0, worstQc = 0, worstRain = 0, rainScale = 0;
  let plumes = 0, plumeFlips = 0, topsDiffer = 0, worstFlux = 0, fluxScale = 0, worstCover = 0, worstWater = 0, waterScale = 0;
  const K0 = K - (layout.PH.CUWATER - layout.PH.CUCOVER) / C;
  const beside = new Set(), switching = new Set();
  for (let i = 0; i < C; i++) {
    const cpu = moist.cumulusBaseFlux[i], gpuFlux = ph.CUMF[i];
    if (moist.cumulusTop[i] > 0 && Math.abs(moist.cumulusTop[i] - ph.CUTOP[i]) <= 1 && Math.abs(cpu - gpuFlux) > 1e-2 * Math.max(cpu, gpuFlux)) beside.add(i);
    if (Math.abs(warmedFlux[i] - cpu) > 1e-3 * cpu) switching.add(i);
    rainScale = Math.max(rainScale, moist.convectivePrecipitation[i]);
  }
  const merged = new Set(), alike = (a, b) => Math.abs(a - b) <= 1e-6 * Math.max(Math.abs(a), Math.abs(b));
  for (let i = 0; i < C; i++) for (let k = 0; k < K - 1; k++) if (alike(model.state[4][k * C + i], model.state[4][(k + 1) * C + i]) !== alike(after[4][k * C + i], after[4][(k + 1) * C + i])) { merged.add(i); break; }
  let rainOver = -Infinity, allowed = 0, fogged = 0, fogFlips = 0, worstFog = 0, fogScale = 0;
  for (let i = 0; i < C; i++) {
    let top = -1;
    for (let k = K - 1; k >= 0; k--) if (moist.trace.convection[k * C + i] !== 0) top = k;
    if (model.radiation.mlmGate[i] >= DECK_CLOSED) decked++;
    else if (top < 0) still++;
    else if (model.radiation.mlmGate[i] > 0.5) opening++;
    else if (pi[i] * sigmaMid[top] > MOIST_DEFAULTS.shallowTop) shallow++;
    else deep++;
    const cpu = moist.cumulusBaseFlux[i], gpuFlux = ph.CUMF[i];
    if (cpu > 0) plumes++;
    if ((cpu > 0) !== (gpuFlux > 0)) plumeFlips++;
    if (moist.cumulusTop[i] > 0 && Math.abs(moist.cumulusTop[i] - ph.CUTOP[i]) > 1) topsDiffer++;
    if (beside.has(i) || merged.has(i) || threshold.has(i) || typedAtEdge.has(i) || downdraftAtEdge.has(i)) continue;
    fluxScale = Math.max(fluxScale, cpu);
    if (!switching.has(i)) worstFlux = Math.max(worstFlux, Math.abs(cpu - gpuFlux));
    for (let k = K0; k < K; k++) {
      worstCover = Math.max(worstCover, Math.abs(moist.cumulusCover[k * C + i] - ph.CUCOVER[(k - K0) * C + i]));
      worstWater = Math.max(worstWater, Math.abs(moist.cumulusWater[k * C + i] - ph.CUWATER[(k - K0) * C + i]));
      waterScale = Math.max(waterScale, moist.cumulusWater[k * C + i]);
    }
    const rainDiff = Math.max(Math.abs(moist.convectivePrecipitation[i] - ph.CONV[i]), Math.abs(moist.largeScalePrecipitation[i] - ph.COND[i]), Math.abs(moist.fogPrecipitation[i] - ph.FOG[i]));
    if (moist.fogPrecipitation[i] > 0) fogged++;
    if ((moist.fogPrecipitation[i] > 0) !== (ph.FOG[i] > 0) && Math.max(moist.fogPrecipitation[i], ph.FOG[i]) > 1e-9) fogFlips++;
    worstFog = Math.max(worstFog, Math.abs(moist.fogPrecipitation[i] - ph.FOG[i])); fogScale = Math.max(fogScale, moist.fogPrecipitation[i]);
    if ((moist.convectivePrecipitation[i] > 0) !== (ph.CONV[i] > 0)) { if (rainDiff > 1e-4 * rainScale + 2 * response[i]) flips++; else residues = Math.max(residues, Math.abs(moist.convectivePrecipitation[i] - ph.CONV[i])); }
    worstRain = Math.max(worstRain, rainDiff);
    rainOver = Math.max(rainOver, rainDiff - (1e-4 * rainScale + 2 * response[i]));
    allowed = Math.max(allowed, 2 * response[i]);
    for (let k = 0; k < K; k++) {
      const x = k * C + i;
      worstTheta = Math.max(worstTheta, Math.abs(model.state[1][x] - after[1][x]));
      thetaOver = Math.max(thetaOver, Math.abs(model.state[1][x] - after[1][x]) - (1e-3 + 2 * thetaResponse[x]));
      thetaAllowed = Math.max(thetaAllowed, 2 * thetaResponse[x]);
      worstQ = Math.max(worstQ, Math.abs(model.state[4][x] - after[4][x]));
      worstQc = Math.max(worstQc, Math.abs(model.state[5][x] - after[5][x]));
    }
  }
  let continuedLeftOut = 0, continuedApart = 0;
  for (const i of continuedColumns) {
    if (beside.has(i) || merged.has(i) || threshold.has(i) || typedAtEdge.has(i) || downdraftAtEdge.has(i)) { continuedLeftOut++; continue; }
    let over = -Infinity;
    for (let k = 0; k < K; k++) over = Math.max(over, Math.abs(model.state[1][k * C + i] - after[1][k * C + i]) - (1e-3 + 2 * thetaResponse[k * C + i]));
    if ((moist.cumulusBaseFlux[i] > 0) !== (ph.CUMF[i] > 0) || Math.abs(moist.cumulusTop[i] - ph.CUTOP[i]) > 1 || Math.abs(moist.cumulusBaseFlux[i] - ph.CUMF[i]) > 2e-3 * fluxScale || over > 0) continuedApart++;
  }
  console.log(`${JSON.stringify(options)}: ${C} random columns: ${ascent ? `${midLevel} of them convect at mid-level, ` : ""}${deep} convect deep, ${shallow} shallow, ${decked} under a deck, ${opening} convect under a deck opening, ${still} still; convective rain differs in sign on ${flips} (beyond the rain's allowance below; within it by at most ${residues.toExponential(1)} kg/m²); engines differ in θ by ${worstTheta.toExponential(1)} K (twice the CPU's response to one f32 ulp of θ at most ${thetaAllowed.toExponential(1)} K), q by ${worstQ.toExponential(1)}, qc by ${worstQc.toExponential(1)}, a step's rain by ${worstRain.toExponential(1)} kg/m² (largest ${rainScale.toFixed(3)}; twice its response to one f32 ulp of the saved subcloud T_v or of θ at most ${allowed.toExponential(1)})`);
  console.log(`  ${plumes} plumes, ${plumeFlips} differ in whether they rise, ${topsDiffer} in their top, ${beside.size} deep in the base flux by more than 1 % (the shallow plume beside it ran on one engine only), ${merged.size} whose dry adjustment merged other layers on the two engines, ${switching.size} whose base flux on the CPU moves by over 1e-3 under one ulp of θ (left out of the flux), ${threshold.size} whose deep plume's condensate at an interface lies within 1e-3 of the conversion threshold, ${typedAtEdge.size} whose test parcel's cloud lies within 1e-3 of 200 hPa or whose w² falls within 1e-3 m²/s² of zero and ${downdraftAtEdge.size} whose downdraft source layer has a moist static energy within 1e-5 of the next lowest (left out); base mass flux differs by ${worstFlux.toExponential(1)} kg/m²/s (largest ${fluxScale.toFixed(3)}), the cumulus fraction by ${worstCover.toExponential(1)}, the plume's condensate by ${worstWater.toExponential(1)} kg/kg (largest ${waterScale.toExponential(1)})`);
  assert.ok(deep > C / 40 && plumes > C / 10 && decked > C / 40 && (ascent || still > C / 40), `${deep} deep, ${plumes} plumes, ${decked} decked, ${still} still`);
  assert.equal(plumeFlips, 0);
  assert.equal(topsDiffer, 0);
  if (nocturnal && options.plumeMomentum !== false) {
    let rooted = 0, apart = 0, fired = 0;
    for (let i = 0; i < C; i++) {
      if (departure[i] >= 0) rooted++;
      if (beside.has(i) || merged.has(i) || threshold.has(i) || typedAtEdge.has(i) || downdraftAtEdge.has(i)) continue;
      if (departure[i] >= 0 && moist.momentum.source[i] < K) fired++;
      if (Math.round(ph.MOMS[i]) !== moist.momentum.source[i]) apart++;
    }
    console.log(`  ${rooted} columns find a departure above the 50 hPa stack on the CPU, ${fired} of them (outside the columns left out) fire a plume there; the engines root the plume in another layer in ${apart} columns`);
    assert.ok(fired >= 50, `${fired} elevated plumes compared`);
    assert.equal(apart, 0, `${apart} columns' plume is rooted in another layer on the GPU`);
  }
  assert.ok(beside.size <= C / 100, `${beside.size} columns' shallow plume beside the deep one ran on one engine only`);
  assert.ok(merged.size <= C / 100, `${merged.size} columns' dry adjustment merged other layers`);
  assert.ok(switching.size <= 0.02 * C, `${switching.size} columns' base flux switches under one ulp of θ`);
  assert.ok(threshold.size <= C / 100 && typedAtEdge.size <= C / 100 && downdraftAtEdge.size <= C / 100, `${threshold.size} columns' plume condensate at the conversion threshold, ${typedAtEdge.size} typed at the edge, ${downdraftAtEdge.size} with a downdraft source at a tie`);
  assert.ok(worstFlux < 2e-3 * fluxScale && worstCover < 1e-4, `base mass flux ${worstFlux}, cover ${worstCover}`);
  assert.ok(waterScale > 0 && worstWater < 1e-3 * waterScale, `plume condensate ${worstWater} against ${waterScale}`);
  assert.equal(flips, 0);
  assert.ok(thetaOver < 0 && worstQ < 2e-6 && worstQc < 5e-7, `θ ${worstTheta} (${thetaOver} above 1e-3 K and twice the response to one ulp), q ${worstQ}, qc ${worstQc}`);
  assert.ok(rainOver < 0, `rain ${worstRain} against ${rainScale}, ${rainOver} above 1e-4 of it and twice the response to one ulp of the saved T_v or of θ`);
  const fogOn = (options.autoconversionFloor ?? MOIST_DEFAULTS.autoconversionFloor) !== 'none' && (('fogDroplets' in options ? options.fogDroplets : MOIST_DEFAULTS.fogDroplets) !== null || (options.fogDeposition ?? MOIST_DEFAULTS.fogDeposition) > 0);
  console.log(`  fog reaches the ground from ${fogged} columns (largest ${(1e3 * fogScale).toFixed(3)} g/m²), the engines differ by ${worstFog.toExponential(1)} kg/m² and on ${fogFlips} in whether more than 1e-9 kg/m² does`);
  assert.equal(fogFlips, 0);
  assert.ok(fogOn ? fogged > C / 10 : fogged === 0, `${fogged} fog columns`);
  const { geopotential, g } = core.diagnostics, ice = model.state[6], shared = new Uint8Array(K * C);
  for (let i = 0; i < C; i++) for (let k = 0; k < K; k++) {
    const below = geopotential[k * C + i] / g < mixingTop[i], iced = ice[i] > 0 ? (concentration[i] > 0 ? concentration[i] : 1) : 0;
    const plumed = moist.cumulusBaseFlux[i] > 0 && pi[i] * sigmaMid[k] >= moist.cumulusTop[i];
    if (!plumed && Math.max(below ? (regime[i] === COUPLED_REGIME ? 1 : 0) : stratiform[i], iced) > 0) shared[k * C + i] = 1;
  }
  return { midLevel, deepened, continued, continuedLeftOut, continuedApart, sentBack, cloud: Float64Array.from(model.state[5]), gpuCloud: Float64Array.from(after[5]), shared, largeScale: Float64Array.from(moist.largeScalePrecipitation), gpuLargeScale: Float64Array.from(ph.COND.subarray(0, C)) };
}

test('the shallow and deep plume and the rain they leave match between the engines on a random set of columns, with the plume under each closure and with its F from the buoyant layers alone, from either source, with either CAPE parcel and its downdraft, carrying momentum with or without a downdraft, the column momentum of each edge exact, with the shallow plume from the lowest layer, raining and overshooting by half without virtual buoyancy, carried through non-buoyant layers by its w² and entraining at 1/z, shut where the EIS of its column passes stratocumulusInversion, under either autoconversion floor and with a shorter lifetime for the cloud above the shallow top', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  await parity({});
  await parity({ plumeClosure: 'maximum', plumeSource: 'lowest', plumeCapeParcel: 'undilute', plumeMassGrowth: 2e-4 });
  await parity({ plumeClosure: 'cape', downdraftShare: 0.5, downdraftEntrainment: 0, plumeRainThreshold: 5e-4, autoconversionFloor: 'boundaryLayer', upperCloudLifetime: 1800 });
  await parity({ plumeConsumption: 'buoyant', cloudLifetime: 7200 });
  await parity({ plumeMomentum: true }, { momentum: true });
  await parity({ plumeMomentum: true, downdraftShare: 0 }, { momentum: true });
  await parity({ autoconversionFloor: 'boundaryLayer' });
  await parity({ autoconversionFloor: 'mixingTop' });
  await parity({ cumulusOvershoot: 1, cumulusAscent: 'buoyant' });
  await parity({ cumulusSource: 'lowest', cumulusRain: 5e-4, cumulusOvershoot: 0.5, cumulusAscent: 'buoyant', virtualBuoyancy: false });
  await parity({ cumulusAscent: 'kinetic', cumulusEntrainmentHeight: 1 });
  await parity({ cumulusAscent: 'kinetic', cumulusSource: 'lowest', cumulusRain: 5e-4 });
  await parity({ stratocumulusInversion: 2 });
  await parity({ stratocumulusInversion: -1e3 });
  await parity({ condensation: 'saturation', iceSaturation: false, iceFall: null });
  await parity({ iceNucleation: true, iceFall: 3.29 });
  await parity({ capeClosure: 'threshold' });
  await parity({ pcapeBoundary: 'signed' });
  await parity({ pcapeMinimum: 0 });
  await parity({ pcapeMinimumOver: 'land', pcapeTimescale: 'sea' }, { land: true });
  await parity({ pcapeTimescale: 'sea' }, { land: true });
  await parity({ excessVelocity: 'convective', ...PREVIOUS_PCAPE_MINIMUM });
  await parity({ plumeSourceDepth: 'boundaryLayer' });
  await parity({ convectionType: 'top' });
  await parity({ convectionType: 'cloudDepth' });
  await parity({ plumePhase: 'liquid' });
  await parity({ plumeConversion: 'zhangMcFarlane' });
  await parity({ plumeEntrainmentLaw: 'gregory' });
});

test("the departure search matches between the engines on random columns with a nocturnal inversion under a conditionally unstable layer: the same layer roots the plume, under 'elevated' and 'elevatedStable', carrying momentum, and with the test parcel's entrainment scaled by the departure layer's own humidity", { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { stratocumulusInversion: null };
  await parity({ ...options, plumeDeparture: 'elevated' }, { nocturnal: true });
  await parity({ ...options, plumeDeparture: 'elevatedStable' }, { nocturnal: true });
  await parity({ ...options, plumeDeparture: 'elevated', plumeMomentum: true }, { momentum: true, nocturnal: true });
  const kept = DEPARTURE.reference;
  DEPARTURE.reference = 'departure';
  try { await parity({ ...options, plumeDeparture: 'elevated' }, { nocturnal: true }); } finally { DEPARTURE.reference = kept; }
});

test('the fog scheme of the floored layers matches between the engines: settling, deposition on the implicit drag and drizzle, under other droplet numbers, under the boundary-layer floor, off, without a floor and as liquid alone without falling ice', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  await parity({ fogDroplets: null, fogDeposition: 0 });
  await parity({ fogDroplets: [50, 300] });
  await parity({ autoconversionFloor: 'none' });
  await parity({ iceFall: null, fogDroplets: [150, 60] });
});

test('the stratiform lifetime matches between the engines on random columns of every regime, mixing top, EIS share and sea-ice cover, and keeps cloud the short lifetime would rain out in the layers its rule gives a long share', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const long = await parity({}), short = await parity({ stratiformLifetime: null });
  let cloudy = 0, kept = 0, gpuKept = 0, shared = 0, keptShared = 0;
  for (let x = 0; x < long.cloud.length; x++) {
    if (!(short.cloud[x] > 0)) continue;
    cloudy++;
    const keeps = long.cloud[x] > short.cloud[x] * (1 + 1e-6);
    if (keeps) kept++;
    if (long.gpuCloud[x] > short.gpuCloud[x] * (1 + 1e-6)) gpuKept++;
    if (long.shared[x]) { shared++; if (keeps) keptShared++; }
  }
  console.log(`after one step the 3 h stratiform lifetime keeps more cloud than the 1 h lifetime alone in ${kept} of ${cloudy} cloudy layers (GPU ${gpuKept}), ${keptShared} of the ${shared} whose long share is positive (outside a plume's layers, below the mixing top of a coupled or ice-covered column or above it under an EIS share)`);
  assert.ok(keptShared > shared / 5 && kept < cloudy && Math.abs(gpuKept - kept) <= cloudy / 100, `${kept} and ${gpuKept} of ${cloudy}, ${keptShared} of the ${shared} with a long share`);
});

test("the Kessler rate of the stratiform rain's evaporation, the IFS code's semi-implicit form of it, the overlap of the rain's precipitating share and the relative-humidity ceiling match between the engines on random columns, under every law, with the share overlapped or held, and the rate lets more of the large-scale rain reach the ground on both", { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const instant = await parity({}), kessler = await parity({ rainEvaporationLaw: 'kessler' });
  const ceiled = await parity({ rainEvaporationLaw: 'kessler', rainCover: 0.5, rainEvaporationCeiling: 0.7 });
  const cloudsc = await parity({ rainEvaporationLaw: 'cloudsc', rainEvaporationCeiling: 0.7 });
  await parity({ rainEvaporationCeiling: 0.7 });
  await parity({ rainEvaporationLaw: 'cloudsc', rainCover: 0.1 });
  const sum = (a) => a.reduce((x, y) => x + y, 0);
  let wetter = 0, gpuWetter = 0;
  for (let i = 0; i < instant.largeScale.length; i++) {
    if (kessler.largeScale[i] > instant.largeScale[i] * (1 + 1e-6) + 1e-9) wetter++;
    if (kessler.gpuLargeScale[i] > instant.gpuLargeScale[i] * (1 + 1e-6) + 1e-9) gpuWetter++;
  }
  console.log(`large-scale rain reaching the ground over the random columns: ${sum(instant.largeScale).toFixed(3)} kg/m² summed under the instant law (GPU ${sum(instant.gpuLargeScale).toFixed(3)}), ${sum(kessler.largeScale).toFixed(3)} under the Kessler rate (GPU ${sum(kessler.gpuLargeScale).toFixed(3)}), more in ${wetter} columns (GPU ${gpuWetter}); ${sum(ceiled.largeScale).toFixed(3)} (GPU ${sum(ceiled.gpuLargeScale).toFixed(3)}) with a_P 0.5 under the 0.85 ceiling; ${sum(cloudsc.largeScale).toFixed(3)} (GPU ${sum(cloudsc.gpuLargeScale).toFixed(3)}) under the semi-implicit form, the share overlapped, and the 0.7 ceiling`);
  assert.ok(wetter > 0 && Math.abs(gpuWetter - wetter) <= Math.max(2, wetter / 100), `${wetter} columns wetter on the CPU, ${gpuWetter} on the GPU`);
  assert.ok(sum(kessler.largeScale) > sum(instant.largeScale) && sum(kessler.gpuLargeScale) > sum(instant.gpuLargeScale));
  assert.ok(sum(cloudsc.largeScale) > sum(instant.largeScale) && sum(cloudsc.gpuLargeScale) > sum(instant.gpuLargeScale));
});

test("with convectionRetype a column its test parcel types deep runs the shallow plume where the deep plume's own cloud, from its base to the last interface it reaches with w² > 0, is no deeper than 200 hPa, as the IFS's code re-types after the ascent, and any other value than true or false is refused", () => {
  const counts = {};
  for (const retype of [false, true]) {
    const model = build({ convectionRetype: retype }), { moist, core } = model, C = model.mesh.nCells, { levels } = core.diagnostics, pi = model.state[0];
    randomColumns(model);
    let deep = 0, shallowCloud = 0;
    for (let i = 0; i < C; i++) {
      model.phases.adjust(i, i + 1, 900);
      if (moist.deep.deep) { deep++; if (!(pi[i] * (levels[moist.deep.base] - levels[moist.deep.top + 1]) > DEEP_CLOUD_DEPTH)) shallowCloud++; }
    }
    counts[retype] = { deep, shallowCloud };
  }
  console.log(`random columns: ${counts.false.deep} run the deep plume as typed by the test parcel, ${counts.false.shallowCloud} of them with a plume cloud no deeper than 200 hPa; re-typed after the ascent ${counts.true.deep} do, ${counts.true.shallowCloud} with such a cloud`);
  assert.ok(counts.false.shallowCloud > 0 && counts.true.shallowCloud === 0, JSON.stringify(counts));
  assert.equal(counts.true.deep, counts.false.deep - counts.false.shallowCloud);
  assert.throws(() => build({ convectionRetype: 'yes' }), /convectionRetype/);
});

test('with convectionRetype the engines agree on random columns, and the GPU refuses a value other than true or false', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  await parity({ convectionRetype: true });
  const { physicsConstants } = await import('../js/gpu/physics.gpu.js');
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  assert.throws(() => physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, convectionRetype: 1 }), /convectionRetype/);
});

test('with midLevelConvection a column that no test parcel types deep and whose shallow plume does not rise lifts the lowest layer above the lowest below 10 km in large-scale ascent with q above 0.8 q_s whose parcel condenses at its top interface and is less than 0.5 K colder in virtual temperature than the half-level air CUININ builds there, at a base flux of −π σ̇/g or (H + L E)/(g z_base), with enthalpy and water exact; sinking, drier air or the option off leave the column as without it', () => {
  const dt = 600;
  const run = (options, { ascent = -0.2, humid = 0.95, sensible = 0, evaporation = 0, lapse = 6.5e-3, low = 1000, high = 3000, profile = null, column = null, buoyancy = -1e-4, depth = 500, ground = null, grid = levels } = {}) => {
    const model = createModel(new Grid(16), { ...(ground === null ? {} : { topography: syntheticTopography(90, 180, () => ground) }), ocean: false, levels: grid, moist: options }), { moist, core } = model, { K, C, dSigma, geopotential, g, piSigmaDot } = core.diagnostics, [pi] = model.state;
    place(model, 0, 101500, column ?? ((z, p) => { const T = 290 - lapse * Math.min(z, 11000); return { T, q: (z > low && z < high ? humid : 0.3) * saturationHumidity(T, p) }; }));
    setDepth(model, 0, depth);
    model.boundaryLayer.buoyancyFlux[0] = buoyancy; model.boundaryLayer.friction[0] = 0.25; model.radiation.mlmGate[0] = 0.3;
    model.radiation.sensibleHeat[0] = sensible; model.radiation.evaporation[0] = evaporation;
    for (let k = 1; k < K; k++) piSigmaDot[k * C] = ascent;
    if (core.diagnostics.omega) for (let k = 0; k < K; k++) core.diagnostics.omega[k * C] = profile ? profile(k) : ascent;
    const { cp, exnerLayer, exnerLower } = core.diagnostics, thetaV = core.arrays.thetaV, b = (K - 1) * C;
    const interfaces = Float64Array.from({ length: K }, (_, k) => (k > 0 ? (geopotential[k * C] + cp * thetaV[k * C] * (exnerLayer[k * C] - exnerLower[(k - 1) * C])) / g : Infinity)), groundHeight = (geopotential[b] - cp * thetaV[b] * (exnerLower[b] - exnerLayer[b])) / g;
    const before = budget(model, 0);
    moist.adjust(model.state, 0, 1, dt);
    const after = budget(model, 0), source = Array.from({ length: K }, (_, k) => k).find((k) => Math.abs(pi[0] * dSigma[k] - moist.deep.sourceMass) < 1e-12 * moist.deep.sourceMass);
    return { deep: { ...moist.deep }, rain: moist.rain[0], entrained: Float64Array.from(moist.plumeEntrained), detrained: Float64Array.from(moist.plumeDetrained), humidity: Float64Array.from(moist.envHumidity), flux: Float64Array.from(moist.cumulusFlux), before, after, state: model.state.map((a) => Float64Array.from(a.subarray(0, a.length))), height: source === undefined ? NaN : geopotential[source * C] / g - geopotential[(K - 1) * C] / g, source, interfaces, groundHeight, z: Float64Array.from({ length: K }, (_, k) => geopotential[k * C] / g), g };
  };
  const same = (a, b) => a.state.every((x, n) => x.every((v, j) => v === b.state[n][j]));
  const off = run({}), on = run({ midLevelConvection: true }), g = on.g;
  console.log(`mid-level: base flux ${on.deep.baseFlux.toExponential(6)} kg/m²/s against 0.2/g ${(0.2 / g).toExponential(6)}, from the layer ${on.height.toFixed(0)} m above the lowest, top ${on.deep.top}, base ${on.deep.base}, rain ${on.rain.toExponential(3)} kg/m²`);
  assert.ok(!off.deep.deep && !off.deep.midLevel && off.deep.baseFlux === 0, 'nothing convects without the option');
  assert.ok(on.deep.deep && on.deep.midLevel && Math.abs(on.deep.baseFlux * g / 0.2 - 1) < 1e-14, `base flux ${on.deep.baseFlux}`);
  assert.ok(on.height > 800 && on.height < 3000 && on.deep.cape > 0 && on.rain > 0, `the plume departs ${on.height} m above the lowest layer with CAPE ${on.deep.cape} J/kg and rain ${on.rain} kg/m²`);
  assert.ok(on.flux[on.deep.base] === 1 && on.flux[on.deep.base - 1] === 1, `the mid-level plume keeps its whole base flux through its first cloud layer: ${on.flux[on.deep.base - 1]}`);
  const lifted = run({ midLevelConvection: true }, { profile: (k) => (k === on.deep.base - 1 ? 0.45 : -0.2) });
  assert.ok(lifted.deep.midLevel && lifted.source === on.source && Math.abs(lifted.deep.baseFlux * g / 0.2 - 1) < 1e-14, `descent at 0.45 Pa/s in the layer above the cloud base leaves the plume in its source layer at -omega/g of that layer: base flux ${lifted.deep.baseFlux}, from ${lifted.height} m against ${on.height} m`);
  const faint = { profile: (k) => (k === on.source ? -1e-7 : -0.2) }, weak = run({ midLevelConvection: true }, faint), weakOff = run({}, faint);
  assert.ok(!weak.deep.midLevel && same(weak, weakOff), `a source layer ascending at 1e-7 Pa/s, which closes the plume below CUMULUS_FLOOR, gets no plume and ends the search, though the layer above ascends at 0.2 Pa/s: base flux ${weak.deep.baseFlux}, from ${weak.height} m`);
  const layered = run({ midLevelConvection: true }, { profile: (k) => (k === on.source ? -0.2 : Math.abs(k - on.source) === 1 ? 0.2 : -0.3) });
  assert.ok(layered.deep.midLevel && layered.source === on.source && Math.abs(layered.deep.baseFlux * g / 0.2 - 1) < 1e-14, `with the source layer ascending at 0.2 Pa/s between layers that descend at 0.2 Pa/s, in ascent of 0.3 Pa/s elsewhere, the plume departs from the source layer (${layered.height} m against ${on.height} m) at -omega/g of that layer: ${layered.deep.baseFlux}`);
  const lee = run({ midLevelConvection: true }, { profile: () => 0.2 }), windward = run({ midLevelConvection: true }, { ascent: 0.2, profile: () => -0.2 });
  assert.ok(!lee.deep.midLevel && same(lee, off), 'no mid-level plume where omega > 0 though pi sigma-dot < 0');
  const single = (rh, lapse = 8e-3) => ({ column: (z, p, memo) => { const T = 290 - lapse * Math.min(z, 11000), picked = !memo.picked && z > 2800; if (picked) memo.picked = true; return { T, q: (picked ? rh : 0.3) * saturationHumidity(T, p) }; }, grid: sigmaInterfaces('bl48') });
  const shallowLcl = run({ midLevelConvection: true }, single(0.805)), shallowLclOff = run({}, single(0.805)), wet = run({ midLevelConvection: true }, single(0.97));
  assert.ok(!shallowLcl.deep.midLevel && same(shallowLcl, shallowLclOff), `a bl48 layer near 3 km at RH 0.805, whose parcel condenses above its top interface though below the next layer's midpoint, gets no mid-level plume: base flux ${shallowLcl.deep.baseFlux} from ${shallowLcl.height} m`);
  assert.ok(wet.deep.midLevel && wet.height > 2000, `the same layer at RH 0.97 departs: base flux ${wet.deep.baseFlux} from ${wet.height} m`);
  const inversion = { column: (z, p) => { const T = z < 800 ? 290 - 6.5e-3 * z : z < 1400 ? 284.8 + 5 * (z - 800) / 600 : 289.8 - 6.5e-3 * (z - 1400); return { T, q: (z > 800 && z < 2400 ? 0.95 : 0.3) * saturationHumidity(T, p) }; } };
  const retried = run({ midLevelConvection: true }, inversion);
  assert.ok(retried.deep.midLevel && retried.height > 1300 && retried.deep.cape > 0, `under a 5 K inversion at 0.8-1.4 km, moist from 0.8 to 2.4 km, the inversion's layers give way to the first one above it: from ${retried.height} m with CAPE ${retried.deep.cape} J/kg`);
  const tooCold = run({ midLevelConvection: true }, single(0.9, 4e-3)), tooColdOff = run({}, single(0.9, 4e-3)), justWarm = run({ midLevelConvection: true }, single(0.9, 5e-3));
  assert.ok(!tooCold.deep.midLevel && same(tooCold, tooColdOff) && justWarm.deep.midLevel, `a bl48 layer at RH 0.9 near 3 km whose parcel is 0.70 K colder at its top interface under 4 K/km gets no mid-level plume, and under 5 K/km, 0.45 K colder, departs: from ${justWarm.height} m`);
  const virtual = run({ midLevelConvection: true }, single(0.99, 3.5e-3)), loaded = run({ midLevelConvection: true }, single(0.97, 1.2e-3));
  assert.ok(virtual.deep.midLevel && loaded.deep.midLevel, `a parcel 0.07 K warmer in virtual temperature though 0.81 K colder in temperature departs (from ${virtual.height} m), and so does one 0.455 K colder that holds 0.10 K of condensate, unloaded as cuascn's first level is (from ${loaded.height} m)`);
  const lowestOnly = { column: (z, p) => { const T = z < 300 ? 290 - 6.5e-3 * z : z < 700 ? 288.05 + 4e-3 * (z - 300) : 289.65 - 6.5e-3 * (z - 700); return { T, q: (z < 60 ? 0.995 : 0.3) * saturationHumidity(T, p) }; } };
  const lowest = run({ midLevelConvection: true }, lowestOnly), lowestOff = run({}, lowestOnly);
  assert.ok(!lowest.deep.midLevel && same(lowest, lowestOff), `a column moist only in its lowest layer, under an inversion at 300-700 m, gets no mid-level plume: the search starts at the second-lowest layer (base flux ${lowest.deep.baseFlux})`);
  assert.ok(windward.deep.midLevel && Math.abs(windward.deep.baseFlux * g / 0.2 - 1) < 1e-14, `where omega < 0 though pi sigma-dot > 0 the plume departs at -omega/g: ${windward.deep.baseFlux}`);
  const mixing = Array.from({ length: on.deep.base - on.deep.top }, (_, n) => on.deep.top + n).filter((k) => on.entrained[k] > 0);
  const spread = Math.max(...mixing.map((k) => on.humidity[k])) - Math.min(...mixing.map((k) => on.humidity[k]));
  assert.ok(mixing.length >= 3 && spread > 0.05 && mixing.every((k) => Math.abs(on.detrained[k] / (on.entrained[k] * (1.6 - on.humidity[k])) - 1) < 1e-12), `the mid-level plume detrains (1.6 - RH) times what it entrains in each of its ${mixing.length} entraining layers: ratios ${mixing.map((k) => (on.detrained[k] / on.entrained[k]).toFixed(4)).join(' ')} against 1.6 - RH ${mixing.map((k) => (1.6 - on.humidity[k]).toFixed(4)).join(' ')}`);
  const deckColumn = (z, p) => {
    const T = (z < 900 ? 290 : 304 + 5e-3 * (z - 900)) * Math.pow(p / 1e5, 0.2857);
    return { T, q: z < 900 ? Math.min(0.0105, saturationHumidity(T, p)) : Math.min(0.003, 0.3 * saturationHumidity(T, p)) };
  };
  const deck = { column: deckColumn, buoyancy: 2e-4, depth: 900, ascent: -0.05, sensible: 10, evaporation: 4e-5 };
  const capped = run({ midLevelConvection: true }, deck), cappedOff = run({}, deck);
  assert.ok(!capped.deep.midLevel && same(capped, cappedOff), `a column under a deck that the stratocumulus cap closes gets no mid-level plume: base flux ${capped.deep.baseFlux} from ${capped.height} m`);
  const noCap = run({ midLevelConvection: true, stratocumulusInversion: null }, deck), noCapOff = run({ stratocumulusInversion: null }, deck);
  assert.ok(!noCap.deep.midLevel && same(noCap, noCapOff), `without the stratocumulus cap the deck column's shallow plume rises and the mid-level type is not tried: base flux ${noCap.deep.baseFlux}`);
  const lidColumn = (z, p) => {
    const T = z < 2200 ? 290 - 6.5e-3 * z : z < 2700 ? 275.7 + 0.01 * (z - 2200) : 280.7 - 6.5e-3 * (Math.min(z, 11000) - 2700);
    return { T, q: (z > 1000 && z < 3000 ? 0.95 : 0.3) * saturationHumidity(T, p) };
  };
  const lid = run({ midLevelConvection: true }, { column: lidColumn }), lidRetyped = run({ midLevelConvection: true, convectionRetype: true }, { column: lidColumn });
  const lidDepth = lid.state[0][0] * (levels[lid.deep.base] - levels[lid.deep.top + 1]);
  assert.ok(lid.deep.midLevel && lidDepth < DEEP_CLOUD_DEPTH && lidRetyped.deep.midLevel && same(lidRetyped, lid), `convectionRetype leaves a mid-level plume ${lidDepth} Pa deep under an inversion as it is`);
  const gradedOmega = run({ midLevelConvection: true }, { profile: (k) => -0.05 - 0.005 * k });
  assert.ok(gradedOmega.deep.midLevel && Math.abs(gradedOmega.deep.baseFlux * g / (0.05 + 0.005 * gradedOmega.source) - 1) < 1e-14, `under omega graded in the vertical the plume departs at -omega/g of its source layer ${gradedOmega.source}: ${gradedOmega.deep.baseFlux * g} Pa/s`);
  const rateFactor = MID_LEVEL.entrainment;
  let singleRate;
  try { MID_LEVEL.entrainment = 1; singleRate = run({ midLevelConvection: true }); } finally { MID_LEVEL.entrainment = rateFactor; }
  const bothMixing = mixing.filter((k) => singleRate.entrained[k] > 0);
  assert.ok(bothMixing.length >= 3 && bothMixing.every((k) => on.entrained[k] === 2 * singleRate.entrained[k]), `the mid-level plume entrains twice eq. 6.7's rate in each of its ${bothMixing.length} entraining layers: ${bothMixing.map((k) => on.entrained[k] / singleRate.entrained[k]).join(' ')}`);
  const pressureColumn = (z, p) => { const T = 290 * Math.pow(p / 101500, 0.19); return { T, q: (p < 90000 && p > 75000 ? 0.95 : 0.3) * saturationHumidity(T, p) }; };
  const driven = { column: pressureColumn, ascent: -1e-4, sensible: 20, evaporation: 2e-5 }, capHeight = MID_LEVEL.height;
  const groundLow = run({ midLevelConvection: true }, { ...driven, ground: 300 }), departed = groundLow.source, aboveGround = (k) => groundLow.z[k] - 300;
  let capAbove, capBelow, groundHigh;
  try {
    MID_LEVEL.height = 0.5 * (aboveGround(departed) + aboveGround(departed - 1));
    capAbove = run({ midLevelConvection: true }, { ...driven, ground: 300 });
    groundHigh = run({ midLevelConvection: true }, { ...driven, ground: 1800 });
    MID_LEVEL.height = 0.5 * (aboveGround(departed) + aboveGround(departed + 1));
    capBelow = run({ midLevelConvection: true }, { ...driven, ground: 300 });
  } finally { MID_LEVEL.height = capHeight; }
  assert.ok(groundLow.deep.midLevel && capAbove.deep.midLevel && same(capAbove, groundLow) && !capBelow.deep.midLevel, `on ground 300 m high the search stops at MID_LEVEL.height above the ground: set between the departure layer (${aboveGround(departed)} m) and the next one up the plume departs as before, set between it and the next one down there is none`);
  assert.ok(groundHigh.deep.midLevel && groundHigh.source === departed && Math.abs(groundHigh.deep.baseFlux / groundLow.deep.baseFlux - 1) < 1e-12, `on ground 1500 m higher, under that bound, the same column departs from the same layer at the same surface-driven flux: ${groundHigh.deep.baseFlux} against ${groundLow.deep.baseFlux}`);
  assert.ok(!same(on, off), 'the mid-level plume moves the column');
  assert.ok(Math.abs(on.after.enthalpy - on.before.enthalpy) < 1e-12 * on.before.enthalpy, `enthalpy ${on.before.enthalpy} → ${on.after.enthalpy}`);
  assert.ok(Math.abs(on.after.water + on.rain - on.before.water) < 1e-12 * on.before.water, `water ${on.before.water} → ${on.after.water} + ${on.rain}`);
  const sinking = run({ midLevelConvection: true }, { ascent: 0.2 }), drier = run({ midLevelConvection: true }, { humid: 0.75 }), drierOff = run({}, { humid: 0.75 });
  assert.ok(!sinking.deep.midLevel && same(sinking, off), 'no mid-level plume under large-scale sinking');
  assert.ok(!drier.deep.midLevel && same(drier, drierOff), 'no mid-level plume where q stays below 0.8 q_s');
  for (const [name, column] of [['a moist layer at 3.5-4.5 km under 4.5 K/km', { lapse: 4.5e-3, low: 3500, high: 4500 }], ['a moist layer below 300 m under 4.5 K/km', { lapse: 4.5e-3, low: -Infinity, high: 300 }]]) {
    const cold = run({ midLevelConvection: true }, column), coldOff = run({}, column);
    assert.ok(!cold.deep.midLevel && same(cold, coldOff), `${name}: no layer's parcel both condenses at the next layer and is less than 0.5 K colder there`);
  }
  const surface = [1, 2].map((n) => run({ midLevelConvection: true }, { ascent: -1e-4, sensible: 10 * n, evaporation: 1e-5 * n }));
  console.log(`mid-level under 1e-4 Pa/s of ascent: base flux ${surface[0].deep.baseFlux.toExponential(6)} at H + L E = ${(10 + 1e-5 * LATENT_HEAT).toFixed(0)} W/m², ${surface[1].deep.baseFlux.toExponential(6)} at twice that`);
  assert.ok(surface[0].deep.midLevel && surface[0].deep.baseFlux > 100 * 1e-4 / g, `base flux ${surface[0].deep.baseFlux}`);
  assert.ok(Math.abs(surface[1].deep.baseFlux / surface[0].deep.baseFlux - 2) < 1e-12, `${surface[1].deep.baseFlux} against twice ${surface[0].deep.baseFlux}`);
  const zBase = surface[0].interfaces[surface[0].source] - surface[0].groundHeight;
  assert.ok(Math.abs(surface[0].deep.baseFlux * g * zBase / (10 + 1e-5 * LATENT_HEAT) - 1) < 1e-12, `the surface closure divides H + L E by g times the height of the source layer's upper interface, ${zBase} m: base flux ${surface[0].deep.baseFlux}`);
  assert.throws(() => build({ midLevelConvection: 'yes' }), /midLevelConvection/);
});

test('with midLevelConvection the engines agree on random columns under a random large-scale ascent, and the GPU refuses a value other than true or false', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { midLevel } = await parity({ midLevelConvection: true }, { ascent: true });
  assert.ok(midLevel >= 50, `${midLevel} columns convect at mid-level on the CPU`);
  const { physicsConstants } = await import('../js/gpu/physics.gpu.js');
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  assert.throws(() => physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, midLevelConvection: 1 }), /midLevelConvection/);
});
test('with midLevelConvection and convectionRetype the engines agree on random columns under a random large-scale ascent', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { midLevel } = await parity({ midLevelConvection: true, convectionRetype: true }, { ascent: true });
  assert.ok(midLevel >= 50, `${midLevel} columns convect at mid-level on the CPU`);
});
test('with midLevelConvection the engines agree on random columns with a deck that the stratocumulus cap closes in every fifth', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const deck = (z, p) => {
    const T = (z < 900 ? 290 : 304 + 5e-3 * (z - 900)) * Math.pow(p / 1e5, 0.2857);
    return { T, q: z < 900 ? Math.min(0.0105, saturationHumidity(T, p)) : Math.min(0.003, 0.3 * saturationHumidity(T, p)) };
  };
  const decks = (model, fields) => {
    for (let i = 0; i < model.mesh.nCells; i += 5) {
      place(model, i, 101500, deck); setDepth(model, i, 900);
      model.boundaryLayer.depth[i] = Math.fround(model.boundaryLayer.depth[i]); model.radiation.mlmGate[i] = Math.fround(0.3);
      fields.buoyancy[i] = model.boundaryLayer.buoyancyFlux[i] = Math.fround(2e-4);
      fields.friction[i] = model.boundaryLayer.friction[i] = 0.25;
      fields.sensible[i] = model.radiation.sensibleHeat[i] = 10;
      fields.evaporation[i] = model.radiation.evaporation[i] = Math.fround(4e-5);
    }
  };
  const { midLevel } = await parity({ midLevelConvection: true }, { ascent: true, columns: decks });
  console.log('mid-level with decks: ' + midLevel + ' columns');
  assert.ok(midLevel >= 30, `${midLevel} columns convect at mid-level on the CPU`);
});

test("a rainEvaporationLaw other than 'deficit', 'kessler' or 'cloudsc', a rainCover other than 'overlap' or a share in (0, 1] and a rainEvaporationCeiling other than null or a relative humidity from 0 to 1 are refused on both engines", async () => {
  const refused = [{ rainEvaporationLaw: 'ifs' }, { rainCover: 'maximum' }, { rainCover: 0 }, { rainCover: 1.5 }, { rainCover: '0.1' }, { rainEvaporationCeiling: -0.1 }, { rainEvaporationCeiling: 1.2 }, { rainEvaporationCeiling: '0.7' }];
  for (const options of refused) assert.throws(() => build(options), /rainEvaporationLaw|rainCover|rainEvaporationCeiling/, `CPU: ${JSON.stringify(options)}`);
  if (!gpuAvailable) return;
  const { physicsConstants } = await import('../js/gpu/physics.gpu.js');
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  for (const options of refused) assert.throws(() => physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, ...options }), /rainEvaporationLaw|rainCover|rainEvaporationCeiling/, `GPU: ${JSON.stringify(options)}`);
});

test('a cumulusMemory other than a finite time of 0 s or more is refused on both engines', async () => {
  const refused = [-1, null, Infinity, '1800'];
  for (const value of refused) assert.throws(() => build({ cumulusMemory: value }), /cumulusMemory/, `CPU: ${value}`);
  if (!gpuAvailable) return;
  const { physicsConstants } = await import('../js/gpu/physics.gpu.js');
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  for (const value of refused) assert.throws(() => physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, cumulusMemory: value }), /cumulusMemory/, `GPU: ${value}`);
});

test("a cumulusAscent other than 'buoyant' or 'kinetic', the kinetic ascent without the kinetic overshoot and a cumulusEntrainmentHeight other than null or a factor above 0 are refused on both engines", async () => {
  const refused = [{ cumulusAscent: 'carried' }, { cumulusAscent: 'kinetic', cumulusOvershoot: 1 }, { cumulusEntrainmentHeight: 0 }, { cumulusEntrainmentHeight: -1 }, { cumulusEntrainmentHeight: Infinity }, { cumulusEntrainmentHeight: '1' }];
  for (const options of refused) assert.throws(() => build(options), /cumulusAscent|cumulusEntrainmentHeight/, `CPU: ${JSON.stringify(options)}`);
  if (!gpuAvailable) return;
  const { physicsConstants } = await import('../js/gpu/physics.gpu.js');
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  for (const options of refused) assert.throws(() => physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, ...options }), /cumulusAscent|cumulusEntrainmentHeight/, `GPU: ${JSON.stringify(options)}`);
});

test("stratocumulusInversion shuts the shallow plume, never the deep one, where the boundary layer reaches the lowest layer's LCL and the column's estimated inversion strength exceeds it: above every column's EIS the step is the default's; below every one, or at 2 K, every deep plume runs as by default and every column whose shallow plume runs steps as by default; and below every one with every boundary layer under its lowest layer's midpoint the step is the default's; all of it on a sea-level world and on one raised 1500 m", () => {
  const step = (options, { below = false, topography = null } = {}) => {
    const model = createModel(new Grid(16), { ocean: false, levels, moist: options, topography }), C = model.mesh.nCells, { deep } = model.moist, deepFlux = new Float64Array(C), deepTop = new Int32Array(C);
    randomColumns(model);
    if (below) {
      const [pi, theta, , , q, qc] = model.state, { geopotential, g, K } = model.core.diagnostics;
      for (let i = 0; i < C; i++) { model.core.diagnoseColumn(i, pi, theta, q, qc); model.boundaryLayer.depth[i] = geopotential[(K - 1) * C + i] / g - 1; }
    }
    for (let i = 0; i < C; i++) {
      deep.deep = false;
      model.phases.adjust(i, i + 1, 900);
      deepFlux[i] = deep.deep ? deep.baseFlux : 0; deepTop[i] = deep.deep ? deep.top : -1;
    }
    const flux = Float64Array.from(model.moist.cumulusBaseFlux);
    return { C, deepFlux, deepTop, flux, shallow: flux.map((f, i) => f - deepFlux[i]), rain: Float64Array.from(model.moist.convectivePrecipitation), state: model.state.map((a) => Float64Array.from(a)) };
  };
  const same = (run, base, i, label) => {
    assert.equal(run.flux[i], base.flux[i], `${label}: column ${i}'s base flux`);
    assert.equal(run.rain[i], base.rain[i], `${label}: column ${i}'s convective rain`);
    run.state.forEach((field, n) => { if (n !== 2 && field.length % run.C === 0) for (let x = i; x < field.length; x += run.C) assert.equal(field[x], base.state[n][x], `${label}: column ${i}, field ${n}`); });
  };
  const RAISED = syntheticTopography(90, 180, () => 1500);
  const base = step({ stratocumulusInversion: null }), every = step({ stratocumulusInversion: 1e3 }), C = base.C;
  assert.deepEqual(every.flux, base.flux);
  for (const [n, field] of every.state.entries()) assert.deepEqual(field, base.state[n], `field ${n} under the highest threshold`);
  const plumes = (values) => values.reduce((n, x) => n + (x > 0 ? 1 : 0), 0), counts = [];
  const raised = step({ stratocumulusInversion: null }, { topography: RAISED });
  for (const [label, threshold, world, control] of [['-1000 K', -1e3, null, base], ['2 K', 2, null, base], ['-1000 K, 1500 m up', -1e3, RAISED, raised]]) {
    const run = step({ stratocumulusInversion: threshold }, { topography: world });
    assert.deepEqual(run.deepFlux, control.deepFlux, `the deep base flux under ${label}`);
    assert.deepEqual(run.deepTop, control.deepTop, `the deep plume's top under ${label}`);
    let shut = 0;
    for (let i = 0; i < C; i++) {
      if (run.shallow[i] === control.shallow[i]) same(run, control, i, label);
      else { assert.equal(run.shallow[i], 0, `${label}: column ${i}'s shallow flux`); shut++; }
    }
    counts.push([label, shut, plumes(run.shallow), world ? C / 100 : C / 40]);
  }
  const low = step({ stratocumulusInversion: null }, { below: true }), lowShut = step({ stratocumulusInversion: -1e3 }, { below: true });
  for (let i = 0; i < C; i++) same(lowShut, low, i, 'below every LCL');
  const raisedLow = step({ stratocumulusInversion: null }, { below: true, topography: RAISED }), raisedLowShut = step({ stratocumulusInversion: -1e3 }, { below: true, topography: RAISED });
  for (let i = 0; i < C; i++) same(raisedLowShut, raisedLow, i, 'below every LCL, 1500 m up');
  console.log(`${C} random columns: ${plumes(base.deepFlux)} deep and ${plumes(base.shallow)} shallow plumes by default; ${counts.map(([label, shut, kept]) => `under ${label} ${shut} shallow plumes shut and ${kept} kept`).join('; ')}; ${plumes(low.shallow)} shallow plumes with every boundary layer under its lowest layer's midpoint, ${plumes(raisedLow.shallow)} 1500 m up`);
  assert.ok(plumes(base.deepFlux) > C / 40 && plumes(base.shallow) > C / 10, 'deep and shallow plumes by default');
  for (const [label, shut, kept, least] of counts) assert.ok(shut > least && kept > least, `${shut} shut, ${kept} kept under ${label}`);
  assert.ok(plumes(low.shallow) > C / 40 && plumes(raisedLow.shallow) > C / 40, "shallow plumes with every boundary layer under its lowest layer's midpoint");
});

test('a stratocumulusInversion other than null or a finite EIS is refused on both engines', async () => {
  const refused = [{ stratocumulusInversion: NaN }, { stratocumulusInversion: Infinity }, { stratocumulusInversion: '7' }];
  for (const options of refused) assert.throws(() => build(options), /stratocumulusInversion/, `CPU: ${JSON.stringify(options)}`);
  if (!gpuAvailable) return;
  const { physicsConstants } = await import('../js/gpu/physics.gpu.js');
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  for (const options of refused) assert.throws(() => physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, ...options }), /stratocumulusInversion/, `GPU: ${JSON.stringify(options)}`);
});

test('the retired Betts–Miller options are refused on both engines', async () => {
  for (const options of [{ convection: 'bettsMiller' }, { convection: 'plume' }, { shallowScheme: 'massFlux' }, { capeThreshold: 100 }]) assert.throws(() => build(options), /retired Betts–Miller/);
  if (!gpuAvailable) return;
  const { physicsConstants } = await import('../js/gpu/physics.gpu.js');
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  assert.throws(() => physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, convection: 'bettsMiller' }), /retired Betts–Miller/);
});

test("with cumulusOvershoot 'kinetic' the plume's reach into its first non-buoyant layer sets the share it detrains there: a share within (0, 1) for the trade-wind inversion that falls as the inversion strengthens, the layer below taking the rest, and the budgets stay exact", () => {
  const column = (options, jump) => {
    const model = build({ ...PREVIOUS_MOIST, ...options });
    place(model, 0, 101500, (z, p) => { const air = tradeWind(z, p); return z >= 1300 ? { T: air.T + jump, q: air.q } : air; });
    setDepth(model, 0, 600);
    model.boundaryLayer.buoyancyFlux[0] = 4e-4;
    model.boundaryLayer.friction[0] = 0.25;
    model.radiation.mlmGate[0] = 0.3;
    return model;
  };
  const run = (options, jump) => {
    const model = column(options, jump), { moist } = model, [pi, theta, , , q, qc] = model.state;
    const before = budget(model, 0), snap = snapshot(model, 0);
    moist.cumulusColumn(0, pi, theta, q, qc, 600);
    const after = budget(model, 0), { top } = moist.cumulus, C = model.core.diagnostics.C;
    assert.ok(top > 0 && moist.cumulusBaseFlux[0] > 0, `a plume topping in layer ${top}`);
    assert.ok(Math.abs(after.enthalpy - before.enthalpy) < 1e-12 * before.enthalpy, `enthalpy ${before.enthalpy} → ${after.enthalpy}`);
    assert.ok(Math.abs(after.water - before.water) < 1e-12 * before.water, `water ${before.water} → ${after.water}`);
    return { share: moist.cumulus.overshoot, top, topWater: q[top * C] + qc[top * C] - snap.q[top] - snap.qc[top], belowWater: q[(top + 1) * C] + qc[(top + 1) * C] - snap.q[top + 1] - snap.qc[top + 1] };
  };
  const full = run({ cumulusOvershoot: 1 }, 0), kinetic = run({}, 0), stronger = run({}, 4);
  console.log(`kinetic overshoot: share ${kinetic.share.toFixed(3)} into layer ${kinetic.top} (the full overshoot's share ${full.share}); with the inversion 4 K stronger ${stronger.share.toFixed(3)}; the top layer's moistening ${kinetic.topWater.toExponential(3)} against ${full.topWater.toExponential(3)}, the layer below's ${kinetic.belowWater.toExponential(3)} against ${full.belowWater.toExponential(3)}`);
  assert.equal(full.share, 1);
  assert.ok(kinetic.share > 0 && kinetic.share < 1, `share ${kinetic.share}`);
  assert.ok(stronger.share < kinetic.share, `${stronger.share} under a stronger inversion against ${kinetic.share}`);
  assert.equal(kinetic.top, full.top);
  assert.ok(Math.abs(kinetic.topWater - kinetic.share * full.topWater) < 1e-3 * Math.abs(full.topWater), `the top layer takes the share: ${kinetic.topWater} against ${kinetic.share * full.topWater}`);
  assert.ok(kinetic.belowWater > full.belowWater, 'the layer below takes the rest');
});

test("on BOMEX's profile (Siebesma et al. 2003, Table B1) the shallow plume carried through non-buoyant layers by its w² and entraining at 1/z tops in the LES's cloud layer, above the plume that ends in its first non-buoyant layer, enthalpy and water conserved", () => {
  const profile = (table) => (z) => { let n = 1; while (n < table.length - 1 && z > table[n][0]) n++; const [z0, v0] = table[n - 1], [z1, v1] = table[n]; return v0 + (v1 - v0) * (z - z0) / (z1 - z0); };
  const thetaL = profile([[0, 298.7], [520, 298.7], [1480, 302.4], [2000, 308.2], [3000, 311.85]]);
  const water = profile([[0, 17e-3], [520, 16.3e-3], [1480, 10.7e-3], [2000, 4.2e-3], [3000, 3e-3]]);
  const run = (options, set) => {
    const model = createModel(new Grid(16), { ocean: false, levels: sigmaInterfaces(set), moist: options }), { moist } = model, [pi, theta, , , q, qc] = model.state, { C, cp, kappa, geopotential, exnerLayer, exnerLower, g } = model.core.diagnostics;
    place(model, 0, 101500, (z, p) => { const T = thetaL(z) * (p / 1e5) ** kappa; return { T, q: Math.max(5e-4, Math.min(water(z), saturationHumidity(T, p))) }; });
    setDepth(model, 0, 500);
    model.boundaryLayer.buoyancyFlux[0] = 5.7e-4;
    model.boundaryLayer.friction[0] = 0.28;
    model.radiation.mlmGate[0] = 0.3;
    const before = budget(model, 0);
    moist.cumulusColumn(0, pi, theta, q, qc, 600);
    const after = budget(model, 0), { top } = moist.cumulus;
    assert.ok(top > 0 && moist.cumulusBaseFlux[0] > 0, `a plume topping in layer ${top}`);
    assert.ok(Math.abs(after.enthalpy - before.enthalpy) < 1e-12 * before.enthalpy, `enthalpy ${before.enthalpy} → ${after.enthalpy}`);
    assert.ok(Math.abs(after.water - before.water) < 1e-12 * before.water, `water ${before.water} → ${after.water}`);
    const thetaV = model.core.arrays.thetaV, upper = (k) => (geopotential[k * C] + cp * thetaV[k * C] * (exnerLayer[k * C] - exnerLower[(k - 1) * C])) / g;
    const base = upper(top + 1);
    return base + moist.cumulus.overshoot * (upper(top) - base);
  };
  const arms = { defaults: { cumulusAscent: 'buoyant', cumulusEntrainmentHeight: null }, carried: { cumulusAscent: 'kinetic', cumulusEntrainmentHeight: null }, perHeight: { cumulusAscent: 'buoyant', cumulusEntrainmentHeight: 1 }, both: { cumulusAscent: 'kinetic', cumulusEntrainmentHeight: 1 }, lowest: { cumulusAscent: 'kinetic', cumulusEntrainmentHeight: 1, cumulusSource: 'lowest' } };
  for (const set of ['bl34', 'bl36', 'bl42']) {
    const tops = Object.fromEntries(Object.entries(arms).map(([name, options]) => [name, run(options, set)]));
    console.log(`BOMEX cumulus top on ${set}, the top layer's base + its share of the layer (km): ${Object.entries(tops).map(([name, z]) => `${name} ${(z / 1000).toFixed(2)}`).join(', ')}`);
    assert.ok(tops.defaults < 1100, `${set}: the plume ending in its first non-buoyant layer tops at ${tops.defaults} m`);
    assert.ok(tops.carried >= tops.defaults, `${set}: carried ${tops.carried} m against ${tops.defaults} m`);
    for (const name of ['perHeight', 'both', 'lowest']) assert.ok(tops[name] > 1600 && tops[name] < 2100, `${set}: ${name} tops at ${tops[name]} m, the LES's core top 1.75–1.8 km and cloud top about 2 km`);
  }
});

const LIGHT_CUMULUS = { cumulusEntrainmentHeight: null, cumulusEntrainment: 1e-4, cumulusRain: 0 };

test("with shallowRetype a column the test parcel leaves shallow runs the deep plume where the shallow plume's own ascent, continued past the shallow top without tendencies, reaches its first cloudy layer with w² > 0, carried from the shallow closure's velocity scale, and has a cloud deeper than 200 hPa from that layer's lower interface to the last interface it passes, as random columns under a light shallow plume do; every other column steps exactly as without the option, the trade-wind cumulus among them, and a value other than true or false is refused", () => {
  const dt = 900, { shallowTop } = MOIST_DEFAULTS;
  const retyping = (shallow, grid = levels) => {
    const stacks = [];
    const [off, on] = [false, true].map((deepen) => {
      const model = createModel(new Grid(16), { ocean: false, levels: grid, moist: { ...shallow, shallowRetype: deepen } }), C = model.mesh.nCells, { moist, core } = model, { levels, sigmaMid } = core.diagnostics;
      randomColumns(model);
      const [pi, theta, , , q, qc] = model.state, deepened = new Set();
      let worst = 0, fired = 0, repeated = 0, based = 0, continued = 0;
      for (let i = 0; i < C; i++) {
        core.diagnoseColumn(i, pi, theta, q, qc);
        const before = budget(model, i), rain = moist.plumeColumn(i, pi, theta, q, qc, dt), after = budget(model, i), d = moist.deep, cu = moist.cumulus;
        if (!deepen) stacks[i] = { S: cu.sourceS, Q: cu.sourceQ, mass: cu.mass, source: cu.source, lcl: cu.lclPressure };
        worst = Math.max(worst, Math.abs(after.enthalpy - before.enthalpy) / before.enthalpy, Math.abs(after.water + rain - before.water) / before.water);
        if (d.deepened) {
          deepened.add(i);
          assert.ok(deepen && !d.parcelDeep && d.rooted && d.shallowDepth > DEEP_CLOUD_DEPTH, `column ${i}: re-typed with a shallow cloud of ${d.shallowDepth} Pa`);
          const stack = stacks[i];
          assert.ok(d.sourceS === stack.S && d.sourceQ === stack.Q && d.sourceMass === stack.mass && d.excessT === 0 && d.excessQ === 0, `column ${i}: the re-typed deep plume departs from the shallow plume's own stack: s_l ${d.sourceS} and ${stack.S} J/kg, q_t ${d.sourceQ} and ${stack.Q}, ${d.sourceMass} and ${stack.mass} Pa, excess ${d.excessT} K and ${d.excessQ}`);
          if (needsContinuation(d, pi[i], levels, sigmaMid)) continued++;
          if (d.deep) {
            fired++;
            assert.equal(moist.cumulusBaseFlux[i], d.baseFlux, `column ${i}: the deep plume runs alone`);
            let base = stack.source;
            while (base > 1 && !(pi[i] * levels[base] <= stack.lcl)) base--;
            assert.equal(d.base, base, `column ${i}: the re-typed deep plume's cloud base is the interface at the shallow plume's condensation level`);
            assert.ok(pi[i] * (levels[d.base] - levels[d.top + 1]) > DEEP_CLOUD_DEPTH, `column ${i}: the re-typed deep plume's cloud is ${pi[i] * (levels[d.base] - levels[d.top + 1])} Pa deep`);
          }
        } else if (d.shallowDepth > 0 && cu.top >= 1 && pi[i] * levels[cu.top - 1] > shallowTop) {
          assert.equal(d.shallowReach, cu.top + 1, `column ${i}: the diagnostic ascent ends where the shallow plume's does`);
          if (moist.cumulusBaseFlux[i] > 0) {
            let cloudy = cu.source - 1;
            while (cloudy > cu.top && !(moist.cumulusWater[cloudy * C + i] > 0)) cloudy--;
            assert.equal(d.shallowBase, cloudy + 1, `column ${i}: the diagnostic cloud starts at the lower interface of the shallow plume's first cloudy layer`);
            based++;
          }
          repeated++;
        }
      }
      return { model, deepened, worst, fired, repeated, based, continued };
    });
    let moved = 0;
    for (let i = 0; i < on.model.mesh.nCells; i++) {
      const a = snapshot(off.model, i), b = snapshot(on.model, i), equal = ['theta', 'q', 'qc'].every((key) => a[key].every((v, k) => v === b[key][k]));
      if (on.deepened.has(i)) { if (!equal) moved++; } else assert.ok(equal, `column ${i} is not re-typed, yet it steps differently`);
    }
    return { deepened: on.deepened.size, fired: on.fired, continued: on.continued, moved, repeated: on.repeated, based: on.based, worst: Math.max(off.worst, on.worst), offDeepened: off.deepened.size };
  };
  const plain = retyping({}), light = retyping(LIGHT_CUMULUS), fine = retyping(LIGHT_CUMULUS, sigmaInterfaces('bl48'));
  console.log(`bl48 random columns under the light shallow plume: ${fine.deepened} re-typed deep, ${fine.continued} of them by the ascent continued past shallowTop from a cloud no deeper than DEEP_CLOUD_DEPTH below it`);
  console.log(`random columns: ${plain.deepened} re-typed deep under the default shallow plume, ${light.deepened} under one entraining 1e-4 /m that keeps no condensate (${light.fired} fire the deep plume, ${light.moved} step differently); ${plain.repeated} and ${light.repeated} shallow plumes that stop below the shallow top end the diagnostic ascent at their own top, ${plain.based} and ${light.based} of them running and starting its cloud at their own first cloudy layer; every other column steps exactly as without the option; enthalpy and water within ${Math.max(plain.worst, light.worst).toExponential(1)}`);
  assert.ok(light.deepened >= 20 && Math.abs(fine.continued - 3) <= 1 && light.fired > 0 && light.moved > 0 && plain.repeated >= 20 && light.repeated >= 20 && plain.based > 0 && light.based > 0 && [plain, light, fine].every((r) => r.moved === r.fired), JSON.stringify({ plain, light, fine }));
  assert.ok(plain.offDeepened === 0 && light.offDeepened === 0 && fine.offDeepened === 0 && Math.max(plain.worst, light.worst, fine.worst) < 1e-15);
  const trade = [false, true].map((deepen) => {
    const model = tradeWindColumn({ shallowRetype: deepen }), [pi, theta, , , q, qc] = model.state;
    model.core.diagnoseColumn(0, pi, theta, q, qc);
    model.moist.plumeColumn(0, pi, theta, q, qc, 600);
    return { snap: snapshot(model, 0), d: { ...model.moist.deep }, flux: model.moist.cumulusBaseFlux[0] };
  });
  console.log(`trade-wind column: the shallow plume's cloud is ${(trade[1].d.shallowDepth / 100).toFixed(1)} hPa deep, base flux ${trade[1].flux.toExponential(3)} kg/m²/s with the option and ${trade[0].flux.toExponential(3)} without`);
  assert.ok(trade[1].d.shallowDepth > 0 && trade[1].d.shallowDepth < DEEP_CLOUD_DEPTH && !trade[1].d.deepened && trade[1].flux > 0, JSON.stringify(trade[1].d));
  for (const key of ['theta', 'q', 'qc']) assert.deepEqual(trade[1].snap[key], trade[0].snap[key]);
  assert.throws(() => build({ shallowRetype: 'yes' }), /shallowRetype/);
});

const needsContinuation = (d, ps, levels, sigmaMid) => {
  let ceiling = d.shallowBase - 1;
  while (ceiling > 1 && ps * sigmaMid[ceiling] > MOIST_DEFAULTS.shallowTop) ceiling--;
  return ps * levels[d.shallowReach] < MOIST_DEFAULTS.shallowTop && ps * (levels[d.shallowBase] - levels[Math.max(d.shallowReach, ceiling + 1)]) <= DEEP_CLOUD_DEPTH;
};

const PAGE_COLUMN_13676 = {"dt":337.5,"pi":101494.40113024111,"theta":[3054.8886076870785,2110.9413649520975,1820.1921609629132,1668.3065751122206,1483.0884561461764,1298.6364960778455,1163.885616732911,1013.5028676910305,874.8338179738957,765.0030420564274,669.7593927653216,593.8539532046864,530.6715428284349,479.14471254294887,443.17516557885824,417.9224745649334,399.0292911879083,388.0540775219285,374.200921511677,354.07304435921634,338.59922374210277,328.57798368631126,324.6353538958846,320.6419441498382,317.41504052946897,315.05046284636757,310.42418753003966,309.2426625953133,306.90000083542407,304.1396705843355,301.9113565318508,299.9234083059886,297.7695582225899,296.99878076720324,295.0792027868827,294.80438402875257,292.8648375547794,287.8811165042776,286.18628885345737,286.1012124480151,286.0485258051846,286.01457996996,285.99385517786254,285.98342543624153,285.9800157296597,285.98201600592904,285.99222645564043,286.04102861618844],"q":[0.000004420723201075017,0.000004379845367709267,0.00000451711231569177,0.000004641158623261906,0.0000043808327147769074,0.000004102350071107816,0.0000031171606670207203,0.0000036988160018668285,0.0000040569189342878315,0.0000029872904051668265,0.000002456148478663317,6.651212904623757e-7,0.0000013952885461895664,0.000002310178851262623,0.000003354896636465968,0.000004194903527820176,0.0000013663247194122255,0.000003963132836540237,0.000006033709453996757,0.000021907684372351807,0.000023571580366914115,0.0000612085636522033,0.00014476234528108182,0.0003354868961888411,0.0008351307543146802,0.0014669471066531825,0.0005634933873287046,0.0002764959016803373,0.00009219694890869269,-0.000001009581714749424,-2.8848519499736504e-7,-0.000001125085710560899,0.000006406012106246861,0.000010330015009707616,0.0004970793442477482,0.000220454153345012,0.0003997442045174921,0.00368551594180363,0.005673311892906608,0.005804869610221586,0.005909699594257936,0.0059993113195499406,0.006078471117161322,0.006150788910067775,0.0061942450391528755,0.006226437483185782,0.006262178789386516,0.006362772421041086],"qc":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000009211694330450904,0.000010248246418038136,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"depth":805.877766933246,"regime":3,"mixingTop":805.877766933246,"cloudLayer":38,"buoyancy":0.0006600150690247001,"friction":0.17775174086504517,"surfaceDrag":0.006910410630024459,"mlmGate":0.04205867755834928,"stratiform":0,"sensible":17.6087223882052,"evaporation":0.00003551381754490235,"cumulusK0":29,"cumulusCover":[0,0,0,0,0,0,0,0.0001757078851368149,0,0,0,0,0,0,0,0,0,0,0],"cumulusWater":[0,0,0,0,0,0,0,0.00014711110022179344,0,0,0,0,0,0,0,0,0,0,0]};
const PAGE_COLUMN_3555 = {"dt":337.5,"pi":99033.06746108439,"theta":[3102.3002279761645,2117.5990433905467,1852.8427965266644,1705.8564912562001,1525.3006180293805,1354.5111878740731,1168.3160652922004,1018.272717508987,887.1122902854764,776.8063226312967,689.5915142349497,605.2248865183666,534.7005464646481,475.8929672806757,434.46506127266554,400.9630054023617,381.1445687182407,364.8127286495149,353.14853864267235,348.28877351690886,344.40942673553366,342.792827396337,340.6640236019349,336.65709774319816,331.12486466450264,328.7493375020934,326.0121645161613,322.0675348159161,317.889984868103,314.33388034657025,312.09108636577326,310.91877589972887,309.33915365199397,307.888893256398,306.6547822488465,305.59767104902915,304.5162680978994,303.84790276064575,303.35933160693145,302.93148945339055,302.9183612865265,302.91501682849344,302.91431122291493,302.9161616670703,302.92089863468817,302.93031276302594,302.95084518136355,303.03020751049576],"q":[0.000004425619126059593,0.000004342598056354839,0.000004423958929910972,0.000004205286229110553,0.0000041679370239080995,0.0000033430861145152167,0.000002475740934458977,0.0000016660833508865741,0.0000010888131190650675,7.392344141963997e-7,7.722073266518058e-8,9.571011600521492e-7,0.0000021287316061971877,0.0000031066628359052366,0.000002914556874113797,0.0000018191188769516828,0.000002661308234261021,0.000005920695395969138,0.00001078485746030204,0.000009394507630896911,0.000058951161863561513,0.0001165329032180707,0.00010438280867699178,0.00005136603369991584,0.00025972494655360044,0.0002304913435511416,0.00045919332298096293,0.0010976971854892328,0.0015319380612758846,0.005255817368282934,0.00926844556527987,0.010728357958316987,0.011613512728149944,0.012264463174956402,0.01297368031152491,0.013760140964062681,0.014500611222856808,0.015505838190287863,0.016113016786198878,0.01629556279053427,0.01630082886116683,0.016303165623252076,0.01630530514945669,0.016307456822123423,0.016309668375235082,0.01631229910986106,0.016316318148494623,0.01632910408734579],"qc":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.0597240715426335e-7,0.0000024029565526814046,0,0,0,0,0,0,0,0,0,0,0,0.00016817669548748387,0.00013340437396978203,0.00007355482465833066,0.000025807765435322942,0.000009093681307620125,0.000004946826467257722,0.000008483349041467395,0.000020649094778766102,0.000005035842738743802,0,0,0,0,0,0,0,0,0],"depth":767.1659843495459,"regime":1,"mixingTop":767.1659843495459,"cloudLayer":-1,"buoyancy":0.0007226390468329481,"friction":0.13536864088306785,"surfaceDrag":0.012051972602621234,"mlmGate":0.5,"stratiform":0,"sensible":24.645496045631113,"evaporation":0.000004382302221062775,"cumulusK0":29,"cumulusCover":[0.002260689077812534,0.0027060352877132155,0.0033487233595176804,0.003782176657685643,0.004183040333503172,0.004620965472531737,0.00489292228236364,0.006031432779177265,0.010444214518071024,0,0,0,0,0,0,0,0,0,0],"cumulusWater":[0.0023216751193501895,0.0019576912729000326,0.00160475700107698,0.0013248772001631611,0.0010992639595943797,0.0009260106961909632,0.0007914371239626211,0.0006010073452279923,0.00021224869705916693,0,0,0,0,0,0,0,0,0,0]};
const replayPageColumn = (c, deepen) => {
  const model = createModel(new Grid(16), { ocean: false, levels: sigmaInterfaces('bl48'), moist: { shallowRetype: deepen } });
  const { moist, core, boundaryLayer: bl, radiation } = model, { K, C, levels } = core.diagnostics, [pi, theta, , , q, qc] = model.state;
  pi[0] = c.pi;
  for (let k = 0; k < K; k++) { theta[k * C] = c.theta[k]; q[k * C] = c.q[k]; qc[k * C] = c.qc[k]; }
  bl.depth[0] = c.depth; bl.regime[0] = c.regime; bl.mixingTop[0] = c.mixingTop; bl.cloudLayer[0] = c.cloudLayer; bl.buoyancyFlux[0] = c.buoyancy; bl.friction[0] = c.friction; bl.surfaceDrag[0] = c.surfaceDrag;
  radiation.mlmGate[0] = c.mlmGate; radiation.stratiform[0] = c.stratiform; radiation.sensibleHeat[0] = c.sensible; radiation.evaporation[0] = c.evaporation;
  c.cumulusCover.forEach((v, n) => { moist.cumulusCover[(c.cumulusK0 + n) * C] = v; moist.cumulusWater[(c.cumulusK0 + n) * C] = c.cumulusWater[n]; });
  core.diagnoseColumn(0, pi, theta, q, qc);
  moist.plumeColumn(0, pi, theta, q, qc, c.dt);
  const d = { ...moist.deep };
  return { d, cu: { ...moist.cumulus }, flux: moist.cumulusBaseFlux[0], snap: snapshot(model, 0), depth: d.deep ? pi[0] * (levels[d.base] - levels[d.top + 1]) : 0, base: pi[0] * levels[d.shallowBase], ps: pi[0], levels };
};
const steppedAlike = (a, b) => ['theta', 'q', 'qc'].every((key) => a.snap[key].every((v, k) => v === b.snap[key][k]));
test("with shallowRetype a re-typed deep plume feeds its stack below the source by cumulusSource, as the shallow plume does: under cumulusSource 'lowest' a re-typed column steps alike under either plumeSource", () => {
  const runs = ['mean', 'lowest'].map((plumeSource) => {
    const model = createModel(new Grid(16), { ocean: false, levels, moist: { ...LIGHT_CUMULUS, cumulusSource: 'lowest', plumeSource, shallowRetype: true } }), C = model.mesh.nCells, { moist, core } = model, { K, levels: edges } = core.diagnostics;
    randomColumns(model);
    const [pi, theta, , , q, qc] = model.state, retyped = new Map();
    for (let i = 0; i < C; i++) {
      core.diagnoseColumn(i, pi, theta, q, qc);
      moist.plumeColumn(i, pi, theta, q, qc, 900);
      if (moist.deep.deepened) retyped.set(i, moist.deep.deep && moist.deep.sourceMass > 1.000001 * pi[i] * (1 - edges[K - 1]));
    }
    return { model, retyped };
  });
  let compared = 0, fed = 0;
  for (const [i, stacked] of runs[0].retyped) {
    if (!runs[1].retyped.has(i)) continue;
    const a = snapshot(runs[0].model, i), b = snapshot(runs[1].model, i);
    for (const key of ['theta', 'q', 'qc']) assert.deepEqual(b[key], a[key], `column ${i}: ${key} differs between plumeSource 'mean' and 'lowest'`);
    compared++;
    if (stacked) fed++;
  }
  console.log(`  ${compared} columns re-type under both plumeSource values, ${fed} of them firing from a stack of more than one layer`);
  assert.ok(fed > 0, `${compared} columns re-type under both, ${fed} of them firing from a stack of more than one layer`);
});
test("with shallowRetype and midLevelConvection a column re-typed deep whose own deep cloud the depth rule finds shallow tries the mid-level type, stepping as without shallowRetype where a shallow plume closed at 1e-9 does not rise", () => {
  const dt = 900;
  const step = (options) => {
    const model = createModel(new Grid(16), { ocean: false, levels, moist: { ...LIGHT_CUMULUS, cumulusClosure: 1e-9, midLevelConvection: true, ...options } }), C = model.mesh.nCells, { moist, core } = model, { K } = core.diagnostics;
    randomColumns(model);
    let seed = 777;
    const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
    for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) core.diagnostics.omega[k * C + i] = 0.6 * random() - 0.4;
    const [pi, theta, , , q, qc] = model.state, columns = [];
    for (let i = 0; i < C; i++) {
      core.diagnoseColumn(i, pi, theta, q, qc);
      moist.plumeColumn(i, pi, theta, q, qc, dt);
      columns.push({ ...moist.deep });
    }
    return { state: model.state, columns, C };
  };
  const both = step({ shallowRetype: true }), alone = step({}), C = both.C;
  let fired = 0, sentBack = 0, apart = 0;
  for (let i = 0; i < C; i++) {
    const d = both.columns[i];
    if (d.deepened && d.deep && !d.midLevel) { fired++; continue; }
    if (d.deepened && d.midLevel) sentBack++;
    if (!both.state.every((a, n) => { for (let x = i; x < a.length; x += C) if (a[x] !== alone.state[n][x]) return false; return true; })) apart++;
  }
  console.log(`random columns under a random omega with midLevelConvection: with shallowRetype too ${fired} run the re-typed deep plume and ${sentBack} are sent back by the depth rule and convect at mid-level; ${apart} of the other columns step otherwise than without shallowRetype`);
  assert.ok(sentBack > 0 && apart === 0, `${fired} run the re-typed deep plume, ${sentBack} sent back convect at mid-level, ${apart} of the others step otherwise than without shallowRetype`);
});
test("with midLevelConvection and shallowRetype the engines agree on random columns under a random large-scale ascent, where columns sent back by the depth rule try the mid-level type under a shallow plume closed at 1e-9", { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { midLevel, deepened, sentBack } = await parity({ midLevelConvection: true, shallowRetype: true, ...LIGHT_CUMULUS, cumulusClosure: 1e-9 }, { ascent: true });
  console.log(`under midLevelConvection and shallowRetype: ${midLevel} columns convect at mid-level on the CPU, ${deepened} re-type deep, ${sentBack} of them sent back by the depth rule to convect at mid-level`);
  assert.ok(midLevel >= 50 && sentBack > 0, `${midLevel} columns convect at mid-level, ${sentBack} of them sent back by the depth rule`);
});
test("with shallowRetype the column of the N=64 page state at 30N (seventeen64 day 365, column 13676, after two steps), whose shallow parcel first holds condensate at 539 hPa, is not re-typed: its w² from the closure's velocity scale vanishes below that layer", () => {
  const [off, on] = [false, true].map((deepen) => replayPageColumn(PAGE_COLUMN_13676, deepen));
  assert.ok(!on.d.deepened && on.d.shallowDepth === 0 && steppedAlike(on, off), `the column is not re-typed: shallow cloud ${on.d.shallowDepth} Pa, re-typed ${on.d.deepened}`);
});
test('with shallowRetype the column of the N=64 page state at 6N (seventeen64 day 365, column 3555, land 121 m up, after two steps, replayed at sea level with its boundary layer as deep above the ground), whose shallow plume holds condensate from 899 hPa, re-types on a cloud deeper than 200 hPa and runs a deep plume deeper than 200 hPa', () => {
  const [off, on] = [false, true].map((deepen) => replayPageColumn(PAGE_COLUMN_3555, deepen));
  assert.ok(on.d.deepened && !on.d.parcelDeep && on.d.shallowDepth > DEEP_CLOUD_DEPTH && Math.abs(on.base - 899e2) < 100, `the column re-types on a shallow cloud ${on.d.shallowDepth} Pa deep from ${on.base} Pa`);
  assert.ok(off.flux > 0 && on.d.deep && on.depth > DEEP_CLOUD_DEPTH && !steppedAlike(on, off), `the re-typed column fires a deep plume ${on.depth} Pa deep (CAPE ${on.d.cape} J/kg) in place of a shallow plume of base flux ${off.flux} kg/m²/s`);
  assert.ok(on.d.sourceS === off.cu.sourceS && on.d.sourceQ === off.cu.sourceQ && on.d.sourceMass === off.cu.mass && on.d.excessT === 0 && on.d.excessQ === 0, `the deep plume departs from the shallow plume's own stack: s_l ${on.d.sourceS} and ${off.cu.sourceS} J/kg, q_t ${on.d.sourceQ} and ${off.cu.sourceQ}, ${on.d.sourceMass} and ${off.cu.mass} Pa, excess ${on.d.excessT} K and ${on.d.excessQ}`);
  let base = off.cu.source;
  while (base > 1 && !(off.ps * off.levels[base] <= off.cu.lclPressure)) base--;
  assert.equal(on.d.base, base, `the deep plume's cloud base is the interface at the shallow plume's condensation level, ${off.cu.lclPressure} Pa`);
});
test("with shallowRetype under a light shallow plume the engines agree on random columns, under cumulusSource 'lowest' too, and the GPU refuses a value other than true or false", { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const lowest = await parity({ shallowRetype: true, cumulusSource: 'lowest', ...LIGHT_CUMULUS });
  assert.ok(lowest.deepened > 0, `${lowest.deepened} columns re-type deep under cumulusSource 'lowest'`);
  const { deepened } = await parity({ shallowRetype: true, ...LIGHT_CUMULUS }), fine = await parity({ shallowRetype: true, ...LIGHT_CUMULUS }, { grid: sigmaInterfaces('bl48') });
  const withRetype = await parity({ shallowRetype: true, convectionRetype: true, ...LIGHT_CUMULUS });
  assert.ok(withRetype.deepened > 0, `${withRetype.deepened} columns re-type deep with convectionRetype on as well`);
  console.log(`  ${deepened} columns re-typed deep on the CPU`);
  console.log(`  on bl48 ${fine.deepened} columns re-typed deep on the CPU, ${fine.continued} of them by the continued ascent; the engines part on ${fine.continuedApart} of those, ${fine.continuedLeftOut} left out`);
  assert.ok(fine.continuedApart === 0 && fine.continuedLeftOut === 0 && Math.abs(fine.continued - 3) <= 1, `on bl48 the engines part on ${fine.continuedApart} of the ${fine.continued} columns re-typed by the continued ascent, ${fine.continuedLeftOut} left out`);
  assert.ok(deepened >= 20 && fine.continued > 0, `${deepened} columns re-typed deep on the CPU; on bl48 ${fine.deepened}, ${fine.continued} of them by the ascent continued past shallowTop from a cloud no deeper than DEEP_CLOUD_DEPTH below it`);
  const { physicsConstants } = await import('../js/gpu/physics.gpu.js');
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  assert.throws(() => physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, shallowRetype: 1 }), /shallowRetype/);
});
