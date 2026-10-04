import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { saturationHumidity, cloudSaturation, liquidFraction, LATENT_HEAT, FUSION_HEAT, MOIST_DEFAULTS, COUPLED_REGIME, DECK_CLOSED, BECHTOLD, SUBCLOUD_LAYERS, DEEP_CLOUD_DEPTH, IFS_ENTRAINMENT, TEST_PARCEL, IFS_PRECIPITATION } from '../js/physics/moist.module.js';
import { VIRTUAL_FACTOR } from '../js/dynamics/sigmaCore.module.js';
import { REGIME } from '../js/physics/boundaryLayer.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }

const levels = sigmaInterfaces('bl34');
const build = (moist = {}, N = 2) => createModel(new Grid(N), { ocean: false, levels, moist });

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

test('autoconversion stays out of the lowest two layers, or with autoconversionFloor: boundaryLayer out of the boundary layer, and rain evaporates only into cloud-free layers', () => {
  const lowest = build(), [pl, tl, , , ql, qcl] = lowest.state, KL = lowest.core.K, CL = lowest.core.diagnostics.C;
  jordanColumn(lowest, 0);
  for (let k = 0; k < KL; k++) qcl[k * CL] = 0;
  qcl[(KL - 1) * CL] = 2e-3; qcl[(KL - 2) * CL] = 2e-3;
  assert.equal(lowest.moist.autoconvertColumn(0, pl, tl, ql, qcl, 600), 0, 'cloud in the lowest two layers does not rain');
  qcl[(KL - 3) * CL] = 2e-3;
  assert.ok(lowest.moist.autoconvertColumn(0, pl, tl, ql, qcl, 600) > 0 || qcl[(KL - 3) * CL] < 2e-3, 'the third layer converts');
  assert.throws(() => build({ autoconversionFloor: 'surface' }), /autoconversionFloor/);
  const model = build({ autoconversionFloor: 'boundaryLayer' }), { moist } = model, [pi, theta, , , q, qc] = model.state;
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

test('with upperCloudLifetime cloud water in the layers above the shallow top converts over that lifetime and the cloud below over cloudLifetime (with no ice falling, all of it as liquid)', () => {
  const dt = 600, model = build({ upperCloudLifetime: 1800, iceFall: null }), plain = build({ iceFall: null }), { K, C, sigmaMid } = model.core.diagnostics;
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

function plumeColumn(options = {}, profile = null, N = 2) {
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

test('under the Bechtold closure a warming subcloud layer over the sea takes its PCAPE_bl over the advective time z_base / ū_bl, and over land over the turnover time, which holds a plume under strong heating back', () => {
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
  const fire = (heating) => {
    const land = plumeColumn({ land: Uint8Array.from({ length: 16 * 16 * 10 + 2 }, (_, i) => (i === 0 ? 1 : 0)) }, null, 16), [lp, lt, , , lq, lqc] = land.state;
    subcloudVirtual(land, 0, heating, dt);
    land.moist.plumeColumn(0, lp, lt, lq, lqc, dt);
    return { ...land.moist.deep };
  };
  const heated = fire(10), calm = fire(0);
  console.log(`land: under +10 K/d below cloud base τ_bl ${heated.boundaryTime.toFixed(0)} s (H / w̄) gives PCAPE_bl ${heated.pcapeBoundary.toFixed(0)} Pa against PCAPE ${heated.pcape.toFixed(0)}: deep ${heated.deep}; with no heating base flux ${calm.baseFlux.toFixed(4)} kg/m²/s`);
  assert.ok(!heated.deep && heated.pcapeBoundary > heated.pcape && Math.abs(heated.boundaryTime - heated.depth / heated.speed) < 1e-9 * heated.boundaryTime, 'strong heating holds the land plume back over the turnover time');
  assert.ok(calm.deep && calm.baseFlux > 0 && calm.pcapeBoundary === 0, 'the calm land column fires');
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
 * Random columns of five kinds (Jordan's sounding warmed or cooled, trade
 * cumulus, stratocumulus, trade wind, a cold dry Jordan) with random
 * boundary layers, deck gates, regimes, mixing tops, sea ice, saved
 * subcloud T_v and surface fluxes, all in single precision, on an N=6 model.
 */
function randomColumns(model) {
  const { moist, core, mesh } = model, C = mesh.nCells, { K } = core.diagnostics;
  const [pi, theta, , surfaceT, q, qc] = model.state;
  let seed = 12345;
  const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  const kinds = new Int32Array(C);
  const buoyancy = new Float64Array(C), friction = new Float64Array(C);
  for (let i = 0; i < C; i++) {
    const kind = Math.floor(random() * 5), warm = 3 * (random() - 0.5), wet = 1 + 0.15 * (random() - 0.3), surface = 100500 + 1500 * random();
    kinds[i] = kind;
    buoyancy[i] = random() < 0.15 ? -2e-4 * random() : 6e-4 * random();
    friction[i] = 0.1 + 0.3 * random();
    if (kind === 0) place(model, i, surface, (z, p) => { const air = jordan(p); const T = air.T + warm * (z < 800 ? 1 : 0.5); return { T, q: Math.min(wet * air.q, 0.95 * saturationHumidity(T, p)) }; });
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
    const model = build(options, 6), C = model.mesh.nCells;
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
  assert.ok(deep > 10 && same > 10 && alone > 10 && both > 0, `${deep} deep, ${same} alike, ${alone} alone, ${both} with both under top`);
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
function handParcel(model) {
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
  const twin = build({}, 2); twin.state.forEach((a, n) => a.set(model.state[n])); twin.boundaryLayer.depth.set(model.boundaryLayer.depth); twin.radiation.mlmGate.set(model.radiation.mlmGate);
  twin.core.diagnoseColumn(0, twin.state[0], twin.state[1], twin.state[4], twin.state[5]);
  twin.moist.plumeColumn(0, twin.state[0], twin.state[1], twin.state[4], twin.state[5], 600);
  let sl = twin.moist.deep.sourceS, qt = twin.moist.deep.sourceQ, w2 = 1, base = 0, deepAt = NaN;
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
  const dt = 900, model = build({}, 6), C = model.mesh.nCells;
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
  assert.ok(typed > 10 && fired > 10 && shallowTyped > 10 && shallowRun > 10, `${typed}, ${fired}, ${shallowTyped}, ${shallowRun}`);
  assert.ok(worst < 1e-15, `enthalpy and water ${worst}`);
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
  const dt = 900, model = build({}, 6), C = model.mesh.nCells, { K, exnerLayer } = model.core.diagnostics, { moist } = model;
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
  pass.setPipeline(kernels.adjust);
  pass.setBindGroup(0, device.createBindGroup({ layout: kernels.adjust.getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) }));
  pass.dispatchWorkgroups(Math.ceil(C / 64));
  pass.end();
  device.queue.submit([encoder.finish()]);
  const gpuState = await gpu.download(), ph = await gpu.downloadPhysics(), gpuBefore = model.state.map((a) => Float64Array.from(a));
  const twin = build({}, 6);
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

async function parity(options, { momentum = false } = {}) {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js');
  const model = build(options, 6), { moist, core, mesh } = model, C = mesh.nCells, { K } = core.diagnostics, dt = 900;
  const [pi, theta, , surfaceT, q, qc] = model.state;
  const { buoyancy, friction, regime, mixingTop, stratiform, concentration, saved, sensible, evaporation } = randomColumns(model);
  let cloudSeed = 4242;
  const cloudRandom = () => { cloudSeed = (cloudSeed * 1103515245 + 12345) % 2147483648; return cloudSeed / 2147483648; };
  for (let x = moist.cumulusK0 * C; x < K * C; x++) if (cloudRandom() < 0.3) { moist.cumulusCover[x] = Math.fround(0.05 * cloudRandom()); moist.cumulusWater[x] = Math.fround(1e-3 * cloudRandom()); }
  const priorCover = Float64Array.from(moist.cumulusCover), priorWater = Float64Array.from(moist.cumulusWater);
  const gpu = await createGpuCore(mesh, { levels, physics: options });
  const { device, buffers, kernels, layout } = gpu;
  gpu.upload(model.state);
  gpu.uploadPhysics({ mlmGate: model.radiation.mlmGate, concentration, cumulusCover: priorCover.subarray(moist.cumulusK0 * C), cumulusWater: priorWater.subarray(moist.cumulusK0 * C) });
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.REGIME, Float32Array.from(regime));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.MIXTOP, Float32Array.from(mixingTop));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.STRAT, Float32Array.from(stratiform));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.DEPTH, Float32Array.from(model.boundaryLayer.depth));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.BUOY, Float32Array.from(buoyancy));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.USTAR, Float32Array.from(friction));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.SUBTV, Float32Array.from(saved));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.SH, Float32Array.from(sensible));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.EVAP, Float32Array.from(evaporation));
  device.queue.writeBuffer(buffers.P, 0, Float32Array.from([dt, 0, 1, 0, 0, 0, 0, 0]));
  const group = device.createBindGroup({ layout: kernels.adjust.getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) });
  const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
  pass.setPipeline(kernels.adjust);
  pass.setBindGroup(0, group);
  pass.dispatchWorkgroups(Math.ceil(C / 64));
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
    for (const a of [moist.precipitation, moist.convectivePrecipitation, moist.largeScalePrecipitation]) a.fill(0);
    for (let x = 0; x < theta.length; x++) theta[x] *= 1 + 2 ** -23;
    model.phases.adjust(0, C, dt);
    for (let i = 0; i < C; i++) { warmed[i] = moist.convectivePrecipitation[i] + moist.largeScalePrecipitation[i]; warmedFlux[i] = moist.cumulusBaseFlux[i]; }
    thetaResponse.set(theta);
    if (momentum) { model.phases.mixMomentum(0, mesh.nEdges, dt); windResponse = Float64Array.from(model.state[2]); }
    model.state.forEach((a, n) => a.set(kept[n])); saved.set(keptSaved); model.seaIce.snow.set(keptSnow); moist.cumulusCover.set(priorCover); moist.cumulusWater.set(priorWater);
    for (const a of [moist.precipitation, moist.convectivePrecipitation, moist.largeScalePrecipitation]) a.fill(0);
  }
  moist.trace.convection = new Float64Array(K * C);
  const threshold = new Set();
  const typedAtEdge = new Set();
  for (let i = 0; i < C; i++) { model.phases.adjust(i, i + 1, dt); if (moist.deep.conversionMargin < 1e-3) threshold.add(i); if (moist.deep.typeMargin < 1e-3) typedAtEdge.add(i); }
  for (let i = 0; i < C; i++) { const now = moist.convectivePrecipitation[i] + moist.largeScalePrecipitation[i]; response[i] = Math.max(Math.abs(response[i] - now), Math.abs(warmed[i] - now)); }
  for (let x = 0; x < K * C; x++) thetaResponse[x] = Math.abs(thetaResponse[x] - theta[x] * (1 + 2 ** -23));
  if (momentum) {
    const E = mesh.nEdges, { dSigma, g } = core.diagnostics, u = model.state[2], before = Float64Array.from(u);
    model.phases.mixMomentum(0, E, dt);
    let moved = 0, worstU = 0, worstColumn = 0, scale = 0, windOver = -Infinity, windAllowed = 0;
    for (let e = 0; e < E; e++) {
      const columnMass = 0.5 * (pi[mesh.cellsOnEdge[2 * e]] + pi[mesh.cellsOnEdge[2 * e + 1]]);
      let was = 0, now = 0, changed = false;
      for (let k = 0; k < K; k++) {
        const x = k * E + e, m = columnMass * dSigma[k] / g;
        was += m * before[x]; now += m * u[x]; scale = Math.max(scale, Math.abs(m * before[x]));
        if (u[x] !== before[x]) changed = true;
        worstU = Math.max(worstU, Math.abs(u[x] - after[2][x]));
        windOver = Math.max(windOver, Math.abs(u[x] - after[2][x]) - (1e-3 + 2 * Math.abs(windResponse[x] - u[x])));
        windAllowed = Math.max(windAllowed, 2 * Math.abs(windResponse[x] - u[x]));
      }
      if (changed) moved++;
      worstColumn = Math.max(worstColumn, Math.abs(now - was));
    }
    console.log(`  momentum: the plume moves the wind on ${moved} of ${E} edges; each edge's column momentum changes by at most ${worstColumn.toExponential(1)} kg/m/s against layer momenta up to ${scale.toFixed(0)}; the engines' winds differ by ${worstU.toExponential(1)} m/s (twice the CPU's response to one f32 ulp of θ at most ${windAllowed.toExponential(1)} m/s)`);
    assert.ok(moved > E / 10, `${moved} edges moved`);
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
  let rainOver = -Infinity, allowed = 0;
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
    if (beside.has(i) || merged.has(i) || threshold.has(i) || typedAtEdge.has(i)) continue;
    fluxScale = Math.max(fluxScale, cpu);
    if (!switching.has(i)) worstFlux = Math.max(worstFlux, Math.abs(cpu - gpuFlux));
    for (let k = K0; k < K; k++) {
      worstCover = Math.max(worstCover, Math.abs(moist.cumulusCover[k * C + i] - ph.CUCOVER[(k - K0) * C + i]));
      worstWater = Math.max(worstWater, Math.abs(moist.cumulusWater[k * C + i] - ph.CUWATER[(k - K0) * C + i]));
      waterScale = Math.max(waterScale, moist.cumulusWater[k * C + i]);
    }
    const rainDiff = Math.max(Math.abs(moist.convectivePrecipitation[i] - ph.CONV[i]), Math.abs(moist.largeScalePrecipitation[i] - ph.COND[i]));
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
  console.log(`${JSON.stringify(options)}: ${C} random columns: ${deep} convect deep, ${shallow} shallow, ${decked} under a deck, ${opening} convect under a deck opening, ${still} still; convective rain differs in sign on ${flips} (beyond the rain's allowance below; within it by at most ${residues.toExponential(1)} kg/m²); engines differ in θ by ${worstTheta.toExponential(1)} K (twice the CPU's response to one f32 ulp of θ at most ${thetaAllowed.toExponential(1)} K), q by ${worstQ.toExponential(1)}, qc by ${worstQc.toExponential(1)}, a step's rain by ${worstRain.toExponential(1)} kg/m² (largest ${rainScale.toFixed(3)}; twice its response to one f32 ulp of the saved subcloud T_v or of θ at most ${allowed.toExponential(1)})`);
  console.log(`  ${plumes} plumes, ${plumeFlips} differ in whether they rise, ${topsDiffer} in their top, ${beside.size} deep in the base flux by more than 1 % (the shallow plume beside it ran on one engine only), ${merged.size} whose dry adjustment merged other layers on the two engines, ${switching.size} whose base flux on the CPU moves by over 1e-3 under one ulp of θ (left out of the flux), ${threshold.size} whose deep plume's condensate at an interface lies within 1e-3 of the conversion threshold and ${typedAtEdge.size} whose test parcel's cloud lies within 1e-3 of 200 hPa or whose w² falls within 1e-3 m²/s² of zero (left out); base mass flux differs by ${worstFlux.toExponential(1)} kg/m²/s (largest ${fluxScale.toFixed(3)}), the cumulus fraction by ${worstCover.toExponential(1)}, the plume's condensate by ${worstWater.toExponential(1)} kg/kg (largest ${waterScale.toExponential(1)})`);
  assert.ok(deep > C / 40 && plumes > C / 10 && decked > C / 40 && still > C / 40, `${deep} deep, ${plumes} plumes, ${decked} decked, ${still} still`);
  assert.equal(plumeFlips, 0);
  assert.equal(topsDiffer, 0);
  assert.ok(beside.size <= C / 100, `${beside.size} columns' shallow plume beside the deep one ran on one engine only`);
  assert.ok(merged.size <= C / 100, `${merged.size} columns' dry adjustment merged other layers`);
  assert.ok(switching.size <= 0.02 * C, `${switching.size} columns' base flux switches under one ulp of θ`);
  assert.ok(threshold.size <= C / 100 && typedAtEdge.size <= C / 100, `${threshold.size} columns' plume condensate at the conversion threshold, ${typedAtEdge.size} typed at the edge`);
  assert.ok(worstFlux < 2e-3 * fluxScale && worstCover < 1e-4, `base mass flux ${worstFlux}, cover ${worstCover}`);
  assert.ok(waterScale > 0 && worstWater < 1e-3 * waterScale, `plume condensate ${worstWater} against ${waterScale}`);
  assert.equal(flips, 0);
  assert.ok(thetaOver < 0 && worstQ < 1e-6 && worstQc < 1e-7, `θ ${worstTheta} (${thetaOver} above 1e-3 K and twice the response to one ulp), q ${worstQ}, qc ${worstQc}`);
  assert.ok(rainOver < 0, `rain ${worstRain} against ${rainScale}, ${rainOver} above 1e-4 of it and twice the response to one ulp of the saved T_v or of θ`);
  const { geopotential, g } = core.diagnostics, ice = model.state[6], shared = new Uint8Array(K * C);
  for (let i = 0; i < C; i++) for (let k = 0; k < K; k++) {
    const below = geopotential[k * C + i] / g < mixingTop[i], iced = ice[i] > 0 ? (concentration[i] > 0 ? concentration[i] : 1) : 0;
    const plumed = moist.cumulusBaseFlux[i] > 0 && pi[i] * sigmaMid[k] >= moist.cumulusTop[i];
    if (!plumed && Math.max(below ? (regime[i] === COUPLED_REGIME ? 1 : 0) : stratiform[i], iced) > 0) shared[k * C + i] = 1;
  }
  return { cloud: Float64Array.from(model.state[5]), gpuCloud: Float64Array.from(after[5]), shared };
}

test('the shallow and deep plume and the rain they leave match between the engines on a random set of columns, with the plume under each closure and with its F from the buoyant layers alone, from either source, with either CAPE parcel and its downdraft, carrying momentum with or without a downdraft, the column momentum of each edge exact, with the shallow plume from the lowest layer, raining and overshooting by half without virtual buoyancy, under either autoconversion floor and with a shorter lifetime for the cloud above the shallow top', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  await parity({});
  await parity({ plumeClosure: 'maximum', plumeSource: 'lowest', plumeCapeParcel: 'undilute', plumeMassGrowth: 2e-4 });
  await parity({ plumeClosure: 'cape', downdraftShare: 0.5, downdraftEntrainment: 0, plumeRainThreshold: 5e-4, autoconversionFloor: 'boundaryLayer', upperCloudLifetime: 1800 });
  await parity({ plumeConsumption: 'buoyant', cloudLifetime: 7200 });
  await parity({ plumeMomentum: true }, { momentum: true });
  await parity({ plumeMomentum: true, downdraftShare: 0 }, { momentum: true });
  await parity({ autoconversionFloor: 'boundaryLayer' });
  await parity({ cumulusSource: 'lowest', cumulusRain: 5e-4, cumulusOvershoot: 0.5, virtualBuoyancy: false });
  await parity({ condensation: 'saturation', iceSaturation: false, iceFall: null });
  await parity({ iceNucleation: true, iceFall: 3.29 });
  await parity({ capeClosure: 'threshold' });
  await parity({ pcapeBoundary: 'signed' });
  await parity({ excessVelocity: 'convective' });
  await parity({ plumeSourceDepth: 'boundaryLayer' });
  await parity({ convectionType: 'top' });
  await parity({ convectionType: 'cloudDepth' });
  await parity({ plumePhase: 'liquid' });
  await parity({ plumeConversion: 'zhangMcFarlane' });
  await parity({ plumeEntrainmentLaw: 'gregory' });
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

test('a cumulusMemory other than a finite time of 0 s or more is refused on both engines', async () => {
  const refused = [-1, null, Infinity, '1800'];
  for (const value of refused) assert.throws(() => build({ cumulusMemory: value }), /cumulusMemory/, `CPU: ${value}`);
  if (!gpuAvailable) return;
  const { physicsConstants } = await import('../js/gpu/physics.gpu.js');
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  for (const value of refused) assert.throws(() => physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, cumulusMemory: value }), /cumulusMemory/, `GPU: ${value}`);
});

test('the retired Betts–Miller options are refused on both engines', async () => {
  for (const options of [{ convection: 'bettsMiller' }, { convection: 'plume' }, { shallowScheme: 'massFlux' }, { capeThreshold: 100 }]) assert.throws(() => build(options), /retired Betts–Miller/);
  if (!gpuAvailable) return;
  const { physicsConstants } = await import('../js/gpu/physics.gpu.js');
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  assert.throws(() => physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, convection: 'bettsMiller' }), /retired Betts–Miller/);
});
