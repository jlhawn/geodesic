import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { saturationHumidity, LATENT_HEAT, MOIST_DEFAULTS, DECK_CLOSED } from '../js/physics/moist.module.js';

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

test('a stratocumulus-topped column under a 1.5 K inversion with dry air above never convects', () => {
  const model = build(), { moist } = model, [pi, theta, , , q, qc] = model.state;
  place(model, 0, 101500, stratocumulus);
  setDepth(model, 0, 1300);
  const { K, C, sigmaMid, geopotential, g } = model.core.diagnostics;
  let cloudy = 0;
  for (let k = 0; k < K; k++) if (qc[k * C] > 0) cloudy++;
  const top = moist.diagnoseParcel(0, pi, theta, q);
  console.log(`stratocumulus: cloud in ${cloudy} layers, parcel CAPE ${moist.parcel.cape.toFixed(1)} J/kg, inhibition ${moist.parcel.inhibition.toFixed(1)} J/kg, top ${top >= 0 ? (pi[0] * sigmaMid[top] / 100).toFixed(0) + ' hPa' : 'none'}, cloud base at ${(geopotential[moist.parcel.base * C] / g).toFixed(0)} m`);
  assert.ok(cloudy >= 2, `a deck of ${cloudy} layers`);
  assert.ok(moist.parcel.cape < MOIST_DEFAULTS.capeThreshold, `CAPE ${moist.parcel.cape}`);
  const before = snapshot(model, 0);
  for (let n = 0; n < 40; n++) assert.equal(moist.convectColumn(0, pi, theta, q, 600, qc), 0, `step ${n} rains`);
  assert.deepEqual(snapshot(model, 0), before, 'the column is untouched');
  assert.ok(Math.abs(moist.activity[0] - 0.5 * Math.exp(-40 * 600 / MOIST_DEFAULTS.activityMemory)) < 1e-12, `activity ${moist.activity[0]}`);
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

test('a deep tropical column convects from cloud base up, heating most between 400 and 500 hPa and not at all below cloud base', () => {
  const model = build(), { moist } = model, [pi, theta, , , q, qc] = model.state;
  const { K, C, sigmaMid, exnerLayer } = model.core.diagnostics, dt = 600;
  jordanColumn(model, 0);
  setDepth(model, 0, 500);
  const before = snapshot(model, 0);
  const rain = moist.convectColumn(0, pi, theta, q, dt, qc);
  const { top, base, cape, inhibition, lclPressure } = moist.parcel;
  const heating = Float64Array.from({ length: K }, (_, k) => (theta[k * C] - before.theta[k]) * exnerLayer[k * C] / dt * 86400);
  let peak = top;
  for (let k = top; k <= base; k++) if (heating[k] > heating[peak]) peak = k;
  const hPa = (k) => pi[0] * sigmaMid[k] / 100;
  console.log(`Jordan (1958): CAPE ${cape.toFixed(0)} J/kg, inhibition ${inhibition.toFixed(1)} J/kg, condensation level ${(lclPressure / 100).toFixed(0)} hPa, cloud base ${hPa(base).toFixed(0)} hPa, top ${hPa(top).toFixed(0)} hPa; rain ${(rain * 86400 / dt).toFixed(1)} mm/d; heating ${Array.from({ length: base - top + 1 }, (_, n) => `${hPa(top + n).toFixed(0)}: ${heating[top + n].toFixed(1)}`).join(', ')} K/d`);
  assert.ok(rain > 0 && cape > 500 && inhibition < 10, `rain ${rain}, CAPE ${cape}, inhibition ${inhibition}`);
  assert.ok(hPa(top) < 300, `top at ${hPa(top)} hPa`);
  assert.ok(hPa(peak) > 400 && hPa(peak) < 500, `heating peaks at ${hPa(peak)} hPa`);
  assert.ok(pi[0] * levels[base + 1] >= lclPressure && pi[0] * levels[base] < lclPressure, 'cloud base holds the condensation level');
  for (let k = base + 1; k < K; k++) assert.ok(theta[k * C] === before.theta[k] && q[k * C] === before.q[k], `layer ${k} below cloud base changed`);
  for (let k = 0; k < top; k++) assert.ok(theta[k * C] === before.theta[k] && q[k * C] === before.q[k], `layer ${k} above the top changed`);
});

test('deep convection keeps column enthalpy and water exact through its rain, its anvil and its downdraft', () => {
  const model = build(), { moist } = model, [pi, theta, , , q, qc] = model.state;
  const dt = 600;
  jordanColumn(model, 0);
  setDepth(model, 0, 500);
  const before = budget(model, 0), snap = snapshot(model, 0);
  model.moist.adjust(model.state, 0, 1, dt);
  const after = budget(model, 0), rain = moist.rain[0], convective = moist.convectivePrecipitation[0], evaporated = moist.falling.evaporated;
  const { base } = moist.parcel, { K, C, exnerLayer } = model.core.diagnostics;
  let cooled = 0;
  for (let k = base + 1; k < K; k++) if (theta[k * C] < snap.theta[k]) cooled++;
  console.log(`deep: rain ${(rain * 86400 / dt).toFixed(1)} mm/d, convective ${(convective * 86400 / dt).toFixed(1)}, downdraft evaporation ${(evaporated * 86400 / dt).toFixed(1)} mm/d into ${cooled} subcloud layers; enthalpy off by ${((after.enthalpy - before.enthalpy) / before.enthalpy).toExponential(1)}, water by ${((after.water + rain - before.water) / before.water).toExponential(1)}`);
  assert.ok(convective > 0 && evaporated > 0 && cooled > 0, `convective rain ${convective}, downdraft ${evaporated}`);
  assert.ok(evaporated <= MOIST_DEFAULTS.downdraftEvaporation * (convective + evaporated) * (1 + 1e-12), 'the downdraft takes at most its share');
  assert.ok(Math.abs(after.enthalpy - before.enthalpy) < 1e-12 * before.enthalpy, `enthalpy ${before.enthalpy} → ${after.enthalpy}`);
  assert.ok(Math.abs(after.water + rain - before.water) < 1e-12 * before.water, `water ${before.water} → ${after.water} + ${rain}`);
  for (let k = base + 1; k < K; k++) assert.ok(theta[k * C] <= snap.theta[k], `subcloud layer ${k} warmed by ${(theta[k * C] - snap.theta[k]) * exnerLayer[k * C]} K`);
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

test('with the Betts–Miller shallow branch, shallowReference: mixingLine and no shallow rain, shallow convection mixes its cloud layer toward the mixing line and never rains, keeping heat and water exact', () => {
  const model = build({ shallowScheme: 'bettsMiller', capeThreshold: 10, shallowReference: 'mixingLine', shallowRain: false }), { moist } = model, [pi, theta, , , q, qc] = model.state;
  const { K, C, sigmaMid, exnerLayer, dSigma } = model.core.diagnostics, dt = 600;
  place(model, 0, 101500, tradeCumulus);
  setDepth(model, 0, 600);
  const before = budget(model, 0), snap = snapshot(model, 0);
  const rain = moist.convectColumn(0, pi, theta, q, dt, qc);
  const { top, base, cape } = moist.parcel, after = budget(model, 0);
  const hPa = (k) => pi[0] * sigmaMid[k] / 100;
  const dq = Array.from({ length: base - top + 1 }, (_, n) => (q[(top + n) * C] - snap.q[top + n]) / dt * 86400 * 1000);
  console.log(`trade cumulus: CAPE ${cape.toFixed(0)} J/kg, cloud base ${hPa(base).toFixed(0)} hPa, top ${hPa(top).toFixed(0)} hPa; moistening ${dq.map((x, n) => `${hPa(top + n).toFixed(0)}: ${x.toFixed(1)}`).join(', ')} g/kg/d; heat off by ${((after.heat - before.heat) / before.heat).toExponential(1)}, water by ${((after.vapour - before.vapour) / before.vapour).toExponential(1)}`);
  assert.ok(top >= 0 && hPa(top) > MOIST_DEFAULTS.shallowTop / 100 && cape > 10, `top ${hPa(top)} hPa, CAPE ${cape}`);
  assert.equal(rain, 0);
  assert.ok(dq.some((x) => x !== 0), 'the cloud layer changes');
  assert.ok(Math.abs(after.heat - before.heat) < 1e-12 * before.heat && Math.abs(after.vapour - before.vapour) < 1e-12 * before.vapour, 'heat and water unchanged');
  const rate = dt / MOIST_DEFAULTS.relaxationTime, reference = moist.reference;
  for (let k = top; k <= base; k++) {
    assert.ok(Math.abs(q[k * C] - snap.q[k] - rate * (reference.q[k] - snap.q[k])) < 1e-15, `layer ${k} moves toward the reference's humidity`);
    assert.ok(Math.abs((theta[k * C] - snap.theta[k]) * exnerLayer[k * C] - rate * (reference.T[k] - snap.theta[k] * exnerLayer[k * C])) < 1e-9, `layer ${k} moves toward the reference's temperature`);
  }
  const above = (top - 1) * C, spread = (k) => reference.q[k] - reference.q[top];
  assert.ok(spread(base) > 0 && reference.q[top] - reference.q[base] < 0, `the reference dries upward from cloud base (${reference.q[base]}) toward the air above the top (${q[above]}) at ${reference.q[top]}`);
  for (let k = 0; k < K; k++) if (k < top || k > base) assert.ok(theta[k * C] === snap.theta[k] && q[k * C] === snap.q[k], `layer ${k} outside the cloud layer changed`);
});

test('with the Betts–Miller shallow branch a trade-cumulus column vents at once on its shallow trigger, at a rate ramped by its CAPE, raining with enthalpy and water exact, but not without the trigger, below half its CAPE threshold, above its stability bound or under an active deck', () => {
  const dt = 600;
  const vent = (options = {}, gate = 0.3) => {
    const model = build({ shallowScheme: 'bettsMiller', ...options }), { moist } = model, [pi, theta, , , q, qc] = model.state;
    place(model, 0, 101500, tradeCumulus);
    setDepth(model, 0, 600);
    model.radiation.mlmGate[0] = gate;
    moist.activity[0] = 0;
    const before = budget(model, 0), snap = snapshot(model, 0);
    moist.diagnoseParcel(0, pi, theta, q);
    const stability = moist.inversionStrength(0, theta, q), { cape, inhibition } = moist.parcel;
    const rain = moist.convectColumn(0, pi, theta, q, dt, qc);
    return { model, moist, rain, before, snap, stability, cape, inhibition };
  };
  const fired = vent();
  const after = budget(fired.model, 0), stored = fired.before;
  let detrained = 0;
  const { K, C, dSigma, g } = fired.model.core.diagnostics, [pi, , , , , qc] = fired.model.state;
  for (let k = 0; k < K; k++) detrained += (qc[k * C] - fired.snap.qc[k]) * pi[0] * dSigma[k] / g;
  console.log(`trade cumulus: CAPE ${fired.cape.toFixed(0)} J/kg, inhibition ${fired.inhibition.toFixed(1)} J/kg, EIS ${fired.stability.toFixed(2)} K; vents at activity 0 with ${(fired.rain * 86400 / dt).toFixed(1)} mm/d of rain and ${(detrained * 86400 / dt).toFixed(1)} mm/d detrained`);
  assert.ok(fired.rain > 0 && fired.moist.activity[0] < 0.5, `rain ${fired.rain}, activity ${fired.moist.activity[0]}`);
  assert.ok(Math.abs(after.enthalpy - stored.enthalpy) < 1e-12 * stored.enthalpy, `enthalpy ${stored.enthalpy} → ${after.enthalpy}`);
  assert.ok(Math.abs(after.water + fired.rain - stored.water) < 1e-12 * stored.water, `water ${stored.water} → ${after.water} + ${fired.rain}`);
  for (const [options, gate, why] of [[{ shallowCape: null }, 0.3, 'without the trigger'], [{ shallowCape: 2 * fired.cape }, 0.3, 'at half its CAPE threshold'],
    [{ shallowStability: fired.stability - 0.1 }, 0.3, 'above its stability bound'], [{}, 0.8, 'under a deck']]) {
    const still = vent(options, gate);
    assert.equal(still.rain, 0, why);
    assert.deepEqual(snapshot(still.model, 0), still.snap, `${why}: the column is untouched`);
  }
  assert.ok(vent({ shallowStability: fired.stability + 0.1 }).rain > 0, 'below its stability bound it vents');
  const half = vent({ shallowCape: fired.cape }).rain, full = vent({ shallowCape: 1e-9 }).rain;
  assert.ok(Math.abs(half - 0.5 * full) < 1e-12 * full, `at its CAPE threshold it vents at half the rate: ${half} against ${full}`);
});

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
  model.moist.activity[0] = 0;
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

test('the deep branch convects a deep tropical column exactly as it did beside the Betts–Miller shallow branch', () => {
  const run = (shallowScheme) => {
    const model = build({ shallowScheme }), { moist } = model, [pi, theta, , , q, qc] = model.state;
    jordanColumn(model, 0);
    setDepth(model, 0, 500);
    moist.activity[0] = 1;
    const rain = moist.convectColumn(0, pi, theta, q, 600, qc);
    return { rain, snap: snapshot(model, 0), fired: moist.parcel.fired };
  };
  const massFlux = run('massFlux'), bettsMiller = run('bettsMiller');
  assert.ok(massFlux.rain > 0 && massFlux.fired, 'it fires');
  assert.equal(massFlux.rain, bettsMiller.rain);
  assert.deepEqual(massFlux.snap, bettsMiller.snap);
});

test('no column convects under an active deck, and its activity decays', () => {
  const model = build(), { moist } = model, [pi, theta, , , q, qc] = model.state;
  jordanColumn(model, 0);
  setDepth(model, 0, 500);
  model.radiation.mlmGate[0] = 0.8;
  const before = snapshot(model, 0);
  assert.equal(moist.convectColumn(0, pi, theta, q, 600, qc), 0);
  assert.deepEqual(snapshot(model, 0), before);
  assert.ok(moist.activity[0] < 0.5, `activity ${moist.activity[0]}`);
  model.radiation.mlmGate[0] = 0.5;
  moist.activity[0] = 0.5;
  assert.ok(moist.convectColumn(0, pi, theta, q, 600, qc) > 0, 'an undecided gate lets it fire');
});

test('the activity switches a column on only after it has passed for a while, and keeps it on for a while after', () => {
  const dt = 600, memory = MOIST_DEFAULTS.activityMemory, expected = Math.ceil(memory * Math.LN2 / dt);
  const model = build(), { moist } = model, [pi, theta, , , q, qc] = model.state;
  jordanColumn(model, 0);
  setDepth(model, 0, 500);
  const fresh = snapshot(model, 0);
  moist.activity[0] = 0;
  let on = -1;
  for (let n = 1; n <= 3 * expected && on < 0; n++) {
    const reset = () => { for (let k = 0; k < model.core.K; k++) { theta[k * model.mesh.nCells] = fresh.theta[k]; q[k * model.mesh.nCells] = fresh.q[k]; } };
    reset();
    if (moist.convectColumn(0, pi, theta, q, dt, qc) > 0) on = n;
  }
  assert.equal(on, expected, `on after ${on} steps`);
  const weak = build({ capeThreshold: 1e6 }), [pw, tw, , , qw, qcw] = weak.state;
  jordanColumn(weak, 0);
  setDepth(weak, 0, 500);
  weak.moist.activity[0] = 1;
  let off = -1;
  for (let n = 1; n <= 3 * expected && off < 0; n++) if (!(weak.moist.convectColumn(0, pw, tw, qw, dt, qcw) > 0)) off = n;
  assert.equal(off, expected, `off after ${off} steps once the column stops passing`);
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

async function parity(options) {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js');
  const model = build(options, 6), { moist, core, mesh } = model, C = mesh.nCells, { K } = core.diagnostics, dt = 900;
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
    moist.activity[i] = [0, 0.5, 1, random()][Math.floor(random() * 4)];
    for (let k = 0; k < K; k++) if (random() < 0.08) qc[k * C + i] = 1e-3 * random();
  }
  for (const a of model.state) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  for (const a of [model.boundaryLayer.depth, model.radiation.mlmGate, moist.activity, buoyancy, friction]) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  model.boundaryLayer.buoyancyFlux.set(buoyancy);
  model.boundaryLayer.friction.set(friction);
  const gpu = await createGpuCore(mesh, { levels, physics: options });
  const { device, buffers, kernels, layout } = gpu;
  gpu.upload(model.state);
  gpu.uploadPhysics({ mlmGate: model.radiation.mlmGate, convectiveActivity: moist.activity });
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.DEPTH, Float32Array.from(model.boundaryLayer.depth));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.BUOY, Float32Array.from(buoyancy));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.USTAR, Float32Array.from(friction));
  device.queue.writeBuffer(buffers.P, 0, Float32Array.from([dt, 0, 1, 0, 0, 0, 0, 0]));
  const group = device.createBindGroup({ layout: kernels.adjust.getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) });
  const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
  pass.setPipeline(kernels.adjust);
  pass.setBindGroup(0, group);
  pass.dispatchWorkgroups(Math.ceil(C / 64));
  pass.end();
  device.queue.submit([encoder.finish()]);
  const after = await gpu.download(), ph = await gpu.downloadPhysics();
  const activityBefore = Float64Array.from(moist.activity);
  moist.trace.convection = new Float64Array(K * C);
  model.phases.adjust(0, C, dt);
  const { sigmaMid } = core.diagnostics;
  const massFlux = (options.shallowScheme ?? MOIST_DEFAULTS.shallowScheme) === 'massFlux';
  let deep = 0, shallow = 0, decked = 0, still = 0, opening = 0, flips = 0, worstTheta = 0, worstQ = 0, worstQc = 0, worstActivity = 0, worstRain = 0, rainScale = 0;
  let plumes = 0, plumeFlips = 0, topsDiffer = 0, worstFlux = 0, fluxScale = 0, worstCover = 0;
  const K0 = K - (layout.PH.CUWATER - layout.PH.CUCOVER) / C;
  for (let i = 0; i < C; i++) {
    let top = -1;
    for (let k = K - 1; k >= 0; k--) if (moist.trace.convection[k * C + i] !== 0) top = k;
    if (model.radiation.mlmGate[i] >= DECK_CLOSED) decked++;
    else if (top < 0) still++;
    else if (model.radiation.mlmGate[i] > 0.5) opening++;
    else if (pi[i] * sigmaMid[top] > MOIST_DEFAULTS.shallowTop) shallow++;
    else deep++;
    if (massFlux) {
      const cpu = moist.cumulusBaseFlux[i], gpuFlux = ph.CUMF[i];
      if (cpu > 0) plumes++;
      if ((cpu > 0) !== (gpuFlux > 0)) plumeFlips++;
      if (moist.cumulusTop[i] > 0 && Math.abs(moist.cumulusTop[i] - ph.CUTOP[i]) > 1) topsDiffer++;
      worstFlux = Math.max(worstFlux, Math.abs(cpu - gpuFlux)); fluxScale = Math.max(fluxScale, cpu);
      for (let k = K0; k < K; k++) worstCover = Math.max(worstCover, Math.abs(moist.cumulusCover[k * C + i] - ph.CUCOVER[(k - K0) * C + i]));
    }
    if ((moist.convectivePrecipitation[i] > 0) !== (ph.CONV[i] > 0)) flips++;
    worstActivity = Math.max(worstActivity, Math.abs(moist.activity[i] - ph.CONVACT[i]));
    worstRain = Math.max(worstRain, Math.abs(moist.convectivePrecipitation[i] - ph.CONV[i]), Math.abs(moist.largeScalePrecipitation[i] - ph.COND[i]));
    rainScale = Math.max(rainScale, moist.convectivePrecipitation[i]);
    for (let k = 0; k < K; k++) {
      const x = k * C + i;
      worstTheta = Math.max(worstTheta, Math.abs(model.state[1][x] - after[1][x]));
      worstQ = Math.max(worstQ, Math.abs(model.state[4][x] - after[4][x]));
      worstQc = Math.max(worstQc, Math.abs(model.state[5][x] - after[5][x]));
    }
  }
  console.log(`${JSON.stringify(options)}: ${C} random columns: ${deep} convect deep, ${shallow} shallow, ${decked} under a deck, ${opening} convect under a deck opening, ${still} still; convective rain differs in sign on ${flips}; engines differ in θ by ${worstTheta.toExponential(1)} K, q by ${worstQ.toExponential(1)}, qc by ${worstQc.toExponential(1)}, the activity by ${worstActivity.toExponential(1)}, a step's rain by ${worstRain.toExponential(1)} kg/m² (largest ${rainScale.toFixed(3)})`);
  if (massFlux) {
    console.log(`  ${plumes} plumes, ${plumeFlips} differ in whether they rise, ${topsDiffer} in their top; base mass flux differs by ${worstFlux.toExponential(1)} kg/m²/s (largest ${fluxScale.toFixed(3)}), the cumulus fraction by ${worstCover.toExponential(1)}`);
    assert.ok(deep > C / 40 && plumes > C / 10 && decked > C / 40 && still > C / 40, `${deep} deep, ${plumes} plumes, ${decked} decked, ${still} still`);
    assert.equal(plumeFlips, 0);
    assert.equal(topsDiffer, 0);
    assert.ok(worstFlux < 2e-3 * fluxScale && worstCover < 1e-4, `base mass flux ${worstFlux}, cover ${worstCover}`);
  } else assert.ok(deep > C / 40 && shallow > C / 40 && decked > C / 40 && still > C / 40 && opening > C / 100, `${deep} deep, ${shallow} shallow, ${decked} decked, ${opening} opening, ${still} still`);
  assert.equal(flips, 0);
  assert.ok(worstTheta < 1e-3 && worstQ < 1e-6 && worstQc < 1e-7, `θ ${worstTheta}, q ${worstQ}, qc ${worstQc}`);
  assert.ok(worstActivity < 1e-5, `activity ${worstActivity}`);
  assert.ok(worstRain < 1e-4 * rainScale, `rain ${worstRain} against ${rainScale}`);
}

test('the triggered convection, the cumulus mass flux and the rain they leave match between the engines on a random set of columns, with the plume on its defaults, from the lowest layer, raining, overshooting by half or kept out of deep columns, under either autoconversion floor, and with the Betts–Miller shallow branch under either shallow reference with or without shallow rain and the shallow stability veto', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  await parity({});
  await parity({ autoconversionFloor: 'boundaryLayer' });
  await parity({ cumulusSource: 'lowest', cumulusRain: 5e-4, cumulusOvershoot: 0.5, cumulusWithDeep: false, virtualBuoyancy: false });
  const bettsMiller = { shallowScheme: 'bettsMiller' };
  await parity(bettsMiller);
  await parity({ ...bettsMiller, shallowReference: 'mixingLine', shallowRain: false, boundaryParcel: true, parcelDepth: 50e2, downdraftEvaporation: 0.25, downdraftSpread: 'fall' });
  await parity({ ...bettsMiller, shallowReference: 'mixingLine', shallowStability: 2 });
  await parity({ ...bettsMiller, shallowRain: false });
});
