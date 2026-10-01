import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { saturationHumidity, LATENT_HEAT, MOIST_DEFAULTS, DECK_CLOSED } from '../js/physics/moist.module.js';
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

test('with upperCloudLifetime cloud water in the layers above the shallow top converts over that lifetime and the cloud below over cloudLifetime', () => {
  const dt = 600, model = build({ upperCloudLifetime: 1800 }), plain = build(), { K, C, sigmaMid } = model.core.diagnostics;
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

function plumeColumn(options = {}, profile = null) {
  const model = build(options);
  if (profile) place(model, 0, 101500, profile); else jordanColumn(model, 0);
  setDepth(model, 0, 500);
  model.boundaryLayer.buoyancyFlux[0] = 4e-4;
  model.boundaryLayer.friction[0] = 0.25;
  model.radiation.mlmGate[0] = 0.3;
  return model;
}

test('the Jordan sounding lifts a deep plume that rains, heats most between 400 and 500 hPa above its cloud-base layer, cools the subcloud layer through its downdraft, and keeps column enthalpy and water exact', () => {
  const model = plumeColumn(), { moist } = model, [pi, theta, , , q, qc] = model.state;
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
    const model = plumeColumn({ plumeCape: 70, plumeEntrainment: 0.1 }, (z, p) => { const air = jordan(p); return { T: air.T, q: Math.min(air.q * (p < 850e2 ? factor : 1), saturationHumidity(air.T, p)) }; });
    model.moist.adjust(model.state, 0, 1, 600);
    return { top: model.moist.cumulusTop[0], deep: model.moist.deep.deep, cape: model.moist.deep.cape, rain: model.moist.rain[0] };
  };
  const moist = run(1), drier = run(0.8), dry = run(0.6), driest = run(0.4);
  console.log(`free-tropospheric humidity × 1, 0.8, 0.6, 0.4: plume tops ${[moist, drier, dry, driest].map((r) => (r.top / 100).toFixed(0)).join(', ')} hPa, CAPE ${[moist, drier, dry, driest].map((r) => r.cape.toFixed(0)).join(', ')} J/kg`);
  assert.ok(moist.deep && drier.deep && dry.deep && !driest.deep);
  assert.ok(moist.top < drier.top && drier.top < dry.top && dry.top < driest.top, 'the top sinks as the air dries');
  assert.ok(moist.rain > drier.rain && drier.rain > dry.rain && driest.rain === 0);
});

test('the deep base flux relaxes the CAPE toward plumeCape over plumeRelaxation, and plumeClosure: maximum gives it at least the shallow closure', () => {
  const flux = (options) => { const model = plumeColumn(options), [pi, theta, , , q, qc] = model.state; model.moist.plumeColumn(0, pi, theta, q, qc, 600); return { ...model.moist.deep }; };
  const hour = flux({}), twoHours = flux({ plumeRelaxation: 7200 }), lower = flux({ plumeCape: 0 });
  const cape0 = MOIST_DEFAULTS.plumeCape;
  assert.ok(Math.abs(hour.baseFlux - (hour.cape - cape0) / (3600 * hour.consumption)) < 1e-12 * hour.baseFlux, `base flux ${hour.baseFlux}`);
  assert.ok(Math.abs(twoHours.baseFlux - 0.5 * hour.baseFlux) < 1e-12 * hour.baseFlux, 'twice the relaxation time, half the flux');
  assert.ok(Math.abs(lower.baseFlux / hour.baseFlux - hour.cape / (hour.cape - cape0)) < 1e-9, 'the flux scales with the CAPE above plumeCape');
  const shallowOnly = flux({ plumeCape: 1e6 }), maximum = flux({ plumeCape: 1e6, plumeClosure: 'maximum' });
  assert.ok(!shallowOnly.deep && maximum.deep && maximum.baseFlux > 0, 'with CAPE below plumeCape the deep plume runs only under plumeClosure: maximum');
});

async function parity(options, { momentum = false } = {}) {
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
    for (let k = 0; k < K; k++) if (random() < 0.08) qc[k * C + i] = 1e-3 * random();
  }
  if (momentum) for (let x = 0; x < model.state[2].length; x++) model.state[2][x] = 20 * (random() - 0.5) + 10 * Math.sin(x / mesh.nEdges);
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
  for (const a of [model.boundaryLayer.depth, model.radiation.mlmGate, buoyancy, friction]) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  model.boundaryLayer.buoyancyFlux.set(buoyancy);
  model.boundaryLayer.friction.set(friction);
  const gpu = await createGpuCore(mesh, { levels, physics: options });
  const { device, buffers, kernels, layout } = gpu;
  gpu.upload(model.state);
  gpu.uploadPhysics({ mlmGate: model.radiation.mlmGate, concentration });
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.REGIME, Float32Array.from(regime));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.MIXTOP, Float32Array.from(mixingTop));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.STRAT, Float32Array.from(stratiform));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.DEPTH, Float32Array.from(model.boundaryLayer.depth));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.BUOY, Float32Array.from(buoyancy));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.USTAR, Float32Array.from(friction));
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
  moist.trace.convection = new Float64Array(K * C);
  model.phases.adjust(0, C, dt);
  if (momentum) {
    const E = mesh.nEdges, { dSigma, g } = core.diagnostics, u = model.state[2], before = Float64Array.from(u);
    model.phases.mixMomentum(0, E, dt);
    let moved = 0, worstU = 0, worstColumn = 0, scale = 0;
    for (let e = 0; e < E; e++) {
      const columnMass = 0.5 * (pi[mesh.cellsOnEdge[2 * e]] + pi[mesh.cellsOnEdge[2 * e + 1]]);
      let was = 0, now = 0, changed = false;
      for (let k = 0; k < K; k++) {
        const x = k * E + e, m = columnMass * dSigma[k] / g;
        was += m * before[x]; now += m * u[x]; scale = Math.max(scale, Math.abs(m * before[x]));
        if (u[x] !== before[x]) changed = true;
        worstU = Math.max(worstU, Math.abs(u[x] - after[2][x]));
      }
      if (changed) moved++;
      worstColumn = Math.max(worstColumn, Math.abs(now - was));
    }
    console.log(`  momentum: the plume moves the wind on ${moved} of ${E} edges; each edge's column momentum changes by at most ${worstColumn.toExponential(1)} kg/m/s against layer momenta up to ${scale.toFixed(0)}; the engines' winds differ by ${worstU.toExponential(1)} m/s`);
    assert.ok(moved > E / 10, `${moved} edges moved`);
    assert.ok(worstColumn < 1e-12 * scale, `column momentum changed by ${worstColumn}`);
    assert.ok(worstU < 1e-3, `winds differ by ${worstU} m/s`);
  }
  const { sigmaMid } = core.diagnostics;
  let deep = 0, shallow = 0, decked = 0, still = 0, opening = 0, flips = 0, worstTheta = 0, worstQ = 0, worstQc = 0, worstRain = 0, rainScale = 0;
  let plumes = 0, plumeFlips = 0, topsDiffer = 0, worstFlux = 0, fluxScale = 0, worstCover = 0, worstWater = 0, waterScale = 0;
  const K0 = K - (layout.PH.CUWATER - layout.PH.CUCOVER) / C;
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
    worstFlux = Math.max(worstFlux, Math.abs(cpu - gpuFlux)); fluxScale = Math.max(fluxScale, cpu);
    for (let k = K0; k < K; k++) {
      worstCover = Math.max(worstCover, Math.abs(moist.cumulusCover[k * C + i] - ph.CUCOVER[(k - K0) * C + i]));
      worstWater = Math.max(worstWater, Math.abs(moist.cumulusWater[k * C + i] - ph.CUWATER[(k - K0) * C + i]));
      waterScale = Math.max(waterScale, moist.cumulusWater[k * C + i]);
    }
    if ((moist.convectivePrecipitation[i] > 0) !== (ph.CONV[i] > 0)) flips++;
    worstRain = Math.max(worstRain, Math.abs(moist.convectivePrecipitation[i] - ph.CONV[i]), Math.abs(moist.largeScalePrecipitation[i] - ph.COND[i]));
    rainScale = Math.max(rainScale, moist.convectivePrecipitation[i]);
    for (let k = 0; k < K; k++) {
      const x = k * C + i;
      worstTheta = Math.max(worstTheta, Math.abs(model.state[1][x] - after[1][x]));
      worstQ = Math.max(worstQ, Math.abs(model.state[4][x] - after[4][x]));
      worstQc = Math.max(worstQc, Math.abs(model.state[5][x] - after[5][x]));
    }
  }
  console.log(`${JSON.stringify(options)}: ${C} random columns: ${deep} convect deep, ${shallow} shallow, ${decked} under a deck, ${opening} convect under a deck opening, ${still} still; convective rain differs in sign on ${flips}; engines differ in θ by ${worstTheta.toExponential(1)} K, q by ${worstQ.toExponential(1)}, qc by ${worstQc.toExponential(1)}, a step's rain by ${worstRain.toExponential(1)} kg/m² (largest ${rainScale.toFixed(3)})`);
  console.log(`  ${plumes} plumes, ${plumeFlips} differ in whether they rise, ${topsDiffer} in their top; base mass flux differs by ${worstFlux.toExponential(1)} kg/m²/s (largest ${fluxScale.toFixed(3)}), the cumulus fraction by ${worstCover.toExponential(1)}, the plume's condensate by ${worstWater.toExponential(1)} kg/kg (largest ${waterScale.toExponential(1)})`);
  assert.ok(deep > C / 40 && plumes > C / 10 && decked > C / 40 && still > C / 40, `${deep} deep, ${plumes} plumes, ${decked} decked, ${still} still`);
  assert.equal(plumeFlips, 0);
  assert.equal(topsDiffer, 0);
  assert.ok(worstFlux < 2e-3 * fluxScale && worstCover < 1e-4, `base mass flux ${worstFlux}, cover ${worstCover}`);
  assert.ok(waterScale > 0 && worstWater < 1e-3 * waterScale, `plume condensate ${worstWater} against ${waterScale}`);
  assert.equal(flips, 0);
  assert.ok(worstTheta < 1e-3 && worstQ < 1e-6 && worstQc < 1e-7, `θ ${worstTheta}, q ${worstQ}, qc ${worstQc}`);
  assert.ok(worstRain < 1e-4 * rainScale, `rain ${worstRain} against ${rainScale}`);
  return { cloud: Float64Array.from(model.state[5]), gpuCloud: Float64Array.from(after[5]) };
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
});

test('the stratiform lifetime matches between the engines on random columns of every regime, mixing top, EIS share and sea-ice cover, and keeps cloud the short lifetime would rain out', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const long = await parity({ cloudLifetime: 3600, stratiformLifetime: 3 * 3600 }), short = await parity({ cloudLifetime: 3600, stratiformLifetime: null });
  let cloudy = 0, kept = 0, gpuKept = 0;
  for (let x = 0; x < long.cloud.length; x++) {
    if (!(short.cloud[x] > 0)) continue;
    cloudy++;
    if (long.cloud[x] > short.cloud[x] * (1 + 1e-6)) kept++;
    if (long.gpuCloud[x] > short.gpuCloud[x] * (1 + 1e-6)) gpuKept++;
  }
  console.log(`after one step the 3 h stratiform lifetime keeps more cloud than the 1 h lifetime alone in ${kept} of ${cloudy} cloudy layers (GPU ${gpuKept})`);
  assert.ok(kept > cloudy / 10 && kept < cloudy && Math.abs(gpuKept - kept) <= cloudy / 100, `${kept} and ${gpuKept} of ${cloudy}`);
});

test('the retired Betts–Miller options are refused on both engines', async () => {
  for (const options of [{ convection: 'bettsMiller' }, { convection: 'plume' }, { shallowScheme: 'massFlux' }, { capeThreshold: 100 }]) assert.throws(() => build(options), /retired Betts–Miller/);
  if (!gpuAvailable) return;
  const { physicsConstants } = await import('../js/gpu/physics.gpu.js');
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  assert.throws(() => physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, convection: 'bettsMiller' }), /retired Betts–Miller/);
});
