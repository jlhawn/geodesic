import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { saturationHumidity, LATENT_HEAT, createMoistPhysics, cloudSaturation, criticalHumidityAt, LIQUID_TEMPERATURE, ICE_TEMPERATURE, KESSLER_EVAPORATION } from '../js/physics/moist.module.js';

const model = createModel(new Grid(16));
const { core, moist } = model;
const { K, C, dSigma, sigmaMid, cp, g, exnerLayer } = core.diagnostics;

function column(surfaceT, humidity, lapse = 6.5e-3) {
  const pi = new Float64Array(C).fill(101325);
  const theta = new Float64Array(K * C);
  const q = new Float64Array(K * C);
  const qc = new Float64Array(K * C);
  core.diagnose(pi, theta);
  for (let i = 0; i < C; i++) {
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      const p = pi[i] * sigmaMid[k];
      const height = 7000 * Math.log(pi[i] / p);
      const temperature = Math.max(200, surfaceT - lapse * height);
      theta[idx] = temperature / exnerLayer[idx];
      q[idx] = humidity * saturationHumidity(temperature, p);
    }
  }
  return [pi, theta, q, qc];
}

const moistEnthalpy = (pi, theta, q, i) => {
  let sum = 0;
  for (let k = 0; k < K; k++) { const idx = k * C + i; sum += (cp * theta[idx] * exnerLayer[idx] + LATENT_HEAT * q[idx]) * pi[i] * dSigma[k] / g; }
  return sum;
};

test('saturation humidity is near 22 g/kg at 300 K and 1000 hPa and falls with cold and pressure', () => {
  assert.ok(Math.abs(saturationHumidity(300, 1e5) - 0.0223) < 5e-4);
  assert.ok(saturationHumidity(273.15, 1e5) < 0.004 && saturationHumidity(273.15, 1e5) > 0.0035);
  assert.ok(saturationHumidity(280, 5e4) > saturationHumidity(280, 1e5));
});

test('saturation adjustment moves exactly the supersaturation into cloud water and conserves moist enthalpy and water', () => {
  const [pi, theta, q, qc] = column(300, 1.3);
  const before = moistEnthalpy(pi, theta, q, 0), waterBefore = moist.columnWater(pi, q, 0);
  core.diagnoseColumn(0, pi, theta, q, qc);
  const formed = moist.condenseColumn(0, pi, theta, q, qc);
  assert.ok(formed > 0);
  for (let k = 0; k < K; k++) {
    const idx = k * C;
    const qs = saturationHumidity(theta[idx] * exnerLayer[idx], pi[0] * sigmaMid[k]);
    assert.ok(q[idx] <= qs * (1 + 1e-3), `layer ${k} still supersaturated: ${q[idx]} > ${qs}`);
    assert.ok(qc[idx] >= 0);
  }
  assert.ok(Math.abs(moist.columnWater(pi, qc, 0) - formed) < 1e-12 * waterBefore);
  assert.ok(Math.abs(moist.columnWater(pi, q, 0) + moist.columnWater(pi, qc, 0) - waterBefore) < 1e-12 * waterBefore);
  assert.ok(Math.abs(moistEnthalpy(pi, theta, q, 0) - before) < 1e-9 * before);
  // dry the air again: the cloud evaporates back, cooling the layers
  for (let k = 0; k < K; k++) q[k * C] *= 0.5;
  const thetaBefore = Float64Array.from(theta), cloudBefore = moist.columnWater(pi, qc, 0);
  const evaporated = -moist.condenseColumn(0, pi, theta, q, qc);
  assert.ok(evaporated > 0 && evaporated <= cloudBefore * (1 + 1e-12));
  for (let k = 0; k < K; k++) assert.ok(theta[k * C] <= thetaBefore[k * C]);
});

test('the cloud layer atop a cloud-topped boundary layer holds the uniform distribution\'s condensate, and the mixed layers below it adjust to saturation', () => {
  const [pi, theta, q, qc] = column(290, 0.5);
  const k0 = sigmaMid.findIndex((s) => s > 0.85);
  core.diagnoseColumn(0, pi, theta, q, qc);
  const humid = (k) => cloudSaturation(theta[k * C] * exnerLayer[k * C], pi[0] * sigmaMid[k], true, LIQUID_TEMPERATURE, ICE_TEMPERATURE, { qs: 0, slope: 0, liquid: 1 }).qs * 0.96;
  q[k0 * C] = humid(k0); q[(k0 + 1) * C] = humid(k0 + 1);
  const T = theta[k0 * C] * exnerLayer[k0 * C], p = pi[0] * sigmaMid[k0];
  const qs = 611.2 * Math.exp(17.67 * (T - 273.15) / (T - 29.65)) * 0.622 / (p - 0.378 * 611.2 * Math.exp(17.67 * (T - 273.15) / (T - 29.65)));
  const a = 1 / (1 + LATENT_HEAT / cp * qs * LATENT_HEAT / (287.06 / 0.622 * T * T));
  const rhc = 0.75 + 0.225 * Math.exp(1 - (pi[0] / p) ** 2), b = a * (1 - rhc) * qs, Q = a * (q[k0 * C] - qs);
  const expected = (Q + b) * (Q + b) / (4 * b);
  const top = new Float64Array(C).fill(1e5), cloudLayer = new Float64Array(C).fill(k0), none = new Float64Array(C).fill(-1);
  const condensed = (boundaryCondensation, layer) => {
    const physics = createMoistPhysics(model.mesh, core, { stabilityLayer: model.radiation.stabilityLayer, boundaryTop: top, boundaryCloudLayer: layer, boundaryCondensation });
    const t = Float64Array.from(theta), v = Float64Array.from(q), c = Float64Array.from(qc);
    physics.condenseColumn(0, pi, t, v, c);
    for (let k = 0; k < K; k++) {
      const idx = k * C, ex = exnerLayer[idx];
      assert.ok(Math.abs(t[idx] - LATENT_HEAT * c[idx] / (cp * ex) - theta[idx]) <= 1e-13 * theta[idx] && Math.abs(v[idx] + c[idx] - q[idx]) <= 1e-15, `layer ${k} keeps θ_l and q_t`);
    }
    return [c[k0 * C], c[(k0 + 1) * C]];
  };
  assert.equal(criticalHumidityAt(p, pi[0], 0.975, 0.75, 2), rhc);
  const [cloud, below] = condensed('cloudLayer', cloudLayer), [old, oldBelow] = condensed('saturation', cloudLayer), [lost] = condensed('cloudLayer', none), [whole, wholeBelow] = condensed('uniform', none);
  console.log(`at ${(p / 100).toFixed(0)} hPa and ${T.toFixed(2)} K, humidity 0.96 of q_s ${(1000 * qs).toFixed(4)} g/kg under RH_c ${rhc.toFixed(4)}: a ${a.toFixed(4)}, Q ${(1000 * Q).toFixed(4)}, b ${(1000 * b).toFixed(4)} g/kg, the cloud layer holds ${(1000 * expected).toFixed(5)} g/kg (model ${(1000 * cloud).toFixed(5)})`);
  assert.ok(Math.abs(cloud - expected) <= 1e-12 * expected && expected > 1e-5, `the cloud layer holds (Q + b)²/(4b): ${cloud} against ${expected}`);
  assert.ok(below === 0 && old === 0 && oldBelow === 0 && lost === 0, 'a mixed layer outside the cloud layer, and every mixed layer under saturation, holds no cloud below saturation');
  assert.ok(Math.abs(whole - expected) <= 1e-12 * expected && wholeBelow > 0, 'with \'uniform\' every mixed layer holds the distribution\'s condensate');
});

test('autoconversion rains out cloud water above the threshold and conserves water', () => {
  const [pi, theta, q, qc] = column(290, 0.5);
  qc.fill(0);
  qc[(K - 3) * C] = 1e-3;
  qc[(K - 5) * C] = 1e-4;
  core.diagnoseColumn(0, pi, theta, q, qc);
  const cloudBefore = moist.columnWater(pi, qc, 0), before = cloudBefore + moist.columnWater(pi, q, 0);
  const rain = moist.autoconvertColumn(0, pi, theta, q, qc, 900);
  assert.ok(rain >= 0 && moist.columnWater(pi, qc, 0) < cloudBefore);
  assert.ok(Math.abs(before - moist.columnWater(pi, qc, 0) - moist.columnWater(pi, q, 0) - rain) < 1e-12 * before);
  assert.ok(qc[(K - 3) * C] < 1e-3 && qc[(K - 3) * C] > 1.5e-4, 'the thick cloud converts toward the threshold');
  assert.ok(qc[(K - 5) * C] < 1e-4 && qc[(K - 5) * C] > 0.7e-4, 'the thin cloud only decays over the cloud lifetime');
});

test('rain evaporates into the dry layers it falls through, conserving water and moist enthalpy, and never oversaturates them', () => {
  const cloudAt = K - 12, converted = (m, mCore = core) => {
    const [pi, theta, q, qc] = column(300, 0.3);
    qc[cloudAt * C] = 2e-3;
    mCore.diagnoseColumn(0, pi, theta, q, qc);
    const water = moist.columnWater(pi, q, 0) + moist.columnWater(pi, qc, 0), enthalpy = moistEnthalpy(pi, theta, q, 0), qBefore = Float64Array.from(q);
    const rain = m.autoconvertColumn(0, pi, theta, q, qc, 900);
    return { pi, theta, q, qc, qBefore, rain, water, enthalpy };
  };
  const { moist: noEvaporation, core: noEvaporationCore } = createModel(new Grid(16), { moist: { rainEvaporation: 0 } });
  const reference = converted(noEvaporation, noEvaporationCore);
  const r = converted(moist);
  assert.ok(reference.rain > 1e-3, `the cloud converts ${reference.rain} kg/m²`);
  assert.ok(r.rain < 0.5 * reference.rain, `rain reaching the ground ${r.rain} of ${reference.rain}`);
  assert.ok(Math.abs(moist.columnWater(r.pi, r.q, 0) + moist.columnWater(r.pi, r.qc, 0) + r.rain - r.water) < 1e-12 * r.water);
  assert.ok(Math.abs(moistEnthalpy(r.pi, r.theta, r.q, 0) - r.enthalpy) < 1e-12 * r.enthalpy);
  for (let k = 0; k < K; k++) {
    const idx = k * C, T = r.theta[idx] * exnerLayer[idx], qs = saturationHumidity(T, r.pi[0] * sigmaMid[k]);
    if (k <= cloudAt) assert.equal(r.q[idx], r.qBefore[idx], `no evaporation at or above the cloud (layer ${k})`);
    assert.ok(r.q[idx] <= qs * (1 + 1e-9), `layer ${k} oversaturated`);
  }
});

test("under rainEvaporationLaw 'kessler' each cloud-free layer takes the IFS's Kessler rate of the rain falling into it, integrated exactly over the step, conserving water and moist enthalpy: rain the instant law evaporates whole reaches the ground, in nearly the same share at half the step, and rainEvaporationCeiling holds every layer it reaches below c + (1 − c) a_P under either law", () => {
  const cloudAt = K - 12, cover = 0.3;
  const run = (options, humidity, dt = 900) => {
    const m = createModel(new Grid(16), { moist: options }), d = m.core.diagnostics;
    const [pi, theta, q, qc] = column(300, humidity);
    q[(K - 2) * C] = 0.9 * saturationHumidity(theta[(K - 2) * C] * exnerLayer[(K - 2) * C], pi[0] * sigmaMid[K - 2]);
    qc[cloudAt * C] = 2e-3;
    m.core.diagnoseColumn(0, pi, theta, q, qc);
    const before = { theta: Float64Array.from(theta), q: Float64Array.from(q) };
    const water = m.moist.columnWater(pi, q, 0) + m.moist.columnWater(pi, qc, 0), enthalpy = moistEnthalpy(pi, theta, q, 0);
    const rain = m.moist.autoconvertColumn(0, pi, theta, q, qc, dt);
    return { pi, theta, q, qc, before, rain, water, enthalpy, exner: d.exnerLayer };
  };
  const law = { rainEvaporationLaw: 'kessler', rainCover: cover }, dt = 900;
  const made = run({ rainEvaporation: 0 }, 0.8), instant = run({}, 0.8), kessler = run(law, 0.8);
  assert.ok(made.before.theta[cloudAt * C] * made.exner[cloudAt * C] > LIQUID_TEMPERATURE && made.rain > 1e-3, `the cloud layer is warm and converts ${made.rain} kg/m²`);
  let rain = made.rain, worst = 0, taken = 0;
  const saturated = { qs: 0, slope: 0, liquid: 1 };
  for (let k = cloudAt + 1; k < K; k++) {
    const idx = k * C, mass = kessler.pi[0] * dSigma[k] / g, p = kessler.pi[0] * sigmaMid[k];
    const { qs, slope } = cloudSaturation(kessler.before.theta[idx] * kessler.exner[idx], p, true, LIQUID_TEMPERATURE, ICE_TEMPERATURE, saturated);
    const gamma = LATENT_HEAT * slope / cp, deficit = Math.max(0, (qs - kessler.before.q[idx]) / (1 + gamma)) * mass;
    const rate = cover * KESSLER_EVAPORATION.rate * Math.pow(Math.sqrt(p / KESSLER_EVAPORATION.pressure) * rain / dt / (cover * KESSLER_EVAPORATION.flux), KESSLER_EVAPORATION.exponent);
    const evaporated = Math.min(rain, deficit, -deficit * Math.expm1(-rate * (1 + gamma) * dt));
    rain = Math.max(0, rain - evaporated);
    worst = Math.max(worst, Math.abs(kessler.q[idx] - kessler.before.q[idx] - evaporated / mass) / (evaporated / mass || 1));
    if (deficit > 0) taken = Math.max(taken, evaporated / deficit);
  }
  const half = { made: run({ rainEvaporation: 0 }, 0.8, dt / 2), kessler: run(law, 0.8, dt / 2) };
  const share = kessler.rain / made.rain, halfShare = half.kessler.rain / half.made.rain;
  console.log(`a ${(1e3 * made.rain).toFixed(1)} g/m² shower over ${dt} s through air at 80 % relative humidity: ${(1e3 * instant.rain).toFixed(1)} g/m² reaches the ground under the instant law, ${(1e3 * kessler.rain).toFixed(1)} (${share.toFixed(3)} of it) under the Kessler rate with a_P ${cover}, its humidity within ${worst.toExponential(1)} of the formula's, each layer taking at most ${taken.toFixed(3)} of its deficit; over ${dt / 2} s ${halfShare.toFixed(3)} of a ${(1e3 * half.made.rain).toFixed(1)} g/m² shower`);
  assert.ok(worst < 1e-9, `the layers' humidity departs from the Kessler formula by ${worst} of what evaporated`);
  assert.ok(Math.abs(kessler.rain - rain) < 1e-12 * made.rain, `rain reaching the ground ${kessler.rain} against ${rain}`);
  assert.ok(taken > 0 && taken < 0.5, `the rate takes ${taken} of a layer's deficit in one step`);
  assert.ok(instant.rain < 1e-3 * made.rain && share > 0.2, `${kessler.rain} reaches the ground under the Kessler rate, ${instant.rain} under the instant law, of ${made.rain}`);
  assert.ok(Math.abs(share - halfShare) < 0.05, `the share reaching the ground is ${share} over ${dt} s and ${halfShare} over ${dt / 2} s`);
  for (const r of [instant, kessler]) {
    assert.ok(Math.abs(r.water - (moist.columnWater(r.pi, r.q, 0) + moist.columnWater(r.pi, r.qc, 0) + r.rain)) < 1e-12 * r.water);
    assert.ok(Math.abs(moistEnthalpy(r.pi, r.theta, r.q, 0) - r.enthalpy) < 1e-12 * r.enthalpy);
  }
  const ceiling = 0.7, limit = ceiling + (1 - ceiling) * cover;
  for (const kind of ['deficit', 'kessler']) {
    const r = run({ rainEvaporationLaw: kind, rainCover: cover, rainEvaporationCeiling: ceiling }, 0.72);
    let reached = 0, highest = 0, above = 0;
    for (let k = cloudAt + 1; k < K; k++) {
      const idx = k * C, p = r.pi[0] * sigmaMid[k];
      const humidity = r.q[idx] / cloudSaturation(r.theta[idx] * r.exner[idx], p, true, LIQUID_TEMPERATURE, ICE_TEMPERATURE, saturated).qs;
      const start = r.before.q[idx] / cloudSaturation(r.before.theta[idx] * r.exner[idx], p, true, LIQUID_TEMPERATURE, ICE_TEMPERATURE, saturated).qs;
      if (start >= limit) { above++; assert.equal(r.q[idx], r.before.q[idx], `${kind}: layer ${k} starts at relative humidity ${start}, above the ceiling ${limit}, and takes no rain`); }
      else { assert.ok(humidity <= limit * (1 + 1e-9), `${kind}: layer ${k} reaches relative humidity ${humidity}, past the ceiling ${limit}`); highest = Math.max(highest, humidity); if (humidity > limit - 5e-3) reached++; }
    }
    console.log(`  ${kind} under rainEvaporationCeiling ${ceiling} through air at 72 %: the moistest layer below the cloud reaches relative humidity ${highest.toFixed(4)} against the ceiling ${limit.toFixed(2)}, ${reached} within 5e-3 of it, ${above} starting above it untouched; ${(1e3 * r.rain).toFixed(1)} g/m² reaches the ground`);
    assert.ok(above >= 1, `${above} layers start above the ceiling`);
    if (kind === 'deficit') assert.ok(reached >= 1, `the instant law takes ${reached} layers to the ceiling`);
  }
});

test('rain falling through saturated air reaches the ground whole', () => {
  const [pi, theta, q, qc] = column(295, 1);
  qc[(K - 12) * C] = 2e-3;
  core.diagnoseColumn(0, pi, theta, q, qc);
  const before = moist.columnWater(pi, qc, 0), qBefore = Float64Array.from(q);
  const rain = moist.autoconvertColumn(0, pi, theta, q, qc, 900);
  assert.ok(Math.abs(before - moist.columnWater(pi, qc, 0) - rain) < 1e-12 * before);
  assert.deepEqual(q, qBefore);
});

test('the plume leaves a stable dry column alone', () => {
  const [pi, theta, q, qc] = column(280, 0.2, 5e-3);
  const copy = Float64Array.from(theta), qCopy = Float64Array.from(q);
  core.diagnoseColumn(0, pi, theta, q, qc);
  assert.equal(moist.plumeColumn(0, pi, theta, q, qc, 450), 0);
  assert.deepEqual(theta, copy);
  assert.deepEqual(q, qCopy);
});

test('transport alone conserves total water', () => {
  const dry = createModel(new Grid(16), { physics: false });
  const init = initializeState(dry, {});
  for (let a = 0; a < init.length; a++) dry.state[a].set(init[a]);
  dry.state[5].set(dry.state[4].map((x) => 0.05 * x));
  const total = () => { const d = dry.diagnostics(); return d.columnWater + d.columnCloud; };
  const water0 = total();
  for (let n = 0; n < 20; n++) dry.step(900);
  const water1 = total();
  assert.ok(Math.abs(water1 - water0) < 1e-10 * water0, `water drifted ${water0} → ${water1}`);
  assert.equal(dry.moist.budget.lost, 0);
});

test('with sources on, the water column change over a day matches evaporation minus precipitation', () => {
  const wet = createModel(new Grid(16));
  const init = initializeState(wet, {});
  for (let a = 0; a < init.length; a++) wet.state[a].set(init[a]);
  const dt = 1350, steps = 64;
  const total = (d) => d.columnWater + d.columnCloud;
  const water0 = total(wet.diagnostics());
  let evaporated = 0, precipitated = 0, diag = null;
  for (let n = 0; n < steps; n++) { wet.step(dt); diag = wet.diagnostics(); evaporated += diag.evaporation * dt; precipitated += diag.precipitation * dt; }
  const water1 = total(diag);
  let area = 0; for (let i = 0; i < wet.mesh.nCells; i++) area += wet.mesh.areaCell[i];
  const rained = (wet.moist.budget.condensation + wet.moist.budget.convection) / area;
  const residual = (water1 - water0) - (evaporated - precipitated);
  console.log(`one day at N=16: evaporation ${evaporated.toFixed(3)} mm, precipitation ${precipitated.toFixed(3)} mm (budget ${rained.toFixed(3)}), water ${water0.toFixed(2)} → ${water1.toFixed(2)} mm of which cloud ${(1000 * diag.columnCloud).toFixed(0)} g/m², residual ${residual.toExponential(2)} mm, lost ${(wet.moist.budget.lost / area).toExponential(1)}`);
  assert.ok(evaporated > 0.5 && evaporated < 15, `evaporation ${evaporated} mm/day out of range`);
  assert.ok(Math.abs(rained - precipitated) < 1e-9 * (1 + rained));
  assert.ok(wet.moist.budget.fog >= 0 && wet.moist.budget.fog <= wet.moist.budget.condensation, `fog ${wet.moist.budget.fog} of the large-scale ${wet.moist.budget.condensation} kg`);
  assert.ok(Math.abs(residual) < 1e-5 * evaporated, `residual ${residual} vs evaporation ${evaporated}`);
  assert.ok(Number.isFinite(diag.maxWind) && diag.maxWind < 80);
});
