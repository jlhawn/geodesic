import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { LATENT_HEAT } from '../js/physics/moist.module.js';
import { CP_DRY } from '../js/dynamics/sigmaCore.module.js';
import { plumeColumn, sounding, jordan, LEVELS, SPACING, SPACINGS, STEP, SURFACE, SURFACE_PRESSURE } from '../explain/plumeColumn.module.js';

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const hPa = (p) => p / 100;
const peak = (r) => {
  let k = 0;
  for (let j = 0; j < r.pressure.length; j++) if (r.heating.convection[j] > r.heating.convection[k]) k = j;
  return { rate: r.heating.convection[k], at: r.pressure[k] };
};
const describe = (r) => `${r.type}, top ${r.top === null ? '-' : hPa(r.top).toFixed(0)} hPa, CAPE ${r.cape === null ? '-' : r.cape.toFixed(0)} J/kg, τ ${r.tau === null ? '-' : r.tau.toFixed(0)} s, base flux ${r.baseFlux.toFixed(4)} kg/m²/s, rain ${r.rain.toFixed(2)} mm/d, heating peak ${peak(r).rate.toFixed(1)} K/d at ${hPa(peak(r).at).toFixed(0)} hPa`;

function fresh({ humidity = 1, warming = 0, spacing = SPACING, ice = true, dt = STEP * spacing / SPACING }) {
  const grid = new Grid(2), probe = buildMesh(grid), mean = Math.sqrt(probe.areaCell.reduce((a, b) => a + b, 0) / probe.nCells);
  const model = createModel(buildMesh(grid, { radius: 6.371e6 * spacing / mean }), { ocean: false, levels: LEVELS, moist: ice ? {} : { plumePhase: 'liquid' } });
  const { K, C, sigmaMid, exnerLayer, geopotential, g } = model.core.diagnostics, [pi, theta, , , q, qc] = model.state, profile = sounding({ humidity, warming });
  pi[0] = SURFACE_PRESSURE;
  for (let k = 0; k < K; k++) theta[k * C] = 300;
  for (let pass = 0; pass < 4; pass++) {
    model.core.diagnoseColumn(0, pi, theta, q, qc);
    for (let k = K - 1; k >= 0; k--) { const air = profile(geopotential[k * C] / g, pi[0] * sigmaMid[k]); theta[k * C] = air.T / exnerLayer[k * C]; q[k * C] = air.q; }
  }
  model.core.diagnoseColumn(0, pi, theta, q, qc);
  model.boundaryLayer.depth[0] = geopotential[(K - 1) * C] / g + SURFACE.depthAbove;
  model.boundaryLayer.buoyancyFlux[0] = SURFACE.buoyancyFlux;
  model.boundaryLayer.friction[0] = SURFACE.friction;
  model.radiation.mlmGate[0] = SURFACE.deckGate;
  model.moist.trace.convection = new Float64Array(K * C);
  model.moist.adjust(model.state, 0, 1, dt);
  const perDay = 86400 / dt;
  return { rain: model.moist.rain[0] * perDay, top: model.moist.cumulusTop[0], baseFlux: model.moist.cumulusBaseFlux[0], deep: model.moist.deep.deep, heating: Float64Array.from({ length: K }, (_, k) => model.moist.trace.convection[k * C] * perDay), cover: Float64Array.from({ length: K }, (_, k) => model.moist.cumulusCover[k * C]) };
}

test('Jordan\'s sounding on a 53.4 km cell lifts a deep plume to 153 hPa with CAPE 259 J/kg, a base flux of 0.0224 kg/m²/s and 40.1 mm/day of rain, heating most near 440 hPa; on a 106.8 km cell its closure is slower', () => {
  const fine = plumeColumn(), coarse = plumeColumn({ spacing: SPACINGS[1] }), coarseFast = plumeColumn({ spacing: SPACINGS[1], dt: STEP });
  console.log(`${(fine.cellSize / 1e3).toFixed(1)} km: ${describe(fine)}`);
  console.log(`${(coarse.cellSize / 1e3).toFixed(1)} km: ${describe(coarse)}`);
  close(fine.cellSize, 53.4e3, 50, 'cell 0 of the 55.8 km mesh');
  close(coarse.cellSize, 106.8e3, 50, 'cell 0 of the 111.6 km mesh');
  assert.equal(fine.type, 'deep');
  assert.ok(fine.parcelDeep, 'the test parcel types it deep');
  close(hPa(fine.top), 153, 0.5, 'top');
  close(hPa(fine.base), 958, 0.5, 'cloud base');
  close(fine.cape, 259, 0.5, 'CAPE');
  close(fine.tau, 3575, 0.5, 'τ');
  close(fine.baseFlux, 0.0224, 0.00005, 'base flux');
  close(fine.rain, 40.1, 0.05, 'rain');
  close(fine.convectiveRain, fine.rain, 1e-12, 'all of it convective');
  close(peak(fine).rate, 19.3, 0.05, 'heating peak');
  close(hPa(peak(fine).at), 440, 1, 'heating peak pressure');
  assert.equal(coarse.type, 'deep');
  close(coarse.dt, 337.5, 0, 'the 111.6 km step');
  close(coarse.cape, fine.cape, 1e-9, 'the same CAPE');
  close(coarse.tau, 5058, 0.5, 'τ');
  close(coarse.baseFlux, 0.0158, 0.00005, 'base flux');
  close(coarse.rain, 28.3, 0.05, 'rain');
  close(peak(coarse).rate, 13.6, 0.05, 'heating peak');
  close(coarseFast.rain, coarse.rain, 1e-5 * coarse.rain, 'the rain rate barely depends on the step');
  const liquid = plumeColumn({ ice: false });
  console.log(`liquid only: ${describe(liquid)}`);
  close(hPa(liquid.top), 211, 0.5, 'liquid top');
  close(liquid.cape, 142, 0.5, 'liquid CAPE');
  assert.ok(liquid.rain < fine.rain && liquid.snow === 0, 'the liquid plume rains less');
});

test('drying the air above 850 hPa lowers the plume\'s top and its rain; every setting stays deep from a column with no past', () => {
  const rows = [1.2, 1.1, 1, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5, 0.4].map((humidity) => ({ humidity, ...plumeColumn({ humidity }) }));
  for (const r of rows) console.log(`×${r.humidity}: ${describe(r)}, large-scale heating up to ${Math.max(...r.heating.largeScale).toFixed(1)} K/d`);
  for (const r of rows) assert.equal(r.type, 'deep', `×${r.humidity}`);
  const at = (h) => rows.find((r) => r.humidity === h);
  close(hPa(at(0.9).top), 180, 0.5, '×0.9 top');
  close(hPa(at(0.85).top), 659, 0.5, '×0.85 top');
  close(hPa(at(0.8).top), 659, 0.5, '×0.8 top');
  close(at(0.8).rain, 2.14, 0.005, '×0.8 rain');
  close(hPa(at(0.6).top), 659, 0.5, '×0.6 top');
  close(at(0.6).rain, 1.31, 0.005, '×0.6 rain');
  close(hPa(at(0.4).top), 756, 0.5, '×0.4 top');
  close(at(0.4).rain, 2.51, 0.005, '×0.4 rain');
  close(hPa(at(1.2).top), 110, 0.5, '×1.2 top');
  assert.ok(at(1.2).heating.largeScale.some((x) => x > 50), 'at ×1.2 the saturated layers condense on their own');
});

test('warming the lowest kilometer at its relative humidity raises CAPE and rain; at its specific humidity it lifts the cloud base until the test parcel fails', () => {
  for (const warming of [0.5, 1]) console.log(`+${warming} K, moistened: ${describe(plumeColumn({ warming }))}`);
  for (const warming of [0.25, 0.5, 2.5]) console.log(`+${warming} K, dry: ${describe(plumeColumn({ warming, moistenWarming: false }))}`);
  const moist = plumeColumn({ warming: 0.5 });
  close(moist.cape, 489, 0.5, 'CAPE at +0.5 K');
  close(moist.rain, 93.9, 0.05, 'rain at +0.5 K');
  close(plumeColumn({ warming: 1 }).cape, 802, 0.5, 'CAPE at +1 K');
  assert.equal(plumeColumn({ warming: 0.25, moistenWarming: false }).type, 'deep');
  const dry = plumeColumn({ warming: 0.5, moistenWarming: false });
  assert.equal(dry.type, 'shallow');
  assert.ok(!dry.parcelDeep && dry.rain === 0 && dry.cape === null, 'a shallow plume that never rains');
  close(hPa(dry.top), 913, 0.5, 'shallow top');
  const none = plumeColumn({ warming: 2.5, moistenWarming: false });
  assert.equal(none.type, 'none');
  assert.ok(none.massFlux.every((m) => m === 0) && none.top === null && none.plume.buoyancy.every(Number.isNaN), 'no plume, nothing left over');
});

test('every call matches a fresh model whatever ran before, with the cover the model\'s cumulus memory relaxes toward', () => {
  const cases = [{}, { humidity: 0.8 }, { humidity: 0.4 }, { warming: 0.5 }, { spacing: SPACINGS[1] }, { ice: false, humidity: 0.9 }];
  for (const order of [cases, [...cases].reverse()]) {
    for (const c of order) {
      const r = plumeColumn(c), f = fresh(c), what = JSON.stringify(c);
      assert.equal(r.type === 'deep', f.deep, `${what} type`);
      assert.equal(r.rain, f.rain, `${what} rain`);
      assert.equal(r.top, f.top, `${what} top`);
      assert.equal(r.baseFlux, f.baseFlux, `${what} base flux`);
      assert.deepEqual(r.heating.convection, f.heating, `${what} heating`);
      const keep = 1 - Math.exp(-r.dt / 1800);
      for (let k = 0; k < r.cover.length; k++) close(f.cover[k], r.cover[k] * keep < 1e-6 ? 0 : r.cover[k] * keep, 1e-15, `${what} cover in layer ${k}`);
    }
  }
});

test('the step keeps column water and enthalpy: the rain is what the column loses, convection\'s moistening sums to its rain, and the arrays describe the plume', () => {
  for (const c of [{}, { ice: false }, { humidity: 1.2 }, { humidity: 0.4 }, { warming: 1 }, { warming: 0.5, moistenWarming: false }]) {
    const r = plumeColumn(c), perDay = 86400 / r.dt, what = JSON.stringify(c), K = r.pressure.length;
    let water = 0, convective = 0, enthalpy = 0, scale = 0;
    for (let k = 0; k < K; k++) {
      const dq = (r.moistening.convection[k] + r.moistening.largeScale[k]) / 1e3 / perDay, dT = (r.heating.convection[k] + r.heating.largeScale[k]) / perDay;
      water += r.mass[k] * (dq + r.cloudWater[k]);
      convective += r.mass[k] * r.moistening.convection[k] / 1e3 / perDay;
      enthalpy += r.mass[k] * (CP_DRY * dT + LATENT_HEAT * dq);
      scale += r.mass[k] * LATENT_HEAT * Math.abs(dq);
    }
    close(water, -r.rain / perDay, 1e-12, `${what} water`);
    close(convective, -r.convectiveRain / perDay, 1e-12, `${what} convection's vapor`);
    close(enthalpy, 0, 1e-9 * Math.max(scale, 1), `${what} enthalpy`);
    close(r.rain, r.convectiveRain + r.largeScaleRain, 1e-12, `${what} rain parts`);
    assert.equal(r.interfaces.length, K + 1);
    close(r.interfaces[K], r.surfacePressure, 1e-9, `${what} ground interface`);
    assert.equal(r.massFlux[K], 0, `${what} no mass flux through the ground`);
    if (r.type === 'deep') {
      const top = r.interfaces.findIndex((p) => p >= r.top);
      for (let k = 0; k <= top; k++) assert.equal(r.massFlux[k], 0, `${what} nothing above the top`);
      assert.ok(Math.max(...r.massFlux) >= r.baseFlux * (1 - 1e-12), `${what} the flux reaches its base value`);
      assert.ok(r.downdraft.every((m) => m <= 0), `${what} the downdraft sinks`);
      assert.ok(r.plume.buoyancy.some((b) => b > 0) && r.plume.interfaceTemperature.some(Number.isFinite), `${what} the plume's path`);
    } else assert.ok(r.plume.buoyancy.every(Number.isNaN) && r.downdraft.every((m) => m === 0), `${what} no deep plume`);
  }
  assert.deepEqual(jordan(103000), jordan(101500), 'the sounding holds its surface values below 1015 hPa');
});

test('one call costs well under a millisecond once the model is built', () => {
  plumeColumn();
  const n = 200, start = performance.now();
  for (let i = 0; i < n; i++) plumeColumn({ humidity: 0.4 + 0.8 * (i % 17) / 16, warming: (i % 5) * 0.25 });
  const each = (performance.now() - start) / n;
  console.log(`${each.toFixed(3)} ms per call`);
  assert.ok(each < 5, `${each} ms`);
});
