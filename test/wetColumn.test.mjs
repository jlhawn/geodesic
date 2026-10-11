import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore } from '../js/dynamics/sigmaCore.module.js';
import { createRadiation } from '../js/physics/radiation.module.js';
import { saturationHumidity, LATENT_HEAT } from '../js/physics/moist.module.js';
import { SEA_DRAG } from '../js/physics/surface.module.js';
import { makeColumn, readColumn, manabeWetherald, DRY_STRATOSPHERE, LEVELS, K } from '../explain/radiativeColumn.module.js';

const MLS = readColumn(JSON.parse(readFileSync(new URL('../explain/data/mls_bl36.json', import.meta.url))));
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const sum = (values) => values.reduce((a, b) => a + b, 0);
const STEP = 3 * 3600;
const layers = (terms) => sum(terms.sunlight) + sum(terms.infrared) + sum(terms.ground) + sum(terms.convection);

function settle(column, convection = 'moist', days = 3000) {
  let calm = 0;
  for (let n = 0; n < days * 8 && calm < 30 * 8; n++) {
    column.step(STEP, { convection });
    calm = Math.abs(column.budget.imbalance) < 0.005 && Math.abs(column.budget.groundNet) < 0.005 ? calm + 1 : 0;
  }
  assert.ok(calm >= 30 * 8, 'the column settles');
  return column;
}

const wet = settle(makeColumn({ albedo: 0.3, humidity: 'relative' }));
const vapor = Float64Array.from(wet.q);

test('fixed relative humidity follows Manabe and Wetherald\'s profile under 100 hPa and is capped at 3 ppm by mass above', () => {
  const column = makeColumn({ temperature: MLS.temperature, humidity: 'relative' });
  for (let k = 0; k < K; k++) {
    const p = column.pressure[k], full = manabeWetherald(p / column.settings.surfacePressure) * saturationHumidity(column.T[k], p);
    close(column.q[k], p < DRY_STRATOSPHERE.pressure ? Math.min(full, 3e-6) : full, 1e-15, `layer ${k}`);
  }
  close(manabeWetherald(1), 0.77, 1e-15, 'at the ground');
  close(manabeWetherald(0.02), 0, 0, 'at σ = 0.02');
  for (let n = 0; n < 8; n++) column.step(STEP, { convection: 'moist' });
  for (let k = 0; k < K; k++) if (column.pressure[k] >= 100e2) close(column.q[k], manabeWetherald(column.pressure[k] / column.settings.surfacePressure) * saturationHumidity(column.T[k], column.pressure[k]), 1e-15, `layer ${k} follows the new temperature`);
  column.evaluate();
  close(column.budget.precipitableWater, sum(Array.from(column.q, (q, k) => q * (column.interfacePressure[k + 1] - column.interfacePressure[k]) / 9.806)), 1e-12, 'precipitable water');
});

test('a fixed q profile is a copy, held while the temperature changes', () => {
  const q = Float64Array.from(MLS.humidity), column = makeColumn({ profile: MLS, humidity: q });
  q.fill(0);
  for (let n = 0; n < 8; n++) column.step(STEP, { convection: 'moist' });
  for (let k = 0; k < K; k++) assert.equal(column.q[k], MLS.humidity[k]);
  close(column.budget.precipitableWater, 28.85, 0.005, 'the midlatitude summer\'s precipitable water, kg/m²');
});

test('under moist convection the ground evaporates by the model\'s bulk formula and pays L·E', () => {
  const column = makeColumn({ profile: MLS, humidity: 'relative' }), b = K - 1;
  column.evaluate({ latent: false });
  const dry = { ...column.budget };
  column.evaluate({ latent: true });
  const { budget } = column, density = column.pressure[b] / (287.06 * column.T[b]);
  close(budget.evaporation, density * SEA_DRAG * 5 * (saturationHumidity(column.Ts, column.settings.surfacePressure) - column.q[b]), 1e-12, 'E');
  close(budget.latent, LATENT_HEAT * budget.evaporation, 1e-9, 'L·E');
  close(dry.groundNet - budget.groundNet, budget.latent, 1e-9, 'the ground pays L·E');
  close(dry.evaporation, 0, 0, 'no evaporation without moist convection');
  close(budget.absorbed - budget.olr, layers(column.terms) + budget.groundNet + budget.latent, 1e-9, 'evaluate: L·E is in transit');
  column.set({ evaporation: false });
  close(column.budget.evaporation, 0, 0, 'evaporation can be switched off');
  for (const settings of [{ humidity: null }, { humidity: 'relative', evaporation: true }]) {
    const c = makeColumn({ profile: MLS, ...settings });
    c.step(STEP, { convection: 'moist' });
    close(c.budget.absorbed - c.budget.olr, layers(c.terms) + c.budget.groundNet, 1e-6, `a moist step closes with ${settings.humidity}`);
    close(sum(c.terms.convection), c.budget.latent, 1e-6, 'convection moves heat and releases L·E');
    if (!settings.humidity) close(c.budget.evaporation, 0, 0, 'dry air does not evaporate');
  }
});

test('the wet column at ground albedo 0.3 settles at 286.6 K: G 133.2 W/m², 12.5 kg/m² of vapor, L·E 59.9 and sensible heat 17.0', () => {
  const { budget } = wet;
  close(wet.Ts, 286.6, 0.3, 'Ts');
  close(budget.greenhouse, 133.2, 0.5, 'G');
  close(budget.precipitableWater, 12.5, 0.05, 'precipitable water, kg/m²');
  close(budget.latent, 59.9, 0.5, 'L·E');
  close(budget.sensible, 17.0, 0.5, 'sensible heat');
  close(budget.evaporation * 86400, 2.07, 0.02, 'evaporation, mm/day');
  close(budget.absorbed, 249.6, 0.5, 'absorbed sunlight');
  close(budget.albedo, 0.267, 0.002, 'planetary albedo');
  close(budget.air, 65.2, 0.5, 'sunlight absorbed in the air');
  close(wet.pressure[wet.convectionTop] / 100, 318, 1, 'convection reaches 318 hPa');
  close(wet.height[wet.convectionTop], 8500, 100, 'about 8.5 km');
  const cold = wet.coldPoint();
  close(cold.pressure / 100, 102, 1, 'cold point pressure');
  close(cold.temperature, 202.7, 0.3, 'cold point temperature');
  close(layers(wet.terms) + budget.groundNet, budget.absorbed - budget.olr, 1e-6, 'a moist step closes');
});

test('the same vapor held fixed with radiation alone leaves the ground at 308.5 K', () => {
  const column = settle(makeColumn({ albedo: 0.3, humidity: vapor }), false);
  close(column.Ts, 308.5, 0.3, 'Ts');
  close(column.Ts - wet.Ts, 21.9, 0.3, 'warmer than with moist convection');
  close(column.budget.precipitableWater, 12.5, 0.05, 'the same vapor');
});

test('doubling CO₂ warms the wet column by 2.36 K with fixed relative humidity and by 0.73 K with its vapor frozen', () => {
  const free = makeColumn({ albedo: 0.3, humidity: 'relative' }), frozen = makeColumn({ albedo: 0.3, humidity: vapor });
  for (const column of [free, frozen]) {
    column.set({ temperature: Float64Array.from(wet.T), groundTemperature: wet.Ts });
    column.set({ carbonDioxide: 2 });
    settle(column);
  }
  close(free.Ts - wet.Ts, 2.36, 0.05, 'fixed relative humidity');
  close(frozen.Ts - wet.Ts, 0.73, 0.05, 'vapor frozen');
  close(free.budget.greenhouse, 144.7, 0.5, 'G, fixed relative humidity');
  close(free.budget.precipitableWater, 15.5, 0.05, 'precipitable water, fixed relative humidity');
  close(frozen.budget.precipitableWater, 12.5, 0.05, 'precipitable water, frozen');
});

test('one cloud layer: the overcast and clear calls blended by the cover reproduce the model\'s own partial cover', () => {
  const mesh = buildMesh(new Grid(2)), C = mesh.nCells, core = createSigmaCore(mesh, { levels: LEVELS });
  const { exnerLayer, dSigma, g, sigmaMid } = core.diagnostics, ps = MLS.surfacePressure;
  const pi = new Float64Array(C).fill(ps), theta = new Float64Array(K * C), q = new Float64Array(K * C), qc = new Float64Array(K * C);
  core.diagnoseColumn(0, pi, theta);
  for (let k = 0; k < K; k++) theta[k * C] = MLS.temperature[k] / exnerLayer[k * C];
  core.diagnoseColumn(0, pi, theta);
  const pdf = createRadiation(mesh, core, { solarConstant: 340.5 * 1.5, ozoneProfile: MLS.ozone });
  for (const [k, water] of [[22, 0.03], [26, 0.05]]) {
    const humidity = Float64Array.from(MLS.humidity);
    humidity[k] = 0.85 * saturationHumidity(MLS.temperature[k], ps * sigmaMid[k]);
    for (let j = 0; j < K; j++) { q[j * C] = humidity[j]; qc[j * C] = 0; }
    qc[k * C] = water / (ps * dSigma[k] / g);
    pdf.column(0, ps, theta, MLS.groundTemperature, 5, undefined, 340.5, null, q, qc, 0.07, 0.07, 1, SEA_DRAG, 0, 0, 0, 0, 'radiation');
    const f = pdf.layerCover[k], column = makeColumn({ profile: MLS, humidity, cloud: { pressure: ps * sigmaMid[k], path: water / f, cover: f } });
    assert.ok(f > 0.1 && f < 0.2, `a partial cover, ${f}`);
    assert.equal(column.cloudLayer, k);
    close(column.budget.absorbed, pdf.budget.absorbedSolar, 1e-9, `absorbed sunlight, layer ${k}`);
    close(column.budget.olr, pdf.budget.outgoingLongwave, 1e-9, `OLR, layer ${k}`);
    close(column.budget.dlr, pdf.budget.downwardLongwave, 1e-9, `DLR, layer ${k}`);
    close(column.terms.infrared[k], pdf.longwave[k * C], 1e-9, `the cloud layer's infrared, layer ${k}`);
  }
});

test('cloud radiative effects scale with the cover and close with the clear sky', () => {
  const full = makeColumn({ profile: MLS, humidity: MLS.humidity, cloud: { pressure: 894e2, path: 0.1, cover: 1 } }), clear = makeColumn({ profile: MLS, humidity: MLS.humidity });
  const part = makeColumn({ profile: MLS, humidity: MLS.humidity, cloud: { pressure: 894e2, path: 0.1, cover: 0.3 } });
  close(part.budget.shortwaveCloudEffect, 0.3 * full.budget.shortwaveCloudEffect, 1e-9, 'SWCRE');
  close(part.budget.longwaveCloudEffect, 0.3 * full.budget.longwaveCloudEffect, 1e-9, 'LWCRE');
  close(full.budget.clearAbsorbed, clear.budget.absorbed, 1e-9, 'clear-sky absorbed sunlight');
  close(full.budget.clearOlr, clear.budget.olr, 1e-9, 'clear-sky OLR');
  close(clear.budget.shortwaveCloudEffect, 0, 0, 'no cloud, no effect');
  assert.equal(clear.cloudLayer, -1);
  for (const column of [full, part]) {
    const { budget } = column;
    close(budget.ozone + budget.oxygen + budget.carbonDioxide + budget.vapor + budget.aerosol + budget.upwardNearInfrared + budget.cloudSolar, budget.air, 1e-9, 'the air\'s sunlight by absorber');
    close(budget.insolation, budget.reflected + budget.absorbed, 1e-9, 'sunlight is reflected or absorbed');
    column.step(STEP, { convection: true });
    close(column.budget.absorbed - column.budget.olr, layers(column.terms) + column.budget.groundNet, 1e-9, `a step with cover ${column.settings.cloud.cover} closes`);
  }
});

test('clouds in the midlatitude summer: a low deck cools by 110 W/m², a thin high cloud warms by 41, and vapor hides the low cloud\'s infrared', () => {
  const cloud = (pressure, path, humidity = MLS.humidity) => makeColumn({ profile: MLS, humidity, cloud: { pressure, path, cover: 1 } });
  const low = cloud(894e2, 0.1), high = cloud(270e2, 0.02), mid = cloud(608e2, 0.1), dryLow = cloud(894e2, 0.1, null), dryHigh = cloud(270e2, 0.02, null);
  close(low.pressure[low.cloudLayer] / 100, 894, 0.5, 'the low cloud\'s layer');
  close(high.pressure[high.cloudLayer] / 100, 270, 0.5, 'the high cloud\'s layer');
  close(low.budget.shortwaveCloudEffect, -116.1, 0.05, 'low SWCRE');
  close(low.budget.longwaveCloudEffect, 6.0, 0.05, 'low LWCRE');
  close(mid.budget.shortwaveCloudEffect, -114.7, 0.1, 'mid SWCRE');
  close(mid.budget.longwaveCloudEffect, 36.6, 0.05, 'mid LWCRE');
  close(high.budget.shortwaveCloudEffect, -38.1, 0.1, 'high SWCRE');
  close(high.budget.longwaveCloudEffect, 79.4, 0.05, 'high LWCRE');
  close(dryLow.budget.longwaveCloudEffect, 18.1, 0.05, 'low LWCRE in dry air');
  close(dryHigh.budget.longwaveCloudEffect, 129.2, 0.1, 'high LWCRE in dry air');
  close(low.heating.infrared[low.cloudLayer], -21.4, 0.05, 'the low cloud layer cools, K/day');
});

test('vapor in the midlatitude summer: the air takes 256.7 W/m² of an overhead sun and 73.2 of the global-mean sun', () => {
  const overhead = makeColumn({ profile: MLS, humidity: MLS.humidity, beam: 1362, mu: 1, albedo: 'sea' }).budget;
  close(overhead.air, 256.7, 0.05, 'the air, sun overhead');
  close(overhead.vapor, 206.4, 0.05, 'vapor, sun overhead');
  const global = makeColumn({ profile: MLS, humidity: MLS.humidity }).budget, dry = makeColumn({ profile: MLS }).budget;
  close(global.air, 73.2, 0.05, 'the air, global-mean sun');
  close(global.vapor, 57.3, 0.05, 'vapor, global-mean sun');
  close(global.ground, 230.2, 0.05, 'the ground, global-mean sun');
  close(dry.air, 15.2, 0.05, 'dry air, global-mean sun');
  close(dry.ground, 282.9, 0.05, 'the ground under dry air');
});
