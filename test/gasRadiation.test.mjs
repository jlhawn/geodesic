import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { clearLongwave, longwaveTableFor, LONGWAVE_TABLE } from '../js/physics/longwave.module.js';
import { sigmaInterfaces, GRAVITY, CP_DRY } from '../js/dynamics/sigmaCore.module.js';
import { BENCHMARK, modelColumn, referenceAt, interfaceAt, layerHeating, referenceHeating, MOLAR } from '../scripts/standardAtmospheres.mjs';
import { PREVIOUS_SHORTWAVE, PREVIOUS_AEROSOL } from './helpers/previousDefaults.mjs';
import { runColumn } from '../scripts/radiationBenchmark.mjs';
import { SOLAR_CONSTANT, GREENHOUSE_GASES, VISIBLE_FRACTION } from '../js/physics/radiation.module.js';
import { ozoneWeights, ozoneAbove } from '../js/physics/ozone.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { ozoneAbsorptivity, visibleVaporAbsorptivity, nearInfraredVaporAbsorptivity, oxygenAbsorptivity, carbonDioxideAbsorptivity, pressureScaling, vaporScaling, OXYGEN, STP_DEPTH } from '../js/physics/shortwaveGases.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const SPECTRAL = { longwaveScheme: 'correlated', solarGases: 'clirad' };
const DT = 900;

test('the longwave g-points follow RRTMG over the standard atmospheres and LBLRTM for doubled CO2', () => {
  const lines = [], upper = [];
  for (const a of ['TROP', 'MLS', 'MLW', 'SAW']) {
    const column = modelColumn(BENCHMARK.atmospheres[a]), ref = BENCHMARK.rrtmgLongwave[a].levels, r = clearLongwave(column), K = column.T.length;
    const olr = r.up[0] - ref[ref.length - 1].up, dlr = r.down[K] - ref[0].down, net = interfaceAt(column, r.net, 20000) - referenceAt(ref, 20000);
    lines.push(`${a} ${olr.toFixed(2)}/${dlr.toFixed(2)}/${net.toFixed(2)}`);
    for (const x of [olr, dlr, net]) assert.ok(Math.abs(x) < 3, `${a}: OLR, DLR and 200 hPa misses ${olr}, ${dlr}, ${net} W/m2`);
    const heat = layerHeating(column, r.net), refHeat = referenceHeating(column, ref);
    upper.push(`${a} ${heat[0].toFixed(2)} (${refHeat[0].toFixed(2)})`);
    assert.ok(Math.abs(heat[0] / refHeat[0] - 1) < 0.2, `${a}: the top layer cools ${heat[0]} against ${refHeat[0]} K/day`);
    for (let k = 1; 0.5 * (column.levels[k] + column.levels[k + 1]) * column.ps < 3000; k++) assert.ok(Math.abs(heat[k] - refHeat[k]) < 0.3, `${a}: layer ${k} cools ${heat[k]} against ${refHeat[k]} K/day`);
  }
  const mls = modelColumn(BENCHMARK.atmospheres.MLS), K = mls.T.length;
  const withCo2 = (vmr) => ({ ...mls, co2: mls.co2.map((_, k) => vmr * MOLAR.co2 / MOLAR.air * (1 - mls.q[k])), ch4: mls.ch4.map((_, k) => 806e-9 * MOLAR.ch4 / MOLAR.air * (1 - mls.q[k])), n2o: mls.n2o.map((_, k) => 275e-9 * MOLAR.n2o / MOLAR.air * (1 - mls.q[k])) });
  const one = clearLongwave(withCo2(287e-6)), two = clearLongwave(withCo2(574e-6)), ref = BENCHMARK.iacono.co2Doubling.longwave;
  const forcing = { toa: one.net[0] - two.net[0], p20000: interfaceAt(mls, one.net, 20000) - interfaceAt(mls, two.net, 20000), surface: two.down[K] - one.down[K] };
  console.log(`top layer's cooling (RRTMG): ${upper.join(', ')} K/day; OLR/DLR/net-200-hPa misses against RRTMG: ${lines.join(', ')} W/m2; doubled CO2 TOA ${forcing.toa.toFixed(2)} 200 hPa ${forcing.p20000.toFixed(2)} surface ${forcing.surface.toFixed(2)} W/m2`);
  for (const level of ['toa', 'p20000', 'surface']) assert.ok(Math.abs(forcing[level] / ref[level] - 1) < 0.1, `doubled CO2 at ${level}: ${forcing[level]} against ${ref[level]}`);
});

test('on bl36 the g-points fitted on bl36 columns follow RRTMG in every layer above 3 hPa the reference covers, the 3-30 hPa layers, the fluxes and doubled CO2', () => {
  const levels = sigmaInterfaces('bl36'), table = longwaveTableFor(levels), lines = [];
  assert.notEqual(table, LONGWAVE_TABLE);
  assert.equal(longwaveTableFor(sigmaInterfaces('bl34')), LONGWAVE_TABLE);
  assert.equal(longwaveTableFor(sigmaInterfaces('cam26')), LONGWAVE_TABLE);
  for (const a of ['TROP', 'MLS', 'MLW', 'SAW']) {
    const column = modelColumn(BENCHMARK.atmospheres[a], levels), ref = BENCHMARK.rrtmgLongwave[a].levels, r = clearLongwave(column, { table }), K = column.T.length;
    const olr = r.up[0] - ref[ref.length - 1].up, dlr = r.down[K] - ref[0].down, net = interfaceAt(column, r.net, 20000) - referenceAt(ref, 20000);
    for (const x of [olr, dlr, net]) assert.ok(Math.abs(x) < 3, `${a}: OLR, DLR and 200 hPa misses ${olr}, ${dlr}, ${net} W/m2`);
    const heat = layerHeating(column, r.net), refHeat = referenceHeating(column, ref), top = ref[ref.length - 1].p, above = [];
    for (let k = 0; 0.5 * (levels[k] + levels[k + 1]) * column.ps < 3000; k++) {
      const p = 0.5 * (levels[k] + levels[k + 1]) * column.ps, covered = (levels[k + 1] * column.ps - Math.max(levels[k] * column.ps, top)) / ((levels[k + 1] - levels[k]) * column.ps);
      if (p < 300) { above.push(`${heat[k].toFixed(2)} (${refHeat[k].toFixed(2)})`); if (covered >= 0.9) assert.ok(Math.abs(heat[k] / refHeat[k] - 1) < 0.15, `${a}: layer ${k} cools ${heat[k]} against ${refHeat[k]} K/day`); }
      else assert.ok(Math.abs(heat[k] - refHeat[k]) < 0.25, `${a}: layer ${k} cools ${heat[k]} against ${refHeat[k]} K/day`);
    }
    lines.push(`${a} ${above.join(' / ')}`);
  }
  const mls = modelColumn(BENCHMARK.atmospheres.MLS, levels), K = mls.T.length;
  const withCo2 = (vmr) => ({ ...mls, co2: mls.co2.map((_, k) => vmr * MOLAR.co2 / MOLAR.air * (1 - mls.q[k])), ch4: mls.ch4.map((_, k) => 806e-9 * MOLAR.ch4 / MOLAR.air * (1 - mls.q[k])), n2o: mls.n2o.map((_, k) => 275e-9 * MOLAR.n2o / MOLAR.air * (1 - mls.q[k])) });
  const one = clearLongwave(withCo2(287e-6), { table }), two = clearLongwave(withCo2(574e-6), { table }), ref = BENCHMARK.iacono.co2Doubling.longwave;
  const forcing = { toa: one.net[0] - two.net[0], p20000: interfaceAt(mls, one.net, 20000) - interfaceAt(mls, two.net, 20000), surface: two.down[K] - one.down[K] };
  console.log(`bl36, the cooling of the layers above 3 hPa (RRTMG): ${lines.join(', ')} K/day; doubled CO2 TOA ${forcing.toa.toFixed(2)} 200 hPa ${forcing.p20000.toFixed(2)} surface ${forcing.surface.toFixed(2)} W/m2`);
  for (const level of ['toa', 'p20000', 'surface']) assert.ok(Math.abs(forcing[level] / ref[level] - 1) < 0.1, `doubled CO2 at ${level}: ${forcing[level]} against ${ref[level]}`);
});

test('on bl60 the g-points fitted on bl60 columns follow RRTMG over each of bl42\'s layers between 1 and 93.6 hPa, the fluxes and doubled CO2', () => {
  const levels = sigmaInterfaces('bl60'), coarse = sigmaInterfaces('bl42'), table = longwaveTableFor(levels), lines = [];
  assert.notEqual(table, longwaveTableFor(coarse));
  for (const a of ['TROP', 'MLS', 'MLW', 'SAW']) {
    const column = modelColumn(BENCHMARK.atmospheres[a], levels), ref = BENCHMARK.rrtmgLongwave[a].levels, r = clearLongwave(column, { table }), K = column.T.length;
    const olr = r.up[0] - ref[ref.length - 1].up, dlr = r.down[K] - ref[0].down, net = interfaceAt(column, r.net, 20000) - referenceAt(ref, 20000);
    for (const x of [olr, dlr, net]) assert.ok(Math.abs(x) < 3, `${a}: OLR, DLR and 200 hPa misses ${olr}, ${dlr}, ${net} W/m2`);
    const spans = [];
    for (let k = 2; k < 11; k++) {
      const top = coarse[k] * column.ps, bottom = coarse[k + 1] * column.ps, per = 86400 * GRAVITY / CP_DRY / (bottom - top);
      const heat = per * (interfaceAt(column, r.net, bottom) - interfaceAt(column, r.net, top)), refHeat = per * (referenceAt(ref, bottom) - referenceAt(ref, top));
      spans.push(`${heat.toFixed(2)} (${refHeat.toFixed(2)})`);
      if (top < 1000) assert.ok(Math.abs(heat / refHeat - 1) < 0.05, `${a}: ${top} to ${bottom} Pa cools ${heat} against ${refHeat} K/day`);
      else assert.ok(Math.abs(heat - refHeat) < 0.15, `${a}: ${top} to ${bottom} Pa cools ${heat} against ${refHeat} K/day`);
    }
    lines.push(`${a} ${spans.join(' / ')}`);
  }
  const mls = modelColumn(BENCHMARK.atmospheres.MLS, levels), K = mls.T.length;
  const withCo2 = (vmr) => ({ ...mls, co2: mls.co2.map((_, k) => vmr * MOLAR.co2 / MOLAR.air * (1 - mls.q[k])), ch4: mls.ch4.map((_, k) => 806e-9 * MOLAR.ch4 / MOLAR.air * (1 - mls.q[k])), n2o: mls.n2o.map((_, k) => 275e-9 * MOLAR.n2o / MOLAR.air * (1 - mls.q[k])) });
  const one = clearLongwave(withCo2(287e-6), { table }), two = clearLongwave(withCo2(574e-6), { table }), ref = BENCHMARK.iacono.co2Doubling.longwave;
  const forcing = { toa: one.net[0] - two.net[0], p20000: interfaceAt(mls, one.net, 20000) - interfaceAt(mls, two.net, 20000), surface: two.down[K] - one.down[K] };
  console.log(`bl60, the cooling over bl42's layers from 1 to 93.6 hPa (RRTMG): ${lines.join(', ')} K/day; doubled CO2 TOA ${forcing.toa.toFixed(2)} 200 hPa ${forcing.p20000.toFixed(2)} surface ${forcing.surface.toFixed(2)} W/m2`);
  for (const level of ['toa', 'p20000', 'surface']) assert.ok(Math.abs(forcing[level] / ref[level] - 1) < 0.1, `doubled CO2 at ${level}: ${forcing[level]} against ${ref[level]}`);
});

test('on bl48 the g-points fitted on bl48 columns follow RRTMG in each of its layers between 1 and 93.6 hPa that RRTMG\'s levels resolve, the fluxes and doubled CO2', () => {
  const levels = sigmaInterfaces('bl48'), table = longwaveTableFor(levels), lines = [];
  assert.notEqual(table, longwaveTableFor(sigmaInterfaces('bl60')));
  const linearAt = (ref, p) => { const j = ref.findIndex((l, i) => i < ref.length - 1 && p <= l.p && p >= ref[i + 1].p); return j < 0 ? (p > ref[0].p ? ref[0].net : ref[ref.length - 1].net) : ref[j].net + (ref[j + 1].net - ref[j].net) * (ref[j].p - p) / (ref[j].p - ref[j + 1].p); };
  for (const a of ['TROP', 'MLS', 'MLW', 'SAW']) {
    const column = modelColumn(BENCHMARK.atmospheres[a], levels), ref = BENCHMARK.rrtmgLongwave[a].levels, r = clearLongwave(column, { table }), K = column.T.length;
    const olr = r.up[0] - ref[ref.length - 1].up, dlr = r.down[K] - ref[0].down, net = interfaceAt(column, r.net, 20000) - referenceAt(ref, 20000);
    for (const x of [olr, dlr, net]) assert.ok(Math.abs(x) < 3, `${a}: OLR, DLR and 200 hPa misses ${olr}, ${dlr}, ${net} W/m2`);
    const heat = layerHeating(column, r.net), refHeat = referenceHeating(column, ref), linear = layerHeating(column, Array.from(column.levels, (s) => linearAt(ref, s * column.ps))), spans = [];
    for (let k = 2; k < 16; k++) {
      const resolved = Math.abs(refHeat[k] - linear[k]) < 0.25, top = levels[k] * column.ps;
      spans.push(`${heat[k].toFixed(2)} (${refHeat[k].toFixed(2)})${resolved ? '' : '*'}`);
      if (resolved) assert.ok(Math.abs(heat[k] - refHeat[k]) < (top < 1000 ? 0.45 : 0.2), `${a}: ${top} to ${levels[k + 1] * column.ps} Pa cools ${heat[k]} against ${refHeat[k]} K/day`);
    }
    lines.push(`${a} ${spans.join(' / ')}`);
  }
  const mls = modelColumn(BENCHMARK.atmospheres.MLS, levels), K = mls.T.length;
  const withCo2 = (vmr) => ({ ...mls, co2: mls.co2.map((_, k) => vmr * MOLAR.co2 / MOLAR.air * (1 - mls.q[k])), ch4: mls.ch4.map((_, k) => 806e-9 * MOLAR.ch4 / MOLAR.air * (1 - mls.q[k])), n2o: mls.n2o.map((_, k) => 275e-9 * MOLAR.n2o / MOLAR.air * (1 - mls.q[k])) });
  const one = clearLongwave(withCo2(287e-6), { table }), two = clearLongwave(withCo2(574e-6), { table }), ref = BENCHMARK.iacono.co2Doubling.longwave;
  const forcing = { toa: one.net[0] - two.net[0], p20000: interfaceAt(mls, one.net, 20000) - interfaceAt(mls, two.net, 20000), surface: two.down[K] - one.down[K] };
  console.log(`bl48, the cooling of its layers from 1 to 93.6 hPa (RRTMG; * where its levels do not resolve the layer): ${lines.join(', ')} K/day; doubled CO2 TOA ${forcing.toa.toFixed(2)} 200 hPa ${forcing.p20000.toFixed(2)} surface ${forcing.surface.toFixed(2)} W/m2`);
  for (const level of ['toa', 'p20000', 'surface']) assert.ok(Math.abs(forcing[level] / ref[level] - 1) < 0.1, `doubled CO2 at ${level}: ${forcing[level]} against ${ref[level]}`);
});

test('the radiation column carries the g-points: a cloudless column at night has the clear-sky longwave of the same layers', () => {
  const column = modelColumn(BENCHMARK.atmospheres.MLW);
  for (const longwaveLines of ['scaled', 'voigt']) {
    const r = runColumn({ ...SPECTRAL, longwaveLines }, column), clear = clearLongwave(column, { table: longwaveTableFor(column.levels, longwaveLines) });
    assert.ok(Math.abs(r.olr - clear.up[0]) < 0.2 && Math.abs(r.dlr - clear.down[column.T.length]) < 0.2, `${longwaveLines}: OLR ${r.olr} against ${clear.up[0]}, DLR ${r.dlr} against ${clear.down[column.T.length]}`);
    assert.ok(Math.abs(r.budget.clearOutgoingLongwave - r.olr) <= 1e-12 * r.olr, 'no cloud: the clear-sky OLR is the OLR');
  }
});

test('the solar gases absorb what RRTMG gives within 2 per cent at an overhead sun, and close the column budget', () => {
  for (const a of ['TROP', 'MLS', 'MLW', 'SAW']) {
    const column = modelColumn(BENCHMARK.atmospheres[a]), ref = BENCHMARK.rrtmgShortwave[a].levels, top = ref[ref.length - 1], sfc = ref[0];
    const r = runColumn(SPECTRAL, column, { beam: top.down, albedo: 0.2, solarConstant: top.down });
    const atmosphere = top.net - sfc.net, miss = r.budget.atmosphereSolar / atmosphere - 1;
    assert.ok(Math.abs(miss) < 0.02, `${a}: atmosphere absorbs ${r.budget.atmosphereSolar} against ${atmosphere}`);
    const layers = r.sw.reduce((s, x) => s + x, 0);
    assert.ok(Math.abs(layers - r.budget.atmosphereSolar) < 1e-9 * top.down, `${a}: the layers take the atmosphere's absorption, ${layers} against ${r.budget.atmosphereSolar}`);
    assert.ok(r.budget.oxygenSolar > 0 && r.budget.carbonDioxideSolar > 0, `${a}: O2 and CO2 absorb`);
  }
});

test('with the near infrared in bands the solar gases absorb what RRTMG gives within 1 percent at an overhead sun and heat each layer between 5 and 150 hPa within 9 percent of it, where CLIRAD\'s near infrared leaves a layer over 25 percent short; the layers take the column\'s absorption and every gas absorbs', () => {
  for (const a of ['TROP', 'MLS', 'MLW', 'SAW']) {
    const column = modelColumn(BENCHMARK.atmospheres[a]), ref = BENCHMARK.rrtmgShortwave[a].levels, top = ref[ref.length - 1], sfc = ref[0], K = column.T.length;
    const refHeat = referenceHeating(column, ref).map((x) => -x), atmosphere = top.net - sfc.net;
    const [clirad, bands] = [{ ...SPECTRAL, nearInfrared: 'clirad' }, { ...SPECTRAL, nearInfrared: 'bands' }].map((o) => {
      const r = runColumn(o, column, { beam: top.down, albedo: 0.2, solarConstant: top.down }), net = [0];
      for (let k = 0; k < K; k++) net.push(net[k] + r.sw[k]);
      const heat = layerHeating(column, net), misses = [];
      for (let k = 0; k < K; k++) { const p = 0.5 * (column.levels[k] + column.levels[k + 1]) * column.ps; if (p > 500 && p < 15000) misses.push(heat[k] / refHeat[k] - 1); }
      return { r, low: Math.min(...misses), high: Math.max(...misses) };
    });
    const miss = bands.r.budget.atmosphereSolar / atmosphere - 1;
    assert.ok(Math.abs(miss) < 0.01, `${a}: the atmosphere absorbs ${bands.r.budget.atmosphereSolar} against ${atmosphere}`);
    assert.ok(bands.low > -0.09 && bands.high < 0.09, `${a}: the layers between 5 and 150 hPa heat ${bands.low} to ${bands.high} off RRTMG`);
    assert.ok(clirad.low < -0.25, `${a}: CLIRAD's near infrared leaves the layers between 5 and 150 hPa at worst ${clirad.low} off RRTMG`);
    const layers = bands.r.sw.reduce((s, x) => s + x, 0);
    assert.ok(Math.abs(layers - bands.r.budget.atmosphereSolar) < 1e-9 * top.down, `${a}: the layers take the atmosphere's absorption, ${layers} against ${bands.r.budget.atmosphereSolar}`);
    assert.ok(bands.r.budget.oxygenSolar > 0 && bands.r.budget.carbonDioxideSolar > 0 && bands.r.budget.methaneSolar > 0 && bands.r.budget.vaporSolar > 0, `${a}: water vapor, O2, CO2 and CH4 absorb`);
  }
});

test('without scattering the solar gases take the CLIRAD absorptivities of their paths down, and of the path down plus 5/3 of the column from the light the surface reflects', () => {
  const grid = new Grid(16), plain = { ...PREVIOUS_SHORTWAVE, ...PREVIOUS_AEROSOL, rayleighDepth: 0, nearInfraredRayleigh: 0, landAerosol: 0, seaAerosol: 0, clearSkyPass: true, nearInfrared: 'clirad' };
  const make = (options) => { const m = createModel(grid, { ocean: false, radiation: { ...plain, ...options } }); initializeState(m, {}).forEach((v, a) => m.state[a].set(v)); return m; };
  const on = make({}), off = make({ upwardAbsorption: false });
  const { core, mesh } = on, { K, C, dSigma, sigmaMid, g, exnerLayer } = core.diagnostics, [pi, theta] = on.state, bottom = (K - 1) * C;
  const q = new Float64Array(K * C);
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) if (sigmaMid[k] > 0.4) q[k * C + i] = 0.6 * saturationHumidity(theta[k * C + i] * (pi[i] * sigmaMid[k] / 1e5) ** 0.2857, pi[i] * sigmaMid[k]);
  for (const m of [on, off]) for (let i = 0; i < C; i++) m.core.diagnoseColumn(i, pi, theta, q, null);
  let worst = 0;
  for (const i of [0, 7]) for (const mu of [1, 0.4, 0.1]) for (const albedo of [0, 0.3]) {
    const beam = SOLAR_CONSTANT * mu, m = 35 / Math.sqrt(1224 * mu * mu + 1);
    let vapour = 0, oxygen = 0, co2 = 0;
    for (let k = 0; k < K; k++) {
      const idx = k * C + i, p = pi[i] * sigmaMid[k], mass = pi[i] * dSigma[k] / g, dry = mass * (1 - q[idx]);
      vapour += q[idx] * mass * 0.1 * vaporScaling(p, theta[idx] * exnerLayer[idx]); oxygen += OXYGEN.mixingRatio * STP_DEPTH * dry * pressureScaling(p); co2 += GREENHOUSE_GASES.carbonDioxide * STP_DEPTH * dry * pressureScaling(p);
    }
    const ozone = ozoneAbove(pi[i], ozoneWeights(mesh.latCell[i], 0));
    const nir = (f) => nearInfraredVaporAbsorptivity(f * vapour) + oxygenAbsorptivity(f * oxygen) + carbonDioxideAbsorptivity(f * co2);
    const down = beam * (ozoneAbsorptivity(m * ozone) + visibleVaporAbsorptivity(m * vapour) + nir(m));
    const visible = VISIBLE_FRACTION * beam - beam * (ozoneAbsorptivity(m * ozone) + visibleVaporAbsorptivity(m * vapour));
    const up = albedo * beam * (nir(m + 5 / 3) - nir(m)) + albedo * visible * -Math.expm1(-0.054209 * ozone * 5 / 3);
    const run = (model) => { model.radiation.column(i, pi[i], theta, 290, 5, undefined, { beam, mu }, q[bottom + i], q, null, albedo, albedo); return { ...model.radiation.budget }; };
    const a = run(on), b = run(off);
    const error = Math.max(Math.abs(b.atmosphereSolar - down) / down, albedo > 0 ? Math.abs(a.atmosphereSolar - b.atmosphereSolar - up) / up : 0);
    worst = Math.max(worst, error);
    assert.ok(error < 1e-4, `cell ${i}, μ ${mu}, albedo ${albedo}: down ${b.atmosphereSolar} against ${down}, up ${a.atmosphereSolar - b.atmosphereSolar} against ${up}`);
    assert.ok(Math.abs(a.absorbedSolar + a.reflectedSolar - beam) < 1e-12 * beam && Math.abs(a.clearAbsorbedSolar - a.absorbedSolar) < 1e-12 * beam, `cell ${i}, μ ${mu}: closes and is its own clear sky`);
  }
  console.log(`12 humid clear columns without scattering: the gases take their absorptivities down and up to ${worst.toExponential(1)}`);
});

async function engines(radiation, levels = undefined) {
  const grid = levels ? { levels } : {};
  const cpu = createModel(new Grid(16), { ocean: false, radiation, ...grid }), gpu = await createGpuModel(new Grid(16), { ocean: false, radiation, ...grid });
  const init = initializeState(cpu, {});
  for (let a = 0; a < init.length; a++) { cpu.state[a].set(init[a]); gpu.state[a].set(init[a]); }
  gpu.load();
  return { cpu, gpu, C: cpu.mesh.nCells };
}

function stats(cpu, gpu) {
  let maxDiff = 0, sumSq = 0, sumRef = 0;
  for (let x = 0; x < cpu.length; x++) { const d = Math.abs(cpu[x] - gpu[x]); maxDiff = Math.max(maxDiff, d); sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
  return { maxDiff, rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)) };
}

for (const grid of [null, 'bl36']) test(`both engines carry the spectral gases alike${grid ? ` on ${grid}, each with its table` : ''}: per-cell sums of absorbed and outgoing radiation, clear and all-sky, at N=16`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpu, gpu, C } = await engines({ ...SPECTRAL, clearSkyPass: true }, grid ? sigmaInterfaces(grid) : undefined);
  const slots = { absorbedSolar: 'ABSSUM', atmosphereSolar: 'ATMSUM', outgoingLongwave: 'OLRSUM', clearAbsorbedSolar: 'ABSCLRSUM', clearOutgoingLongwave: 'OLRCLRSUM' };
  let kept = Array.from({ length: C }, (_, x) => x);
  const compare = async (names, limit) => {
    const device = await gpu.gpu.downloadPhysics();
    const apart = (name, x) => Math.abs(cpu.radiation.summed[name][x] - device[slots[name]][x]) > 1e-3 * Math.abs(cpu.radiation.summed[name][x]);
    const flipped = names.every((name) => name.startsWith('clear')) ? [] : Array.from({ length: C }, (_, x) => x).filter((x) => (apart('outgoingLongwave', x) || apart('absorbedSolar', x)) && !apart('clearOutgoingLongwave', x) && !apart('clearAbsorbedSolar', x));
    assert.ok(flipped.length <= 2, `${flipped.length} cells whose all-sky sums parted with their clear-sky sums alike: ${flipped}`);
    if (flipped.length) kept = kept.filter((x) => !flipped.includes(x));
    return names.map((name) => {
      const s = stats(kept.map((x) => cpu.radiation.summed[name][x]), kept.map((x) => device[slots[name]][x]));
      assert.ok(s.rmsRel < limit, `${name}: per-cell rms ${s.rmsRel}`);
      return `${name} rms ${s.rmsRel.toExponential(1)} (max ${s.maxDiff.toExponential(1)})`;
    });
  };
  for (let s = 0; s < 4; s++) { cpu.step(DT); await gpu.step(DT); }
  const early = await compare(Object.keys(slots), 1e-5);
  const K = cpu.core.diagnostics.K, device = await gpu.gpu.downloadPhysics(), top = [];
  for (let k = 0; 0.5 * (cpu.core.levels[k] + cpu.core.levels[k + 1]) < 0.03; k++) {
    const s = stats(kept.map((x) => cpu.radiation.longwave[k * C + x]), kept.map((x) => device.LWH[k * C + x]));
    assert.ok(s.rmsRel < 1e-4, `layer ${k}'s longwave heating: per-cell rms ${s.rmsRel}`);
    top.push(s.rmsRel.toExponential(1));
  }
  early.push(`the ${top.length} layers above 30 hPa' longwave heating rms ${top.join(', ')}`);
  for (let s = 4; s < 24; s++) { cpu.step(DT); await gpu.step(DT); }
  const later = await compare(['clearAbsorbedSolar', 'clearOutgoingLongwave'], 1e-4);
  console.log(`spectral gases${grid ? ` on ${grid}` : ''}, CPU against GPU per-cell sums at N=16, 4 steps: ${early.join(', ')}; 24 steps: ${later.join(', ')}`);
  gpu.destroy();
});
