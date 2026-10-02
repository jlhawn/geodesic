import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { topographyFromInt16, syntheticTopography } from '../js/geography.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { REFERENCE_PRESSURE, RAYLEIGH_BANDS, NEAR_INFRARED_RAYLEIGH, VISIBLE_FRACTION, LAND_AEROSOL, SEA_AEROSOL, SOLAR_CONSTANT } from '../js/physics/radiation.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

function random(seed) {
  let s = seed >>> 0;
  return () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
}

test('over a black surface the clear column reflects the two-stream value of its Rayleigh sub-bands and aerosol, the weighted τ/(τ + 2μ) of the visible beam less ozone and aerosol absorption, at three zenith angles over sea and over land (with the fixed ozone share)', () => {
  const grid = new Grid(4), C = createModel(grid, { ocean: false }).mesh.nCells, land = Uint8Array.from({ length: C }, (_, i) => i % 2);
  const model = createModel(grid, { ocean: false, radiation: { land, clearSkyPass: true, solarGases: 'lacisHansen' } });
  const init = initializeState(model, {});
  init.forEach((values, a) => model.state[a].set(values));
  const { radiation, core } = model, { K } = core.diagnostics, [pi, theta] = model.state;
  const dry = new Float64Array(K * C);
  for (let i = 0; i < C; i++) core.diagnoseColumn(i, pi, theta, dry, dry);
  let checked = 0, worst = 0;
  for (const i of [0, 1]) for (const mu of [1, 0.5, 0.15]) {
    const beam = SOLAR_CONSTANT * mu;
    radiation.column(i, pi[i], theta, 290, 5, undefined, beam, null, dry, null, 0, 0);
    const b = radiation.budget, aerosol = land[i] ? LAND_AEROSOL : SEA_AEROSOL;
    const ozone = 0.03 * beam, visible = VISIBLE_FRACTION * beam - ozone, nearInfrared = NEAR_INFRARED_RAYLEIGH * pi[i] / REFERENCE_PRESSURE;
    const taken = visible * (1 - Math.exp(-0.05 * aerosol * 35 / Math.sqrt(1224 * mu * mu + 1)));
    const depths = RAYLEIGH_BANDS.map(([w, tau]) => [w, tau * pi[i] / REFERENCE_PRESSURE + 0.3 * 0.95 * aerosol]);
    const depth = depths.reduce((s, [w, d]) => s + w * d, 0);
    const reflected = (visible - taken) * depths.reduce((s, [w, d]) => s + w * d / (d + 2 * mu), 0) + (1 - VISIBLE_FRACTION) * beam * nearInfrared / (nearInfrared + 2 * mu);
    const direct = (visible - taken) * depths.reduce((s, [w, d]) => s + w * Math.exp(-d / mu), 0) + (1 - VISIBLE_FRACTION) * beam * Math.exp(-nearInfrared / mu);
    const error = Math.max(Math.abs(b.reflectedSolar - reflected) / reflected, Math.abs(b.aerosolSolar - taken) / taken, Math.abs(b.surfaceDirect - direct) / direct);
    worst = Math.max(worst, error);
    assert.ok(error < 1e-12, `cell ${i} (${land[i] ? 'land' : 'sea'}), μ ${mu}: reflected ${b.reflectedSolar} against ${reflected}, aerosol ${b.aerosolSolar} against ${taken}, direct ${b.surfaceDirect} against ${direct}`);
    assert.ok(Math.abs(b.absorbedSolar + b.reflectedSolar - beam) < 1e-12 * beam && Math.abs(b.clearAbsorbedSolar - b.absorbedSolar) < 1e-12 * beam, `cell ${i}, μ ${mu}: the column closes and is its own clear sky`);
    checked++;
    if (mu === 0.5) console.log(`${land[i] ? 'land' : 'sea'} (aerosol ${aerosol}, mean scattering depth ${depth.toFixed(4)}) at μ 0.5: reflects ${(b.reflectedSolar / beam).toFixed(4)} of the beam over a black surface, aerosol absorbs ${(taken / beam).toFixed(4)}`);
  }
  console.log(`${checked} clear columns over a black surface match the two-stream to ${worst.toExponential(1)}`);
});

test('in every column of a state with real geography, cloud, deck and cumulus, the absorbed and reflected sunlight add up to the incident beam, the layers take what the atmosphere absorbs, and aerosol heats the lowest kilometres', () => {
  const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
  const model = createModel(new Grid(12), { topography, ocean: false, radiation: { clearSkyPass: true } });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  model.land.initialize();
  for (let n = 0; n < 6; n++) model.step(1800);
  const { radiation, core } = model, { K, C, geopotential, g } = core.diagnostics, [pi, theta, , surfaceT, q, qc, ice] = model.state;
  for (let i = 0; i < C; i++) core.diagnoseColumn(i, pi, theta, q, qc);
  const bottom = (K - 1) * C;
  let lit = 0, cloudy = 0, worstClosure = 0, worstLayers = 0, aerosolTotal = 0;
  for (let i = 0; i < C; i++) {
    const mu = radiation.cosZenith(i), beam = radiation.insolation(i);
    const onLand = model.geography.land[i], area = onLand ? 0 : model.seaIce.cover(i, ice[i]);
    const direct = onLand ? model.land.albedo(i) : model.seaIce.albedo(ice[i], mu, model.seaIce.snow[i], area), diffuse = onLand ? direct : model.seaIce.albedo(ice[i], null, model.seaIce.snow[i], area);
    radiation.column(i, pi[i], theta, surfaceT[i], 5, undefined, beam, q[bottom + i], q, qc, direct, diffuse, 1, undefined, onLand ? 0 : 1 - area, model.boundaryLayer.depth[i] - geopotential[bottom + i] / g, 0, model.boundaryLayer.mixingTop[i] - geopotential[bottom + i] / g);
    const b = radiation.budget;
    if (!(beam > 0)) continue;
    lit++;
    if (b.cloudCover > 0 || b.stratus > 0) cloudy++;
    worstClosure = Math.max(worstClosure, Math.abs(b.absorbedSolar + b.reflectedSolar - beam) / beam, Math.abs(b.clearAbsorbedSolar + (beam - b.clearAbsorbedSolar) - beam) / beam);
    assert.ok(b.reflectedSolar >= 0 && b.aerosolSolar >= 0 && b.atmosphereSolar >= b.aerosolSolar && b.clearAbsorbedSolar <= beam, `cell ${i}: every term is a share of the beam`);
    let layers = 0;
    for (let k = 0; k < K; k++) layers += radiation.layerFlux[k] - radiation.longwave[k * C + i] - (k === K - 1 ? b.sensibleHeat : 0);
    worstLayers = Math.max(worstLayers, Math.abs(layers - b.atmosphereSolar) / beam);
    aerosolTotal += b.aerosolSolar;
  }
  console.log(`${lit} sunlit columns of ${C} (${cloudy} with cloud): absorbed + reflected equals the beam to ${worstClosure.toExponential(1)} of it, the layers' shortwave heating equals the atmosphere's absorption to ${worstLayers.toExponential(1)}; aerosol absorbs ${(aerosolTotal / lit).toFixed(2)} W/m² on the sunlit mean`);
  assert.ok(lit > 0.4 * C && cloudy > 0.2 * lit, `${lit} sunlit, ${cloudy} cloudy`);
  assert.ok(worstClosure < 1e-12, `closure ${worstClosure}`);
  assert.ok(worstLayers < 1e-12, `layers ${worstLayers}`);
  const { levels } = core.diagnostics, exponent = 7e3 / 2000;
  let below = 0;
  for (let k = 0; k < K; k++) if (0.5 * (levels[k] + levels[k + 1]) > Math.exp(-2000 / 7e3)) below += levels[k + 1] ** exponent - levels[k] ** exponent;
  console.log(`the aerosol's mass profile puts ${below.toFixed(3)} of its absorption in the layers whose midpoints lie below σ = exp(−2 km / 7 km)`);
  assert.ok(below > 0.55 && aerosolTotal > 0, `${below} of the aerosol in the lowest 2 km`);
});

async function enginePair(radiation, seed) {
  const topography = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 : -4000));
  const cpu = createModel(new Grid(6), { topography, ocean: false, radiation }), gpu = await createGpuModel(new Grid(6), { topography, ocean: false, radiation });
  const init = initializeState(cpu, {});
  const { K, sigmaMid } = cpu.core, C = cpu.mesh.nCells, rnd = random(seed);
  const [pi, theta] = init, q = init[4], qc = init[5];
  for (let i = 0; i < C; i++) {
    const cloudy = rnd() < 0.7, top = Math.floor(K * (0.3 + 0.6 * rnd()));
    for (let k = 0; k < K; k++) {
      const x = k * C + i, T = theta[x] * (pi[i] * sigmaMid[k] / 1e5) ** 0.2857;
      if (sigmaMid[k] > 0.3) q[x] = (0.2 + 0.7 * rnd()) * saturationHumidity(T, pi[i] * sigmaMid[k]);
      qc[x] = cloudy && k >= top && k < top + 3 && sigmaMid[k] > 0.3 ? 3e-4 * rnd() : 0;
    }
  }
  for (let i = 0; i < C; i++) if (cpu.geography.land[i]) init[6][i] = 0;
  for (let a = 0; a < init.length; a++) { cpu.state[a].set(init[a]); gpu.state[a].set(init[a]); }
  for (const m of [cpu, gpu]) { m.time = 0.37 * 86400; m.land.initialize(); }
  gpu.load();
  await cpu.step(900); await gpu.step(900);
  const device = await gpu.gpu.downloadPhysics();
  const pick = (slot) => Float64Array.from(device[slot].subarray(0, C));
  return {
    C, land: cpu.geography.land,
    cpu: { absorbed: Float64Array.from(cpu.radiation.summed.absorbedSolar), atmosphere: Float64Array.from(cpu.radiation.summed.atmosphereSolar), reflected: Float64Array.from(cpu.radiation.summed.reflectedSolar), clear: Float64Array.from(cpu.radiation.summed.clearAbsorbedSolar), down: Float64Array.from(cpu.radiation.surfaceShortwave), insolation: Float64Array.from(cpu.radiation.summed.insolation) },
    gpu: { absorbed: pick('ABSSUM'), atmosphere: pick('ATMSUM'), reflected: pick('REFLSUM'), clear: pick('ABSCLRSUM'), down: pick('SWDN') },
    done: () => gpu.destroy(),
  };
}

test('both engines give a random set of sunlit columns over sea and land, clear and cloudy, the same absorbed, atmospheric, reflected, clear-sky and surface sunlight with the Rayleigh and aerosol scattering, and the same change from turning it off', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const on = await enginePair({ clearSkyPass: true }, 11), off = await enginePair({ clearSkyPass: true, rayleighDepth: 0, nearInfraredRayleigh: 0, landAerosol: 0, seaAerosol: 0, upwardAbsorption: false }, 11);
  const { C, land } = on;
  let lit = 0, landLit = 0, worst = 0, worstEffect = 0, largest = 0, name = '';
  for (let i = 0; i < C; i++) {
    if (!(on.cpu.insolation[i] > 0)) continue;
    lit++;
    if (land[i]) landLit++;
    for (const key of ['absorbed', 'atmosphere', 'reflected', 'clear', 'down']) {
      const scale = on.cpu.insolation[i];
      const d = Math.abs(on.cpu[key][i] - on.gpu[key][i]) / scale;
      if (d > worst) { worst = d; name = `${key} at cell ${i}`; }
      const effect = on.cpu[key][i] - off.cpu[key][i], gpuEffect = on.gpu[key][i] - off.gpu[key][i];
      largest = Math.max(largest, Math.abs(effect) / scale);
      worstEffect = Math.max(worstEffect, Math.abs(effect - gpuEffect) / scale);
    }
  }
  console.log(`one step at N=6 from a random humid, cloudy state: ${lit} sunlit columns (${landLit} over land); the engines' fluxes differ by at most ${worst.toExponential(1)} of the beam (${name}), the scattering moves them by up to ${largest.toFixed(3)} of it and the engines' change differs by ${worstEffect.toExponential(1)}`);
  assert.ok(lit > 100 && landLit > 20, `${lit} sunlit, ${landLit} over land`);
  assert.ok(worst < 2e-5, `engines apart by ${worst} of the beam at ${name}`);
  assert.ok(largest > 0.02 && worstEffect < 2e-5, `effect ${largest}, engines' effect apart by ${worstEffect}`);
  on.done(); off.done();
});

test('the light the surface reflects loses to vapour what the path it crossed coming down plus 5/3 of the column\'s absorbs beyond the first, and in the visible what the aerosol absorbs over 5/3 of its depth (Lacis-Hansen vapour)', async () => {
  const { waterVaporAbsorptivity } = await import('../js/physics/radiation.module.js');
  const grid = new Grid(4), C = createModel(grid, { ocean: false }).mesh.nCells, land = Uint8Array.from({ length: C }, (_, i) => i % 2);
  const make = (options) => {
    const model = createModel(grid, { ocean: false, radiation: { land, clearSkyPass: true, ...options } });
    initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
    return model;
  };
  const plain = { rayleighDepth: 0, nearInfraredRayleigh: 0, landAerosol: 0, seaAerosol: 0, solarGases: 'lacisHansen' };
  const on = make(plain), off = make({ ...plain, upwardAbsorption: false }), hazy = make({ rayleighDepth: 0, nearInfraredRayleigh: 0, aerosolAsymmetry: 1, solarGases: 'lacisHansen' }), hazyOff = make({ rayleighDepth: 0, nearInfraredRayleigh: 0, aerosolAsymmetry: 1, upwardAbsorption: false, solarGases: 'lacisHansen' });
  const { core } = on, { K, dSigma, sigmaMid, g } = core.diagnostics, [pi, theta] = on.state, bottom = (K - 1) * C;
  const q = new Float64Array(K * C);
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) if (sigmaMid[k] > 0.4) q[k * C + i] = 0.6 * saturationHumidity(theta[k * C + i] * (pi[i] * sigmaMid[k] / 1e5) ** 0.2857, pi[i] * sigmaMid[k]);
  for (const m of [on, off, hazy, hazyOff]) for (let i = 0; i < C; i++) m.core.diagnoseColumn(i, pi, theta, q, null);
  let worst = 0, worstAerosol = 0, largest = 0;
  for (const i of [0, 1]) for (const mu of [1, 0.4, 0.1]) for (const albedo of [0.3, 0.8]) {
    const beam = SOLAR_CONSTANT * mu, run = (m) => { m.radiation.column(i, pi[i], theta, 290, 5, undefined, beam, q[bottom + i], q, null, albedo, albedo); return { ...m.radiation.budget }; };
    const a = run(on), b = run(off);
    let path = 0;
    for (let k = 0; k < K; k++) path += q[k * C + i] * pi[i] * dSigma[k] / g * Math.sqrt(sigmaMid[k]) * 0.1;
    const down = path * 35 / Math.sqrt(1224 * mu * mu + 1);
    const expected = albedo * 0.97 * beam * (waterVaporAbsorptivity(down + 5 / 3 * path) - waterVaporAbsorptivity(down));
    const error = Math.abs(b.reflectedSolar - a.reflectedSolar - expected) / expected;
    worst = Math.max(worst, error);
    largest = Math.max(largest, expected / beam);
    assert.ok(error < 1e-12, `cell ${i}, μ ${mu}, albedo ${albedo}: vapour takes ${b.reflectedSolar - a.reflectedSolar} of the reflected light, expected ${expected}`);
    assert.ok(Math.abs(a.absorbedSolar + a.reflectedSolar - beam) < 1e-12 * beam && Math.abs(a.clearAbsorbedSolar - a.absorbedSolar) < 1e-12 * beam && Math.abs(a.atmosphereSolar - b.atmosphereSolar - expected) < 1e-12 * beam, `cell ${i}, μ ${mu}: closes, is its own clear sky, and heats the air by what it takes`);
    const h = run(hazy), hb = run(hazyOff), aerosol = land[i] ? LAND_AEROSOL : SEA_AEROSOL;
    const reflectedVisible = (VISIBLE_FRACTION * beam - 0.03 * beam - hb.aerosolSolar) * albedo;
    const aerosolUp = h.aerosolSolar - hb.aerosolSolar;
    worstAerosol = Math.max(worstAerosol, Math.abs(aerosolUp / reflectedVisible - -Math.expm1(-0.05 * aerosol * 5 / 3)) / -Math.expm1(-0.05 * aerosol * 5 / 3));
  }
  console.log(`12 humid clear columns: the vapour takes up to ${largest.toFixed(4)} of the beam from the reflected light, as the up path's absorptivity gives to ${worst.toExponential(1)}; an absorbing aerosol takes 1 - exp(-(1 - ω) τ 5/3) of the visible light the surface reflects to ${worstAerosol.toExponential(1)}`);
  assert.ok(largest > 0.005, `the up path takes ${largest} of the beam`);
  assert.ok(worstAerosol < 1e-12, `aerosol share off by ${worstAerosol}`);
});

test('both engines take the same light from the surface\'s reflection on its way up, and give it the same change', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const on = await enginePair({ clearSkyPass: true }, 11), off = await enginePair({ clearSkyPass: true, upwardAbsorption: false }, 11);
  const { C } = on;
  let lit = 0, worst = 0, worstEffect = 0, largest = 0, name = '';
  for (let i = 0; i < C; i++) {
    if (!(on.cpu.insolation[i] > 0)) continue;
    lit++;
    for (const key of ['absorbed', 'atmosphere', 'reflected', 'clear', 'down']) {
      const scale = on.cpu.insolation[i];
      const d = Math.abs(on.cpu[key][i] - on.gpu[key][i]) / scale;
      if (d > worst) { worst = d; name = `${key} at cell ${i}`; }
      const effect = on.cpu[key][i] - off.cpu[key][i], gpuEffect = on.gpu[key][i] - off.gpu[key][i];
      largest = Math.max(largest, Math.abs(effect) / scale);
      worstEffect = Math.max(worstEffect, Math.abs(effect - gpuEffect) / scale);
    }
  }
  console.log(`one step at N=6: ${lit} sunlit columns; engines apart by ${worst.toExponential(1)} of the beam (${name}); the upward absorption moves the fluxes by up to ${largest.toFixed(4)} of it and the engines' change differs by ${worstEffect.toExponential(1)}`);
  assert.ok(worst < 2e-5, `engines apart by ${worst} at ${name}`);
  assert.ok(largest > 0.002 && worstEffect < 2e-5, `effect ${largest}, engines' effect apart by ${worstEffect}`);
  on.done(); off.done();
});
