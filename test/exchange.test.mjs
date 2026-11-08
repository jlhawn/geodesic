import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { SEA_DRAG, LAND_DRAG } from '../js/physics/surface.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { VIRTUAL_FACTOR, sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { transfer, psiMomentum, psiHeat, andreasScalar, seaIceRoughness, charnockRoughness, createBlend, referenceCoefficient, exchangeMode, airViscosity, KARMAN, heatProfile, profileTemperature } from '../js/physics/exchange.module.js';
import { PREVIOUS_LEVELS, PREVIOUS_SURFACE } from './helpers/previousDefaults.mjs';
import { cpuFrameFields } from '../js/cpuFrame.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { readRanges } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};

const close = (actual, expected, tolerance, what) => assert.ok(Math.abs(actual - expected) <= tolerance * Math.abs(expected), `${what}: ${actual} against ${expected}`);

test('neutral coefficients follow the logarithmic law at the level\'s own height', () => {
  const c = transfer(0, 20, 0.1, 1e-3);
  close(c.drag, KARMAN ** 2 / Math.log(20.1 / 0.1) ** 2, 1e-12, 'C_D');
  close(c.drag, 5.68888e-3, 1e-5, 'C_D by hand');
  close(c.heat, 3.04485e-3, 1e-5, 'C_H by hand');
  close(transfer(0, 20, 2, 2).drag, 0.16 / Math.log(11) ** 2, 1e-12, 'forest C_D');
});

test('the stability functions are the IFS\'s Paulson and Holtslag–De Bruin forms', () => {
  close(psiMomentum(-1), 1.1162322497683264, 1e-12, 'Ψm(−1)');
  close(psiHeat(-1), 1.8812272842144175, 1e-12, 'Ψh(−1)');
  close(psiMomentum(-0.1), 0.28361371121278045, 1e-12, 'Ψm(−0.1)');
  close(psiHeat(-0.1), 0.5342837819484251, 1e-12, 'Ψh(−0.1)');
  close(psiMomentum(0.5), -2.308799761502047, 1e-12, 'Ψm(0.5)');
  close(psiHeat(0.5), -2.3484004793410485, 1e-12, 'Ψh(0.5)');
  close(psiMomentum(2), -7.456539416565598, 1e-12, 'Ψm(2)');
  close(psiHeat(2), -8.020764957086806, 1e-12, 'Ψh(2)');
  close((psiMomentum(1e-6) - psiMomentum(0)) / 1e-6, -5, 1e-4, 'stable slope at neutral');
  assert.ok(Math.abs(psiMomentum(0)) < 1e-12 && Math.abs(psiHeat(0)) < 1e-12);
});

test('the 2 m temperature\'s heat profile is the exchange\'s: at the lowest layer\'s height it is κ√C_D/C_H of transfer() at the same Richardson number and roughness, unstable to z/L past 30', () => {
  let largest = 0;
  for (const [z0m, z0h] of [[0.1, 1e-3], [2, 2], [0.013, 1.3e-4], [2e-4, 3e-5]]) for (const ri of [-2, -0.3, 0, 0.1, 0.6, 1.2, 1.8, 3]) {
    const c = transfer(ri, 21, z0m, z0h);
    largest = Math.max(largest, c.zeta);
    close(heatProfile(21, c.zeta, 21, z0m, z0h), KARMAN * Math.sqrt(c.drag) / c.heat, 1e-12, `Ri ${ri}, z0m ${z0m}, z0h ${z0h}, ζ ${c.zeta}`);
  }
  assert.ok(largest > 30, `the largest ζ is ${largest}`);
});

test('the 2 m temperature lies along that profile with the WMO\'s open-terrain roughness wherever z0m exceeds 0.03 m and the profile uniform above z/L = 5, as the IFS interpolates its own (Cy41r2 Part IV eq. 3.95, section 3.2)', () => {
  const g = 9.80616, cp = 1004.64, z = 21, surfaceT = 280, airT = 284;
  const share = (zeta, z0m, z0h, top = z) => (profileTemperature(2, surfaceT, airT, top, zeta, z0m, z0h, g, cp) + 2 * g / cp - surfaceT) / (airT + g * top / cp - surfaceT);
  for (const zeta of [-2, -0.5, 0, 1, 3, 20, 400]) {
    assert.equal(share(zeta, 2, 2), share(zeta, 0.1, 1e-3), `ζ ${zeta}: forest and grass both take the WMO roughness`);
    assert.equal(share(zeta, 0.1, 1e-3), share(zeta, 0.03, 0.003), `ζ ${zeta}: the WMO roughness is z0m 0.03 m and z0h 0.003 m`);
    assert.ok(zeta >= 5 * z / 2 || share(zeta, 0.013, 1.3e-4) !== share(zeta, 0.03, 0.003), `ζ ${zeta}: a cell under 0.03 m keeps its own roughness`);
    close(profileTemperature(z, surfaceT, airT, z, zeta, 0.1, 1e-3, g, cp), airT, 1e-12, `ζ ${zeta}: the lowest layer's air at its height`);
  }
  close(share(0, 0.013, 1.3e-4), Math.log(2.013 / 1.3e-4) / Math.log(21.013 / 1.3e-4), 1e-12, 'neutral: the cell\'s own log law');
  close(share(0, 2, 2), Math.log(2.03 / 0.003) / Math.log(21.03 / 0.003), 1e-12, 'neutral: the WMO roughness\'s log law');
  for (const zeta of [6, 20, 50]) close(share(zeta, 0.1, 1e-3), share(2 * zeta, 0.1, 1e-3, 2 * z), 1e-12, `ζ ${zeta}: uniform above h = 5 L, a lowest layer twice as high at the same L gives the same share`);
  assert.ok(Math.abs(share(3, 0.1, 1e-3) - share(6, 0.1, 1e-3, 2 * z)) > 1e-3, 'below h = 5 L the lowest layer\'s height matters');
  for (const zeta of [5 * z / 2, 400]) close(share(zeta, 0.1, 1e-3), 1, 1e-12, `ζ ${zeta}: with 5 L at or under 2 m the 2 m air is the lowest layer's`);
  for (const [zeta, weight] of [[0, 0.736], [-0.5, 0.832], [1, 0.526], [10, 0.438], [20, 0.596]]) assert.ok(Math.abs(share(zeta, 2, 2) - weight) < 5e-4, `ζ ${zeta}: ${share(zeta, 2, 2)} against the review's IFS construction, ${weight}`);
});

test('the bulk Richardson number gives back the Obukhov length it came from, and its coefficients', () => {
  for (const [zeta, ri, drag, heat] of [[-1, -0.45398989837073, 0.009052788535721315, 0.004743155717988523], [0.5, 0.1059101692383713, 0.002762622390942275, 0.001713778310570432]]) {
    const c = transfer(ri, 20, 0.1, 1e-3);
    close(c.zeta, zeta, 1e-3, `ζ from Ri ${ri}`);
    close(c.drag, drag, 1e-3, `C_D at ζ ${zeta}`);
    close(c.heat, heat, 1e-3, `C_H at ζ ${zeta}`);
  }
  let previous = Infinity;
  for (let ri = -2; ri <= 2; ri += 0.05) {
    const c = transfer(ri, 20, 0.1, 1e-3);
    assert.ok(c.drag > 0 && c.heat > 0 && c.drag < previous, `C_D falls as the layer stabilises (Ri ${ri.toFixed(2)})`);
    previous = c.drag;
  }
});

test('roughness by surface: COARE 3.5 over the sea, the IFS over sea ice, Andreas over snow and ice, and the blend of tiles', () => {
  const nu = airViscosity(20);
  close(nu, 1.326e-5 * (1 + 0.13084 + 0.0033204 - 0.00003872), 1e-12, 'viscosity at 20 °C');
  const sea = charnockRoughness(10, 10, 1.5e-5, undefined, 40);
  const ustar = KARMAN * 10 / Math.log((10 + sea.momentum) / sea.momentum), alpha = 0.0017 * (ustar / KARMAN * Math.log((10 + sea.momentum) / sea.momentum)) - 0.005;
  close(sea.momentum, alpha * ustar * ustar / 9.81 + 0.11 * 1.5e-5 / ustar, 1e-6, 'Charnock z0 solves its relation');
  close(sea.heat, Math.min(1.6e-4, 5.8e-5 / Math.pow(sea.momentum * ustar / 1.5e-5, 0.72)), 1e-6, 'COARE scalar roughness');
  close(seaIceRoughness(1), 1e-3, 1e-12, 'compact ice');
  close(seaIceRoughness(0.5), 0.93e-3 * 0.5 + 6.05e-3, 1e-12, 'half cover');
  close(seaIceRoughness(0), 0.93e-3 + 6.05e-3 * Math.exp(-4.25), 1e-12, 'open');
  close(andreasScalar(1e-3, 100 * 1.5e-5 / 1e-3, 1.5e-5), 1e-3 * 0.002099805469839084, 1e-9, 'Andreas rough');
  close(andreasScalar(1e-4, 1 * 1.5e-5 / 1e-4, 1.5e-5), 1e-4 * Math.exp(0.149), 1e-9, 'Andreas transition');
  close(andreasScalar(1e-5, 0.1 * 1.5e-5 / 1e-5, 1.5e-5), 1e-5 * Math.exp(1.25), 1e-9, 'Andreas smooth');
  const single = createBlend(10).clear().add(1, 0.1, 1e-3).finish();
  close(single.momentum, 0.1, 1e-12, 'one tile keeps its z0m');
  close(single.heat, 1e-3, 1e-12, 'one tile keeps its z0h');
  const mixed = createBlend(10).clear().add(0.5, 2, 2).add(0.5, 0.1, 1e-3).finish();
  close(mixed.momentum, 1.1664147664251205, 1e-12, 'forest and grass z0m');
  close(mixed.heat, 1.030731926828355, 1e-12, 'forest and grass z0h');
  close(referenceCoefficient(2), 0.41 ** 2 / (Math.log(1.92 / 0.01476) * Math.log(1.92 / 0.001476)), 1e-12, 'FAO-56 reference at 2 m');
  close(1 / referenceCoefficient(2), 208, 3e-3, 'FAO-56 eq. 4: r_a 208/u₂');
});

test('the exchange is physical by default and fixed where an option set names a coefficient', () => {
  assert.equal(exchangeMode({}, {}), 'roughness');
  assert.equal(exchangeMode({ dragCoefficient: 1.3e-3 }, {}), 'fixed');
  assert.equal(exchangeMode({}, { dragCoefficient: 2e-3 }), 'fixed');
  assert.equal(exchangeMode({ exchange: 'fixed' }, {}), 'fixed');
  assert.throws(() => exchangeMode({ exchange: 'roughness', dragCoefficient: 1e-3 }, {}));
  const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
  const fixed = createModel(new Grid(16), { topography, surface: { dragCoefficient: 1.1e-3 } });
  for (let i = 0; i < fixed.mesh.nCells; i++) assert.equal(fixed.exchange.drag[i], fixed.geography.land[i] ? LAND_DRAG : 1.1e-3);
  assert.equal(fixed.exchange.heat, fixed.exchange.drag);
  assert.equal(createModel(new Grid(16), {}).exchange.mode, 'roughness');
  assert.equal(createModel(new Grid(16), { surface: { exchange: 'fixed' } }).exchange, null);
  assert.ok(SEA_DRAG === 1.2e-3);
});

const topography = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 : -4000));

function prepare(model) {
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const { mesh, geography } = model, C = mesh.nCells;
  let seed = 7;
  const random = () => { seed = (seed * 16807) % 2147483647; return seed / 2147483647; };
  const concentration = new Float64Array(C);
  for (let i = 0; i < C; i++) {
    if (geography.land[i]) { model.state[6][i] = 0; continue; }
    if (Math.abs(mesh.latCell[i]) > 0.9) { model.state[6][i] = 1 + random(); concentration[i] = random() < 0.3 ? 1 : 0.1 + 0.85 * random(); model.state[3][i] = 255 + 15 * random(); }
  }
  model.seaIce.load(model.state[6], concentration);
  if (model.load) model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  const soil = new Float64Array(C), snow = new Float64Array(C), vegetation = new Float64Array(C), canopy = new Float64Array(C);
  for (let i = 0; i < C; i++) if (geography.land[i]) {
    soil[i] = 60 + 150 * random(); vegetation[i] = random(); canopy[i] = vegetation[i] * random();
    snow[i] = random() < 0.3 ? 20 + 40 * random() : 0;
  }
  model.land.load({ soil, snow, vegetation, canopy, seasonLength: new Float64Array(C).fill(0.6), seasonWarmth: new Float64Array(C).fill(8), rainMean: new Float64Array(C).fill(2), demandMean: new Float64Array(C).fill(3) }, model.state[6]);
  return model;
}

const bottomLayer = (levels) => levels[levels.length - 1] - levels[levels.length - 2];
for (const [name, levels] of [['cam26', PREVIOUS_LEVELS], ['bl48', sigmaInterfaces('bl48')]]) test(`on ${name} both engines give every surface the same transfer coefficients and fluxes, to an f32 bound that grows as the lowest layer thins`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const bound = 2e-4 * bottomLayer(PREVIOUS_LEVELS) / bottomLayer(levels);
  const cpu = prepare(createModel(new Grid(16), { topography, levels }));
  const gpu = prepare(await createGpuModel(new Grid(16), { topography, levels }));
  cpu.step(900); await gpu.step(900); await gpu.settle();
  const PH = gpu.gpu.layout.PH, C = cpu.mesh.nCells;
  const [drag, heat, reference, evaporation, surfaceFlux, t2m] = await readRanges(gpu.gpu.device, gpu.gpu.buffers.PH, ['DRAG', 'HEATX', 'REFX', 'EVAP', 'SFLUX', 'T2M'].map((name) => ({ offset: PH[name], length: C })));
  const { land, iceSheet } = cpu.geography;
  const classOf = (i) => (land[i] ? (iceSheet[i] ? 'ice sheet' : cpu.land.snow[i] > 30 ? 'snow' : cpu.land.canopy[i] > 0.5 ? 'forest' : 'other land') : cpu.state[6][i] > 0 ? 'sea ice' : 'open sea');
  const worst = {};
  for (let i = 0; i < C; i++) {
    const k = classOf(i), w = (worst[k] ??= { cells: 0, drag: 0, heat: 0, reference: 0, evaporation: 0, flux: 0, cd: 0, t2m: 0, finite: 0 });
    w.cells++;
    w.cd += cpu.exchange.drag[i];
    w.drag = Math.max(w.drag, Math.abs(drag[i] - cpu.exchange.drag[i]) / Math.max(cpu.exchange.drag[i], 1e-5));
    w.heat = Math.max(w.heat, Math.abs(heat[i] - cpu.exchange.heat[i]) / Math.max(cpu.exchange.heat[i], 1e-5));
    w.reference = Math.max(w.reference, Math.abs(reference[i] / cpu.exchange.reference[i] - 1));
    w.t2m = Math.max(w.t2m, Math.abs(t2m[i] - cpu.exchange.temperature2m[i]));
    if (Number.isFinite(t2m[i]) && Number.isFinite(cpu.exchange.temperature2m[i])) w.finite++;
    w.evaporation = Math.max(w.evaporation, (Math.abs(evaporation[i] - cpu.radiation.evaporation[i]) * 2.5e6 - 1e-5 * Math.abs(cpu.radiation.evaporation[i]) * 2.5e6));
    w.flux = Math.max(w.flux, Math.abs(surfaceFlux[i] - cpu.radiation.surfaceFlux[i]) - 1e-5 * Math.abs(cpu.radiation.surfaceFlux[i]));
  }
  for (const [k, w] of Object.entries(worst)) {
    console.log(`${k}: ${w.cells} cells, mean C_D ${(w.cd / w.cells * 1e3).toFixed(3)}e-3; engines apart: C_D ${w.drag.toExponential(1)} and C_H ${w.heat.toExponential(1)} relative to the coefficient or 1e-5 if larger, reference ${w.reference.toExponential(1)} relative; latent heat ${w.evaporation.toExponential(1)}, net surface flux ${w.flux.toExponential(1)} W/m² beyond 1e-5 of each cell's own`);
    assert.ok(w.drag < bound && w.heat < bound && w.reference < bound / 20, `${k} coefficients against ${bound.toExponential(1)}`);
    assert.ok(w.evaporation < 0.05 && w.flux < 0.1, `${k} fluxes`);
    console.log(`${k}: the 2 m temperature ${w.t2m.toExponential(1)} K apart, ${w.finite} of ${w.cells} finite on both`);
    assert.ok(w.finite === w.cells && w.t2m < 2e-3, `${k} 2 m temperature`);
  }
  for (const k of ['ice sheet', 'snow', 'forest', 'other land', 'sea ice', 'open sea']) assert.ok(worst[k] && worst[k].cells > 0, `the grid has ${k}`);
});

test('snow smooths grass and bare soil but leaves the trees, and the implicit drag applies the stress it stores', () => {
  const model = prepare(createModel(new Grid(16), { topography, orography: false, moist: { plumeMomentum: false } }));
  model.step(900);
  const { land, iceSheet } = model.geography, C = model.mesh.nCells;
  const i = [...Array(C).keys()].find((n) => land[n] && !iceSheet[n]);
  const args = (snow, cover, trees) => [i, model.state[0], model.state[1], model.state[4], model.state[5], 285, 6, 0, snow, cover, trees];
  const forest = model.exchange.cell(...args(0, 1, 1)).momentum, forestSnow = model.exchange.cell(...args(60, 1, 1)).momentum;
  const grass = model.exchange.cell(...args(0, 1, 0)).momentum, grassSnow = model.exchange.cell(...args(60, 1, 0)).momentum, halfSnow = model.exchange.cell(...args(15, 1, 0)).momentum;
  close(forestSnow, forest, 1e-12, 'trees under snow');
  close(forest, 2, 1e-9, 'forest');
  close(grass, 0.1, 1e-9, 'grass');
  close(grassSnow, 1.3e-3, 1e-9, 'buried grass');
  assert.ok(halfSnow < grass && halfSnow > grassSnow);
  const { surfaceStress, surfaceDrag } = model.boundaryLayer, E = model.mesh.nEdges, bottom = (model.core.K - 1) * E;
  for (let e = 0; e < E; e += 7) {
    const a = model.mesh.cellsOnEdge[2 * e], b = model.mesh.cellsOnEdge[2 * e + 1];
    close(surfaceStress[e], 0.5 * (surfaceDrag[a] + surfaceDrag[b]) * model.state[2][bottom + e], 1e-12, `edge ${e}'s stress`);
  }
  const stress = model.surface.stress(model.state);
  for (let e = 0; e < E; e++) assert.equal(stress[e], surfaceStress[e]);
});

test('the gustiness is the free-convection velocity of the step\'s own flux, and the land\'s surface humidity is the one its evaporation implies', () => {
  const model = prepare(createModel(new Grid(16), { topography }));
  model.step(900);
  const { land, iceSheet } = model.geography, C = model.mesh.nCells, { exchange } = model;
  const i = [...Array(C).keys()].find((n) => land[n] && !iceSheet[n]), sea = [...Array(C).keys()].find((n) => !land[n] && model.state[6][n] <= 0);
  const call = (cell, skin, lowest, wetness, depth) => exchange.cell(cell, model.state[0], model.state[1], model.state[4], model.state[5], skin, lowest, 0, 0, 1, 1, wetness, depth);
  const air = model.state[1][(model.core.K - 1) * C + i] * model.core.diagnostics.exnerLayer[(model.core.K - 1) * C + i];
  for (const [cell, beta, skin] of [[i, 1, air + 8], [sea, 1.2, model.state[3][sea] + 3]]) {
    const { z, ri } = call(cell, skin, 0.5, 1, 800);
    const speed = exchange.wind[cell], buoyancy = -ri * speed ** 3 * exchange.heat[cell] / z;
    assert.ok(buoyancy > 0, 'a heated surface');
    close(speed, Math.hypot(0.5, Math.max(0.2, beta * Math.cbrt(buoyancy * 800))), 2e-3, `U² = |v|² + max(0.2, β w*)² with β ${beta}`);
    let stable = skin - 20, unstable = skin;
    for (let n = 0; n < 60; n++) { const mid = 0.5 * (stable + unstable); if (call(cell, mid, 0.5, 1, 800).ri > 0) stable = mid; else unstable = mid; }
    call(cell, stable, 0.5, 1, 800); const below = exchange.wind[cell];
    call(cell, unstable, 0.5, 1, 800); const above = exchange.wind[cell];
    assert.ok(Math.abs(above - below) < 1e-6, `the gust is continuous where the buoyancy flux changes sign (${below} against ${above})`);
  }
  call(i, air - 5, 0.5, 1, 800);
  close(exchange.wind[i], Math.hypot(0.5, 0.2), 1e-12, 'a stable surface keeps COARE\'s 0.2 m/s');
  call(i, air + 8, 0, 1, 0);
  assert.ok(exchange.wind[i] > 0.5, `with no wind and no depth yet the surface layer's own depth sets w* (${exchange.wind[i]})`);
  const wet = call(i, air + 2, 3, 1, 800).ri, dry = call(i, air + 2, 3, 0, 800).ri;
  assert.ok(wet < dry, `a transpiring surface is the more unstable (Ri ${wet} against ${dry})`);
  const airHumidity = prepare(createModel(new Grid(16), { topography, surface: { landHumidity: 'air' } }));
  airHumidity.step(900);
  const other = (wetness) => airHumidity.exchange.cell(i, airHumidity.state[0], airHumidity.state[1], airHumidity.state[4], airHumidity.state[5], air + 2, 3, 0, 0, 1, 1, wetness, 800).ri;
  assert.equal(other(1), other(0), '\'air\' takes the lowest layer\'s humidity whatever the wetness');
});

test('a saturated, fog-filled lowest layer over a sea 3 K warmer counts its vapor only, so its Ri_b, C_H and gust stay those of moist air; with its cloud water as loading it would read far more unstable', () => {
  const model = prepare(createModel(new Grid(16), { topography, surface: PREVIOUS_SURFACE }));
  model.step(900);
  const { exchange, state, core, mesh } = model, C = mesh.nCells, K = core.K, b = (K - 1) * C;
  const { exnerLayer, exnerLower, sigmaMid, g } = core.diagnostics;
  const i = [...Array(C).keys()].find((n) => !model.geography.land[n] && !(state[6][n] > 0) && Math.abs(mesh.latCell[n]) < 0.5);
  const pi = state[0][i], airT = state[1][b + i] * exnerLayer[b + i], skin = airT + 3;
  state[4][b + i] = saturationHumidity(airT, pi * sigmaMid[K - 1]);
  state[5][b + i] = 5e-3;
  const { z, ri, momentum, heat } = { ...exchange.cell(i, state[0], state[1], state[4], state[5], skin, 0.5, 0, 0, 0, 0, 0, 800) };
  const speed = exchange.wind[i], theta = state[1][b + i];
  const surfaceV = skin / exnerLower[b + i] * (1 + VIRTUAL_FACTOR * saturationHumidity(skin, pi));
  const richardson = (airV) => g * z * (airV - surfaceV) / (0.5 * (airV + surfaceV) * speed * speed);
  const vapor = richardson(theta * (1 + VIRTUAL_FACTOR * state[4][b + i])), loaded = richardson(theta * (1 + VIRTUAL_FACTOR * state[4][b + i] - state[5][b + i]));
  close(ri, vapor, 1e-12, 'Ri_b with the vapor alone');
  const neutral = transfer(0, z, momentum, heat).heat;
  assert.ok(exchange.heat[i] <= 5 * neutral, `C_H ${exchange.heat[i]} against neutral ${neutral}`);
  assert.ok(Math.sqrt(speed * speed - 0.25) <= 2, `gust ${Math.sqrt(speed * speed - 0.25)} m/s`);
  assert.ok(ri < 0 && loaded < 1.3 * ri, `the loading would read Ri_b ${loaded} against ${ri}`);
});

const surfaceKinds = (model) => Array.from({ length: model.mesh.nCells }, (_, i) => (model.geography.land[i] ? 'land' : !(model.state[6][i] > 0) ? 'open sea' : model.seaIce.concentration[i] >= 1 ? 'full ice' : 'part ice'));

test('seaHumidity takes the open water\'s share of the skin\'s saturation humidity: from one state the first step\'s evaporation, C_H and the fixed exchange\'s buoyancy flux are those of a factor of 1 to the bit over land and full ice, the evaporation falls over open water and partial ice and the fixed exchange\'s buoyancy flux over open water, and Ri_b takes 1 − (1 − seaHumidity)(1 − A) of q_s at concentration A', () => {
  const first = (surface) => { const model = prepare(createModel(new Grid(16), { topography, surface })), kind = surfaceKinds(model); model.step(900); return { model, kind }; };
  let rough = null;
  for (const fixed of [false, true]) {
    const exchange = fixed ? { exchange: 'fixed' } : {}, { model: base, kind } = first({ ...exchange, seaHumidity: 1 }), { model: salty } = first(exchange);
    const C = base.mesh.nCells, fall = {};
    for (let i = 0; i < C; i++) {
      const k = kind[i], e0 = base.radiation.evaporation[i], e1 = salty.radiation.evaporation[i], b0 = base.boundaryLayer.buoyancyFlux[i], b1 = salty.boundaryLayer.buoyancyFlux[i];
      if (k === 'land' || k === 'full ice') assert.ok(e1 === e0 && b1 === b0 && salty.exchange.heat[i] === base.exchange.heat[i], `${k} cell ${i}: evaporation ${e1} against ${e0}, buoyancy flux ${b1} against ${b0}`);
      else {
        assert.ok(e0 > 0 ? e1 < e0 : e1 === 0, `${k} cell ${i}: evaporation ${e1} against ${e0}`);
        if (fixed && k === 'open sea') assert.ok(b1 < b0, `cell ${i}: the fixed exchange's buoyancy flux ${b1} against ${b0}`);
        const a = base.mesh.areaCell[i], f = (fall[k] ??= { area: 0, base: 0, salty: 0 });
        f.area += a; f.base += a * e0; f.salty += a * e1;
      }
    }
    console.log(`${fixed ? 'fixed' : 'roughness'} exchange, seaHumidity 0.98: ${Object.entries(fall).map(([k, f]) => `${k} evaporation ${(2.5e6 * f.base / f.area).toFixed(1)} → ${(2.5e6 * f.salty / f.area).toFixed(1)} W/m² (−${(100 * (1 - f.salty / f.base)).toFixed(1)} %)`).join(', ')}`);
    assert.ok(fall['open sea'].salty < (1 - 0.02) * fall['open sea'].base && fall['part ice'].salty < fall['part ice'].base, `the mean evaporation falls by at least the factor's 2 % over open sea (${fall['open sea'].salty / fall['open sea'].base}) and falls over partial ice (${fall['part ice'].salty / fall['part ice'].base})`);
    if (!fixed) rough = { salty, kind };
  }
  const { salty: { exchange, state, core, mesh }, kind } = rough, C = mesh.nCells, b = (core.K - 1) * C;
  const { exnerLayer, exnerLower, g } = core.diagnostics;
  const i = [...Array(C).keys()].find((n) => kind[n] === 'open sea' && Math.abs(mesh.latCell[n]) < 0.5);
  const pi = state[0][i], skin = state[1][b + i] * exnerLayer[b + i] + 3, airV = state[1][b + i] * (1 + VIRTUAL_FACTOR * state[4][b + i]);
  for (const A of [0, 0.5, 1]) {
    const { z, ri } = { ...exchange.cell(i, state[0], state[1], state[4], state[5], skin, 0.5, A, 0, 0, 0, 0, 800) };
    const speed = exchange.wind[i], surfaceV = skin / exnerLower[b + i] * (1 + VIRTUAL_FACTOR * saturationHumidity(skin, pi) * (1 - (1 - 0.98) * (1 - A)));
    close(ri, g * z * (airV - surfaceV) / (0.5 * (airV + surfaceV) * speed * speed), 1e-12, `Ri_b at concentration ${A}`);
  }
  assert.throws(() => createModel(new Grid(16), { topography, surface: { seaHumidity: 0 } }), /seaHumidity must be in \(0, 1\]/);
});

test('both engines take seaHumidity alike: from one state the first step\'s evaporation over every surface and its fall from a factor of 1 agree, and with the fixed exchange so does the fall of the boundary layer\'s buoyancy flux', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  for (const fixed of [false, true]) {
    const exchange = fixed ? { exchange: 'fixed' } : {}, runs = [];
    let kind = null;
    for (const surface of [{ ...exchange, seaHumidity: 1 }, exchange]) {
      const cpu = prepare(createModel(new Grid(16), { topography, surface })), gpu = prepare(await createGpuModel(new Grid(16), { topography, surface }));
      kind = surfaceKinds(cpu);
      cpu.step(900); await gpu.step(900); await gpu.settle();
      const PH = gpu.gpu.layout.PH, C = cpu.mesh.nCells;
      const [evaporation, buoyancy] = await readRanges(gpu.gpu.device, gpu.gpu.buffers.PH, ['EVAP', 'BUOY'].map((name) => ({ offset: PH[name], length: C })));
      runs.push({ cpuE: Float64Array.from(cpu.radiation.evaporation), gpuE: Float64Array.from(evaporation), cpuB: Float64Array.from(cpu.boundaryLayer.buoyancyFlux), gpuB: Float64Array.from(buoyancy) });
    }
    const [d, s] = runs, worst = {};
    for (let i = 0; i < kind.length; i++) {
      const w = (worst[kind[i]] ??= { cells: 0, evaporation: 0, fall: 0, meanFall: 0, buoyancyFall: 0, buoyancyScale: 0, buoyancyFallScale: 0 });
      const cpuFall = d.cpuE[i] - s.cpuE[i], gpuFall = d.gpuE[i] - s.gpuE[i];
      w.cells++;
      w.meanFall += 2.5e6 * cpuFall;
      w.evaporation = Math.max(w.evaporation, 2.5e6 * (Math.abs(s.gpuE[i] - s.cpuE[i]) - 1e-5 * Math.abs(s.cpuE[i])));
      w.fall = Math.max(w.fall, 2.5e6 * (Math.abs(gpuFall - cpuFall) - 2e-5 * Math.abs(d.cpuE[i])));
      w.buoyancyFall = Math.max(w.buoyancyFall, Math.abs((d.gpuB[i] - s.gpuB[i]) - (d.cpuB[i] - s.cpuB[i])));
      w.buoyancyScale = Math.max(w.buoyancyScale, Math.abs(d.cpuB[i]));
      w.buoyancyFallScale = Math.max(w.buoyancyFallScale, Math.abs(d.cpuB[i] - s.cpuB[i]));
    }
    for (const [k, w] of Object.entries(worst)) {
      console.log(`${fixed ? 'fixed' : 'roughness'} exchange, ${k}: ${w.cells} cells, the evaporation falls ${(w.meanFall / w.cells).toFixed(2)} W/m² in the mean; engines apart by ${w.evaporation.toExponential(1)} W/m² in the evaporation and ${w.fall.toExponential(1)} W/m² in its fall beyond 1e-5 of it${fixed ? `, ${w.buoyancyFall.toExponential(1)} m²/s³ in the buoyancy flux's fall (largest fall ${w.buoyancyFallScale.toExponential(1)})` : ''}`);
      assert.ok(w.evaporation < 0.05 && w.fall < 0.05, `${k} evaporation ${w.evaporation}, its fall ${w.fall} W/m²`);
      if (fixed) assert.ok(w.buoyancyFall <= 1e-3 * w.buoyancyScale + 1e-2 * w.buoyancyFallScale, `${k} buoyancy flux's fall apart by ${w.buoyancyFall}`);
    }
    for (const k of ['open sea', 'part ice']) assert.ok(worst[k].meanFall > 0, `${k}: the evaporation falls`);
  }
});

test('surface.virtualDensity lowers the evaporation alike on both engines: the first step, and with the sea skin the second', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  for (const [label, skin, steps] of [['first step', {}, 1], ['sea skin, second step', { warmLayer: true, coolSkin: true }, 2]]) {
    const runs = [];
    for (const surface of [{ virtualDensity: true, ...skin }, skin]) {
      const cpu = prepare(createModel(new Grid(16), { topography, surface })), gpu = prepare(await createGpuModel(new Grid(16), { topography, surface }));
      for (let n = 0; n < steps; n++) { cpu.step(900); await gpu.step(900); }
      await gpu.settle(); await gpu.sync();
      const C = cpu.mesh.nCells, lowest = (cpu.core.K - 1) * C;
      const [evaporation] = await readRanges(gpu.gpu.device, gpu.gpu.buffers.PH, [{ offset: gpu.gpu.layout.PH.EVAP, length: C }]);
      const split = Array.from({ length: C }, (_, i) => Math.abs(gpu.state[4][lowest + i] - cpu.state[4][lowest + i]) > 1e-2 * cpu.state[4][lowest + i]);
      runs.push({ cpu: Float64Array.from(cpu.radiation.evaporation), gpu: Float64Array.from(evaporation), split });
      gpu.destroy();
    }
    const [on, off] = runs, split = on.split.map((s, i) => s || off.split[i]), parted = split.filter(Boolean).length;
    assert.ok(parted <= 0.005 * split.length, `${label}: the lowest layer's humidity parts between the engines on ${parted} of ${split.length} cells`);
    let apart = 0, moveApart = 0, largest = 0, mean = 0;
    for (let i = 0; i < on.cpu.length; i++) {
      if (split[i]) continue;
      const cpuMove = on.cpu[i] - off.cpu[i], gpuMove = on.gpu[i] - off.gpu[i];
      apart = Math.max(apart, 2.5e6 * (Math.abs(on.gpu[i] - on.cpu[i]) - 1e-5 * Math.abs(on.cpu[i])));
      moveApart = Math.max(moveApart, 2.5e6 * (Math.abs(gpuMove - cpuMove) - 2e-5 * Math.abs(on.cpu[i])));
      largest = Math.max(largest, -2.5e6 * cpuMove);
      mean += 2.5e6 * cpuMove / on.cpu.length;
    }
    console.log(`${label}: the evaporation moves by ${mean.toFixed(3)} W/m² in the mean and falls by up to ${largest.toFixed(2)}; the engines are apart by ${apart.toExponential(1)} W/m², their moves by ${moveApart.toExponential(1)}, leaving out ${parted} cells where a moist decision parts them`);
    assert.ok(apart < 0.05 && moveApart < 0.05, `${label}: evaporation ${apart}, its move ${moveApart} W/m²`);
    assert.ok(mean < 0 && largest > 0.1, `${label}: the evaporation moves by ${mean} W/m² in the mean and falls by up to ${largest}`);
  }
});

test('surface.gustiness moves the first step alike on both engines: the evaporation, the sensible heat and the lowest wind, with the fixed exchange and with the roughness one', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  for (const [label, exchange] of [['fixed', { exchange: 'fixed' }], ['roughness', {}]]) {
    const runs = [];
    for (const surface of [{ ...exchange, gustiness: 6 }, exchange]) {
      const cpu = prepare(createModel(new Grid(16), { topography, surface })), gpu = prepare(await createGpuModel(new Grid(16), { topography, surface }));
      cpu.step(900); await gpu.step(900); await gpu.settle(); await gpu.sync();
      const { nCells: C, nEdges: E } = cpu.mesh, PH = gpu.gpu.layout.PH, lowest = (cpu.core.K - 1) * E;
      const [evaporation, heat] = await readRanges(gpu.gpu.device, gpu.gpu.buffers.PH, ['EVAP', 'SFLUX'].map((name) => ({ offset: PH[name], length: C })));
      runs.push({
        cpu: { E: Float64Array.from(cpu.radiation.evaporation, (x) => 2.5e6 * x), H: Float64Array.from(cpu.radiation.surfaceFlux), u: Float64Array.from(cpu.state[2].subarray(lowest, lowest + E)) },
        gpu: { E: Float64Array.from(evaporation, (x) => 2.5e6 * x), H: Float64Array.from(heat), u: Float64Array.from(gpu.state[2].subarray(lowest, lowest + E)) },
      });
      gpu.destroy();
    }
    const [on, off] = runs, flux = {};
    for (const name of ['E', 'H']) {
      const f = (flux[name] = { apart: 0, moveApart: 0, mean: 0, largest: 0 }), n = on.cpu[name].length;
      for (let i = 0; i < n; i++) {
        const cpuMove = on.cpu[name][i] - off.cpu[name][i], gpuMove = on.gpu[name][i] - off.gpu[name][i];
        f.apart = Math.max(f.apart, Math.abs(on.gpu[name][i] - on.cpu[name][i]) - 1e-5 * Math.abs(on.cpu[name][i]));
        f.moveApart = Math.max(f.moveApart, Math.abs(gpuMove - cpuMove) - 2e-5 * Math.abs(on.cpu[name][i]));
        f.mean += cpuMove / n;
        f.largest = Math.max(f.largest, Math.abs(cpuMove));
      }
    }
    let gap = 0, move = 0;
    for (let e = 0; e < on.cpu.u.length; e++) {
      const cpuMove = on.cpu.u[e] - off.cpu.u[e], gpuMove = on.gpu.u[e] - off.gpu.u[e];
      gap += (gpuMove - cpuMove) ** 2;
      move += cpuMove ** 2;
    }
    const windMove = Math.sqrt(move / on.cpu.u.length), windApart = Math.sqrt(gap / move);
    console.log(`${label} exchange, gust 6 against 3 m/s: the evaporation moves ${flux.E.mean.toFixed(3)} W/m² in the mean, ${flux.E.largest.toFixed(2)} at most, the sensible heat ${flux.H.mean.toFixed(3)} and ${flux.H.largest.toFixed(2)}; the engines are apart by ${flux.E.apart.toExponential(1)} and ${flux.H.apart.toExponential(1)} W/m², their moves by ${flux.E.moveApart.toExponential(1)} and ${flux.H.moveApart.toExponential(1)}; the lowest wind moves ${windMove.toExponential(2)} m/s rms, the engines' moves apart by ${windApart.toExponential(1)} of it`);
    for (const name of ['E', 'H']) assert.ok(flux[name].apart < 0.05 && flux[name].moveApart < 0.05, `${label} ${name}: apart ${flux[name].apart}, its move ${flux[name].moveApart} W/m²`);
    assert.ok(windApart < 0.02, `${label}: the lowest wind's moves apart by ${windApart} of them`);
    assert.ok(flux.E.largest > 1 && windMove > 1e-4, `${label}: the gust moves the evaporation by up to ${flux.E.largest} W/m² and the lowest wind by ${windMove} m/s`);
  }
});

test('surface.virtualDensity takes the moist density alike on both engines in the implicit drag, the explicit drag, the fog deposition\'s cap and the boundary layer\'s kinematic surface fluxes', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const boundaryLayer = { countergradient: 7.2 }, moist = { fogDroplets: null, fogDepositionLimit: 1e-3 };
  for (const [label, exchange] of [['roughness', {}], ['fixed', { exchange: 'fixed' }]]) {
    const names = label === 'roughness' ? ['SDRAG', 'CGH', 'CGW', 'FOG'] : ['CGH', 'CGW'], runs = [];
    for (const surface of [{ ...exchange, virtualDensity: true }, exchange]) {
      const cpu = prepare(createModel(new Grid(16), { topography, surface, boundaryLayer, moist })), gpu = prepare(await createGpuModel(new Grid(16), { topography, surface, boundaryLayer, moist }));
      cpu.step(900); await gpu.step(900); await gpu.settle(); await gpu.sync();
      const { nCells: C, nEdges: E } = cpu.mesh, PH = gpu.gpu.layout.PH, lowest = (cpu.core.K - 1) * E, bl = cpu.boundaryLayer;
      const read = await readRanges(gpu.gpu.device, gpu.gpu.buffers.PH, names.map((name) => ({ offset: PH[name], length: C })));
      const mine = { SDRAG: bl.surfaceDrag, CGH: bl.kinematicHeat, CGW: bl.kinematicWater, FOG: cpu.moist.fogPrecipitation }, low = (cpu.core.K - 1) * C;
      runs.push({
        split: Array.from({ length: C }, (_, i) => Math.abs(gpu.state[4][low + i] - cpu.state[4][low + i]) > 1e-2 * cpu.state[4][low + i]), cellsOnEdge: cpu.mesh.cellsOnEdge,
        flux: { CGH: Float64Array.from(cpu.radiation.sensibleHeat), CGW: Float64Array.from(cpu.radiation.evaporation) },
        cpu: { ...Object.fromEntries(names.map((name) => [name, Float64Array.from(mine[name])])), u: Float64Array.from(cpu.state[2].subarray(lowest, lowest + E)) },
        gpu: { ...Object.fromEntries(names.map((name, n) => [name, Float64Array.from(read[n])])), u: Float64Array.from(gpu.state[2].subarray(lowest, lowest + E)) },
      });
      gpu.destroy();
    }
    const [on, off] = runs, report = [], { cellsOnEdge } = on;
    const split = on.split.map((s, i) => s || off.split[i]), parted = split.filter(Boolean).length;
    assert.ok(parted <= 0.005 * split.length, `${label}: the lowest layer's humidity parts between the engines on ${parted} of ${split.length} cells`);
    const left = (name, i) => (name === 'u' ? split[cellsOnEdge[2 * i]] || split[cellsOnEdge[2 * i + 1]] : split[i]);
    report.push(`${parted} cells left out where a moist decision parts the engines`);
    for (const name of [...names, 'u']) {
      let apart = 0, scale = 0, gap = 0, move = 0;
      for (let i = 0; i < on.cpu[name].length; i++) {
        if (left(name, i)) continue;
        const cpuMove = on.cpu[name][i] - off.cpu[name][i], gpuMove = on.gpu[name][i] - off.gpu[name][i];
        apart += (on.gpu[name][i] - on.cpu[name][i]) ** 2; scale += on.cpu[name][i] ** 2;
        gap += (gpuMove - cpuMove) ** 2; move += cpuMove ** 2;
      }
      const relative = Math.sqrt(apart / scale), moved = Math.sqrt(move / scale), moveApart = move > 0 ? Math.sqrt(gap / move) : 0;
      report.push(`${name} apart ${relative.toExponential(1)}, moved ${moved.toExponential(1)}, moves apart ${moveApart.toExponential(1)}`);
      assert.ok(relative < (name === 'FOG' ? 2e-3 : 1e-4), `${label} ${name}: the engines are apart by ${relative} of it`);
      if (name === 'SDRAG' || name === 'u' || name === 'FOG') assert.ok(moved > (name === 'u' ? 1e-5 : 1e-3) && moveApart < (name === 'FOG' ? 1e-2 : 0.02), `${label} ${name}: the option moves it by ${moved} of it, the engines' moves apart by ${moveApart}`);
      else assert.ok(moveApart < 0.05, `${label} ${name}: the engines' moves apart by ${moveApart}`);
      if (label === 'roughness' && (name === 'CGH' || name === 'CGW')) {
        let kinematic = 0, flux = 0, cells = 0;
        for (let i = 0; i < on.flux[name].length; i++) {
          if (!(Math.abs(off.flux[name][i]) > (name === 'CGH' ? 1 : 1e-6)) || !(Math.abs(off.cpu[name][i]) > 0)) continue;
          kinematic += on.cpu[name][i] / off.cpu[name][i] - 1; flux += on.flux[name][i] / off.flux[name][i] - 1; cells++;
        }
        report.push(`${name} over ${cells} cells: the flux moves ${(flux / cells).toExponential(1)}, its kinematic value ${(kinematic / cells).toExponential(1)}`);
        assert.ok(flux < 0 && kinematic > 0.5 * flux, `${label} ${name}: the flux moves ${flux / cells} on average, its kinematic value ${kinematic / cells}, as if divided by the dry density`);
      }
    }
    console.log(`${label} exchange, virtualDensity against none (relative rms): ${report.join('; ')}`);
  }
  assert.throws(() => createModel(new Grid(16), { topography, boundaryLayer: { surfaceExcess: -8.5 } }), /surfaceExcess must be 0 or more/);
  await assert.rejects(createGpuModel(new Grid(16), { topography, boundaryLayer: { surfaceExcess: -8.5 } }), /surfaceExcess must be 0 or more/);
});

test('both engines refuse under boundaryLayer the options they take from the surface or the orography', async () => {
  for (const [key, value, group] of [['virtualDensity', true, 'surface'], ['gustiness', 1, 'surface'], ['implicitDrag', false, 'surface'], ['relativeWind', true, 'surface'], ['seaHumidity', 0.5, 'surface'], ['formDrag', false, 'orography']]) {
    const refused = new RegExp(`boundaryLayer\\.${key} is not a boundary-layer option: set ${group}\\.${key}`);
    assert.throws(() => createModel(new Grid(16), { topography, boundaryLayer: { [key]: value } }), refused);
    if (gpuAvailable) await assert.rejects(createGpuModel(new Grid(16), { topography, boundaryLayer: { [key]: value } }), refused);
  }
});

test('the boundary layer\'s K-profile options, diagnoseFirst, surfaceExcess and the subsidence share, and its search ceiling raised and lowered, move both engines\' mixing alike through two model steps, in stable, convective and entraining columns apart and at their lowest interface, under the roughness and the fixed exchange', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cases = [
    ['defaults', {}, []],
    ['ground, stable profile, Prandtl coefficient 3.6, countergradient', { profileHeight: 'ground', stableProfile: 'holtslagBoville', prandtlCoefficient: 3.6, countergradient: 7.2 }, ['mixing stable', 'mixing convective', 'scalarMixing stable', 'scalarMixing convective']],
    ['ground, capped scale, countergradient', { profileHeight: 'ground', velocityScale: 'capped', countergradient: 7.2 }, ['mixing stable', 'mixing convective']],
    ['stable profile, Prandtl coefficient 3.6', { stableProfile: 'holtslagBoville', prandtlCoefficient: 3.6 }, ['mixing stable', 'scalarMixing stable', 'scalarMixing convective']],
    ['capped scale, countergradient', { velocityScale: 'capped', countergradient: 7.2 }, ['mixing convective']],
    ['diagnoseFirst, surfaceExcess 8.5', { diagnoseFirst: true, surfaceExcess: 8.5 }, ['mixing stable', 'mixing convective', 'scalarMixing convective']],
    ['the subsidence share', { entrainment: { subsidence: true } }, []],
    ['the search ceiling raised to sigma 0.3', { searchTop: 0.3 }, []],
    ['the search ceiling lowered to sigma 0.9', { searchTop: 0.9 }, ['mixing convective']],
  ];
  const apartBound = { stable: 1e-3, convective: 5e-4, entraining: 1e-3 }, failures = [];
  for (const [label, surface] of [['roughness', {}], ['fixed', { exchange: 'fixed' }]]) {
    const runs = {};
    let C = 0, K = 0, kTop = Infinity, g = 0, ground = null;
    for (const [name, boundaryLayer] of cases) {
      const cpu = prepare(createModel(new Grid(16), { topography, surface, boundaryLayer })), gpu = prepare(await createGpuModel(new Grid(16), { topography, surface, boundaryLayer }));
      for (let s = 0; s < 2; s++) { cpu.step(900); await gpu.step(900); }
      await gpu.settle(); await gpu.sync();
      const bl = cpu.boundaryLayer, PH = gpu.gpu.layout.PH, key = { mixing: 'MIX', scalarMixing: 'MIXH', nonlocal: 'CG' };
      C = cpu.mesh.nCells; K = cpu.core.K; kTop = Math.min(kTop, bl.kTop); g = cpu.core.diagnostics.g; ground = Float64Array.from(cpu.surfaceGeopotential, (p) => p / g);
      const low = (K - 1) * C, fields = Object.keys(key).filter((field) => bl[field]);
      const read = await readRanges(gpu.gpu.device, gpu.gpu.buffers.PH, [...fields.map((field) => ({ offset: PH[key[field]], length: bl[field].length })), { offset: PH.DEPTH, length: C }]);
      runs[name] = {
        split: Array.from({ length: C }, (_, i) => Math.abs(gpu.state[4][low + i] - cpu.state[4][low + i]) > 1e-2 * cpu.state[4][low + i]),
        classes: { stable: Array.from(bl.buoyancyFlux, (b) => !(b > 0)), convective: Array.from(bl.buoyancyFlux, (b) => b > 0), entraining: Array.from(bl.entrainment, (w) => w > 0) },
        cpu: Object.fromEntries(fields.map((field) => [field, Float64Array.from(bl[field])])),
        gpu: Object.fromEntries(fields.map((field, n) => [field, Float64Array.from(read[n])])),
        depth: Float64Array.from(bl.depth), gpuDepth: Float64Array.from(read[fields.length]), geopotential: Float64Array.from(cpu.core.diagnostics.geopotential),
        bad: [...fields.flatMap((field, n) => [bl[field], read[n]].map((a) => a.subarray(kTop * C, (K - 1) * C))), bl.depth, read[fields.length], cpu.state[4].subarray(low, low + C), gpu.state[4].subarray(low, low + C)].reduce((n, a) => n + a.filter((v) => !Number.isFinite(v)).length, 0),
      };
      gpu.destroy();
    }
    const control = runs.defaults;
    const columnGap = (a, b, i) => {
      let largest = 0, gap = 0;
      for (let k = kTop; k < K - 1; k++) { largest = Math.max(largest, Math.abs(a[k * C + i]), Math.abs(b[k * C + i])); gap = Math.max(gap, Math.abs(b[k * C + i] - a[k * C + i])); }
      return largest > 0 ? gap / largest : 0;
    };
    const clearance = (run, i) => Math.abs(run.depth[i] - (run.geopotential[(K - 2) * C + i] + run.geopotential[(K - 1) * C + i]) / (2 * g));
    const topsApart = (run, i) => Math.abs(run.depth[i] - run.gpuDepth[i] - ground[i]);
    const flipped = (run, i) => topsApart(run, i) > Math.max(2.5e-3 * (run.depth[i] - run.geopotential[(K - 1) * C + i] / g), 1e-2);
    for (const [name, , acting] of cases) {
      const run = runs[name], where = `${label}, ${name}`, report = [];
      if (run.bad) failures.push(`${where}: ${run.bad} coefficients, depths or lowest-layer humidities are not finite on one engine or the other`);
      const out = run.split.map((s, i) => s || control.split[i]), splits = out.filter(Boolean).length;
      report.push(`${splits} cells out where a moist decision parts the engines`);
      if (!(splits <= 0.005 * C)) failures.push(`${where}: the lowest layer's humidity parts between the engines on ${splits} of ${C} cells`);
      for (const field of Object.keys(run.cpu)) {
        const mine = run.cpu[field], theirs = run.gpu[field], base = name !== 'defaults' && control.cpu[field] ? [control.cpu[field], control.gpu[field]] : null;
        const kept = new Uint8Array(C);
        let parted = 0, tops = 0, flips = 0, topRatio = 0, topRoom = 0, topApart = 0, worst = 0, worstAt = -1;
        const loose = [];
        for (let i = 0; i < C; i++) {
          if (out[i]) continue;
          const own = columnGap(mine, theirs, i), controlGap = base ? columnGap(base[0], base[1], i) : 0;
          if (own > 1e-2 || controlGap > 1e-2) {
            parted++;
            const unexplained = [], atTop = [];
            const why = [own > 1e-2 && [run, own], controlGap > 1e-2 && [control, controlGap]].filter(Boolean).map(([r, gap]) => {
              const apart = topsApart(r, i), room = clearance(r, i);
              if (room < 1 && apart <= 1e-2 && gap <= 4 * apart / room + 1e-3) { atTop.push([gap * room / (2 * apart), room, apart]); return 'top'; }
              if (flipped(r, i)) return 'flip';
              unexplained.push(`column ${i}${r === control ? ' in the defaults' : ''} ${gap.toExponential(1)} apart, its top ${room.toFixed(2)} m from the lowest interface and ${(1e3 * apart).toFixed(1)} mm from the GPU's`);
              return 'loose';
            });
            if (why.includes('loose')) loose.push(unexplained.join(', '));
            else if (why.includes('flip')) flips++;
            else { tops++; for (const [ratio, room, apart] of atTop) { topRatio = Math.max(topRatio, ratio); topRoom = Math.max(topRoom, room); topApart = Math.max(topApart, apart); } }
            continue;
          }
          kept[i] = 1;
          if (own > worst) { worst = own; worstAt = i; }
        }
        if (!(loose.length <= 2 && flips <= 0.005 * C && tops <= 0.005 * C)) failures.push(`${where}, ${field}: of ${parted} columns parted by over 1 % of their largest coefficient, ${tops} where the top meets the lowest interface, ${flips} where the engines' tops lie apart and ${loose.length} that neither explains (${loose.join('; ')})`);
        const classes = [];
        for (const [kind, member] of Object.entries(run.classes)) {
          let apart = 0, scale = 0, lowApart = 0, lowScale = 0, move = 0, gap = 0, before = 0, columns = 0;
          for (let i = 0; i < C; i++) {
            if (!kept[i] || !member[i]) continue;
            columns++;
            const lowest = (K - 2) * C + i;
            lowApart += (theirs[lowest] - mine[lowest]) ** 2; lowScale += Math.max(mine[lowest] ** 2, theirs[lowest] ** 2);
            for (let k = kTop; k < K - 1; k++) {
              const x = k * C + i;
              apart += (theirs[x] - mine[x]) ** 2; scale += Math.max(mine[x] ** 2, theirs[x] ** 2);
              if (!base) continue;
              const cpuMove = mine[x] - base[0][x], gpuMove = theirs[x] - base[1][x];
              move += cpuMove ** 2; gap += (gpuMove - cpuMove) ** 2; before += base[0][x] ** 2;
            }
          }
          const relative = scale > 0 ? Math.sqrt(apart / scale) : 0, lowRelative = lowScale > 0 ? Math.sqrt(lowApart / lowScale) : 0, moved = before > 0 ? Math.sqrt(move / before) : 0, moveGap = move > 0 ? Math.sqrt(gap / move) : gap > 0 ? Infinity : 0;
          classes.push(`${kind} (${columns}) ${relative.toExponential(1)} apart, ${lowRelative.toExponential(1)} at the lowest interface${base ? `, moved ${moved.toExponential(1)}, the moves ${moveGap.toExponential(1)} apart` : ''}`);
          if (!(relative < apartBound[kind] && lowRelative < 2e-4)) failures.push(`${where}, ${field}, ${kind} columns: the engines are ${relative} rms apart, ${lowRelative} at the lowest interface`);
          if (acting.includes(`${field} ${kind}`) && !(columns >= 20 && moved > 1e-3 && moveGap < 0.02)) failures.push(`${where}, ${field}, ${kind} columns: ${columns} of them, the options move it ${moved} rms, the engines' moves ${moveGap} of the move apart`);
        }
        let seat = '';
        if (worstAt >= 0) {
          let at = kTop, top = 0, largest = 0;
          for (let k = kTop; k < K - 1; k++) { const x = k * C + worstAt, d = Math.abs(theirs[x] - mine[x]); largest = Math.max(largest, Math.abs(mine[x])); if (d > top) { top = d; at = k; } }
          const z = (run.geopotential[at * C + worstAt] + run.geopotential[(at + 1) * C + worstAt]) / (2 * g);
          seat = ` (${run.classes.stable[worstAt] ? 'stable' : 'convective'}${run.classes.entraining[worstAt] ? ', entraining' : ''}: its largest gap at ${z.toFixed(0)} m with the top at ${run.depth[worstAt].toFixed(0)} m, where the coefficient is ${(Math.abs(mine[at * C + worstAt]) / largest).toExponential(1)} of the column's largest)`;
        }
        report.push(`${field}: ${parted} columns parted, ${tops} where the top meets the lowest interface (within ${topRoom.toFixed(2)} m of it, their tops up to ${(topApart * 1e3).toFixed(1)} mm apart and their gaps up to ${topRatio.toFixed(2)} of 2|Δh|/(h − z)), ${flips} where the engines' tops lie apart, ${loose.length} unexplained${loose.length ? ` (${loose.join('; ')})` : ''}; ${classes.join('; ')}; the worst kept column ${worst.toExponential(1)} of its largest${seat}`);
        if (name !== 'the subsidence share' || field !== 'mixing') continue;
        let columns = 0, off = 0, worstShare = 0;
        for (let i = 0; i < C; i++) {
          if (!kept[i] || !run.classes.entraining[i]) continue;
          let at = -1, largest = 0, top = 0;
          for (let k = kTop; k < K - 1; k++) { const x = k * C + i, m = Math.abs(mine[x] - base[0][x]); top = Math.max(top, Math.abs(base[0][x])); if (m > largest) { largest = m; at = x; } }
          if (!(largest > 1e-3 * top)) continue;
          columns++;
          const share = Math.abs(theirs[at] - base[1][at] - (mine[at] - base[0][at])) / largest;
          worstShare = Math.max(worstShare, share);
          if (share > 0.05) off++;
        }
        report.push(`the share moves ${columns} entraining columns' exchange by over 1e-3 of their largest coefficient, the engines' moves at that interface part by over 5 % in ${off}, by ${worstShare.toExponential(1)} of it at worst`);
        if (!(columns >= 20 && off <= 1)) failures.push(`${where}: the share moves ${columns} entraining columns, the engines' moves part in ${off}`);
      }
      console.log(`${label} exchange, ${name}: ${report.join('; ')}`);
    }
  }
  assert.deepEqual(failures, []);
});

test('the offline surface exchange diagnoses the state as it stands and takes no step: a diagnoseFirst layer holding deposits as the plain layer does, and after a step', async () => {
  const { diagnoseExchange } = await import('../scripts/surfaceExchange.mjs');
  const [first, plain] = [{ diagnoseFirst: true }, {}].map((boundaryLayer) => prepare(createModel(new Grid(16), { topography, boundaryLayer })));
  first.radiation.sensibleHeat.fill(150); first.radiation.evaporation.fill(6e-5);
  diagnoseExchange(first); diagnoseExchange(plain);
  for (const name of ['depth', 'regime', 'friction', 'buoyancyFlux']) assert.deepEqual(first.boundaryLayer[name], plain.boundaryLayer[name], name);
  first.step(900);
  diagnoseExchange(first);
  assert.ok(first.boundaryLayer.depth.every(Number.isFinite));
});

test('the GPU\'s boundary-layer diagnosis refuses a D_EXL read, which would take the last RK stage\'s exner, and its options leave none', async () => {
  const { PHYSICS_KERNELS, boundaryLayerDiagnosis } = await import('../js/gpu/physics.gpu.js');
  const body = PHYSICS_KERNELS.pblDiagnose;
  assert.equal(boundaryLayerDiagnosis(body, {}), body);
  for (const options of [{ profileHeight: 'ground' }, { profileHeight: 'ground', stableProfile: 'holtslagBoville', prandtlCoefficient: 3.6, diagnoseFirst: true, surfaceExcess: 8.5 }]) assert.ok(boundaryLayerDiagnosis(body, options).includes('blLift'));
  assert.throws(() => boundaryLayerDiagnosis(body.replace('fn blExchange(', 'fn blStale(i: i32) -> f32 { return D[D_EXL + i]; }\nfn blExchange('), {}), /reads D_EXL/);
});

test('surface.virtualDensity\'s GPU rewrite finds each of its densities in every kernel that lays one, and refuses a kernel whose fog cap has moved', async () => {
  const { PHYSICS_KERNELS, virtualDensityKernels } = await import('../js/gpu/physics.gpu.js');
  const rewritten = virtualDensityKernels(PHYSICS_KERNELS);
  for (const name of ['physics', 'physicsSurface', 'pblDiagnose', 'adjustMix', 'adjustPlume', 'adjustRain']) assert.ok(rewritten[name]?.includes('VIRT'), name);
  const cap = 'FOG_VMAX * pressure / (RGAS * temperature))';
  assert.throws(() => virtualDensityKernels({ ...PHYSICS_KERNELS, adjustRain: PHYSICS_KERNELS.adjustRain.replace(cap, 'FOG_VMAX * pressure / (temperature * RGAS))') }), /adjustRain kernel lays no surface density/);
});

test('the CPU model refuses ocean.gustiness: its ocean and ice drift take the gust floor of surface.gustiness', () => {
  assert.throws(() => createModel(new Grid(16), { topography, ocean: { gustiness: 5 } }), /ocean\.gustiness is not an option/);
});

test('both coupled engines hand the air the same sea velocity under surface.relativeWind, with the ice drift and without it, and surface.gustiness and surface.virtualDensity move the ocean and the ice drift alike through their stress under the fixed exchange, and both refuse ocean.gustiness', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const configs = [
    ['drift', { surface: { relativeWind: true } }],
    ['no drift', { surface: { relativeWind: true }, ocean: { iceDrift: false } }],
    ['fixed exchange', { surface: { exchange: 'fixed' } }],
    ['fixed exchange, gustiness 6', { surface: { exchange: 'fixed', gustiness: 6 } }],
    ['fixed exchange, virtualDensity', { surface: { exchange: 'fixed', virtualDensity: true } }],
  ];
  const flows = {}, drifts = {}, failures = [], report = [];
  for (const [name, { surface, ocean = {} }] of configs) {
    const options = { topography, surface, ocean: { everySteps: 2, ...ocean } };
    const cpu = prepare(createModel(new Grid(16), options)), gpu = prepare(await createGpuModel(new Grid(16), options));
    for (let n = 0; n < 8; n++) { cpu.step(900); await gpu.step(900); }
    await gpu.settle(); await gpu.sync();
    if (surface.relativeWind) {
      const [flow] = await readRanges(gpu.gpu.device, gpu.gpu.buffers.PH, [{ offset: gpu.gpu.layout.PH.UO, length: cpu.mesh.nEdges }]);
      flows[name] = { cpu: Float64Array.from(cpu.seaFlow), gpu: Float64Array.from(flow) };
    } else flows[name] = { cpu: Float64Array.from(cpu.ocean.serialize().u), gpu: Float64Array.from((await gpu.ocean.serialize()).u) };
    if (!surface.relativeWind) drifts[name] = { cpu: Float64Array.from(cpu.ocean.iceDrift.velocity), gpu: Float64Array.from((await gpu.oceanEngine.restartArrays()).iceVelocity) };
    gpu.destroy();
  }
  const rms = (a, b = null) => Math.sqrt(a.reduce((sum, v, e) => sum + (v - (b ? b[e] : 0)) ** 2, 0) / a.length);
  for (const [name] of configs) {
    const { cpu, gpu } = flows[name], scale = rms(cpu), apart = rms(gpu, cpu) / scale;
    report.push(`${name}: ${scale.toExponential(1)} m/s rms, the engines ${apart.toExponential(1)} of it apart`);
    if (!(scale > 1e-4 && apart < 2e-4)) failures.push(`${name}: the velocity is ${scale} m/s rms, the engines ${apart} of it apart`);
  }
  for (const [name, from, least, spread] of [['no drift', 'drift', 1e-2, 0.05], ['fixed exchange, gustiness 6', 'fixed exchange', 1e-2, 0.05], ['fixed exchange, virtualDensity', 'fixed exchange', 1e-4, 0.5]]) {
    const cpuMove = flows[name].cpu.map((v, e) => v - flows[from].cpu[e]), gpuMove = flows[name].gpu.map((v, e) => v - flows[from].gpu[e]);
    const move = rms(cpuMove) / rms(flows[from].cpu), gap = rms(gpuMove, cpuMove) / rms(cpuMove);
    report.push(`${name} against ${from}: it moves ${move.toExponential(1)} of it, the engines' moves ${gap.toExponential(1)} apart`);
    if (!(move > least && gap < spread)) failures.push(`${name} against ${from}: the velocity moves ${move} of it, the engines' moves ${gap} of the move apart`);
  }
  const [still, moved] = [drifts['fixed exchange'], drifts['fixed exchange, virtualDensity']];
  const cpuDrift = moved.cpu.map((v, e) => v - still.cpu[e]), gpuDrift = moved.gpu.map((v, e) => v - still.gpu[e]);
  const driftScale = rms(still.cpu), driftApart = rms(still.gpu, still.cpu) / driftScale, driftMove = rms(cpuDrift) / driftScale, driftGap = rms(gpuDrift, cpuDrift) / rms(cpuDrift);
  report.push(`the ice drift: ${driftScale.toExponential(1)} m/s rms, the engines ${driftApart.toExponential(1)} of it apart; virtualDensity moves it ${driftMove.toExponential(1)} of it, the engines' moves ${driftGap.toExponential(1)} apart`);
  if (!(still.cpu.length === still.gpu.length && driftApart < 2e-4 && driftMove > 1e-4 && driftGap < 0.1)) failures.push(`the ice drift: the engines ${driftApart} of it apart, virtualDensity moves it ${driftMove} of it, the engines' moves ${driftGap} of the move apart`);
  console.log(`the sea velocity under the air, then the ocean's: ${report.join('; ')}`);
  assert.deepEqual(failures, []);
  await assert.rejects(createGpuModel(new Grid(16), { topography, ocean: { gustiness: 5 } }), /ocean\.gustiness is not an option/);
});

test('under the fixed exchange the 2 m temperature is NaN at every cell on both engines', async () => {
  const options = { topography, surface: { exchange: 'fixed' } };
  const cpu = prepare(createModel(new Grid(16), options));
  cpu.step(900);
  assert.equal(cpu.exchange.fixed, true);
  assert.ok(cpuFrameFields(cpu, { fields: ['t2m'] }).t2m.every(Number.isNaN), 'the CPU frame');
  if (!gpuAvailable) return;
  const gpu = prepare(await createGpuModel(new Grid(16), options));
  await gpu.step(900);
  await gpu.settle();
  const [t2m] = await readRanges(gpu.gpu.device, gpu.gpu.buffers.PH, [{ offset: gpu.gpu.layout.PH.T2M, length: cpu.mesh.nCells }]);
  assert.ok(t2m.every(Number.isNaN), 'the GPU');
});
