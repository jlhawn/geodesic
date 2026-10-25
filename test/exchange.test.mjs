import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { SEA_DRAG, LAND_DRAG } from '../js/physics/surface.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { VIRTUAL_FACTOR, sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { transfer, psiMomentum, psiHeat, andreasScalar, seaIceRoughness, charnockRoughness, createBlend, referenceCoefficient, exchangeMode, airViscosity, KARMAN } from '../js/physics/exchange.module.js';
import { PREVIOUS_LEVELS } from './helpers/previousDefaults.mjs';

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
  const [drag, heat, reference, evaporation, surfaceFlux] = await readRanges(gpu.gpu.device, gpu.gpu.buffers.PH, ['DRAG', 'HEATX', 'REFX', 'EVAP', 'SFLUX'].map((name) => ({ offset: PH[name], length: C })));
  const { land, iceSheet } = cpu.geography;
  const classOf = (i) => (land[i] ? (iceSheet[i] ? 'ice sheet' : cpu.land.snow[i] > 30 ? 'snow' : cpu.land.canopy[i] > 0.5 ? 'forest' : 'other land') : cpu.state[6][i] > 0 ? 'sea ice' : 'open sea');
  const worst = {};
  for (let i = 0; i < C; i++) {
    const k = classOf(i), w = (worst[k] ??= { cells: 0, drag: 0, heat: 0, reference: 0, evaporation: 0, flux: 0, cd: 0 });
    w.cells++;
    w.cd += cpu.exchange.drag[i];
    w.drag = Math.max(w.drag, Math.abs(drag[i] - cpu.exchange.drag[i]) / Math.max(cpu.exchange.drag[i], 1e-5));
    w.heat = Math.max(w.heat, Math.abs(heat[i] - cpu.exchange.heat[i]) / Math.max(cpu.exchange.heat[i], 1e-5));
    w.reference = Math.max(w.reference, Math.abs(reference[i] / cpu.exchange.reference[i] - 1));
    w.evaporation = Math.max(w.evaporation, (Math.abs(evaporation[i] - cpu.radiation.evaporation[i]) * 2.5e6 - 1e-5 * Math.abs(cpu.radiation.evaporation[i]) * 2.5e6));
    w.flux = Math.max(w.flux, Math.abs(surfaceFlux[i] - cpu.radiation.surfaceFlux[i]) - 1e-5 * Math.abs(cpu.radiation.surfaceFlux[i]));
  }
  for (const [k, w] of Object.entries(worst)) {
    console.log(`${k}: ${w.cells} cells, mean C_D ${(w.cd / w.cells * 1e3).toFixed(3)}e-3; engines apart: C_D ${w.drag.toExponential(1)} and C_H ${w.heat.toExponential(1)} relative to the coefficient or 1e-5 if larger, reference ${w.reference.toExponential(1)} relative; latent heat ${w.evaporation.toExponential(1)}, net surface flux ${w.flux.toExponential(1)} W/m² beyond 1e-5 of each cell's own`);
    assert.ok(w.drag < bound && w.heat < bound && w.reference < bound / 20, `${k} coefficients against ${bound.toExponential(1)}`);
    assert.ok(w.evaporation < 0.05 && w.flux < 0.1, `${k} fluxes`);
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
  const model = prepare(createModel(new Grid(16), { topography }));
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
