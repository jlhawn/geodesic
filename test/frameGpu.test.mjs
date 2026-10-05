import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { levelFields, verticalVelocity, smoothCells, dewPoint, wetBulb, miseryIndex } from '../js/levels.module.js';
import { VERTICAL_MEMORY } from '../js/frames.module.js';
import { depthFields } from '../js/ocean/layered.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { FIELDS, CLOUD_TYPES, CLOUD_LOW_PRESSURE, CLOUD_HIGH_PRESSURE, CLOUD_OPACITY_PATH, CLOUD_SEEN, visibleCloudHeights } from '../js/frames.module.js';
import { createModel } from '../js/model.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 + 2500 * Math.exp(-(((lat - 0.3) / 0.3) ** 2)) : -4000));

function worst(reference, got, { skipNaN = false } = {}) {
  let max = 0, at = -1;
  for (let x = 0; x < reference.length; x++) {
    if (skipNaN && Number.isNaN(reference[x])) { assert.ok(Number.isNaN(got[x]), `NaN expected at ${x}`); continue; }
    const d = Math.abs(reference[x] - got[x]);
    if (!(d <= max)) { max = d; at = x; }
  }
  return { max, at };
}

test('the GPU frame matches the fields and diagnostics computed from the full state in double precision', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const model = await createGpuModel(new Grid(6), { topography });
  const { mesh, core, state } = model;
  const C = mesh.nCells, E = mesh.nEdges;
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < C; i++) if (model.geography.land[i]) state[6][i] = 0;
  model.load();
  model.ocean.initialize(state[3], state[6]);
  model.land.initialize();
  for (let n = 0; n < 12; n++) await model.step(900);
  const runoffPattern = Float32Array.from({ length: C }, (_, i) => (model.geography.land[i] ? 0.25 + 1e-3 * (i % 97) : 0));
  model.gpu.device.queue.writeBuffer(model.gpu.buffers.PH, 4 * model.gpu.layout.PH.RUNOFF, runoffPattern);

  const physics = await model.gpu.downloadPhysics();
  await model.sync();
  const ocean = await model.oceanEngine.download();
  const [pi, theta, u, surfaceT, q, qc, ice] = state;
  core.diagnose(pi, theta, q, qc);

  const all = Object.keys(FIELDS);
  const first = await model.beginFrame({ level: 'surface', fields: all, diagnostics: true });
  assert.deepEqual(Object.keys(first.fields).sort(), all.slice().sort());

  for (const level of ['surface', 850, 500, 250]) {
    const frame = level === 'surface' ? first : await model.beginFrame({ level, fields: all });
    const layerWind = (k) => cellVector(mesh, u.subarray(k * E, (k + 1) * E));
    const reference = levelFields(core, pi, theta, layerWind, level, q);
    const comfort = (fn) => Float64Array.from(reference.temperature, (t, i) => fn(t - 273.15, reference.humidity[i], reference.speed[i]) + 273.15);
    const vertical = smoothCells(mesh, verticalVelocity(mesh, core, pi, u, level, reference.temperature));
    const checks = {
      temperature: [reference.temperature, 2e-3], height: [reference.height, 0.1], humidity: [reference.humidity, 5e-5], speed: [reference.speed, 5e-4], wind: [reference.vector, 5e-4],
      dewPoint: [comfort(dewPoint), 5e-3], wetBulb: [comfort(wetBulb), 5e-3], misery: [comfort(miseryIndex), 5e-3], vertical: [vertical, 2e-3],
    };
    for (const [name, [values, tolerance]] of Object.entries(checks)) {
      const { max, at } = worst(values, frame.fields[name]);
      assert.ok(max < tolerance, `${name} at ${level}: ${max} at ${at} (${values[at]} against ${frame.fields[name][at]})`);
    }
  }

  const f = first.fields;
  const { R, g, exnerLayer } = core.diagnostics, phis = model.surfaceGeopotential, bottom = (core.K - 1) * C;
  const mslp = Float64Array.from(pi, (p, i) => p * Math.exp(phis[i] / (R * (theta[bottom + i] * exnerLayer[bottom + i] + 0.00325 * phis[i] / g))));
  const water = Float64Array.from({ length: C }, (_, i) => model.moist.columnWater(pi, q, i));
  const condensed = Float64Array.from({ length: C }, (_, i) => model.moist.columnWater(pi, qc, i));
  const { dSigma, g: gravity } = core.diagnostics, K0 = model.moist.cumulusK0;
  const cumulus = Float64Array.from({ length: C }, (_, i) => { let path = 0; for (let k = K0; k < core.K; k++) { const slot = (k - K0) * C + i; path += pi[i] * dSigma[k] / gravity * physics.CUCOVER[slot] * physics.CUWATER[slot]; } return path; });
  const cloud = Float64Array.from(condensed, (w, i) => w + cumulus[i] + physics.DECKF[i] * physics.DECK[i]);
  const { sigmaMid } = core.diagnostics;
  const condensedWhere = (inside) => Float64Array.from({ length: C }, (_, i) => { let path = 0; for (let k = 0; k < core.K; k++) if (inside(pi[i] * sigmaMid[k])) path += pi[i] * dSigma[k] / gravity * qc[k * C + i]; return path; });
  const types = {
    cloudLow: condensedWhere((p) => p > CLOUD_LOW_PRESSURE), cloudMid: condensedWhere((p) => p <= CLOUD_LOW_PRESSURE && p > CLOUD_HIGH_PRESSURE), cloudHigh: condensedWhere((p) => p <= CLOUD_HIGH_PRESSURE),
    cloudCumulus: cumulus, cloudDeck: Float64Array.from({ length: C }, (_, i) => physics.DECKF[i] * physics.DECK[i]),
  };
  for (const name of ['cloudLow', 'cloudMid', 'cloudHigh']) assert.ok(types[name].some((v) => v > 1e-4), `some ${name}`);
  const visible = { top: new Float64Array(C), base: new Float64Array(C) }, seen = { top: Float64Array.from(f.cloudTop), base: Float64Array.from(f.cloudBase) };
  {
    const z = new Float64Array(core.K), paths = new Float64Array(core.K), out = [0, 0], geopotential = core.diagnostics.geopotential;
    for (let i = 0; i < C; i++) {
      let total = physics.DECKF[i] * physics.DECK[i];
      for (let k = 0; k < core.K; k++) {
        z[k] = geopotential[k * C + i] / gravity;
        const plume = k >= K0 ? physics.CUCOVER[(k - K0) * C + i] * physics.CUWATER[(k - K0) * C + i] : 0;
        paths[k] = pi[i] * dSigma[k] / gravity * (Math.max(0, qc[k * C + i]) + plume);
        total += paths[k];
      }
      visibleCloudHeights(core.K, z, phis[i] / gravity, paths, physics.DECKF[i] * physics.DECK[i], Math.max(physics.DEPTH[i], physics.MLMTOP[i]) + phis[i] / gravity, out);
      const opacity = 1 - Math.exp(-total / CLOUD_OPACITY_PATH);
      if (Math.abs(opacity - CLOUD_SEEN) < 1e-4) { visible.top[i] = visible.base[i] = seen.top[i] = seen.base[i] = NaN; continue; }
      visible.top[i] = out[0]; visible.base[i] = out[1];
      if (opacity < CLOUD_SEEN) assert.ok(f.cloudTop[i] === 0 && f.cloudBase[i] === 0, `a clear column has no cloud heights at ${i}`);
      else assert.ok(f.cloudBase[i] <= f.cloudTop[i], `cloud base ${f.cloudBase[i]} above the top ${f.cloudTop[i]} at ${i}`);
    }
  }
  assert.ok(f.cloudTop.some((height) => height > 6000), 'some cloud tops above 6 km');
  const sea = (source) => Float64Array.from({ length: C }, (_, i) => (model.geography.land[i] ? NaN : source[i]));
  const currents = cellVector(mesh, ocean.u1);
  for (let i = 0; i < C; i++) if (model.geography.land[i]) currents.fill(0, 3 * i, 3 * i + 3);
  const current = sea(Float64Array.from({ length: C }, (_, i) => Math.hypot(currents[3 * i], currents[3 * i + 1], currents[3 * i + 2])));
  const checks = {
    ps: [pi, 0.1], mslp: [mslp, 0.5], water: [water, 1e-3], cloud: [cloud, 1e-6], rain: [physics.RAIN.subarray(0, C), 1e-6], ice: [ice, 1e-6], concentration: [model.seaIce.concentration, 1e-6],
    albedo: [physics.ADIF.subarray(0, C), 1e-6], shortwave: [physics.SWDN.subarray(0, C), 1e-3], longwave: [physics.OLR.subarray(0, C), 1e-3],
    soil: [physics.SOIL.subarray(0, C), 1e-4], snow: [physics.SNOW.subarray(0, C), 1e-4],
    sst: [sea(ocean.T1), 1e-3], sss: [sea(ocean.S1), 1e-4], layerDepth: [sea(ocean.h1), 1e-3], thermocline: [ocean.thermoclineDepth, 1e-2], ssh: [sea(ocean.eta), 1e-5],
    current: [current, 1e-5], currents: [currents, 1e-5],
    ...Object.fromEntries(Object.entries(types).map(([name, values]) => [name, [values, 1e-6]])),
  };
  const misses = ['top', 'base'].map((name) => {
    const { max, at } = worst(visible[name], seen[name], { skipNaN: true });
    assert.ok(max < 2, `cloud ${name}: ${max} m at ${at} (${visible[name][at]} against ${seen[name][at]})`);
    return max;
  });
  console.log(`cloud tops up to ${Math.max(...f.cloudTop).toFixed(0)} m, the GPU's within ${misses[0].toFixed(2)} m (tops) and ${misses[1].toFixed(2)} m (bases) of double precision`);
  for (const [name, [values, tolerance]] of Object.entries(checks)) {
    const { max, at } = worst(values, f[name], { skipNaN: true });
    assert.ok(max < tolerance, `${name}: ${max} at ${at} (${values[at]} against ${f[name][at]})`);
  }

  const summed = Float64Array.from({ length: C }, (_, i) => CLOUD_TYPES.reduce((total, name) => total + f[name][i], 0));
  { const { max, at } = worst(f.cloud, summed); assert.ok(max <= 1e-6 * Math.max(...f.cloud), `the types sum to the cloud: ${max} at ${at}`); }

  const land = physics.LAND, d = first.diagnostics;
  let area = 0, mass = 0, ts = 0, piMin = Infinity, piMax = -Infinity, maxWind = 0, w = 0, c = 0, rain = 0, iceArea = 0, iceVolume = 0, olr = 0, absorbed = 0, landArea = 0, landT = 0, soil = 0;
  for (let i = 0; i < C; i++) {
    const a = mesh.areaCell[i];
    area += a; mass += a * pi[i]; ts += a * surfaceT[i]; piMin = Math.min(piMin, pi[i]); piMax = Math.max(piMax, pi[i]);
    w += a * water[i]; c += a * condensed[i]; rain += a * physics.RAIN[i]; olr += a * physics.OLR[i]; absorbed += a * physics.ABS[i];
    if (ice[i] > 0) { const cover = model.seaIce.cover(i, ice[i]); iceArea += a * cover; iceVolume += a * cover * ice[i]; }
    if (land[i] > 0.5) { landArea += a; landT += a * surfaceT[i]; soil += a * physics.SOIL[i]; }
  }
  for (const x of u) maxWind = Math.max(maxWind, Math.abs(x));
  const close = (name, expected, got, relative = 1e-5) => assert.ok(Math.abs(got - expected) <= relative * Math.abs(expected) + 1e-12, `${name}: ${got} against ${expected}`);
  close('mass', mass / area, d.mass); close('meanSurfaceT', ts / area, d.meanSurfaceT);
  close('piMin', piMin, d.piMin); close('piMax', piMax, d.piMax); close('maxWind', maxWind, d.maxWind);
  close('columnWater', w / area, d.columnWater); close('columnCloud', c / area, d.columnCloud, 1e-4);
  close('precipitation', rain / area / model.time, d.precipitation, 1e-4);
  close('iceFraction', iceArea / area, d.iceFraction); close('iceThickness', iceArea > 0 ? iceVolume / iceArea : 0, d.iceThickness);
  close('outgoingLongwave', olr / area, d.instantaneous.outgoingLongwave); close('absorbedSolar', absorbed / area, d.instantaneous.absorbedSolar);
  close('landMeanT', landT / landArea, d.landMeanT); close('soilWater', soil / landArea, d.soilWater);
  let runoff = 0;
  for (let i = 0; i < C; i++) { runoff += mesh.areaCell[i] * physics.RUNOFF[i]; assert.ok(Math.abs(model.land.runoff[i] - physics.RUNOFF[i]) <= 1e-6 * Math.abs(physics.RUNOFF[i]) + 1e-12, `runoff at ${i}`); }
  assert.ok(runoff > 0);
  close('runoff', runoff / area, d.runoff);
  let oceanArea = 0, depth = 0, ssh = 0;
  for (let i = 0; i < C; i++) if (!model.geography.land[i]) { const a = mesh.areaCell[i]; oceanArea += a; depth += a * ocean.h1[i]; ssh = Math.max(ssh, Math.abs(ocean.eta[i])); }
  close('oceanUpperDepth', depth / oceanArea, d.oceanUpperDepth); close('oceanSSH', ssh, d.oceanSSH);
  console.log(`GPU frame at N=6: Ts ${d.meanSurfaceT.toFixed(3)} K, OLR ${d.outgoingLongwave.toFixed(2)} W/m², precipitation ${(86400 * d.precipitation).toFixed(3)} mm/day, ocean h1 ${d.oceanUpperDepth.toFixed(2)} m`);

  const column = await model.oceanEngine.serialize();
  const L = model.oceanEngine.layers, cellOcean = model.oceanEngine.cellOcean;
  for (const depth of [137, 903]) {
    const at = await model.beginFrame({ depth, fields: ['sst', 'current', 'currents', 'upwelling'] });
    assert.equal(at.depth, depth);
    const expected = depthFields(mesh, L, { h: column.h, u: column.u, temperature: (k, i) => column.T[k * C + i], cellOcean }, depth);
    const speed = expected.speed;
    let magnitude = 0;
    for (let i = 0; i < C; i++) if (!Number.isNaN(expected.upwelling[i])) magnitude = Math.max(magnitude, Math.abs(expected.upwelling[i]));
    assert.ok(magnitude > 0);
    for (const [name, [values, tolerance]] of Object.entries({ sst: [expected.temperature, 1e-3], current: [speed, 1e-5], currents: [expected.current, 1e-5], upwelling: [expected.upwelling, 1e-3 * magnitude] })) {
      const { max, at: worstAt } = worst(values, at.fields[name], { skipNaN: true });
      assert.ok(max < tolerance, `${name} at ${depth} m: ${max} at ${worstAt} (${values[worstAt]} against ${at.fields[name][worstAt]})`);
    }
  }

  const floor = await model.beginFrame({ depth: 4500, fields: ['sst', 'current', 'currents', 'upwelling'] });
  for (let i = 0; i < C; i++) {
    assert.ok(Number.isNaN(floor.fields.sst[i]) && Number.isNaN(floor.fields.current[i]) && Number.isNaN(floor.fields.upwelling[i]), `the sea floor lies above 4500 m at ${i}`);
    for (let c = 0; c < 3; c++) assert.equal(floor.fields.currents[3 * i + c], 0);
  }

  const unsubscribed = await model.beginFrame({ level: 'surface', fields: ['sst'] });
  assert.deepEqual(Object.keys(unsubscribed.fields), ['sst']);
  assert.equal(unsubscribed.diagnostics, null);
  const paused = await model.beginFrame({ diagnostics: true });
  assert.ok(paused.diagnostics.precipitation > 0, 'a frame with no time elapsed reports the recent rain rate');

  const held = await model.beginFrame({ level: 250, fields: ['vertical'] });
  const again = await model.beginFrame({ level: 250, fields: ['vertical'] });
  assert.deepEqual(again.fields.vertical, held.fields.vertical, 'no time elapsed, the memory holds');
  const dt = 900;
  await model.step(dt);
  const later = await model.beginFrame({ level: 250, fields: ['vertical'] });
  await model.sync();
  core.diagnose(pi, theta, q, qc);
  const fresh = smoothCells(mesh, verticalVelocity(mesh, core, pi, u, 250, levelFields(core, pi, theta, (k) => cellVector(mesh, u.subarray(k * E, (k + 1) * E)), 250, q).temperature));
  const keep = Math.exp(-dt / VERTICAL_MEMORY);
  const blend = Float64Array.from(fresh, (w, i) => keep * held.fields.vertical[i] + (1 - keep) * w);
  { const { max, at } = worst(blend, later.fields.vertical); assert.ok(max < 2e-3, `the two-hour memory blends the new frame in: ${max} at ${at}`); }
});

test('the EIS deck shows in the cloud field, the same in both engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const run = async (stratus) => {
    const radiation = { stratus, mixedLayerDeck: false };
    const cpu = createModel(new Grid(6), { ocean: false, radiation }), gpu = await createGpuModel(new Grid(6), { ocean: false, radiation });
    const C = cpu.mesh.nCells, { K, sigmaMid } = cpu.core;
    const init = initializeState(cpu, {});
    for (let a = 0; a < init.length; a++) { cpu.state[a].set(init[a]); gpu.state[a].set(init[a]); }
    for (let k = 0; k < K; k++) if (sigmaMid[k] < 0.75) for (let i = 0; i < C; i++) { cpu.state[1][k * C + i] += 10; gpu.state[1][k * C + i] += 10; }
    gpu.load();
    for (let n = 0; n < 2; n++) { cpu.step(900); await gpu.step(900); }
    const frame = await gpu.beginFrame({ fields: ['cloud', ...CLOUD_TYPES, 'cloudTop', 'cloudBase'] }), physics = await gpu.gpu.downloadPhysics(), K0 = gpu.moist.cumulusK0, { dSigma, g } = cpu.core.diagnostics, pi = cpu.state[0];
    const gpuCumulus = Float64Array.from({ length: C }, (_, i) => { let path = 0; for (let k = K0; k < K; k++) { const slot = (k - K0) * C + i; path += pi[i] * dSigma[k] / g * physics.CUCOVER[slot] * physics.CUWATER[slot]; } return path; });
    const parted = [];
    for (let i = 0; i < C; i++) if ((cpu.moist.cumulusCloudPath(pi, i) > 0) !== (gpuCumulus[i] > 0)) parted.push(i);
    const kept = (values) => Float64Array.from(values, (x, i) => (parted.includes(i) ? NaN : x));
    const parts = Array.from({ length: C }, (_, i) => cpu.cloudParts(i));
    const types = Object.fromEntries(CLOUD_TYPES.map((name) => [name, { cpu: kept(parts.map((part) => part[name])), gpu: kept(frame.fields[name]) }]));
    const summed = Float64Array.from({ length: C }, (_, i) => CLOUD_TYPES.reduce((total, name) => total + frame.fields[name][i], 0));
    cpu.core.diagnose(cpu.state[0], cpu.state[1], cpu.state[4], cpu.state[5]);
    const heights = cpu.cloudHeights();
    const seen = { top: { cpu: kept(heights.top), gpu: kept(frame.fields.cloudTop) }, base: { cpu: kept(heights.base), gpu: kept(frame.fields.cloudBase) } };
    return { cpu: kept(Array.from({ length: C }, (_, i) => cpu.cloudWater(i))), gpu: kept(frame.fields.cloud), types, summed, all: frame.fields.cloud, deck: Float64Array.from({ length: C }, (_, i) => cpu.radiation.stratusFraction[i] * cpu.radiation.stratus[i]), parted, seen };
  };
  const on = await run(true), off = await run(false);
  for (const [label, r] of Object.entries({ on, off })) {
    for (const [name, { cpu, gpu }] of Object.entries(r.types)) { const { max, at } = worst(cpu, gpu, { skipNaN: true }); assert.ok(max < 1e-4, `deck ${label}: ${name} differs between the engines by ${max} at ${at}`); }
    for (const [name, { cpu, gpu }] of Object.entries(r.seen)) { const { max, at } = worst(cpu, gpu, { skipNaN: true }); assert.ok(max < 5, `deck ${label}: the cloud ${name} differs between the engines by ${max} m at ${at} (${cpu[at]} against ${gpu[at]})`); }
    const { max, at } = worst(r.all, r.summed);
    assert.ok(max <= 1e-6 * Math.max(...r.all), `deck ${label}: the GPU's types sum to its cloud to ${max} at ${at}`);
  }
  assert.ok(on.parted.length + off.parted.length <= 0.02 * on.cpu.length, `the cumulus fired on one engine only in ${on.parted.length} + ${off.parted.length} of ${on.cpu.length} columns`);
  let at = 0;
  for (let i = 0; i < on.deck.length; i++) if (on.deck[i] > on.deck[at] && Number.isFinite(on.cpu[i]) && Number.isFinite(off.cpu[i])) at = i;
  const engines = worst(on.cpu, on.gpu, { skipNaN: true });
  console.log(`under a 10 K inversion the thickest deck adds ${(1000 * on.deck[at]).toFixed(1)} g/m² to cell ${at}'s cloud: ${(1000 * on.cpu[at]).toFixed(1)} against ${(1000 * off.cpu[at]).toFixed(1)} g/m² without it; the engines' cloud fields differ by at most ${engines.max.toExponential(1)} kg/m²`);
  assert.ok(on.deck[at] > 1e-3, `deck ${on.deck[at]} kg/m²`);
  assert.ok(on.cpu[at] > off.cpu[at] + 0.5 * on.deck[at] && on.gpu[at] > off.gpu[at] + 0.5 * on.deck[at], `cloud ${on.cpu[at]} (GPU ${on.gpu[at]}) against ${off.cpu[at]} (GPU ${off.gpu[at]}) without the deck`);
  assert.ok(engines.max < 1e-4, `cloud differs between the engines by ${engines.max} at ${engines.at}`);
  const apart = worst(off.cpu, off.gpu, { skipNaN: true });
  assert.ok(apart.max < 1e-4, `without the deck the engines' cloud fields differ by ${apart.max} kg/m² at ${apart.at} (${off.cpu[apart.at]} against ${off.gpu[apart.at]})`);
});

test('the mixed-layer deck shows in the cloud field, the same in both engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const radiation = {};
  const surface = { exchange: 'fixed' };
  const cpu = createModel(new Grid(6), { ocean: false, radiation, surface }), gpu = await createGpuModel(new Grid(6), { ocean: false, radiation, surface });
  const C = cpu.mesh.nCells, { K, sigmaMid } = cpu.core;
  const init = initializeState(cpu, {});
  for (let a = 0; a < init.length; a++) cpu.state[a].set(init[a]);
  const theta = cpu.state[1], q = cpu.state[4];
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) {
    if (sigmaMid[k] < 0.85) theta[k * C + i] += 10;
    else { theta[k * C + i] = theta[(K - 1) * C + i]; q[k * C + i] = q[(K - 1) * C + i]; }
  }
  for (let a = 0; a < init.length; a++) gpu.state[a].set(cpu.state[a]);
  cpu.radiation.mlmSubsidence.fill(-1e-3);
  gpu.radiation.mlmSubsidence.fill(-1e-3);
  gpu.load();
  for (let n = 0; n < 2; n++) { cpu.step(900); await gpu.step(900); }
  const frame = await gpu.beginFrame({ fields: ['cloud', ...CLOUD_TYPES, 'cloudTop', 'cloudBase'] });
  const reference = Float64Array.from({ length: C }, (_, i) => cpu.cloudWater(i));
  const parts = Array.from({ length: C }, (_, i) => cpu.cloudParts(i));
  for (const name of CLOUD_TYPES) { const { max, at } = worst(parts.map((part) => part[name]), frame.fields[name]); assert.ok(max < 1e-4, `${name} differs between the engines by ${max} at ${at}`); }
  assert.ok(parts.some((part) => part.cloudDeck > 0.05), 'the deck shows in its own field');
  let decked = 0, added = 0;
  for (let i = 0; i < C; i++) if (cpu.radiation.mlmCover[i] > 0) { decked++; added += cpu.radiation.stratusFraction[i] * cpu.radiation.stratus[i]; }
  const engines = worst(reference, frame.fields.cloud);
  console.log(`the mixed-layer deck on ${decked} of ${C} cells adds ${(1000 * added / Math.max(1, decked)).toFixed(1)} g/m² to their cloud; the engines' cloud fields differ by at most ${engines.max.toExponential(1)} kg/m²`);
  assert.ok(decked > C / 2 && added / decked > 0.05, `deck on ${decked} cells adding ${added / decked} kg/m²`);
  assert.ok(engines.max < 1e-4, `cloud differs between the engines by ${engines.max} at ${engines.at}`);

  cpu.core.diagnose(cpu.state[0], cpu.state[1], cpu.state[4], cpu.state[5]);
  const heights = cpu.cloudHeights();
  for (const [name, values] of Object.entries({ cloudTop: heights.top, cloudBase: heights.base })) {
    const { max, at } = worst(values, frame.fields[name]);
    assert.ok(max < 5, `${name} differs between the engines by ${max} m at ${at}`);
    console.log(`over the mixed-layer deck the engines' ${name} differ by at most ${max.toFixed(2)} m`);
  }
});

test('frames leave the model state bit-identical whichever fields they carry', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const run = async (fields) => {
    const model = await createGpuModel(new Grid(6), { topography });
    const init = initializeState(model, {});
    for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
    for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
    model.load();
    model.ocean.initialize(model.state[3], model.state[6]);
    model.land.initialize();
    for (let n = 0; n < 10; n++) { await model.step(900); await model.beginFrame({ fields }); }
    const physics = await model.gpu.downloadPhysics();
    await model.sync();
    return { state: model.state.map((array) => Float64Array.from(array)), physics };
  };
  const plain = await run(['cloud']), heights = await run(['cloud', 'cloudTop', 'cloudBase']);
  plain.state.forEach((array, a) => assert.deepEqual(heights.state[a], array, `state array ${a}`));
  for (const [name, array] of Object.entries(plain.physics)) assert.deepEqual(heights.physics[name], array, `physics ${name}`);
});
