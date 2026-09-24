import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { levelFields, dewPoint, wetBulb, miseryIndex } from '../js/levels.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { FIELDS } from '../js/frames.module.js';

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
    const checks = {
      temperature: [reference.temperature, 2e-3], height: [reference.height, 0.1], humidity: [reference.humidity, 5e-5], speed: [reference.speed, 5e-4], wind: [reference.vector, 5e-4],
      dewPoint: [comfort(dewPoint), 5e-3], wetBulb: [comfort(wetBulb), 5e-3], misery: [comfort(miseryIndex), 5e-3],
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
  const cloud = Float64Array.from({ length: C }, (_, i) => model.moist.columnWater(pi, qc, i));
  const sea = (source) => Float64Array.from({ length: C }, (_, i) => (model.geography.land[i] ? NaN : source[i]));
  const currents = cellVector(mesh, ocean.u1);
  for (let i = 0; i < C; i++) if (model.geography.land[i]) currents.fill(0, 3 * i, 3 * i + 3);
  const current = sea(Float64Array.from({ length: C }, (_, i) => Math.hypot(currents[3 * i], currents[3 * i + 1], currents[3 * i + 2])));
  const checks = {
    ps: [pi, 0.1], mslp: [mslp, 0.5], water: [water, 1e-3], cloud: [cloud, 1e-6], rain: [physics.RAIN.subarray(0, C), 1e-6], ice: [ice, 1e-6],
    albedo: [physics.ADIF.subarray(0, C), 1e-6], shortwave: [physics.SWDN.subarray(0, C), 1e-3], longwave: [physics.OLR.subarray(0, C), 1e-3],
    soil: [physics.SOIL.subarray(0, C), 1e-4], snow: [physics.SNOW.subarray(0, C), 1e-4],
    sst: [sea(ocean.T1), 1e-3], sss: [sea(ocean.S1), 1e-4], layerDepth: [sea(ocean.h1), 1e-3], thermocline: [ocean.thermoclineDepth, 1e-2], ssh: [sea(ocean.eta), 1e-5],
    current: [current, 1e-5], currents: [currents, 1e-5],
  };
  for (const [name, [values, tolerance]] of Object.entries(checks)) {
    const { max, at } = worst(values, f[name], { skipNaN: true });
    assert.ok(max < tolerance, `${name}: ${max} at ${at} (${values[at]} against ${f[name][at]})`);
  }

  const land = physics.LAND, d = first.diagnostics;
  let area = 0, mass = 0, ts = 0, piMin = Infinity, piMax = -Infinity, maxWind = 0, w = 0, c = 0, rain = 0, iceArea = 0, iceVolume = 0, olr = 0, absorbed = 0, landArea = 0, landT = 0, soil = 0;
  for (let i = 0; i < C; i++) {
    const a = mesh.areaCell[i];
    area += a; mass += a * pi[i]; ts += a * surfaceT[i]; piMin = Math.min(piMin, pi[i]); piMax = Math.max(piMax, pi[i]);
    w += a * water[i]; c += a * cloud[i]; rain += a * physics.RAIN[i]; olr += a * physics.OLR[i]; absorbed += a * physics.ABS[i];
    if (ice[i] > 0) { iceArea += a; iceVolume += a * ice[i]; }
    if (land[i] > 0.5) { landArea += a; landT += a * surfaceT[i]; soil += a * physics.SOIL[i]; }
  }
  for (const x of u) maxWind = Math.max(maxWind, Math.abs(x));
  const close = (name, expected, got, relative = 1e-5) => assert.ok(Math.abs(got - expected) <= relative * Math.abs(expected) + 1e-12, `${name}: ${got} against ${expected}`);
  close('mass', mass / area, d.mass); close('meanSurfaceT', ts / area, d.meanSurfaceT);
  close('piMin', piMin, d.piMin); close('piMax', piMax, d.piMax); close('maxWind', maxWind, d.maxWind);
  close('columnWater', w / area, d.columnWater); close('columnCloud', c / area, d.columnCloud, 1e-4);
  close('precipitation', rain / area / model.time, d.precipitation, 1e-4);
  close('iceFraction', iceArea / area, d.iceFraction); close('iceThickness', iceArea > 0 ? iceVolume / iceArea : 0, d.iceThickness);
  close('outgoingLongwave', olr / area, d.outgoingLongwave); close('absorbedSolar', absorbed / area, d.absorbedSolar);
  close('landMeanT', landT / landArea, d.landMeanT); close('soilWater', soil / landArea, d.soilWater);
  let runoff = 0;
  for (let i = 0; i < C; i++) { runoff += mesh.areaCell[i] * physics.RUNOFF[i]; assert.ok(Math.abs(model.land.runoff[i] - physics.RUNOFF[i]) <= 1e-6 * Math.abs(physics.RUNOFF[i]) + 1e-12, `runoff at ${i}`); }
  assert.ok(runoff > 0);
  close('runoff', runoff / area, d.runoff);
  let oceanArea = 0, depth = 0, ssh = 0;
  for (let i = 0; i < C; i++) if (!model.geography.land[i]) { const a = mesh.areaCell[i]; oceanArea += a; depth += a * ocean.h1[i]; ssh = Math.max(ssh, Math.abs(ocean.eta[i])); }
  close('oceanUpperDepth', depth / oceanArea, d.oceanUpperDepth); close('oceanSSH', ssh, d.oceanSSH);
  console.log(`GPU frame at N=6: Ts ${d.meanSurfaceT.toFixed(3)} K, OLR ${d.outgoingLongwave.toFixed(2)} W/m², precipitation ${(86400 * d.precipitation).toFixed(3)} mm/day, ocean h1 ${d.oceanUpperDepth.toFixed(2)} m`);

  const unsubscribed = await model.beginFrame({ level: 'surface', fields: ['sst'] });
  assert.deepEqual(Object.keys(unsubscribed.fields), ['sst']);
  assert.equal(unsubscribed.diagnostics, null);
  const paused = await model.beginFrame({ diagnostics: true });
  assert.ok(paused.diagnostics.precipitation > 0, 'a frame with no time elapsed reports the recent rain rate');
});
