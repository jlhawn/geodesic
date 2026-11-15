import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces, createSigmaCore } from '../js/dynamics/sigmaCore.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { saturationHumidity, placeCloudIce } from '../js/physics/moist.module.js';
import { createRadiation, cloudOptics, liquidShare, CLOUD_OPTICS, SOLAR_CONSTANT, DAY } from '../js/physics/radiation.module.js';
import { physicsFunctions, PHYSICS_FUNCTIONS } from '../js/gpu/physics.gpu.js';
import { PREVIOUS_AEROSOL, PREVIOUS_LEVELS, PREVIOUS_RADIATION } from './helpers/previousDefaults.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};
const N = 16, DT = 1350, bl34 = sigmaInterfaces('bl34'), bl48 = sigmaInterfaces('bl48'), mesh = buildMesh(new Grid(N)), C = mesh.nCells;
const WIDTH = CLOUD_OPTICS.liquidTemperature - CLOUD_OPTICS.iceTemperature;
const close = (a, b, tolerance = 1e-9) => Math.abs(a - b) <= tolerance * Math.max(Math.abs(b), 1e-30);

/*
 * radiation.cloudPhase 'memory': a layer's resolved cloud water takes the liquid share 1 − φ of the
 * moist physics' glaciated fraction; its cumulus water and the deck water keep the ramp's α(T). The
 * optics per kilogram are linear in the liquid share, so a layer's combined water has the
 * water-weighted mean.
 */

test("cloudPhase is 'ramp' by default, 'memory' needs moist.iceDeposition and the optics of cloudOptics on both engines, and anything else is refused", async () => {
  const phase = { ...PREVIOUS_AEROSOL };
  const plain = createModel(new Grid(N), { ocean: false, radiation: phase }), ramp = createModel(new Grid(N), { ocean: false, radiation: { ...phase, cloudPhase: 'ramp' }, moist: { iceDeposition: true } });
  for (const model of [plain, ramp]) assert.doesNotThrow(() => model.radiation);
  assert.throws(() => createModel(new Grid(N), { ocean: false, radiation: { cloudPhase: 'memory' } }), /cloudPhase 'memory' needs moist.iceDeposition/);
  assert.throws(() => createModel(new Grid(N), { ocean: false, radiation: { cloudPhase: 'memory' }, moist: { iceDeposition: false } }), /needs moist.iceDeposition/);
  assert.throws(() => createModel(new Grid(N), { ocean: false, moist: false, radiation: { cloudPhase: 'memory' } }), /needs moist.iceDeposition/);
  assert.throws(() => createModel(new Grid(N), { ocean: false, radiation: { cloudPhase: 'frozen' }, moist: { iceDeposition: true } }), /cloudPhase must be 'ramp' or 'memory', not frozen/);
  assert.doesNotThrow(() => createModel(new Grid(N), { ocean: false, radiation: { cloudPhase: 'memory' }, moist: { iceDeposition: true } }));
  for (const gray of [{ cloudScattering: 150 }, { cloudAbsorption: 100 }]) assert.throws(() => createModel(new Grid(N), { ocean: false, radiation: { cloudPhase: 'memory', ...gray }, moist: { iceDeposition: true } }), /not a gray cloudScattering or cloudAbsorption/, JSON.stringify(gray));
  if (!gpuAvailable) return;
  const { mesh: base } = createModel(new Grid(N), { ocean: false, levels: bl34 });
  await assert.rejects(createGpuCore(base, { levels: bl34, physics: { cloudPhase: 'memory' } }), /cloudPhase 'memory' needs iceDeposition/);
  await assert.rejects(createGpuCore(base, { levels: bl34, physics: { cloudPhase: 'frozen', iceDeposition: true } }), /cloudPhase must be 'ramp' or 'memory', not frozen/);
  for (const gray of [{ cloudScattering: 150 }, { cloudAbsorption: 100 }]) await assert.rejects(createGpuCore(base, { levels: bl34, physics: { cloudPhase: 'memory', iceDeposition: true, ...gray } }), /not a gray cloudScattering or cloudAbsorption/, JSON.stringify(gray));
});

test('cloudOptics takes a liquid share that replaces the ramp, and the optics of a mixture of two shares are the water-weighted mean', () => {
  for (const T of [215, 255, 268, 290]) for (const continental of [false, true]) {
    const ramp = cloudOptics(T, continental), same = cloudOptics(T, continental, CLOUD_OPTICS, undefined, liquidShare(T));
    for (const key of Object.keys(ramp)) assert.equal(same[key], ramp[key], `${T} K ${key}`);
    const a = cloudOptics(T, continental, CLOUD_OPTICS, undefined, 0.8), b = cloudOptics(T, continental, CLOUD_OPTICS, undefined, 0.1), both = cloudOptics(T, continental, CLOUD_OPTICS, undefined, 0.3 * 0.8 + 0.7 * 0.1);
    for (const key of ['liquid', 'visible', 'solar', 'infrared']) assert.ok(close(both[key], 0.3 * a[key] + 0.7 * b[key], 1e-12), `${T} K ${key}`);
  }
  assert.equal(cloudOptics(255, false, CLOUD_OPTICS, undefined, 0).liquid, 0);
  assert.equal(cloudOptics(255, false, CLOUD_OPTICS, undefined, 1).liquid, 1);
});

const SIDES = ['absorbedSolar', 'reflectedSolar', 'outgoingLongwave'];

/*
 * One cloudy layer of column 0 at temperature T under a 0.5 sun, with `water` kg/kg of resolved cloud
 * and the glaciated fraction `phi` stored, with cumulus of cover and water `cumulus`; the model's
 * radiation columns' sunlight, longwave and heating.
 */
function layerRun({ memory = false, iceDeposition = memory, shortwave = 'sequential', T = 255, phi = 0.2, water = 1e-4, cumulus = null, limits = null, share = null }) {
  const radiation = { ...PREVIOUS_AEROSOL, shortwave, mixedLayerDeck: false, ...(memory ? { cloudPhase: 'memory' } : {}), ...(limits ? { iceTemperature: limits[0], liquidTemperature: limits[1] } : {}) };
  const moist = { ...(iceDeposition ? { iceDeposition: true } : {}), ...(limits ? { iceTemperature: CLOUD_OPTICS.iceTemperature, liquidTemperature: CLOUD_OPTICS.liquidTemperature } : {}) };
  const model = createModel(new Grid(N), { ocean: false, radiation, moist });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  const { core, radiation: rad, mesh: m } = model, { K, g, dSigma, sigmaMid, exnerLayer } = core.diagnostics, c = m.nCells, [pi, theta] = model.state;
  const i = 0, k = Math.floor(0.6 * K), idx = k * c + i, dry = new Float64Array(K * c), qc = new Float64Array(K * c), q = new Float64Array(K * c);
  core.diagnoseColumn(i, pi, theta, dry, dry);
  theta[idx] = T / exnerLayer[idx];
  core.diagnoseColumn(i, pi, theta, dry, dry);
  q[idx] = saturationHumidity(T, pi[i] * sigmaMid[k]);
  qc[idx] = water;
  if (model.moist.cloudIceShare) model.moist.cloudIceShare[idx] = phi;
  if (cumulus) { model.moist.cumulusCover[idx] = cumulus.cover; model.moist.cumulusWater[idx] = cumulus.water; }
  const mass = pi[i] * dSigma[k] / g, mu = 0.5, beam = SOLAR_CONSTANT * mu;
  rad.column(i, pi[i], theta, 290, 5, undefined, { beam, mu }, null, q, qc, 0.1, 0.1);
  const out = Object.fromEntries(SIDES.map((name) => [name, rad.budget[name]]));
  out.longwave = Float64Array.from({ length: K }, (_, l) => (rad.longwave ? rad.longwave[l * c + i] : 0));
  out.mass = mass;
  return out;
}

const sameRun = (a, b, tolerance = 1e-9, label = '') => {
  for (const name of SIDES) assert.ok(close(a[name], b[name], tolerance), `${label} ${name}: ${a[name]} against ${b[name]}`);
  for (let k = 0; k < a.longwave.length; k++) assert.ok(Math.abs(a.longwave[k] - b.longwave[k]) <= tolerance * 100, `${label} longwave ${k}: ${a.longwave[k]} against ${b.longwave[k]}`);
};
const tuned = (share, T) => [T - share * WIDTH, T - share * WIDTH + WIDTH];

for (const shortwave of ['sequential', 'layered']) {
  test(`${shortwave} shortwave: a layer at 255 K with φ = 0.2 and only resolved water has the optics of liquid share 0.8, not the ramp's ${liquidShare(255).toFixed(3)}`, () => {
    const memory = layerRun({ memory: true, shortwave, phi: 0.2 }), expected = layerRun({ shortwave, limits: tuned(0.8, 255) }), ramp = layerRun({ shortwave });
    sameRun(memory, expected, 1e-9, 'memory against liquid share 0.8');
    const apart = Math.abs(memory.reflectedSolar - ramp.reflectedSolar) / ramp.reflectedSolar, hot = Math.abs(memory.outgoingLongwave - ramp.outgoingLongwave) / ramp.outgoingLongwave;
    console.log(`${shortwave}: memory against ramp at 255 K, φ 0.2: reflected sunlight differs by ${apart.toExponential(2)}, outgoing longwave by ${hot.toExponential(2)}`);
    assert.ok(apart > 1e-3 && hot > 1e-5, 'the memory moves the radiation');
    for (const [phi, share] of [[0, 1], [1, 0], [0.5, 0.5]]) sameRun(layerRun({ memory: true, shortwave, phi }), layerRun({ shortwave, limits: tuned(share, 255) }), 1e-9, `φ ${phi}`);
  });

  test(`${shortwave} shortwave: with cumulus water in the layer the optics are the water-weighted mean of liquid share 1 − φ for the resolved water and the ramp's for the cumulus`, () => {
    const water = 1e-4, cumulus = { cover: 0.3, water: 2e-4 }, T = 255, phi = 0.2;
    const memory = layerRun({ memory: true, shortwave, T, phi, water, cumulus }), ramp = layerRun({ shortwave, T, water, cumulus });
    const { mass } = memory, resolved = water * mass, extra = cumulus.cover * cumulus.water * mass, share = (resolved * (1 - phi) + extra * liquidShare(T)) / (resolved + extra);
    sameRun(memory, layerRun({ shortwave, T, water, cumulus, limits: tuned(share, T) }), 1e-9, 'memory against the mean share');
    assert.ok(Math.abs(share - (1 - phi)) > 0.01 && Math.abs(share - liquidShare(T)) > 0.01, `the mean share ${share} lies between 1 − φ and the ramp`);
    assert.ok(Math.abs(memory.reflectedSolar - ramp.reflectedSolar) > 1e-4 * ramp.reflectedSolar);
  });

  test(`${shortwave} shortwave: where moist ignores the stored share, at or below iceTemperature and at or above liquidTemperature, 'memory' takes the ramp`, () => {
    for (const [T, phi] of [[220, 0], [220, 0.5], [235, 0], [285, 1], [285, 0.5], [274, 1]]) {
      const cumulus = { cover: 0.3, water: 2e-4 };
      sameRun(layerRun({ memory: true, shortwave, T, phi }), layerRun({ shortwave, T }), 1e-11, `${T} K, φ ${phi}`);
      sameRun(layerRun({ memory: true, shortwave, T, phi, cumulus }), layerRun({ shortwave, T, cumulus }), 1e-11, `${T} K, φ ${phi}, with cumulus`);
    }
    const inside = layerRun({ memory: true, shortwave, T: 236, phi: 0 }), ramp = layerRun({ shortwave, T: 236 });
    assert.ok(Math.abs(inside.reflectedSolar - ramp.reflectedSolar) > 1e-4 * ramp.reflectedSolar, 'just inside the ramp the stored share counts');
  });

  test(`${shortwave} shortwave: under 'ramp' the stored share changes nothing, and under 'memory' a share that is the ramp's 1 − α(T) gives the ramp's radiation`, () => {
    const alpha = liquidShare(255);
    const plain = layerRun({ shortwave, iceDeposition: false }), stored = layerRun({ shortwave, iceDeposition: true, phi: 0.7 });
    sameRun(stored, plain, 0, 'ramp with a stored share');
    sameRun(layerRun({ memory: true, shortwave, phi: 1 - alpha }), plain, 1e-11, 'memory at the ramp share');
    const cumulus = { cover: 0.3, water: 2e-4 };
    sameRun(layerRun({ memory: true, shortwave, phi: 1 - alpha, cumulus }), layerRun({ shortwave, cumulus }), 1e-11, 'memory at the ramp share with cumulus');
  });
}

/*
 * The mixed-layer deck of mixedLayerColumn in physics.test.mjs: a 289 K mixed layer under a 298 K
 * free troposphere over a 292 K sea at the sun's zenith, sinking at 0.4 of the column's rate,
 * with `resolved` kg/kg of cloud in the layer nearest 255 K.
 */
function deckRun({ phase, share, resolved = 0, limits = null, deckCloud = 'sigma' }) {
  const core = createSigmaCore(mesh, { levels: PREVIOUS_LEVELS }), { K, g, geopotential, piSigmaDot, exnerLayer } = core.diagnostics;
  const options = { ...PREVIOUS_RADIATION, mixedLayerDeck: true, subsidenceMemory: 1e-9, prognosticHeight: false, gateMemory: 0, deckCloud, cloudPhase: phase, ...(limits ? { iceTemperature: limits[0], liquidTemperature: limits[1] } : {}) };
  const radiation = createRadiation(mesh, core, options);
  if (phase === 'memory') radiation.useCloudIce(new Float64Array(K * C).fill(share));
  radiation.setTime(0);
  let noon = 0;
  for (let i = 0; i < C; i++) if (radiation.insolation(i) > radiation.insolation(noon)) noon = i;
  const pi = new Float64Array(C).fill(1e5), theta = new Float64Array(K * C), q = new Float64Array(K * C), qc = new Float64Array(K * C);
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) {
    const sigma = core.sigmaMid[k], mixed = sigma > 0.9;
    theta[k * C + i] = mixed ? 289 : 298 + 60 * (0.9 - sigma);
    q[k * C + i] = mixed ? 9e-3 : 1.5e-3 * Math.min(1, sigma / 0.8) ** 3;
  }
  core.diagnose(pi, theta, q, qc);
  for (let k = 1; k < K; k++) piSigmaDot.fill(0.4 * (1 - core.levels[k]), k * C, (k + 1) * C);
  let cloudy = 0;
  for (let k = 0; k < K; k++) if (core.sigmaMid[k] < 0.8 && Math.abs(theta[k * C + noon] * exnerLayer[k * C + noon] - 255) < Math.abs(theta[cloudy * C + noon] * exnerLayer[cloudy * C + noon] - 255)) cloudy = k;
  if (resolved > 0) qc[cloudy * C + noon] = resolved;
  core.diagnose(pi, theta, q, qc);
  const depth = geopotential[(K - 3) * C + noon] / g + 100 - geopotential[(K - 1) * C + noon] / g, bottom = (K - 1) * C + noon;
  const flux = radiation.column(noon, 1e5, theta, 292, 5, radiation.opticalDepth(mesh.latCell[noon]), undefined, q[bottom], q, qc, 0.07, 0.06, 1, 1.5e-3, 1, depth, 900);
  const out = Object.fromEntries([...SIDES, 'mlmCover', 'mlmWater'].map((name) => [name, radiation.budget[name]]));
  return { ...out, flux, layers: Array.from(radiation.layerFlux), T: theta[cloudy * C + noon] * exnerLayer[cloudy * C + noon] };
}

for (const deckCloud of ['sigma', 'top']) {
  test(`with the mixed-layer deck on (deckCloud '${deckCloud}') the deck's water takes the ramp and the resolved water above it takes 1 − φ`, () => {
    const bare = deckRun({ phase: 'memory', share: 0.3, deckCloud }), plain = deckRun({ phase: 'ramp', deckCloud });
    assert.ok(bare.mlmCover > 0 && bare.mlmWater > 0.03, `the deck forms: cover ${bare.mlmCover}, water ${bare.mlmWater}`);
    assert.deepEqual(bare, plain);
    const resolved = 1e-4, probe = deckRun({ phase: 'ramp', resolved, deckCloud }), phi = 0.2;
    const memory = deckRun({ phase: 'memory', share: phi, resolved, deckCloud }), tunedRamp = deckRun({ phase: 'ramp', resolved, deckCloud, limits: tuned(1 - phi, probe.T) });
    assert.ok(memory.mlmCover > 0 && tunedRamp.mlmCover === memory.mlmCover && tunedRamp.mlmWater === memory.mlmWater, 'the deck is the same under the tuned ramp');
    for (const name of [...SIDES, 'flux']) assert.ok(close(memory[name], tunedRamp[name], 1e-9), `${name}: ${memory[name]} against ${tunedRamp[name]}`);
    memory.layers.forEach((value, k) => assert.ok(Math.abs(value - tunedRamp.layers[k]) <= 1e-9 * Math.max(1, Math.abs(value)), `layer ${k}: ${value} against ${tunedRamp.layers[k]}`));
    assert.ok(Math.abs(memory.reflectedSolar - probe.reflectedSolar) > 1e-4 * probe.reflectedSolar, 'the memory moves the radiation');
    console.log(`deck '${deckCloud}': cover ${memory.mlmCover}, water ${(1000 * memory.mlmWater).toFixed(1)} g/m², cloud at ${probe.T.toFixed(1)} K; reflected sunlight ${memory.reflectedSolar.toFixed(3)} against ${probe.reflectedSolar.toFixed(3)} under 'ramp'`);
  });
}

async function cores(levels) {
  const { mesh: base } = createModel(new Grid(N), { ocean: false, levels });
  const make = (physics) => createGpuCore(base, { levels, physics });
  return { base, make };
}

test("with cloudPhase 'ramp' the GPU kernels, the shared functions and the PH layout are those of the default; 'memory' changes the physics, radiation and layered kernels, and nothing else", { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  assert.equal(physicsFunctions({ cloudPhase: 'ramp' }), PHYSICS_FUNCTIONS);
  assert.equal(physicsFunctions({}), PHYSICS_FUNCTIONS);
  assert.notEqual(physicsFunctions({ cloudPhase: 'memory' }), PHYSICS_FUNCTIONS);
  for (const levels of [bl34, bl48]) {
    const { make } = await cores(levels);
    const plain = await make({}), ramp = await make({ cloudPhase: 'ramp' }), ice = await make({ iceDeposition: true }), iceRamp = await make({ iceDeposition: true, cloudPhase: 'ramp' }), memory = await make({ iceDeposition: true, cloudPhase: 'memory' });
    assert.deepEqual(ramp.kernelBodies, plain.kernelBodies);
    assert.deepEqual(ramp.layout.PH, plain.layout.PH);
    assert.deepEqual(iceRamp.kernelBodies, ice.kernelBodies);
    assert.deepEqual(iceRamp.layout.PH, ice.layout.PH);
    for (const bodies of [plain.kernelBodies, ice.kernelBodies]) for (const [name, body] of Object.entries(bodies)) for (const word of ['memoryShare', 'cloudOpticsAt', 'rampOptics']) assert.ok(!body.includes(word), `${name} has ${word} under 'ramp'`);
    for (const [name, body] of Object.entries(plain.kernelBodies)) assert.ok(!body.includes('PH_CICE') || name.startsWith('adjust'), `${name} reads PH_CICE with the option off`);
    assert.deepEqual(Object.keys(memory.kernelBodies), Object.keys(ice.kernelBodies));
    const changed = Object.keys(memory.kernelBodies).filter((name) => memory.kernelBodies[name] !== ice.kernelBodies[name]);
    assert.deepEqual(changed.sort(), ['layered', 'physics', 'radiation']);
    assert.deepEqual(memory.layout.PH, ice.layout.PH);
    for (const name of changed) assert.ok(memory.kernelBodies[name].includes('PH_CICE'), `${name} reads the memory`);
  }
});

const CLOUD_LAYERS = [[0.8, 0.95, 4e-4, (i) => i % 4 === 0], [0.45, 0.6, 2e-4, (i) => i % 7 < 3], [0.2, 0.3, 6e-5, (i) => i % 5 === 1], [0.55, 0.9, 1e-4, (i) => i % 11 === 2], [0.3, 0.45, 1e-4, (i) => i % 3 === 0]];
const UPPER_LAYERS = CLOUD_LAYERS.filter(([from, to]) => to <= 0.6);

function cloudier(model, layers = CLOUD_LAYERS) {
  const { K, sigmaMid, exnerLayer } = model.core.diagnostics, [pi, theta, , , q, qc] = model.state;
  for (const [from, to, water, test] of layers) for (let k = 0; k < K; k++) {
    if (sigmaMid[k] < from || sigmaMid[k] > to) continue;
    for (let i = 0; i < C; i++) {
      if (!test(i)) continue;
      const idx = k * C + i;
      q[idx] = saturationHumidity(theta[idx] * exnerLayer[idx], pi[i] * sigmaMid[k]);
      qc[idx] = water;
    }
  }
}

function stats(cpu, gpu) {
  let sumSq = 0, sumRef = 0;
  const errors = new Float64Array(cpu.length);
  for (let x = 0; x < cpu.length; x++) { const d = Math.abs(cpu[x] - gpu[x]); errors[x] = d; sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
  const scale = Math.sqrt(sumRef / cpu.length);
  if (!Number.isFinite(sumSq) || !Number.isFinite(scale)) return { rmsRel: Infinity, p99Rel: Infinity };
  return { rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)), p99Rel: errors.sort()[Math.floor(0.99 * (cpu.length - 1))] / Math.max(scale, 1e-300) };
}

const HELD = [['absorbed', 'RADABS'], ['back', 'RADBACK'], ['emission', 'RADEMIT'], ['outgoing', 'RADOLR'], ['escape', 'RADT'], ['down', 'RADSWDN'], ['atmosphere', 'RADATM'], ['top', 'RADTOA'], ['reflected', 'RADREFL']];

/*
 * The CPU model and a GPU core of the same options, one full radiation call from model time `time`,
 * with the glaciated fraction drawn at random in the cloudy layers so that it differs from 1 − α(T).
 */
async function pair(phase, levels, shortwave, { time = 200 * DAY, deck = false, cumulus = false, steps = 1, every = 4 } = {}) {
  const radiation = { stratus: true, ...PREVIOUS_AEROSOL, shortwave, radiationEvery: every, mixedLayerDeck: deck, cloudPhase: phase }, moist = { iceDeposition: true, plumeEntrainmentLaw: 'gregory' };
  const model = createModel(new Grid(N), { ocean: false, levels, radiation, moist });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  if (deck) {
    const { sigmaMid } = model.core.diagnostics, [pi, theta, , surfaceT, q, qc] = model.state;
    for (let k = 0; k < model.core.diagnostics.K; k++) for (let i = 0; i < C; i++) {
      const s = sigmaMid[k], mixed = s > 0.9;
      theta[k * C + i] = mixed ? 289 : Math.max(298 + 60 * (0.9 - s), 215 / s ** 0.2857);
      q[k * C + i] = mixed ? 9e-3 : 1.5e-3 * Math.min(1, s / 0.8) ** 3;
      qc[k * C + i] = 0;
    }
    pi.fill(1e5); surfaceT.fill(292);
    model.core.diagnose(pi, theta, q, qc);
    cloudier(model, UPPER_LAYERS);
  } else cloudier(model);
  placeCloudIce(model, null);
  let seed = 4242;
  const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; }, share = model.moist.cloudIceShare, { K } = model.core.diagnostics;
  for (let x = 0; x < K * C; x++) if (model.state[5][x] > 0) share[x] = Math.fround(random());
  for (const a of model.state) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  model.time = time;
  const mean = Float64Array.from({ length: K }, (_, k) => { let s = 0, a = 0; for (let i = 0; i < C; i++) { s += model.mesh.areaCell[i] * model.state[1][k * C + i]; a += model.mesh.areaCell[i]; } return s / a; });
  const gpu = await createGpuCore(model.mesh, { levels, nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, divergenceDamping: model.core.divergenceDamping, referenceTheta: mean, physics: { ...radiation, ...moist } });
  gpu.upload(model.state);
  const upload = { cloudIceShare: Float32Array.from(share) }, mixture = { layers: 0, both: 0 };
  if (cumulus) {
    const layers = (gpu.layout.PH.CUWATER - gpu.layout.PH.CUCOVER) / C, k0 = K - layers, { cumulusCover, cumulusWater } = model.moist, exner = model.core.diagnostics.exnerLayer;
    for (let x = k0 * C; x < K * C; x++) {
      if (random() < 0.7) continue;
      cumulusCover[x] = Math.fround(0.1 + 0.4 * random());
      cumulusWater[x] = Math.fround(2e-5 + 1.8e-4 * random());
      mixture.layers++;
      const T = model.state[1][x] * exner[x], ramped = liquidShare(T);
      if (model.state[5][x] > 1e-7 && Math.abs(share[x] - (1 - ramped)) > 0.05 && ramped > 0 && ramped < 1) mixture.both++;
    }
    upload.cumulusCover = Float32Array.from(cumulusCover.subarray(k0 * C));
    upload.cumulusWater = Float32Array.from(cumulusWater.subarray(k0 * C));
  }
  gpu.uploadPhysics(upload);
  for (let n = 0; n < steps; n++) { const t = model.time; model.step(DT); await gpu.stepModel(DT, t); }
  const physics = await gpu.downloadPhysics(), decks = { cpu: 0, gpu: 0 };
  for (let i = 0; i < C; i++) { if (model.radiation.mlmCover[i] > 0) decks.cpu++; if (physics.MLMCOVER[i] > 0) decks.gpu++; }
  const run = { model, physics, K, mixture, decks };
  for (const buffer of Object.values(gpu.buffers)) buffer.destroy();
  return run;
}

for (const [label, levels, options] of [['bl34', bl34, { cumulus: true }], ['bl48', bl48, { cumulus: true }], ['bl48 with the mixed-layer deck', bl48, { deck: true, steps: 3, every: 2 }]]) {
  for (const shortwave of ['sequential', 'layered']) {
    test(`at N=16 on ${label} with ${shortwave} shortwave the GPU's radiation under cloudPhase 'memory' matches the CPU's by layer and at the surface and the top, and differs from 'ramp' where φ ≠ 1 − α`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
      const ramp = await pair('ramp', levels, shortwave, options), memory = await pair('memory', levels, shortwave, options), { K } = memory, { cp, g, dSigma } = memory.model.core.diagnostics, pi = memory.model.state[0];
      const rate = (k, i) => 86400 * g / (cp * pi[i] * dSigma[k]);
      const lit = (run, i) => run.model.radiation.held.mu[i] > 0;
      const columns = { total: 0, sunlit: 0, parted: 0 };
      let worstSw = 0, worstLw = 0, changeSw = 0, changeLw = 0, worstChangeSw = 0, worstChangeLw = 0, worstTop = 0, worstSurface = 0, where = null, scaledHeat = 0, worstFluxSw = 0, worstFluxLw = 0, baseSwMax = 0, baseLwMax = 0, moved = 0;
      for (let i = 0; i < C; i++) {
        columns.total++;
        let sunlit = lit(memory, i), apartSw = 0, apartLw = 0;
        if (sunlit) columns.sunlit++;
        for (let k = 0; k < K; k++) {
          const x = k * C + i, r = rate(k, i);
          apartSw = Math.max(apartSw, Math.abs(ramp.model.radiation.held.shortwave[x] - ramp.physics.RADSW[x]) * r);
          apartLw = Math.max(apartLw, Math.abs(ramp.model.radiation.longwave[x] - ramp.physics.LWH[x]) * r);
        }
        baseSwMax = Math.max(baseSwMax, apartSw); baseLwMax = Math.max(baseLwMax, apartLw);
        if (apartSw > 1e-3 || apartLw > 1e-3) { columns.parted++; continue; }
        for (let k = 0; k < K; k++) {
          const x = k * C + i, r = rate(k, i);
          const swGap = Math.abs(memory.model.radiation.held.shortwave[x] - memory.physics.RADSW[x]) * r, lwGap = Math.abs(memory.model.radiation.longwave[x] - memory.physics.LWH[x]) * r;
          const swCpu = (memory.model.radiation.held.shortwave[x] - ramp.model.radiation.held.shortwave[x]) * r, swGpu = (memory.physics.RADSW[x] - ramp.physics.RADSW[x]) * r;
          const lwCpu = (memory.model.radiation.longwave[x] - ramp.model.radiation.longwave[x]) * r, lwGpu = (memory.physics.LWH[x] - ramp.physics.LWH[x]) * r;
          worstSw = Math.max(worstSw, swGap); worstLw = Math.max(worstLw, lwGap);
          const heatSw = Math.abs(memory.model.radiation.held.shortwave[x]) * r, heatLw = Math.abs(memory.model.radiation.longwave[x]) * r;
          const baseSw = Math.abs(ramp.model.radiation.held.shortwave[x] - ramp.physics.RADSW[x]) * r, baseLw = Math.abs(ramp.model.radiation.longwave[x] - ramp.physics.LWH[x]) * r;
          scaledHeat = Math.max(scaledHeat, swGap / (1e-3 + 3e-4 * heatSw + baseSw), lwGap / (1e-3 + 3e-4 * heatLw + baseLw));
          changeSw = Math.max(changeSw, Math.abs(swCpu)); changeLw = Math.max(changeLw, Math.abs(lwCpu));
          worstChangeSw = Math.max(worstChangeSw, Math.abs(swCpu - swGpu)); worstChangeLw = Math.max(worstChangeLw, Math.abs(lwCpu - lwGpu));
          if (Math.abs(swCpu - swGpu) > worstChangeSw - 1e-12 && Math.abs(swCpu - swGpu) >= worstChangeSw) where = { i, k, phi: memory.model.moist.cloudIceShare[x], qc: memory.model.state[5][x], T: memory.model.state[1][x] * memory.model.core.diagnostics.exnerLayer[x], swCpu, swGpu, base: (ramp.model.radiation.held.shortwave[x] - ramp.physics.RADSW[x]) * r };
          worstFluxSw = Math.max(worstFluxSw, Math.abs(swCpu - swGpu) / r); worstFluxLw = Math.max(worstFluxLw, Math.abs(lwCpu - lwGpu) / r);
          if (Math.abs(swCpu) > 1e-3 || Math.abs(lwCpu) > 1e-3) moved++;
        }
        const beam = Math.max(1, memory.model.radiation.held.mu[i] * SOLAR_CONSTANT);
        worstTop = Math.max(worstTop, Math.abs(memory.model.radiation.held.top[i] - memory.physics.RADTOA[i]) / beam);
        worstSurface = Math.max(worstSurface, Math.abs(memory.model.radiation.held.absorbed[i] - memory.physics.RADABS[i]) / beam);
      }
      const flux = {};
      for (const [cpuName, gpuName] of HELD) {
        for (const [run, label2] of [[ramp, 'ramp'], [memory, 'memory']]) for (const [who, values] of [['CPU', run.model.radiation.held[cpuName]], ['GPU', run.physics[gpuName].subarray(0, C)]]) { const bad = values.findIndex((v) => !Number.isFinite(v)); assert.ok(bad < 0, `${label2} ${who} ${cpuName} is ${values[bad]} at column ${bad}`); }
        const a = stats(memory.model.radiation.held[cpuName], memory.physics[gpuName].subarray(0, C)), b = stats(ramp.model.radiation.held[cpuName], ramp.physics[gpuName].subarray(0, C));
        flux[cpuName] = [a.p99Rel, b.p99Rel];
        assert.ok(Number.isFinite(a.p99Rel) && a.p99Rel < 3 * b.p99Rel + 2e-4, `${cpuName}: 99th percentile ${a.p99Rel} against ${b.p99Rel} under 'ramp'`);
      }
      console.log(`${label} ${shortwave}: ${columns.sunlit} of ${columns.total} columns sunlit, ${columns.parted} left out where a discrete decision parts the engines under 'ramp' too (baseline shortwave ${baseSwMax.toExponential(1)}, longwave ${baseLwMax.toExponential(1)} K/day); under 'memory' the engines' layer heating differs by shortwave ${worstSw.toExponential(1)}, longwave ${worstLw.toExponential(1)} K/day; scaled heating gap ${scaledHeat.toFixed(2)}, the change's flux gaps ${worstFluxSw.toExponential(1)} and ${worstFluxLw.toExponential(1)} W/m2; the memory changes a layer by at most ${changeSw.toFixed(3)} (shortwave) and ${changeLw.toFixed(3)} (longwave) K/day on ${moved} layers and the engines differ on that change by ${worstChangeSw.toExponential(1)} and ${worstChangeLw.toExponential(1)}; surface and top absorbed sunlight differ by ${worstSurface.toExponential(1)} and ${worstTop.toExponential(1)} of the beam; flux 99th percentiles (memory, ramp): ${Object.entries(flux).map(([name, [a, b]]) => `${name} ${a.toExponential(1)}, ${b.toExponential(1)}`).join('; ')}`);
      if (where) console.log(`worst change gap at ${JSON.stringify(where)}`);
      console.log(`${label} ${shortwave}: ${memory.mixture.layers} cumulus layers, ${memory.mixture.both} of them holding resolved water with φ off the ramp; mixed-layer deck on ${memory.decks.cpu} columns on the CPU and ${memory.decks.gpu} on the GPU; largest change gaps ${(worstChangeSw / changeSw).toExponential(1)} and ${(worstChangeLw / changeLw).toExponential(1)} of the changes`);
      if (options.cumulus) assert.ok(memory.mixture.layers > 50 && memory.mixture.both > 5, `${memory.mixture.layers} cumulus layers, ${memory.mixture.both} mixed`);
      if (options.deck) assert.ok(memory.decks.cpu > 0.3 * C && memory.decks.gpu > 0.3 * C, `decks on ${memory.decks.cpu} and ${memory.decks.gpu} columns`);
      else assert.equal(memory.decks.cpu + memory.decks.gpu, 0);
      assert.ok(worstChangeSw < 0.05 * changeSw && worstChangeLw < 0.05 * changeLw, `the change differs by ${worstChangeSw} of ${changeSw} and ${worstChangeLw} of ${changeLw} K/day`);
      assert.ok(columns.parted < 0.1 * columns.total && moved > 20, `${columns.parted} columns parted, ${moved} layers moved`);
      assert.ok(changeSw > 0.02 && changeLw > 0.02, `the memory moves the CPU's heating by ${changeSw} and ${changeLw} K/day`);
      assert.ok(scaledHeat < 1, `heating differs by ${worstSw} and ${worstLw} K/day, ${scaledHeat} of 1e-3 + 3e-4 of the heating`);
      const allowed = options.deck ? [1e-2, 1e-2] : [1e-3, 3e-3];
      assert.ok(worstFluxSw < allowed[0] && worstFluxLw < allowed[1], `the change in a layer's flux differs by ${worstFluxSw} and ${worstFluxLw} W/m2`);
      const beamGap = options.deck ? 1e-4 : 1e-5;
      assert.ok(worstTop < beamGap && worstSurface < beamGap, `top ${worstTop}, surface ${worstSurface} of the beam`);
    });
  }
}
