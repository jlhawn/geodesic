import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';
import { HEAT_REFERENCE, SALT_REFERENCE } from '../js/gpu/layeredOcean.gpu.js';
import { columnMomentum } from './helpers/layered.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { readRanges } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};

const DT = 1800, SPINUP = 24, DEG = Math.PI / 180, LIGHT = 250;
const mesh = buildMesh(new Grid(16)), C = mesh.nCells, E = mesh.nEdges, AREA = mesh.areaCell.reduce((s, a) => s + a, 0);

/*
 * An aquaplanet's sea a little way into its life: createOcean's start on
 * Grid(16) under a zonal surface from 28 °C at the equator to just above
 * freezing at the poles, with no ice, so that each column holds thick
 * deep classes under a thermocline that rises poleward and outcrops, thin
 * classes at the outcrops and tokens above them; then `steps` steps under
 * easterly trades and westerlies, each reading back the surface it wrote,
 * so the currents, the eddy transport and the mixed layer's exchanges are
 * all under way and no heat crosses the surface. The sea publishes the
 * share of the sunlight that passes its mixed layer, and none is passed
 * to it.
 */
const CARRIED = { carryMomentum: true };
function spunUp(steps = SPINUP, options = {}) {
  const ocean = createOcean(mesh, { everySteps: 1, penetration: 'IB', ...options });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => {
    const polar = Math.max(0, Math.min(1, (Math.abs(lat) / DEG - 50) / 20));
    return Math.max(FREEZING_POINT + 0.3, 273.15 + 30 * Math.cos(lat) ** 2 - 2 - 4 * polar);
  });
  const ice = new Float64Array(C), flux = new Float64Array(C), stress = new Float64Array(E);
  for (let e = 0; e < E; e++) {
    const x = mesh.xEdge[3 * e], y = mesh.xEdge[3 * e + 1], r = Math.hypot(x, y) || 1e-12;
    stress[e] = -0.1 * Math.cos(3 * mesh.latEdge[e]) * ((-y / r) * mesh.nEdge[3 * e] + (x / r) * mesh.nEdge[3 * e + 1]);
  }
  ocean.initialize(surfaceT, ice);
  for (let n = 0; n < steps; n++) ocean.advance(surfaceT, ice, flux, stress, DT);
  return { ocean, surfaceT, ice, flux, stress, options };
}

/*
 * What one phase did to the sea, from cell-by-cell differences of two
 * snapshots of the stores (h, heat as h·(T − HEAT_REFERENCE), salt h·S)
 * and of the sunlight waiting to be deposited (J/m²), summed over every
 * layer and cell in double precision: the heat, the stores' with the
 * waiting light's, as a global-mean flux over the step (W/m²), and
 * `stored`, the stores' alone, the salt and the volume as global means
 * (psu·m, m), `offMean`, the heat and salt a column changed beyond what
 * its change of volume carries at the column's own mean temperature and
 * salinity, and `momentum`, the largest change of an edge's column
 * momentum Σ h_e u relative to its Σ h_e |u|.
 */
function books(before, after, layers, rhoCp, D) {
  let heat = 0, waiting = 0, salt = 0, volume = 0, offHeat = 0, offSalt = 0;
  for (let i = 0; i < C; i++) {
    let dh = 0, dq = 0, dw = 0, h = 0, q = 0, w = 0;
    for (let k = 0; k < layers; k++) {
      const n = k * C + i;
      dh += after.h[n] - before.h[n]; dq += after.q[n] - before.q[n]; dw += after.w[n] - before.w[n];
      h += before.h[n]; q += before.q[n]; w += before.w[n];
    }
    const a = mesh.areaCell[i];
    heat += a * dq; salt += a * dw; volume += a * dh; waiting += a * (after.pen[i] - before.pen[i]);
    if (h > 0) { offHeat += a * (dq - q / h * dh); offSalt += a * (dw - w / h * dh); }
  }
  const watts = (x) => rhoCp * x / (AREA * DT);
  const was = columnMomentum(mesh, layers, before, D), is = columnMomentum(mesh, layers, after, D);
  let momentum = 0;
  for (let e = 0; e < E; e++) momentum = Math.max(momentum, Math.abs(is.momentum[e] - was.momentum[e]) / was.scale[e]);
  return { heat: watts(heat) + waiting / (AREA * DT), stored: watts(heat), salt: salt / AREA, volume: volume / AREA, offMean: { heat: watts(offHeat), salt: offSalt / AREA }, momentum };
}

/*
 * `phases` as [name, run, prepare]: a phase's `prepare`, when it has one,
 * runs outside the books, before the snapshot the phase is measured from.
 */
async function follow(phases, snapshot, steps, layers, rhoCp, D) {
  const record = Object.fromEntries(phases.map(([name]) => [name, []]));
  let before = await snapshot();
  for (let n = 0; n < steps; n++) {
    for (const [name, run, prepare] of phases) {
      if (prepare) { await prepare(); before = await snapshot(); }
      await run();
      const after = await snapshot();
      record[name].push(books(before, after, layers, rhoCp, D));
      before = after;
    }
  }
  return record;
}

/*
 * LIGHT W/m² of absorbed sunlight over the step on every sea cell, the
 * share each passes below its mixed layer as the pending light.
 */
const passing = (transmission) => Float64Array.from(transmission, (t) => LIGHT * DT * t);
const meanOf = (perCell) => perCell.reduce((s, x, i) => s + mesh.areaCell[i] * x, 0) / (AREA * DT);

const largest = (rows, pick) => Math.max(...rows.map((row) => Math.abs(pick(row))));
const mean = (rows, pick) => rows.reduce((s, row) => s + pick(row), 0) / rows.length;
const show = (x) => x.toExponential(1);
const kept = (name) => (name === 'rescale' ? [(r) => r.offMean.heat, (r) => r.offMean.salt] : [(r) => r.heat, (r) => r.salt]);
const volume = (r) => r.volume;
function report(label, record) {
  for (const [name, rows] of Object.entries(record)) {
    const [heat, salt] = kept(name);
    console.log(`${label} ${name}: heat ${name === 'rescale' ? 'off the columns\' means' : 'changed'} by ${show(mean(rows, heat))} W/m² on average, ${show(largest(rows, heat))} at most; salt ${show(mean(rows, salt))}, ${show(largest(rows, salt))} psu·m; volume ${show(mean(rows, volume))}, ${show(largest(rows, volume))} m${name === 'mixedLayer' ? `; an edge's column momentum by ${show(largest(rows, (r) => r.momentum))} of its Σ h|u| at most` : ''}`);
  }
}

test('in double precision every phase of the CPU ocean\'s step keeps the heat and salt it holds, the rescale to the averaged free surface changes them only with its water, at each column\'s mean temperature and salinity, the mixed layer\'s exchanges keep each edge\'s column momentum, and the surface read deposits all the sunlight that passed the mixed layer', async () => {
  const { ocean, surfaceT, ice, flux, stress } = spunUp(SPINUP, CARRIED);
  ocean.setStress(stress, ice);
  const snapshot = async () => ({ h: Float64Array.from(ocean.h), q: Float64Array.from(ocean.Q, (q, n) => q - HEAT_REFERENCE * ocean.h[n]), w: Float64Array.from(ocean.W), u: Float64Array.from(ocean.u), eta: Float64Array.from(ocean.eta), pen: Float64Array.from(ocean.penetrated) });
  const light = [];
  const record = await follow([
    ['readSurface', () => ocean.readSurface(surfaceT, ice)],
    ['dynamics', () => ocean.phases.dynamics(DT)],
    ['rescale', () => ocean.phases.rescale()],
    ['eddy', () => ocean.phases.eddy(DT)],
    ['shift', () => ocean.phases.shift()],
    ['mixedLayer', () => ocean.mixedLayer(DT)],
    ['salt and surface', () => { ocean.salt(); ocean.writeSurface(surfaceT, flux, DT); }],
    ['deposit', () => ocean.readSurface(surfaceT, ice), () => { ocean.penetrated.set(passing(ocean.transmission)); light.push(meanOf(ocean.penetrated)); }],
  ], snapshot, 4, ocean.layers, ocean.rhoCp, ocean.D);
  report('CPU', record);
  for (const [name, rows] of Object.entries(record)) {
    const [heat, salt] = kept(name);
    assert.ok(largest(rows, heat) < 1e-6, `${name}: heat ${largest(rows, heat)} W/m²`);
    assert.ok(largest(rows, salt) < 1e-9, `${name}: salt ${largest(rows, salt)} psu·m`);
    if (name !== 'rescale') assert.ok(largest(rows, volume) < 1e-12, `${name}: volume ${largest(rows, volume)} m`);
  }
  assert.ok(largest(record.rescale, (r) => r.heat) > 1e-3, 'the rescale moves water, and heat with it');
  assert.ok(largest(record.mixedLayer, (r) => r.momentum) < 1e-12, `the mixed layer changes an edge's momentum by ${largest(record.mixedLayer, (r) => r.momentum)}`);
  const deposited = record.deposit.map((r, n) => r.stored - light[n]);
  console.log(`CPU deposit: ${show(mean(light, (x) => x))} W/m² of sunlight passed the mixed layer each step; the stores took it to ${show(largest(deposited, (x) => x))} W/m², and none waits after it`);
  assert.ok(Math.min(...light) > 0.5 && largest(deposited, (x) => x) < 1e-6 && ocean.penetrated.every((x) => x === 0), `the stores took ${record.deposit.map((r) => r.stored)} of ${light} W/m²`);
});

/*
 * A GPU ocean holding a CPU ocean's sea, uploaded from its h·T and h·S
 * as a file from before Qref and Wref would be, the surface it last
 * wrote in the atmosphere's buffers and its wind on the sea; `stores`
 * reads the device's h, heat and salt as stored and the light waiting
 * to be deposited, `snapshot` as books takes them, the salt as h·S.
 */
async function onGpu({ ocean: cpu, surfaceT, ice, stress, options }) {
  const model = await createGpuModel(mesh, { ocean: { everySteps: 1, penetration: 'IB', ...options } });
  const { gpu, oceanEngine: ocean } = model, layers = ocean.layers, { OS } = ocean.layout, PH = gpu.layout.PH;
  ocean.upload({ ...cpu.serialize(), Q: cpu.Q, W: cpu.W, previousT0: cpu.T0, previousIceVolume: ice, fresh: new Float64Array(C), buoyancy: new Float64Array(C), flux: new Float64Array(C), capacity: cpu.capacity }, surfaceT, ice);
  gpu.uploadSurfaceTemperature(surfaceT);
  gpu.device.queue.writeBuffer(gpu.buffers.S, 4 * gpu.layout.S.ICE, new Float32Array(C));
  ocean.setStress(stress, ice);
  const stores = async () => {
    await gpu.device.queue.onSubmittedWorkDone();
    const [[h, q, w], [pen]] = await Promise.all([readRanges(gpu.device, ocean.buffers.S, [{ offset: OS.OH, length: layers * C }, { offset: OS.OQ, length: layers * C }, { offset: OS.OW, length: layers * C }]), readRanges(gpu.device, gpu.buffers.PH, [{ offset: PH.PEN, length: C }])]);
    return [h, q, w, pen];
  };
  const snapshot = async () => {
    const [h, q, w, pen] = await stores();
    const [[u], [eta]] = await Promise.all([readRanges(gpu.device, ocean.buffers.S, [{ offset: OS.OU, length: layers * E }]), readRanges(gpu.device, ocean.buffers.OD, [{ offset: ocean.layout.OD.ETA, length: C }])]);
    return { h: Float64Array.from(h), q: Float64Array.from(q), w: Float64Array.from(w, (x, n) => x + SALT_REFERENCE * h[n]), u: Float64Array.from(u), eta: Float64Array.from(eta), pen: Float64Array.from(pen) };
  };
  const light = [];
  const pass = async () => {
    const [[transmission], [pen]] = await Promise.all([readRanges(gpu.device, gpu.buffers.PH, [{ offset: PH.TRANS, length: C }]), readRanges(gpu.device, gpu.buffers.PH, [{ offset: PH.PEN, length: C }])]);
    const added = Float32Array.from(passing(transmission));
    light.push(meanOf(added));
    gpu.device.queue.writeBuffer(gpu.buffers.PH, 4 * PH.PEN, Float32Array.from(added, (x, i) => x + pen[i]));
  };
  const phases = [
    ['readSurface', () => ocean.readSurfaceFromAtmosphere()],
    ['dynamics', () => ocean.phases.dynamics(DT)],
    ['rescale', () => ocean.phases.rescale()],
    ['eddy', () => ocean.phases.eddy(DT)],
    ['shift', () => ocean.phases.shift()],
    ['mixedLayer', () => ocean.mixedLayer(DT)],
    ['salt and surface', () => { ocean.salt(DT); ocean.writeSurface(DT); }],
    ['deposit', () => ocean.readSurfaceFromAtmosphere(), pass],
  ];
  return { model, ocean, stores, light, follow: (steps) => follow(phases, snapshot, steps, layers, ocean.options.density * ocean.options.specificHeat, ocean.D) };
}

function within(label, record, floors) {
  for (const [name, [meanHeat, stepHeat, meanSalt, stepSalt, meanVolume, stepVolume]] of Object.entries(floors)) {
    const rows = record[name], [heat, salt] = kept(name);
    assert.ok(Math.abs(mean(rows, heat)) < meanHeat && largest(rows, heat) < stepHeat, `${label} ${name}: heat ${mean(rows, heat)} W/m² on average, ${largest(rows, heat)} at most`);
    assert.ok(Math.abs(mean(rows, salt)) < meanSalt && largest(rows, salt) < stepSalt, `${label} ${name}: salt ${mean(rows, salt)} psu·m on average, ${largest(rows, salt)} at most`);
    assert.ok(Math.abs(mean(rows, volume)) < meanVolume && largest(rows, volume) < stepVolume, `${label} ${name}: volume ${mean(rows, volume)} m on average, ${largest(rows, volume)} at most`);
  }
}

/*
 * The most each phase of a GPU step may change the stored heat (W/m²),
 * salt (psu·m) and volume (m) of spunUp's sea, as [mean over the steps,
 * largest step] for each, the rescale's heat and salt beyond what its
 * water carries at its columns' means and its volume by design: two to
 * six times what an M1 Max leaves here. The surface read adds nothing
 * when the surface has not changed, and without freshwater or ice the
 * velocity shift and the salt and surface write leave the stores alone.
 * The mixed layer's exchanges move whole multiples of the coarser store's
 * spacing; the eddy transport moves them across each edge and gives what
 * rounding leaves in a column to its finest stores; the RK4 and the
 * rescale round each cell's update, and the RK4's water carries its salt.
 * The deposit of the passed sunlight changes each class by whole multiples
 * of its store's spacing and keeps the rest waiting, so only the sums of
 * its float32 subtractions round, and the surface read deposits what
 * waited from the step before.
 */
const EXACT = [1e-9, 1e-9, 1e-12, 1e-12, 1e-15, 1e-15];
const GPU_FLOORS = {
  readSurface: EXACT,
  dynamics: [0.02, 0.1, 3e-5, 3e-4, 1e-6, 1e-5],
  rescale: [0.01, 0.2, 2e-8, 1e-5, Infinity, Infinity],
  eddy: [0.01, 0.01, 1e-6, 2e-6, 3e-8, 6e-8],
  shift: EXACT,
  mixedLayer: [1e-4, 2e-3, 1e-9, 1e-8, 1e-12, 1e-11],
  'salt and surface': EXACT,
  deposit: [1e-7, 2e-7, 1e-12, 1e-12, 1e-15, 1e-15],
};

test('in single precision every phase of the GPU ocean\'s step keeps the heat, salt and volume it stores to its rounding, and the rescale changes heat and salt only with its water; the surface read deposits the sunlight that passed the mixed layer, keeping what its stores cannot take waiting', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const sea = await onGpu(spunUp());
  const record = await sea.follow(16);
  sea.model.destroy();
  report('GPU', record);
  within('GPU', record, GPU_FLOORS);
  const taken = record.deposit.map((r, n) => r.stored / sea.light[n]), overall = record.deposit.reduce((s, r) => s + r.stored, 0) / sea.light.reduce((s, x) => s + x, 0);
  console.log(`GPU deposit: ${show(mean(sea.light, (x) => x))} W/m² of sunlight passed the mixed layer each step; the stores took ${Math.min(...taken).toFixed(4)} to ${Math.max(...taken).toFixed(4)} of each step's, with what waited from the step before, and ${overall.toFixed(5)} of all of it`);
  assert.ok(Math.min(...taken) > 0.95 && Math.max(...taken) < 1.05 && Math.abs(overall - 1) < 0.01, `the stores took ${taken} of the light`);
});

/*
 * The float32 spacing of a store holding x, as the kernels' grain() has
 * it: 2^(e − 24) for x = m·2^e with ½ ≤ |m| < 1.
 */
const bits = new Float32Array(1), word = new Uint32Array(bits.buffer);
const spacing = (x) => { bits[0] = x; return x === 0 ? 2 ** -24 : 2 ** Math.max(((word[0] >>> 23) & 255) - 150, -126); };

test('the GPU deposit changes each class\'s heat by whole multiples of its store\'s spacing and keeps the rest waiting, so that in every column the stores\' change and the waiting light add up to the light passed, to the rounding of one subtraction', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const sea = await onGpu(spunUp()), { model, ocean, stores } = sea, { gpu } = model, PH = gpu.layout.PH, layers = ocean.layers;
  const rhoCp = Math.fround(ocean.options.density * ocean.options.specificHeat);
  ocean.readSurfaceFromAtmosphere();
  const [, transmission] = await Promise.all([stores(), readRanges(gpu.device, gpu.buffers.PH, [{ offset: PH.TRANS, length: C }]).then(([t]) => t)]);
  const given = Float32Array.from(passing(transmission));
  gpu.device.queue.writeBuffer(gpu.buffers.PH, 4 * PH.PEN, given);
  const [, before] = await stores();
  ocean.readSurfaceFromAtmosphere();
  const [, after, , waiting] = await stores();
  let classes = 0, off = 0, cells = 0, worst = 0, kept = 0;
  for (let i = 0; i < C; i++) {
    if (!(given[i] > 0)) continue;
    cells++;
    let change = 0;
    for (let k = 0; k < layers; k++) {
      const n = k * C + i, dq = after[n] - before[n];
      if (dq === 0) continue;
      classes++;
      if (dq % spacing(before[n]) !== 0) off++;
      change += dq;
    }
    worst = Math.max(worst, Math.abs(rhoCp * change + waiting[i] - given[i]) / given[i]);
    kept = Math.max(kept, Math.abs(waiting[i]) / given[i]);
  }
  model.destroy();
  console.log(`${cells} columns took the light passed in ${classes} class changes, ${off} of them not whole multiples of their store's spacing; the stores and the waiting light add up to the light passed within ${worst.toExponential(2)} of it, the most left waiting ${(100 * kept).toFixed(2)} % of a column's light`);
  assert.ok(cells > 0 && off === 0, `${off} of ${classes} changes are not whole multiples of their store's spacing`);
  assert.ok(worst < 1e-6, `a column's stores and waiting light differ from the light passed by ${worst} of it`);
});

test('with carryMomentum the GPU mixed layer\'s exchanges keep each edge\'s column momentum to single-precision rounding, and its heat, salt and volume as without', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const sea = await onGpu(spunUp(SPINUP, CARRIED));
  const record = await sea.follow(8);
  sea.model.destroy();
  report('GPU carried', { mixedLayer: record.mixedLayer });
  within('GPU carried', record, { mixedLayer: GPU_FLOORS.mixedLayer });
  assert.ok(largest(record.mixedLayer, (r) => r.momentum) < 1e-5, `the mixed layer changes an edge's momentum by ${largest(record.mixedLayer, (r) => r.momentum)}`);
});

/*
 * The eddy transport straight after createOcean's start, when its
 * interfaces are steepest and the transport strongest; bounds as for
 * GPU_FLOORS.
 */
const FRESH_EDDY = { eddy: [0.01, 0.01, 3e-6, 3e-6, 6e-8, 1e-7] };

test('from a fresh start the GPU ocean\'s eddy transport keeps the volume, heat and salt it moves to its rounding', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const sea = await onGpu(spunUp(0));
  const record = await sea.follow(8);
  sea.model.destroy();
  report('GPU fresh', { eddy: record.eddy });
  within('GPU fresh', record, FRESH_EDDY);
});

test('the GPU ocean loads a saved W of h·S, as files from before Wref carry it, with one rounding of each layer\'s salt, and a saved Wref as it is', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const start = spunUp(0), { ocean: cpu, surfaceT, ice } = start;
  const first = await onGpu(start);
  const [, , w] = await first.stores();
  let rounded = 0;
  for (let n = 0; n < w.length; n++) if (w[n] !== Math.fround(cpu.W[n] - SALT_REFERENCE * cpu.h[n])) rounded++;
  assert.equal(rounded, 0, `${rounded} of ${w.length} stored salts are not h·S − ${SALT_REFERENCE}·h rounded once`);
  const saved = { ...(await first.ocean.serialize()), ...(await first.ocean.restartArrays()) };
  assert.ok(saved.Wref && saved.Qref && !saved.W && !saved.Q, 'the restart arrays carry Wref and Qref');
  const legacy = { ...saved, W: Float64Array.from(saved.Wref, (x, n) => x + SALT_REFERENCE * saved.h[n]) };
  delete legacy.Wref;
  for (const [label, file, salt] of [['Wref', saved, (n) => saved.Wref[n]], ['W', legacy, (n) => Math.fround(legacy.W[n] - SALT_REFERENCE * saved.h[n])]]) {
    const model = await createGpuModel(mesh, { ocean: { everySteps: 1 } }), ocean = model.oceanEngine, { OS } = ocean.layout;
    ocean.upload(file, surfaceT, ice);
    await model.gpu.device.queue.onSubmittedWorkDone();
    const [q, stored] = await readRanges(model.gpu.device, ocean.buffers.S, [{ offset: OS.OQ, length: w.length }, { offset: OS.OW, length: w.length }]);
    let off = 0;
    for (let n = 0; n < w.length; n++) if (stored[n] !== salt(n) || q[n] !== saved.Qref[n]) off++;
    assert.equal(off, 0, `loaded from ${label}, ${off} of ${w.length} layers' heat or salt are not as saved`);
    model.destroy();
  }
  first.model.destroy();
});
