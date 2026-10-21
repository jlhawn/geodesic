import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { gridMesh } from '../js/stateGrid.module.js';
import { bathymetryFrom, createOcean, displacedWeight, inertialStepLimit } from '../js/ocean/layered.module.js';
import { seawaterDensity } from '../js/ocean/seawater.module.js';
import { zonalWindOnEdges, totalHeatSalt } from './helpers/layered.mjs';
import { syntheticTopography, createGeography, encodeLandIce, decodeLandIce, readTopography, geographyOptions, FRESH_GEOGRAPHY, ICE_DENSITY, SEAWATER_DENSITY } from '../js/geography.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const deg = 180 / Math.PI;
const sparse = ({ rows, cols, points, count, ...fields }) => {
  const index = Uint32Array.from(Array.from(count.keys()).filter((p) => count[p] > 0)), take = (Kind, field) => Kind.from(index, (p) => (field ? field[p] : 0));
  return decodeLandIce(encodeLandIce({
    rows, cols, points, index, count: take(Uint8Array, count), floating: take(Uint8Array, fields.floating), thickness: take(Uint16Array, fields.thickness),
    floatThickness: take(Uint16Array, fields.floatThickness), floatSurface: take(Uint16Array, fields.floatSurface), floatBed: take(Int16Array, fields.floatBed),
  }));
};

function shelfBands(N = 8, slope = 0) {
  const mesh = gridMesh(new Grid(N)), rows = 90, cols = 180, points = 225, freeboard = 50, H = 400, draft = H * ICE_DENSITY / SEAWATER_DENSITY, deep = -800, shallow = -Math.round(draft + 30);
  const zone = (lat) => (lat < -70 ? 'grounded' : lat < -54 ? 'shallow' : lat < -38 ? 'deep' : 'sea');
  const raster = syntheticTopography(rows, cols, (lat) => ({ grounded: 500, shallow: freeboard, deep: freeboard, sea: -4000 })[zone(lat * deg)]);
  const pixels = (f) => Array.from({ length: rows * cols }, (_, n) => { const lat = 90 - (Math.floor(n / cols) + 0.5) * 180 / rows; return f(zone(lat), lat); }), afloat = (z) => z === 'shallow' || z === 'deep';
  const shelf = (z, lat) => (z === 'deep' ? Math.round(H + slope * (-46 - lat) / 8) : H);
  const ice = sparse({
    rows, cols, points, count: pixels((z) => (z === 'sea' ? 0 : points)), thickness: pixels((z, lat) => (z === 'grounded' ? 1500 : afloat(z) ? shelf(z, lat) : 0)),
    floating: pixels((z) => (afloat(z) ? points : 0)), floatThickness: pixels((z, lat) => (afloat(z) ? shelf(z, lat) : 0)), floatSurface: pixels((z) => (afloat(z) ? freeboard : 0)),
    floatBed: pixels((z) => (z === 'deep' ? deep : z === 'shallow' ? shallow : 0)),
  });
  return { mesh, raster, zone, draft, deep, options: { landBridges: {}, seaStraits: {}, landIce: ice } };
}

async function prepare(model) {
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  if (model.load) model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  return model;
}

// Under a shelf the water starts on ISOMIP+'s warm profile: in a neutral column rounding decides convection, and f32 and f64 round apart.
function stratify(cpu, gpu) {
  const ocean = cpu.ocean, C = cpu.mesh.nCells, cavity = cpu.geography.cavity ?? new Uint8Array(C), saved = { ...ocean.serialize(), ...ocean.restartArrays(), flux: new Float64Array(C) };
  const isomip = (z) => [271.25 + 2.9 * z / 720, 33.8 + 0.9 * z / 720];
  for (let i = 0; i < C; i++) {
    if (!cavity[i]) continue;
    let z = cpu.geography.iceDraft[i];
    for (let k = 0; k < ocean.layers; k++) {
      const n = k * C + i, hk = ocean.h[n], [t, sal] = isomip(z + 0.5 * hk);
      saved.Q[n] = hk * t;
      saved.W[n] = hk * sal;
      z += hk;
    }
    saved.previousT0[i] = saved.Q[i] / ocean.h[i];
    saved.shelfLoad[i] = displacedWeight(cpu.geography.iceDraft[i], isomip);
  }
  ocean.load(saved, cpu.state[3], cpu.state[6]);
  gpu.ocean.load(saved, gpu.state[3], gpu.state[6]);
}

test('a cavity is land to the atmosphere and sea beneath its shelf, from the draft of the floating ice down to the bed, and stays land where less than cavityDepth of water lies under the shelf', () => {
  const { mesh, raster, zone, draft, deep, options } = shelfBands(), shelf = createGeography(mesh, raster, options), cavity = createGeography(mesh, raster, { ...options, iceShelves: 'cavity' });
  for (const name of ['land', 'landFraction', 'elevation', 'iceSheet', 'iceFraction', 'iceThickness', 'shelfFraction', 'shelfThickness', 'wetDepth', 'edgeOcean', 'edgeSill', 'coastEdges']) {
    assert.deepEqual(Array.from(cavity[name]), Array.from(shelf[name]), `the atmosphere's ${name}`);
  }
  const { ocean } = cavity, D = bathymetryFrom(mesh, ocean), pure = new Int8Array(mesh.nCells).fill(-1), seen = { grounded: 0, shallow: 0, deep: 0, sea: 0 };
  for (let i = 0; i < mesh.nCells; i++) {
    const lat = mesh.latCell[i] * deg, z = zone(lat);
    if (zone(lat - 5) !== z || zone(lat + 5) !== z) continue;
    seen[z]++; pure[i] = 1;
    const wet = { grounded: 0, shallow: 0, deep: 1, sea: 1 }[z];
    assert.deepEqual([cavity.cavity[i], 1 - ocean.land[i]], [z === 'deep' ? 1 : 0, wet], `${z} cell ${i} at ${lat.toFixed(1)}°`);
    if (z === 'deep') {
      assert.deepEqual([cavity.land[i], cavity.iceSheet[i]], [1, 1]);
      assert.ok(Math.abs(cavity.iceDraft[i] - draft) < 1e-9 && Math.abs(ocean.elevation[i] - deep) < 1e-9 && D[i] >= -deep, `cavity ${i}: draft ${cavity.iceDraft[i]} m, bed ${ocean.elevation[i]} m, depth ${D[i]} m`);
    } else assert.equal(cavity.iceDraft[i], 0);
    if (z === 'sea') assert.equal(ocean.elevation[i], -4000);
  }
  assert.ok(Object.values(seen).every((n) => n >= 5), JSON.stringify(seen));
  let inner = 0;
  for (let e = 0; e < mesh.nEdges; e++) {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
    assert.equal(ocean.edgeOcean[e], !ocean.land[a] && !ocean.land[b] ? 1 : 0);
    if (!ocean.edgeOcean[e]) assert.equal(ocean.edgeSill[e], 0);
    if (pure[a] > 0 && pure[b] > 0 && cavity.cavity[a] && cavity.cavity[b]) { inner++; assert.equal(ocean.edgeSill[e], -deep); }
  }
  assert.ok(inner >= 5, `${inner} edges inside the cavity`);
  const thin = createGeography(mesh, raster, { ...options, iceShelves: 'cavity', cavityDepth: 20 });
  for (let i = 0; i < mesh.nCells; i++) if (pure[i] > 0 && zone(mesh.latCell[i] * deg) === 'shallow') assert.equal(thin.cavity[i], 1, `shallow cell ${i} holds 30 m under its shelf`);
});

test('under a shelf the ocean takes nothing from the air: the temperature, ice, rain, evaporation and stress there change nothing in the water, and the ocean writes no surface temperature, flux, light share or fresh water there', () => {
  const { mesh, raster, options } = shelfBands(), geography = createGeography(mesh, raster, { ...options, iceShelves: 'cavity' });
  const C = mesh.nCells, E = mesh.nEdges, { cavity } = geography, cellsOf = (e) => [mesh.cellsOnEdge[2 * e], mesh.cellsOnEdge[2 * e + 1]];
  const inside = (e) => cellsOf(e).every((c) => cavity[c]), front = (e) => geography.ocean.edgeOcean[e] && cavity[cellsOf(e)[0]] !== cavity[cellsOf(e)[1]];
  const volume = (ocean) => { let sum = 0; for (let n = 0; n < ocean.h.length; n++) sum += mesh.areaCell[n % C] * ocean.h[n]; return sum; };
  const run = (over) => {
    const ocean = createOcean(mesh, { geography, everySteps: 1 });
    const surfaceT = Float64Array.from({ length: C }, (_, i) => (cavity[i] ? over.T : 285 - 15 * Math.sin(mesh.latCell[i]) ** 2)), flux = new Float64Array(C);
    const ice = Float64Array.from({ length: C }, (_, i) => (cavity[i] ? over.ice : 0));
    const stress = Float64Array.from({ length: E }, (_, e) => (inside(e) ? over.stress : 0.1 * Math.cos(3 * mesh.latCell[cellsOf(e)[0]])));
    const evaporation = Float64Array.from({ length: C }, (_, i) => (cavity[i] ? over.evaporation : 3e-5)), rain = Float64Array.from({ length: C }, (_, i) => (cavity[i] ? over.rain : 0.05));
    ocean.initialize(surfaceT, ice);
    const before = volume(ocean);
    for (let n = 0; n < 24; n++) { ocean.accumulate(evaporation, rain, 3600); ocean.advance(surfaceT, ice, flux, stress, 3600); }
    return { ocean, surfaceT, flux, before, after: volume(ocean) };
  };
  const a = run({ T: 240, ice: 2, stress: 0.8, evaporation: 1e-4, rain: 3 }), b = run({ T: 300, ice: 0, stress: 0, evaporation: 0, rain: 0 });
  for (const name of ['h', 'u', 'Q', 'W', 'eta']) {
    let differ = 0;
    for (let n = 0; n < a.ocean[name].length; n++) if (a.ocean[name][n] !== b.ocean[name][n]) differ++;
    assert.equal(differ, 0, `${name}: ${differ} values differ`);
  }
  let cells = 0, inner = 0, fronts = 0, moving = 0;
  for (let i = 0; i < C; i++) {
    if (!cavity[i]) continue;
    cells++;
    assert.deepEqual([a.surfaceT[i], b.surfaceT[i], a.flux[i], b.flux[i], a.ocean.transmission[i], a.ocean.fresh[i], a.ocean.calved[i]], [240, 300, 0, 0, 0, 0, 0], `cavity cell ${i}`);
  }
  for (let e = 0; e < E; e++) {
    if (inside(e)) { inner++; assert.equal(a.ocean.stress[e], 0, `edge ${e} inside the cavity`); }
    if (front(e)) { fronts++; moving = Math.max(moving, Math.abs(a.ocean.u[e])); }
  }
  assert.ok(cells >= 10 && inner >= 10 && fronts >= 10, `${cells} cavity cells, ${inner} edges inside, ${fronts} on the front`);
  assert.ok(moving > 1e-4, `the water at the front moves (${moving} m/s), so the cavity is not walled off`);
  assert.ok(Math.abs(a.after - a.before) < 1e-12 * a.before, `the sea and its cavities hold their water: ${a.after - a.before} m³ of ${a.before}`);
});

test('under a shelf the boundary layer is stirred by the shelf\'s drag on the water beneath it, u*² = C_d|u|², and by neither the wind nor the release of convection: still water there has none, moving water has one whatever the stress on the front, and that stress deepens the open water\'s boundary layer beside it', () => {
  const { mesh, raster, options } = shelfBands(), geography = createGeography(mesh, raster, { ...options, iceShelves: 'cavity' });
  const C = mesh.nCells, E = mesh.nEdges, { cavity } = geography, cellsOf = (e) => [mesh.cellsOnEdge[2 * e], mesh.cellsOnEdge[2 * e + 1]];
  const front = (e) => geography.ocean.edgeOcean[e] && cavity[cellsOf(e)[0]] !== cavity[cellsOf(e)[1]];
  const eastward = (e) => { const x = mesh.xEdge[3 * e], y = mesh.xEdge[3 * e + 1]; return (x * mesh.nEdge[3 * e + 1] - y * mesh.nEdge[3 * e]) / Math.hypot(x, y); };
  const mixed = (scale, speed) => {
    const ocean = createOcean(mesh, { geography, everySteps: 1 }), surfaceT = new Float64Array(C).fill(300), ice = new Float64Array(C);
    ocean.initialize(surfaceT, ice);
    for (let e = 0; e < E; e++) if (geography.ocean.edgeOcean[e]) ocean.u[e] = speed * eastward(e);
    ocean.setStress(Float64Array.from({ length: E }, (_, e) => (front(e) ? 0.3 * scale : 0.05)), ice);
    ocean.columnPhysics(3600);
    return ocean;
  };
  const still = mixed(10, 0), calm = mixed(1, 0.1), stormy = mixed(10, 0.1), L = calm.layers, beside = new Set();
  for (let e = 0; e < E; e++) if (front(e)) for (const c of cellsOf(e)) if (!cavity[c]) beside.add(c);
  let stirred = 0, cavities = 0;
  for (let i = 0; i < C; i++) {
    if (!cavity[i]) continue;
    cavities++;
    for (let k = 0; k < L; k++) assert.deepEqual([stormy.Q[k * C + i], stormy.W[k * C + i]], [calm.Q[k * C + i], calm.W[k * C + i]], `cavity cell ${i}, layer ${k}`);
    assert.equal(still.boundaryDepths[i], 0, `cavity cell ${i}'s boundary layer in still water`);
    assert.equal(stormy.boundaryDepths[i], calm.boundaryDepths[i], `cavity cell ${i}'s boundary layer under a stormy front`);
    if (calm.boundaryDepths[i] > 1) stirred++;
  }
  assert.ok(cavities >= 10 && stirred === cavities, `${stirred} of ${cavities} cavity columns have a boundary layer under water moving at 0.1 m/s`);
  let deepened = 0;
  for (const i of beside) if (stormy.boundaryDepths[i] > calm.boundaryDepths[i]) deepened++;
  assert.ok(beside.size >= 10 && deepened > beside.size / 2, `${deepened} of ${beside.size} open cells beside the front deepened their boundary layer`);
});

test('both engines run the cavities alike at N=8: the water under the shelves tracks across engines as the open sea does without them, with no sea ice, ocean flux or light share in a cavity and the land\'s skin over it', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { raster, options } = shelfBands();
  const pair = async (iceShelves) => {
    const setup = { topography: raster, geography: { ...options, iceShelves } };
    const cpu = await prepare(createModel(new Grid(8), setup)), gpu = await prepare(await createGpuModel(new Grid(8), setup));
    stratify(cpu, gpu);
    for (let n = 0; n < 48; n++) { cpu.step(900); await gpu.step(900); }
    await gpu.sync();
    return { cpu, gpu, a: cpu.ocean.serialize(), b: await gpu.ocean.serialize(), end: await gpu.gpu.downloadPhysics() };
  };
  const shelves = await pair('cavity'), bare = await pair('land');
  const { cavity } = shelves.cpu.geography, mesh = shelves.cpu.mesh, C = mesh.nCells, E = mesh.nEdges, L = shelves.cpu.ocean.layers;
  const edgeWhere = (e) => ['sea', 'front', 'inside'][cavity[mesh.cellsOnEdge[2 * e]] + cavity[mesh.cellsOnEdge[2 * e + 1]]];
  const spread = ({ cpu, a, b }) => {
    const sums = {}, add = (where, name, d) => { const s = (sums[`${where} ${name}`] ??= { sq: 0, n: 0, max: 0 }); s.sq += d * d; s.n++; s.max = Math.max(s.max, Math.abs(d)); };
    for (let i = 0; i < C; i++) {
      if (cpu.ocean.D[i] === 0) continue;
      const where = cavity[i] ? 'cavity' : 'sea';
      add(where, 'eta', a.eta[i] - b.eta[i]);
      for (let k = 0; k < L; k++) {
        const n = k * C + i;
        add(where, 'h', a.h[n] - b.h[n]);
        if (a.h[n] > 1 && b.h[n] > 1) { add(where, 'T', a.T[n] - b.T[n]); add(where, 'S', a.S[n] - b.S[n]); }
      }
    }
    for (let e = 0; e < E; e++) if (cpu.ocean.sillEdge[e] > 0) for (let k = 0; k < L; k++) add(edgeWhere(e), 'u', a.u[k * E + e] - b.u[k * E + e]);
    return Object.fromEntries(Object.entries(sums).map(([key, s]) => [key, { rms: Math.sqrt(s.sq / s.n), max: s.max }]));
  };
  const withShelves = spread(shelves), without = spread(bare), show = (o) => Object.entries(o).map(([k, v]) => `${k} rms ${v.rms.toExponential(1)} max ${v.max.toExponential(1)}`).join(', ');
  console.log(`N=8, twelve hours, CPU − GPU with cavities: ${show(withShelves)}`);
  console.log(`N=8, twelve hours, CPU − GPU with the shelves as land: ${show(without)}`);
  for (const name of ['h', 'u', 'T', 'S', 'eta']) {
    const base = without[`sea ${name}`].rms;
    for (const where of name === 'u' ? ['front', 'inside', 'sea'] : ['cavity', 'sea']) assert.ok(withShelves[`${where} ${name}`].rms < 3 * base, `${where} ${name}: rms ${withShelves[`${where} ${name}`].rms} against ${base} in the open sea without cavities`);
  }
  let skin = 0;
  for (let i = 0; i < C; i++) {
    if (!cavity[i]) continue;
    const { cpu, gpu, end } = shelves;
    assert.deepEqual([cpu.state[6][i], gpu.state[6][i], end.OFLUX[i], end.TRANS ? end.TRANS[i] : 0], [0, 0, 0, 0], `cavity cell ${i}: sea ice, the ocean's flux and its light share`);
    skin = Math.max(skin, Math.abs(cpu.state[3][i] - gpu.state[3][i]));
  }
  assert.ok(skin < 0.05, `the skin over the cavities: ${skin} K apart`);
});

test('with the shelves loading the water and its top at their sloping draft, every phase of the ocean step but the rescale keeps its water, heat and salt to rounding, and the rescale corrects the columns under and beside the shelves as little as the open sea\'s', () => {
  const { mesh, raster, options } = shelfBands(8, 250), geography = createGeography(mesh, raster, { ...options, iceShelves: 'cavity' }), DT = 1350;
  const C = mesh.nCells, ocean = createOcean(mesh, { geography, everySteps: 1 }), L = ocean.layers, P = ocean.phases, parts = Math.max(1, Math.ceil(DT / inertialStepLimit(mesh))), part = DT / parts;
  const surfaceT = new Float64Array(C).fill(290), ice = new Float64Array(C), flux = new Float64Array(C), wind = zonalWindOnEdges(mesh, (lat) => 0.1 * Math.cos(2 * lat));
  ocean.initialize(surfaceT, ice);
  ocean.advance(surfaceT, ice, flux, wind, DT);
  const column = (i) => { let sum = 0; for (let k = 0; k < L; k++) sum += ocean.h[k * C + i]; return sum; };
  const totals = () => { const { heat, salt } = totalHeatSalt(ocean, mesh); let volume = 0; for (let i = 0; i < C; i++) if (ocean.cellOcean[i]) volume += mesh.areaCell[i] * column(i); return [volume, heat, salt]; };
  const start = totals(), drift = [0, 0, 0], kept = (phase) => { const a = totals(); phase(); totals().forEach((x, m) => { drift[m] += x - a[m]; }); };
  const beside = new Uint8Array(C);
  for (let e = 0; e < mesh.nEdges; e++) { const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1]; if (ocean.cellOcean[a] && ocean.cellOcean[b] && geography.cavity[a] !== geography.cavity[b]) { beside[a] = !geography.cavity[a]; beside[b] = !geography.cavity[b]; } }
  const group = (i) => (geography.cavity[i] ? 'under' : beside[i] ? 'beside' : 'open'), corrections = { under: [0, 0], beside: [0, 0], open: [0, 0] };
  for (let n = 0; n < 24; n++) {
    kept(() => { ocean.readSurface(surfaceT, ice); ocean.setStress(wind, ice); });
    for (let m = 0; m < parts; m++) {
      kept(() => P.dynamics(part));
      const before = Float64Array.from({ length: C }, (_, i) => (ocean.cellOcean[i] ? column(i) : 0));
      P.rescale();
      for (let i = 0; i < C; i++) if (ocean.cellOcean[i]) { const c = corrections[group(i)]; c[0] += (column(i) - before[i]) ** 2; c[1]++; }
      kept(() => { P.eddy(part); P.shift(); });
    }
    kept(() => { P.regrid(DT); P.column(DT); ocean.salt(); ocean.writeSurface(surfaceT, flux, DT); });
  }
  const end = totals(), rms = Object.fromEntries(Object.entries(corrections).map(([g, [sq, n]]) => [g, Math.sqrt(sq / n)]));
  assert.ok(Math.abs(end[0] - start[0]) < 1e-12 * start[0], `volume ${start[0]} → ${end[0]}`);
  assert.ok(Math.abs(drift[1]) < 1e-12 * start[1] && Math.abs(drift[2]) < 1e-12 * start[2], `outside the rescale, heat ${drift[1] / start[1]} and salt ${drift[2] / start[2]} of their totals`);
  assert.ok(rms.under < 10 * rms.open && rms.beside < 10 * rms.open, `rescale corrections per step, rms: under the shelves ${rms.under} m, beside them ${rms.beside} m, open sea ${rms.open} m`);
  let under = 0;
  for (let e = 0; e < mesh.nEdges; e++) if (geography.cavity[mesh.cellsOnEdge[2 * e]] && geography.cavity[mesh.cellsOnEdge[2 * e + 1]] && ocean.sillEdge[e] > 0) under = Math.max(under, Math.abs(ocean.u[e]));
  assert.ok(under > 1e-4, `the water under the shelves moves: ${under} m/s`);
});

test('a homogeneous ocean under a sloping shelf with a front stays as still for a month as it does with the shelves as land, on each engine', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { raster, options } = shelfBands(8, 250), T = 271.65, S = 34.4, rho = seawaterDensity(T, S);
  const rest = async (iceShelves) => {
    const setup = { topography: raster, geography: { ...options, iceShelves } }, { mesh, geography } = createModel(new Grid(8), setup), C = mesh.nCells, E = mesh.nEdges;
    const cavity = geography.cavity ?? new Uint8Array(C), cpu = createOcean(mesh, { geography, everySteps: 1 }), L = cpu.layers, ice = new Float64Array(C), calm = new Float64Array(E), flux = new Float64Array(C);
    cpu.initialize(new Float64Array(C).fill(T), ice);
    const saved = { ...cpu.serialize(), ...cpu.restartArrays(), flux: new Float64Array(C) };
    for (let i = 0; i < C; i++) {
      let sum = 0;
      for (let k = 0; k < L; k++) sum += saved.h[k * C + i];
      for (let k = 0; k < L; k++) { const n = k * C + i; if (sum > 0) saved.h[n] *= cpu.D[i] / sum; saved.Q[n] = saved.h[n] * T; saved.W[n] = saved.h[n] * S; saved.T[n] = T; saved.S[n] = S; }
      saved.eta[i] = 0; saved.previousT0[i] = T;
      if (saved.shelfLoad) saved.shelfLoad[i] = cavity[i] ? rho * geography.iceDraft[i] : 0;
    }
    saved.u.fill(0);
    const gpu = (await createGpuModel(new Grid(8), { ...setup, ocean: { everySteps: 1 } })).oceanEngine, cpuSkin = new Float64Array(C).fill(T), gpuSkin = new Float64Array(C).fill(T);
    cpu.load(saved, cpuSkin, ice);
    gpu.upload(saved, gpuSkin, ice);
    for (let n = 0; n < 720; n++) { cpu.advance(cpuSkin, ice, flux, calm, 3600); await gpu.advance(gpuSkin, ice, calm, 3600); }
    const b = await gpu.download(), most = (a) => a.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
    return { cpu: most(cpu.u), gpu: most(b.u), cpuEta: most(cpu.eta), gpuEta: most(b.eta) };
  };
  const shelves = await rest('cavity'), bare = await rest('land'), show = (r) => Object.entries(r).map(([k, v]) => `${k} ${v.toExponential(1)}`).join(', ');
  console.log(`a month at rest, homogeneous water: with the shelves ${show(shelves)}; with them as land ${show(bare)}`);
  // Rounding walks both oceans off rest; 10× leaves room for the cavity's added cells and catches a load, draft or pin that pushes the water.
  for (const [u, eta] of [['cpu', 'cpuEta'], ['gpu', 'gpuEta']]) assert.ok(shelves[u] < 10 * bare[u] && shelves[eta] < 10 * bare[eta], `${u}: ${show(shelves)} against ${show(bare)}`);
});

test('on ETOPO 2022 the cavities lie under the Antarctic shelves near their observed extent, Ross and Filchner–Ronne the largest, the Amundsen Sea shelves resolved at N=128', () => {
  const topography = readTopography('etopo2022'), { geography: base } = geographyOptions(FRESH_GEOGRAPHY);
  const regions = {
    Ross: (la, lo) => la < -77.3 && (lo > 158 || lo < -145),
    'Filchner–Ronne': (la, lo) => la < -74.4 && lo > -85 && lo < -28,
    Amundsen: (la, lo) => la > -76.5 && la < -72.5 && lo > -136 && lo < -98,
    Amery: (la, lo) => la > -73.8 && la < -68.5 && lo > 66 && lo < 75,
    Larsen: (la, lo) => la > -70 && la < -64 && lo > -66 && lo < -58,
  };
  const least = { 64: { all: [100, 1.3], Ross: [30, 0.38], 'Filchner–Ronne': [30, 0.38], Amundsen: [3, 0], Amery: [4, 0] }, 128: { all: [400, 1.2], Ross: [120, 0.38], 'Filchner–Ronne': [110, 0.36], Amundsen: [15, 0], Amery: [15, 0] } };
  for (const N of [64, 128]) {
    const mesh = gridMesh(new Grid(N)), shelf = createGeography(mesh, topography, base), g = createGeography(mesh, topography, { ...base, iceShelves: 'cavity' }), D = bathymetryFrom(mesh, g.ocean);
    for (const name of ['land', 'elevation', 'iceSheet']) assert.deepEqual(Array.from(g[name]), Array.from(shelf[name]), `N=${N}: the atmosphere's ${name}`);
    const out = {}, add = (k, i) => { const o = (out[k] ??= { n: 0, area: 0, draft: 0, water: 0 }); o.n++; o.area += mesh.areaCell[i] / 1e12; o.draft += g.iceDraft[i] * mesh.areaCell[i]; o.water += (D[i] - g.iceDraft[i]) * mesh.areaCell[i]; };
    for (let i = 0; i < mesh.nCells; i++) {
      if (!g.cavity[i]) continue;
      const la = mesh.latCell[i] * deg, lo = mesh.lonCell[i] * deg;
      assert.deepEqual([g.land[i], g.ocean.land[i], g.iceSheet[i]], [1, 0, 1], `N=${N}: cavity ${i} at ${la.toFixed(1)}°, ${lo.toFixed(1)}°`);
      assert.ok(la < -60 || la > 75, `N=${N}: cavity at ${la.toFixed(1)}°`);
      assert.ok(g.iceDraft[i] > 0 && D[i] - g.iceDraft[i] >= 50, `N=${N}: cavity ${i}: draft ${g.iceDraft[i].toFixed(0)} m over ${D[i].toFixed(0)} m`);
      add(Object.keys(regions).find((r) => regions[r](la, lo)) ?? (la < 0 ? 'other Antarctic' : 'Greenland'), i); add('all', i);
    }
    for (const [k, o] of Object.entries(out)) console.log(`N=${N} ${k}: ${o.n} cavities, ${o.area.toFixed(3)} million km², mean draft ${(o.draft / o.area / 1e12).toFixed(0)} m, mean water beneath ${(o.water / o.area / 1e12).toFixed(0)} m`);
    for (const [k, [n, area]] of Object.entries(least[N])) assert.ok(out[k] && out[k].n >= n && out[k].area >= area && out[k].area < 1.7, `N=${N} ${k}: ${JSON.stringify(out[k])}`);
  }
});
