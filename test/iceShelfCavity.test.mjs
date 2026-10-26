import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { gridMesh } from '../js/stateGrid.module.js';
import { bathymetryFrom, createOcean, displacedWeight, inertialStepLimit } from '../js/ocean/layered.module.js';
import { N, shelfBands, warmOcean, restingShelf, steppedDraft } from './helpers/iceShelf.mjs';
import { basalMelt, baseFreezing, SHELF_MELT } from '../js/ocean/iceShelfMelt.module.js';
import { decodeClimatology, encodeClimatology } from '../js/ocean/climatology.module.js';
import { zonalWindOnEdges, totalHeatSalt } from './helpers/layered.mjs';
import { syntheticTopography, createGeography, encodeLandIce, decodeLandIce, readTopography, geographyOptions, FRESH_GEOGRAPHY, ICE_DENSITY, SEAWATER_DENSITY } from '../js/geography.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const deg = 180 / Math.PI;

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

test('under a shelf the boundary layer is stirred by the shelf\'s drag on the water beneath it, u*² = C_d|u|², and by neither the wind nor, while its base melts, the release of convection: still water there has none, moving water has one whatever the stress on the front, and that stress deepens the open water\'s boundary layer beside it', () => {
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

test(`both engines run the cavities alike at N=${N}: the water under the shelves tracks across engines as the open sea does without them, with no sea ice, ocean flux or light share in a cavity and the land's skin over it`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { raster, options } = shelfBands();
  const pair = async (iceShelves) => {
    const setup = { topography: raster, geography: { ...options, iceShelves }, ocean: { iceDrift: false } };
    const cpu = await prepare(createModel(new Grid(N), setup)), gpu = await prepare(await createGpuModel(new Grid(N), setup));
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
  console.log(`N=${N}, twelve hours, CPU − GPU with cavities: ${show(withShelves)}`);
  console.log(`N=${N}, twelve hours, CPU − GPU with the shelves as land: ${show(without)}`);
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
  const total = (o) => cavity.reduce((sum, c, i) => sum + (c ? mesh.areaCell[i] * o.shelfMelted[i] : 0), 0), melt = { cpu: total(shelves.a), gpu: total(shelves.b) };
  const pools = await shelves.gpu.oceanEngine.readShelfMelt(), record = shelves.cpu.ocean.shelfBudget;
  console.log(`twelve hours of melt under the shelves: CPU ${melt.cpu.toExponential(4)} kg, GPU ${melt.gpu.toExponential(4)} kg; last step's southern discharge, melt and draw: CPU ${[record.discharge[0], record.melt[0], record.drawn[0]].map((x) => x.toExponential(3))}, GPU ${[pools.pools[0], pools.pools[2], pools.pools[6]].map((x) => x.toExponential(3))}`);
  assert.ok(melt.cpu > 0 && Math.abs(melt.gpu - melt.cpu) < 1e-3 * melt.cpu, `melt: CPU ${melt.cpu} kg, GPU ${melt.gpu} kg`);
  assert.ok(Math.abs(pools.pools[2] - record.melt[0]) < 1e-3 * record.melt[0] && Math.abs(pools.pools[4] + pools.pools[6] - pools.pools[2]) < 1e-6 * pools.pools[2], `the GPU's southern pool ${Array.from(pools.pools)} against the CPU's ${JSON.stringify(record)}`);
});

test('the three-equation melt balances heat, salt and the freezing point of the interface, takes the root that vanishes with the thermal driving, keeps it in f32, and freezes supercooled water onto the base', () => {
  const k = { rho0: 1025, rhoCp: 1025 * 3985, latentHeat: 3.34e5 }, ustar = 5e-4, A = k.rhoCp * SHELF_MELT.heatTransfer * ustar / k.latentHeat, B = k.rho0 * SHELF_MELT.saltTransfer * ustar;
  const freezing = baseFreezing(1027 * 9.81 * 580), liquidus = (S) => SHELF_MELT.liquidusSlope * S + freezing;
  const Tw = 273.15 + 0.44, Sw = 34.53, { melt, Tb, Sb } = basalMelt(Tw, Sw, freezing, ustar, k);
  const balance = [melt - A * (Tw - Tb), melt * Sb - B * (Sw - Sb), Tb - liquidus(Sb)];
  assert.ok(balance.every((r) => Math.abs(r) < 1e-15), `residuals ${balance}`);
  let lo = 0, hi = Sw;
  for (let n = 0; n < 200; n++) { const mid = 0.5 * (lo + hi); if (A * (Tw - liquidus(mid)) * mid - B * (Sw - mid) > 0) hi = mid; else lo = mid; }
  assert.ok(Math.abs(Sb - lo) < 1e-12, `S_b ${Sb} against ${lo} by bisection`);
  assert.ok(Math.abs(melt / 1.2744e-4 - 1) < 1e-4 && Math.abs(Sb - 19.160) < 1e-3 && Math.abs(Tb - 271.6953) < 1e-4, `580 m draft, 0.44 °C, 34.53 psu, u* 5e-4: m ${melt} kg/m²/s, S_b ${Sb}, T_b ${Tb}`);
  const still = basalMelt(liquidus(Sw), Sw, freezing, ustar, k), frozen = basalMelt(liquidus(Sw) - 0.05, Sw, freezing, ustar, k);
  assert.deepEqual([still.melt, still.Sb], [0, Sw]);
  assert.ok(frozen.melt < 0 && frozen.Sb > Sw && frozen.Tb < liquidus(Sw), `0.05 K supercooled: m ${frozen.melt}, S_b ${frozen.Sb}`);
  const f = Math.fround, zero = baseFreezing(1027 * 9.81 * 580, 0), single = (T, S, fr, u) => {
    const a1 = f(f(k.rhoCp * SHELF_MELT.heatTransfer / k.latentHeat) * u), b1 = f(f(k.rho0 * SHELF_MELT.saltTransfer) * u), sl = f(SHELF_MELT.liquidusSlope), d = f(T - f(f(sl * S) + fr));
    const a = f(-a1 * sl), b = f(f(a1 * f(f(sl * S) - d)) - b1), c = f(f(a1 * d) * S), q = f(-0.5 * f(b + (b > 0 ? 1 : -1) * f(Math.sqrt(Math.max(0, f(f(b * b) - f(f(4 * a) * c)))))));
    const fresher = q !== 0 ? f(c / q) : 0;
    return f(f(b1 * fresher) / f(S - fresher));
  };
  for (const driving of [1, 1e-2, 1e-4, 1e-6]) {
    const T = f(SHELF_MELT.liquidusSlope * Sw + zero + driving), exact = basalMelt(T, f(Sw), f(zero), f(ustar), k).melt, low = single(T, f(Sw), f(zero), f(ustar));
    assert.ok(low > 0 && Math.abs(low - exact) < (1e-6 + 4e-7 / driving) * exact, `f32 at ${driving} K above freezing: ${low} against ${exact}`);
  }
});

test('with the shelves loading the water and its top at their sloping draft, every phase of the ocean step but the rescale keeps its water, heat and salt to rounding beyond what the melt and the outlets\' calving take and give, and the rescale corrects the columns under and beside the shelves as little as the open sea\'s', () => {
  const { mesh, raster, options } = shelfBands(250), geography = createGeography(mesh, raster, { ...options, iceShelves: 'cavity' }), DT = 1350;
  const C = mesh.nCells, ocean = createOcean(mesh, { geography, everySteps: 1, shelfMeltLayer: 0 }), L = ocean.layers, P = ocean.phases, parts = Math.max(1, Math.ceil(DT / inertialStepLimit(mesh))), part = DT / parts;
  const surfaceT = new Float64Array(C).fill(290), ice = new Float64Array(C), flux = new Float64Array(C), wind = zonalWindOnEdges(mesh, (lat) => 0.1 * Math.cos(2 * lat));
  ocean.initialize(surfaceT, ice);
  ocean.advance(surfaceT, ice, flux, wind, DT);
  const column = (i) => { let sum = 0; for (let k = 0; k < L; k++) sum += ocean.h[k * C + i]; return sum; };
  const totals = () => { const { heat, salt } = totalHeatSalt(ocean, mesh); let volume = 0; for (let i = 0; i < C; i++) if (ocean.cellOcean[i]) volume += mesh.areaCell[i] * column(i); return [volume, heat, salt]; };
  const start = totals(), drift = [0, 0, 0], kept = (phase) => { const a = totals(); phase(); totals().forEach((x, m) => { drift[m] += x - a[m]; }); };
  const sums = () => { const q = new Float64Array(C), w = new Float64Array(C); for (let k = 0; k < L; k++) for (let i = 0; i < C; i++) { q[i] += ocean.Q[k * C + i]; w[i] += ocean.W[k * C + i]; } return [q, w]; };
  const worst = [0, 0];
  let melting = 0, excess = 0, covered = 0;
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
    kept(() => P.regrid(DT));
    const top = Array.from({ length: C }, (_, i) => [ocean.Q[i] / ocean.h[i], ocean.W[i] / ocean.h[i]]), [q0, w0] = sums();
    P.column(DT);
    const [q1, w1] = sums();
    for (let i = 0; i < C; i++) {
      if (!ocean.cellOcean[i]) continue;
      const m = ocean.melted[i], water = m / 1000, under = geography.cavity[i], heat = under ? -(3.34e5 * m / ocean.rhoCp + water * (top[i][0] - ocean.meltBase[i])) : 0, salt = under ? -top[i][1] * water : 0;
      if (under) melting += mesh.areaCell[i] * m;
      worst[0] = Math.max(worst[0], Math.abs(q1[i] - q0[i] - heat) / q0[i]); worst[1] = Math.max(worst[1], Math.abs(w1[i] - w0[i] - salt) / w0[i]);
    }
    const runoff = Float64Array.from({ length: C }, (_, j) => (!ocean.cellOcean[j] && mesh.latCell[j] < -Math.PI / 3 ? (n % 2 ? 1e-4 : 1) : 0));
    ocean.accumulate(null, null, DT, runoff, runoff);
    const calved = Float64Array.from(ocean.calved), fresh = Float64Array.from(ocean.fresh), drawn = Float64Array.from(ocean.shelfDraw), W0 = Float64Array.from(ocean.W.subarray(0, C));
    ocean.salt();
    const B = ocean.shelfBudget, pool = (i) => (mesh.fCell[i] < 0 ? 0 : 1), sum = [[0, 0], [0, 0], [0, 0], [0, 0]];
    for (let i = 0; i < C; i++) {
      if (!ocean.cellOcean[i]) continue;
      const a = mesh.areaCell[i], p = pool(i);
      if (geography.cavity[i]) { sum[0][p] += a * ocean.melted[i]; sum[1][p] += a * (ocean.shelfDraw[i] - drawn[i]); continue; }
      sum[2][p] += a * calved[i]; sum[3][p] += a * ocean.calved[i];
      const s = W0[i] / ocean.h[i], expected = s * (fresh[i] + (1 - B.keep[p]) * calved[i]) / 1000;
      worst[1] = Math.max(worst[1], Math.abs(ocean.W[i] - W0[i] - expected) / W0[i]);
    }
    for (let p = 0; p < 2; p++) {
      const close = (x, y) => Math.abs(x - y) <= 1e-12 * Math.max(Math.abs(x), Math.abs(y), 1);
      assert.ok(close(B.melt[p], sum[0][p]) && close(B.discharge[p], sum[2][p]) && B.moved[p] === (B.discharge[p] > 0 ? Math.min(B.melt[p], B.discharge[p]) : 0), `pool ${p}: ${JSON.stringify(B)} against ${sum.map((x) => x[p])}`);
      assert.ok(close(sum[3][p], B.discharge[p] - B.moved[p]) && close(sum[1][p], B.melt[p] - B.moved[p]), `pool ${p}: the outlets calve ${sum[3][p]} of ${B.discharge[p]} kg and the store pays ${sum[1][p]} of ${B.melt[p]} kg melted`);
      if (B.melt[p] > B.discharge[p] && B.discharge[p] > 0) excess++; else if (B.discharge[p] > B.melt[p] && B.melt[p] > 0) covered++;
    }
    const Q0 = Float64Array.from(ocean.Q.subarray(0, C)), calving = Float64Array.from(ocean.calved);
    ocean.meltCalved();
    for (let i = 0; i < C; i++) if (ocean.cellOcean[i] && !geography.cavity[i]) worst[0] = Math.max(worst[0], Math.abs(ocean.Q[i] - Q0[i] + 3.34e5 * calving[i] / ocean.rhoCp) / Q0[i]);
    kept(() => ocean.writeSurface(surfaceT, flux, DT));
  }
  assert.ok(worst[0] < 1e-12 && worst[1] < 1e-12 && melting > 0, `the melt's heat and salt and the outlets' share: worst column ${worst} of its heat and salt, ${melting} kg melted`);
  assert.ok(excess > 0 && covered > 0, `steps where the melt outran the discharge ${excess}, where the discharge covered it ${covered}`);
  const end = totals(), rms = Object.fromEntries(Object.entries(corrections).map(([g, [sq, n]]) => [g, Math.sqrt(sq / n)]));
  assert.ok(Math.abs(end[0] - start[0]) < 1e-12 * start[0], `volume ${start[0]} → ${end[0]}`);
  assert.ok(Math.abs(drift[1]) < 1e-12 * start[1] && Math.abs(drift[2]) < 1e-12 * start[2], `outside the rescale, heat ${drift[1] / start[1]} and salt ${drift[2] / start[2]} of their totals`);
  assert.ok(rms.under < 10 * rms.open && rms.beside < 10 * rms.open, `rescale corrections per step, rms: under the shelves ${rms.under} m, beside them ${rms.beside} m, open sea ${rms.open} m`);
  let under = 0;
  for (let e = 0; e < mesh.nEdges; e++) if (geography.cavity[mesh.cellsOnEdge[2 * e]] && geography.cavity[mesh.cellsOnEdge[2 * e + 1]] && ocean.sillEdge[e] > 0) under = Math.max(under, Math.abs(ocean.u[e]));
  assert.ok(under > 1e-4, `the water under the shelves moves: ${under} m/s`);
});

test('with the melt off, a homogeneous ocean under a sloping shelf with a front stays as still for a month as it does with the shelves as land, on each engine', { skip: !gpuAvailable && 'webgpu not installed' }, () => restingShelf(['gpu']));

test('under a shelf whose draft steps 300 m from cell to cell over 200 m of water, where neighbors share no depth, every edge passes the sill floor, and with the melt off the ocean stays as still for a month as with the shelves as land, homogeneous or on ISOMIP+\'s warm profile, on each engine', { skip: !gpuAvailable && 'webgpu not installed' }, () => steppedDraft(['gpu']));

test('from an atlas a cavity starts with its nearest open column\'s water at the same depth, that column\'s bottom water below its floor, no colder than the in-situ freezing point, on each engine', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { raster, options } = shelfBands(340, { water: 200, shelfSea: 700 }), setup = { topography: raster, geography: { ...options, iceShelves: 'cavity' }, ocean: { iceDrift: false } }, { mesh, geography } = createModel(new Grid(N), setup), C = mesh.nCells;
  const open = (z) => [-2.6 + 0.001 * Math.min(z, 1000), 34 + 0.0005 * Math.min(z, 1000)], depths = [0, 10, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 3000, 4000, 5000];
  const grid = { nLon: 180, nLat: 90, lon0: -179, dLon: 2, lat0: -89, dLat: 2, depths, T: [], S: [] };
  for (const z of depths) for (let r = 0; r < grid.nLat; r++) for (let c = 0; c < grid.nLon; c++) { const [t, s] = open(z), under = grid.lat0 + r * grid.dLat < -42; grid.T.push(under ? t - 3 : t); grid.S.push(under ? s - 0.5 : s); }
  const atlas = decodeClimatology(encodeClimatology(grid)), ice = new Float64Array(C), surfaceT = new Float64Array(C).fill(275);
  const cpu = createOcean(mesh, { geography, climatology: atlas }), started = cpu.initialize(surfaceT, ice), L = cpu.layers, gpu = (await createGpuModel(new Grid(N), setup)).oceanEngine;
  gpu.initialize(surfaceT, ice, { climatology: atlas });
  const cells = [...geography.cavity.keys()].filter((i) => geography.cavity[i]), source = new Map(), filled = (c) => cpu.cellOcean[c] && !geography.cavity[c];
  assert.equal(started.analytic, cells.length, 'the atlas fills every open column and no cavity');
  for (const i of cells) {
    let ring = [i], found = -1;
    const seen = new Set(ring);
    for (let r = 0; r < 40 && found < 0 && ring.length; r++) {
      const next = [];
      for (const j of ring) for (let m = 0; m < mesh.nEdgesOnCell[j]; m++) {
        const c = mesh.cellsOnCell[mesh.maxEdges * j + m];
        if (seen.has(c)) continue;
        seen.add(c);
        if (filled(c)) { found = c; break; }
        if (cpu.cellOcean[c]) next.push(c);
      }
      ring = next;
    }
    source.set(i, found);
  }
  const engines = { cpu: cpu.serialize(), gpu: await gpu.serialize() }, freezing = (s, load, depth) => SHELF_MELT.liquidusSlope * s + baseFreezing(9.81 * (load + 1025 * depth));
  const mean = (a, b, f) => { let sum = 0; for (let m = 0; m < 16; m++) sum += f(a + (m + 0.5) * (b - a) / 16); return sum / 16; };
  const bottomWater = (st, s) => { let m = -1; for (let k = 0; k < L; k++) if (st.h[k * C + s] > 1e-9) m = k * C + s; return [st.T[m], st.S[m]]; };
  for (const [name, st] of Object.entries(engines)) {
    const tight = name === 'cpu' ? 1e-9 : 2e-4, bottoms = new Map();
    let cold = -Infinity, above = 0, capped = 0, below = 0, raised = 0, worstT = 0, worstS = 0;
    for (const i of cells) {
      const draft = geography.iceDraft[i], load = st.shelfLoad[i], s = source.get(i), floor = cpu.D[s];
      let z = draft;
      for (let k = 0; k < L; k++) {
        const n = k * C + i, hk = st.h[n], T = st.T[n], S = st.S[n];
        if (hk > 0.1) {
          cold = Math.max(cold, freezing(S, load, z + 0.5 * hk - draft) - T);
          if (z + hk < floor - 150) {
            above++;
            if (open(z + 0.5 * hk)[0] + 273.15 < freezing(open(z + 0.5 * hk)[1], load, z + 0.5 * hk - draft)) capped++;
            worstT = Math.max(worstT, Math.abs(T - mean(z, z + hk, (y) => Math.max(open(y)[0] + 273.15, freezing(open(y)[1], load, y - draft)))));
            worstS = Math.max(worstS, Math.abs(S - mean(z, z + hk, (y) => open(y)[1])));
          } else if (z > floor) {
            below++;
            const [bottomT, bottomS] = bottoms.get(s) ?? bottoms.set(s, bottomWater(st, s)).get(s), expected = Math.max(bottomT, freezing(bottomS, load, z + 0.5 * hk - draft));
            if (expected > bottomT) raised++;
            assert.ok(Math.abs(T - expected) < tight && Math.abs(S - bottomS) < tight, `${name}: cell ${i} below its source's ${floor.toFixed(0)} m floor holds ${T}, ${S} against its bottom water ${bottomT}, ${bottomS}, no colder than ${expected}`);
          }
        }
        z += hk;
      }
    }
    console.log(`${name}: ${above} cavity layers 150 m above their source's floor (${capped} at the in-situ freezing point) within ${worstT.toFixed(3)} K and ${worstS.toFixed(4)} psu of the open water at their depth, ${below} below their source's floor on its bottom water (${raised} raised to the in-situ freezing point); the coldest ${cold.toExponential(1)} K from the in-situ freezing point`);
    assert.ok(capped > 0 && above > capped && below > 0, `${name}: ${above} layers above a floor, ${capped} capped, ${below} below`);
    assert.ok(worstT < 0.1 && worstS < 0.05 && cold < tight, `${name}: ${worstT} K, ${worstS} psu, ${cold} K below freezing`);
  }
  const column = (st, i, field) => { let sum = 0; for (let k = 0; k < L; k++) sum += st.h[k * C + i] * st[field][k * C + i]; return sum; };
  for (const i of cells) for (const field of ['T', 'S']) {
    const a = column(engines.cpu, i, field), b = column(engines.gpu, i, field);
    assert.ok(Math.abs(a - b) < 1e-5 * Math.abs(a), `cell ${i} ${field}: cpu ${a}, gpu ${b}`);
  }
  for (const i of cells) assert.ok(Math.abs(engines.cpu.shelfLoad[i] - engines.gpu.shelfLoad[i]) < 1e-5 * engines.cpu.shelfLoad[i], `cell ${i} load: cpu ${engines.cpu.shelfLoad[i]}, gpu ${engines.gpu.shelfLoad[i]}`);
});

test('under a still shelf the release of convection from an unstable top reaches the boundary layer only while the base freezes: supercooled water deepens it on both engines, and water as far above its freezing point melts the base and has none', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { raster, options } = shelfBands(250), setup = { topography: raster, geography: { ...options, iceShelves: 'cavity' }, ocean: { everySteps: 1, iceDrift: false } };
  const { mesh, geography } = createModel(new Grid(N), setup), C = mesh.nCells, cells = [...geography.cavity.keys()].filter((i) => geography.cavity[i]);
  const columns = async (offset) => {
    const cpu = createOcean(mesh, { geography, everySteps: 1 }), L = cpu.layers, ice = new Float64Array(C), skin = new Float64Array(C).fill(275);
    cpu.initialize(skin, ice);
    const saved = { ...cpu.serialize(), ...cpu.restartArrays(), flux: new Float64Array(C) };
    for (const i of cells) {
      let z = 0;
      for (let k = 0; k < L; k++) {
        const n = k * C + i, S = 34.5 + (k === 0 ? 0.05 : 0) + (z >= 10 ? 0.1 : 0), T = SHELF_MELT.liquidusSlope * S + baseFreezing(9.81 * saved.shelfLoad[i]) + offset;
        saved.Q[n] = saved.h[n] * T; saved.W[n] = saved.h[n] * S; saved.T[n] = T; saved.S[n] = S; z += saved.h[n];
        if (k === 0) saved.previousT0[i] = T;
      }
    }
    saved.u.fill(0);
    const gpu = (await createGpuModel(new Grid(N), setup)).oceanEngine;
    cpu.load(saved, skin, ice);
    gpu.upload(saved, Float64Array.from(skin), ice);
    cpu.phases.column(3600);
    gpu.phases.column(3600);
    const { melted } = await gpu.readShelfMelt(), { boundaryDepths } = await gpu.restartArrays();
    return cells.map((i) => [cpu.melted[i], melted[i], cpu.boundaryDepths[i], boundaryDepths[i]]);
  };
  const frozen = await columns(-0.05), melting = await columns(0.05), least = (rows, m) => Math.min(...rows.map((r) => r[m]));
  console.log(`0.05 K supercooled under a top 0.05 psu saltier than the 10 m below it: boundary layer CPU ${least(frozen, 2).toFixed(1)}–${Math.max(...frozen.map((r) => r[2])).toFixed(1)} m, GPU ${least(frozen, 3).toFixed(1)}–${Math.max(...frozen.map((r) => r[3])).toFixed(1)} m; melt CPU ${frozen[0][0].toExponential(4)}, GPU ${frozen[0][1].toExponential(4)} kg/m² in an hour`);
  for (const [cpuMelt, gpuMelt, cpuDepth, gpuDepth] of frozen) assert.ok(cpuMelt < 0 && gpuMelt < 0 && cpuDepth > 1 && gpuDepth > 1, `supercooled: melt ${cpuMelt} / ${gpuMelt}, boundary layer ${cpuDepth} / ${gpuDepth} m`);
  for (const [cpuMelt, gpuMelt, cpuDepth, gpuDepth] of melting) assert.ok(cpuMelt > 0 && gpuMelt > 0 && cpuDepth === 0 && gpuDepth === 0, `above freezing: melt ${cpuMelt} / ${gpuMelt}, boundary layer ${cpuDepth} / ${gpuDepth} m`);
});

test('with discharge at the outlets both engines take the melt out of the calving alike, when the discharge covers the melt and when the melt outruns it and draws on the store, with T_w from the top layer or the 20 m layer', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { raster, options } = shelfBands(250), geo = { ...options, iceShelves: 'cavity' }, { mesh, geography } = createModel(new Grid(N), { topography: raster, geography: geo });
  const C = mesh.nCells, E = mesh.nEdges, area = mesh.areaCell, cells = [...geography.cavity.keys()].filter((i) => geography.cavity[i]), calm = new Float64Array(E), flux = new Float64Array(C);
  const total = (field) => cells.reduce((sum, i) => sum + area[i] * field[i], 0), near = (x, y, tol) => Math.abs(x - y) <= tol * Math.max(Math.abs(x), Math.abs(y));
  for (const shelfMeltLayer of [0, 20]) {
    const { ocean: cpu, surfaceT, ice } = warmOcean(mesh, geography, { shelfMeltLayer });
    const model = await createGpuModel(new Grid(N), { topography: raster, geography: geo, ocean: { everySteps: 1, shelfMeltLayer, iceDrift: false } }), gpu = model.oceanEngine, core = model.gpu;
    gpu.upload({ ...cpu.serialize(), ...cpu.restartArrays(), flux: new Float64Array(C) }, Float64Array.from(surfaceT), ice);
    const shed = new Float32Array(C), land = Array.from({ length: C }, (_, j) => !cpu.cellOcean[j] && mesh.latCell[j] < -Math.PI / 3);
    for (const amount of [10, 0.3]) {
      const runoff = Float64Array.from(land, (x) => (x ? amount : 0));
      runoff.forEach((x, j) => { shed[j] += x; });
      cpu.accumulate(null, null, 3600, runoff, runoff);
      core.device.queue.writeBuffer(core.buffers.PH, 4 * core.layout.PH.RUNOFF, shed);
      core.device.queue.writeBuffer(core.buffers.PH, 4 * core.layout.PH.DISCHARGE, shed);
      gpu.accumulateFreshwater(3600);
      const outlets = [...cpu.calved.keys()].filter((i) => cpu.calved[i] > 0), cooling = Math.max(...outlets.map((i) => 3.34e5 * cpu.calved[i] / cpu.rhoCp / cpu.h[i])), skin = Float64Array.from(surfaceT);
      cpu.advance(surfaceT, ice, flux, calm, 3600);
      await gpu.advance(skin, ice, calm, 3600);
      const B = cpu.shelfBudget, { pools, melted, shelfDraw } = await gpu.readShelfMelt(), b = await gpu.download();
      const dT = Math.max(...outlets.map((i) => Math.abs(cpu.Q[i] / cpu.h[i] - b.T1[i]))), dS = Math.max(...outlets.map((i) => Math.abs(cpu.W[i] / cpu.h[i] - b.S1[i])));
      const label = `${shelfMeltLayer || 3} m layer, ${amount} kg/m² discharged`;
      console.log(`${label}: southern discharge ${B.discharge[0].toExponential(3)} / ${pools[0].toExponential(3)} kg, melt ${B.melt[0].toExponential(4)} / ${pools[2].toExponential(4)}, calving kept ${B.keep[0].toFixed(4)} / ${pools[8].toFixed(4)}; ${outlets.length} outlets cooled by up to ${cooling.toFixed(3)} K differ by ${dT.toExponential(1)} K and ${dS.toExponential(1)} psu`);
      assert.ok(amount > 1 ? B.discharge[0] > B.melt[0] && B.keep[0] > 0 : B.melt[0] > B.discharge[0] && B.keep[0] === 0 && pools[8] === 0 && B.drawn[0] > 0, `${label}: ${JSON.stringify(B)}, GPU keeps ${pools[8]}`);
      assert.ok(near(pools[0], B.discharge[0], 1e-5) && near(pools[2], B.melt[0], 1e-3) && near(pools[4], B.moved[0], 1e-3) && Math.abs(pools[8] - B.keep[0]) < 1e-3 && near(total(melted), total(cpu.melted), 1e-3), `${label}: GPU pools ${Array.from(pools)} against ${JSON.stringify(B)}`);
      assert.ok(outlets.length > 0 && dT < 1e-3 && dS < 1e-3 && cooling > 0.02, `${label}: outlets differ by ${dT} K and ${dS} psu against a calving cooling of up to ${cooling} K`);
      if (amount < 1) assert.ok(near(total(shelfDraw), total(cpu.shelfDraw), 1e-3) && total(cpu.shelfDraw) > 0, `${label}: the store pays ${total(cpu.shelfDraw)} kg on the CPU, ${total(shelfDraw)} on the GPU`);
    }
  }
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
