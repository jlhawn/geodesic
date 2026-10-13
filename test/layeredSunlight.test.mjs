import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean, transmitted, penetrationShares, PENETRATION_DEPTH, THIN } from '../js/ocean/layered.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';

const J = (z) => 0.67 * Math.exp(-z) + 0.33 * Math.exp(-z / 17) - (0.67 * Math.exp(-PENETRATION_DEPTH) + 0.33 * Math.exp(-PENETRATION_DEPTH / 17));

test('Jerlov type IB passes 18 % of the absorbed sunlight below 10 m, 10 % below 20 m, 1.7 % below 50 m and 0.09 % below 100 m, and none below 300 m', () => {
  for (const [z, share] of [[10, 0.183], [20, 0.102], [50, 0.0174], [100, 9.2e-4], [200, 2.6e-6]]) {
    assert.ok(Math.abs(transmitted(z) - J(z)) < 1e-15, `${z} m: ${transmitted(z)} against the two bands' ${J(z)}`);
    assert.ok(Math.abs(transmitted(z) / share - 1) < 0.02, `${z} m: ${transmitted(z)}, not ${share}`);
  }
  assert.equal(transmitted(PENETRATION_DEPTH), 0);
  assert.equal(transmitted(4000), 0);
  assert.ok(Math.abs(transmitted(0) - J(0)) < 1e-15 && J(0) < 1 && J(0) > 1 - 1e-8);
});

test('each receiver takes the light absorbed over its depths, a thin class passes its part to the receiver below, the last takes all that passes it, and none reaches past the floor or a column without a receiver', () => {
  const L = 8;
  for (const h0 of [10, 50]) {
    const h = Float64Array.from([h0, 5.5, 3, 20, 40, 0.01, 100, 1000]), depth = h.reduce((s, x) => s + x, 0);
    const share = new Float64Array(L), taken = penetrationShares(h, L, depth, 'IB', share);
    assert.equal(taken, transmitted(h0));
    const top = (k) => h.slice(0, k).reduce((s, x) => s + x, 0);
    const expected = [0, J(top(1)) - J(top(2)), 0, J(top(2)) - J(top(4)), J(top(4)) - J(top(5)), 0, J(top(5)) - J(top(7)), J(top(7))].map((x) => x / taken);
    for (let k = 0; k < L; k++) assert.ok(Math.abs(share[k] - expected[k]) < 1e-14, `h0 ${h0} m, class ${k}: ${share[k]} against ${expected[k]}`);
    assert.ok(Math.abs(share.reduce((s, x) => s + x, 0) - 1) < 1e-15, `h0 ${h0} m: the shares sum to ${share.reduce((s, x) => s + x, 0)}`);
    console.log(`h0 ${h0} m: ${(100 * taken).toFixed(2)} % passes; shares ${Array.from(share, (x) => x.toFixed(4)).join(', ')}`);
  }
  const shelf = Float64Array.from([10, 20, 30]), share = new Float64Array(3);
  assert.equal(penetrationShares(shelf, 3, 60, 'IB', share), transmitted(10) - transmitted(60), 'the light that would pass the floor stays in the mixed layer');
  assert.ok(Math.abs(share[1] - (J(10) - J(30)) / (J(10) - J(60))) < 1e-14 && Math.abs(share[2] - (J(30) - J(60)) / (J(10) - J(60))) < 1e-14);
  for (const column of [[50, 3, 0.01, 0.01], [50, THIN, 2]]) {
    const h = Float64Array.from(column), out = new Float64Array(h.length).fill(7);
    assert.equal(penetrationShares(h, h.length, h.reduce((s, x) => s + x, 0), 'IB', out), 0, `${column}: no receiver, no light`);
    assert.ok(out.every((x) => x === 0));
  }
});

/*
 * An aquaplanet's CPU model at N = 8 from initializeState with its ocean
 * started under it, the sun up over part of the planet.
 */
function started(ocean) {
  const model = createModel(new Grid(8), { ocean });
  initializeState(model).forEach((values, s) => model.state[s].set(values));
  model.ocean.initialize(model.state[3], model.state[6]);
  model.time = 3 * 3600;
  return model;
}

test('with penetration the CPU model\'s step leaves the air and the net surface flux as they are without it, and keeps out of each open-water temperature exactly the light it adds to the pending store', () => {
  const plain = started({ penetration: null }), lit = started({ penetration: 'IB' }), C = plain.mesh.nCells, dt = 900;
  assert.ok(plain.ocean.transmission.every((x) => x === 0) && lit.ocean.transmission.some((x) => x > 0.005));
  const surfaceT = Float64Array.from(lit.state[3]), ice = Float64Array.from(lit.state[6]);
  plain.step(dt); lit.step(dt);
  for (const a of [0, 1, 2, 4, 5, 6]) assert.deepEqual(lit.state[a], plain.state[a], `state ${a}`);
  assert.deepEqual(lit.radiation.surfaceFlux, plain.radiation.surfaceFlux);
  assert.deepEqual(lit.totals, plain.totals);
  let open = 0, passed = 0, worst = 0;
  for (let i = 0; i < C; i++) {
    const pen = dt * (lit.radiation.surfaceSolar[i] * lit.ocean.transmission[i]);
    if (ice[i] > 0 || plain.state[6][i] > 0) {
      assert.equal(lit.ocean.penetrated[i], 0, `cell ${i} under ice`);
      assert.equal(lit.state[3][i], plain.state[3][i]);
      continue;
    }
    open++;
    passed += pen;
    assert.equal(lit.ocean.penetrated[i], pen);
    const kept = lit.ocean.capacity[i] * (plain.state[3][i] - lit.state[3][i]);
    worst = Math.max(worst, Math.abs(kept - pen));
    assert.ok(Math.abs(kept - pen) < 4 * Number.EPSILON * lit.ocean.capacity[i] * surfaceT[i], `cell ${i}: the temperature kept out ${kept} J/m², the store took ${pen}`);
  }
  console.log(`${open} open-water cells passed ${(passed / open / dt).toFixed(2)} W/m² on average below the mixed layer; the temperature kept it out to ${worst.toExponential(1)} J/m²`);
  assert.ok(passed > 0);
});

test('over one ocean interval of the CPU model the pending store holds exactly the light that passed the mixed layer, and the next surface read deposits all of it below the layer by depth', () => {
  const model = started({ penetration: 'IB' }), { ocean, mesh, radiation, state } = model, C = mesh.nCells, L = ocean.layers, dt = 900;
  const passed = new Float64Array(C);
  for (let n = 1; n < 4; n++) {
    const open = Uint8Array.from(state[6], (h) => (h > 0 ? 0 : 1));
    model.step(dt);
    for (let i = 0; i < C; i++) if (open[i]) passed[i] += dt * (radiation.surfaceSolar[i] * ocean.transmission[i]);
  }
  assert.deepEqual(ocean.penetrated, passed, 'the pending store holds the light passed');
  const pending = Float64Array.from(ocean.penetrated), Q = Float64Array.from(ocean.Q), share = new Float64Array(L), column = new Float64Array(L);
  ocean.readSurface(state[3], state[6]);
  let cells = 0, total = 0, off = 0;
  for (let i = 0; i < C; i++) {
    if (!(pending[i] > 0)) continue;
    cells++; total += mesh.areaCell[i] * pending[i];
    for (let k = 0; k < L; k++) column[k] = ocean.h[k * C + i];
    penetrationShares(column, L, ocean.D[i] + ocean.eta[i], 'IB', share);
    let deposited = 0, rounding = 0;
    for (let k = 1; k < L; k++) {
      const dq = ocean.Q[k * C + i] - Q[k * C + i], grain = ocean.rhoCp * 2 * Number.EPSILON * Math.abs(Q[k * C + i]);
      deposited += ocean.rhoCp * dq; rounding += grain;
      assert.ok(Math.abs(ocean.rhoCp * dq - pending[i] * share[k]) <= 1e-12 * pending[i] + grain, `cell ${i} class ${k}: ${ocean.rhoCp * dq} J/m² deposited, ${pending[i] * share[k]} its share`);
    }
    off = Math.max(off, Math.abs(deposited - pending[i]) / (1e-12 * pending[i] + rounding));
    assert.equal(ocean.penetrated[i], 0);
  }
  console.log(`${cells} cells deposited ${(total / mesh.areaCell.reduce((s, a) => s + a, 0) / (3 * dt)).toFixed(3)} W/m² of the planet's mean below their mixed layers, each within ${off.toFixed(2)} of the rounding of its stores`);
  assert.ok(cells > 0 && off <= 1, `the deposit is off the pending light by ${off} of its stores' rounding`);
});

test('a CPU ocean restarted from a state without the pending light loads with none, and one from a state with it takes it back', () => {
  const model = started({ penetration: 'IB' }), { ocean, mesh, state } = model, C = mesh.nCells;
  const restart = { ...ocean.serialize(), Q: ocean.Q, W: ocean.W, previousT0: ocean.T0, previousIceVolume: new Float64Array(C), fresh: new Float64Array(C), buoyancy: new Float64Array(C), flux: new Float64Array(C), capacity: ocean.capacity };
  const pending = Float64Array.from({ length: C }, (_, i) => (ocean.cellOcean[i] ? 1000 + i : 0));
  const reloaded = createOcean(mesh, { penetration: 'IB' });
  reloaded.load({ ...restart, penetrated: pending }, Float64Array.from(state[3]), state[6]);
  assert.deepEqual(reloaded.penetrated, pending);
  reloaded.load(restart, Float64Array.from(state[3]), state[6]);
  assert.ok(reloaded.penetrated.every((x) => x === 0));
  for (let i = 0; i < C; i++) assert.ok(Math.abs(reloaded.transmission[i] - ocean.transmission[i]) < 1e-12, `cell ${i}'s transmission after the load`);
});

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { createForcedOcean } = gpuAvailable ? await import('../js/gpu/forcing.gpu.js') : {};
const { readRanges } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};
const LIT = { everySteps: 1, penetration: 'IB' };
const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
const bits = new Float32Array(1), word = new Uint32Array(bits.buffer);
const spacing = (x) => { bits[0] = x; return x === 0 ? 2 ** -24 : 2 ** Math.max(((word[0] >>> 23) & 255) - 150, -126); };
const physicsArrays = (model, names) => readRanges(model.gpu.device, model.gpu.buffers.PH, names.map((name) => ({ offset: model.gpu.layout.PH[name], length: model.mesh.nCells })));
const rms = (a, b, mask = null) => { let d = 0, r = 0; for (let x = 0; x < a.length; x++) if (!mask || mask[x]) { d += (a[x] - b[x]) ** 2; r += a[x] * a[x]; } return Math.sqrt(d / Math.max(r, 1e-300)); };

test('at N=8 the GPU ocean publishes the CPU\'s transmission to 1e-6, deposits the same passed light class by class to its stores\' spacing, and tracks the CPU through one and twenty-one steps that each pass it light', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cpuModel = createModel(new Grid(8), { topography });
  initializeState(cpuModel, {}).forEach((values, a) => cpuModel.state[a].set(values));
  const { mesh, geography } = cpuModel, C = mesh.nCells, dt = 1350;
  for (let i = 0; i < C; i++) if (geography.land[i]) cpuModel.state[6][i] = 0;
  const surfaceT0 = Float64Array.from(cpuModel.state[3]), ice = Float64Array.from(cpuModel.state[6]), stress = new Float64Array(mesh.nEdges);
  for (let e = 0; e < mesh.nEdges; e++) { const lon = Math.atan2(mesh.xEdge[3 * e + 1], mesh.xEdge[3 * e]); stress[e] = 0.08 * Math.sin(2 * mesh.latEdge[e]) * Math.cos(lon) + 0.01 * Math.sin(5 * lon); }
  const cpu = createOcean(mesh, { geography, ...LIT }), gpuModel = await createGpuModel(new Grid(8), { topography, ocean: LIT }), gpu = gpuModel.oceanEngine, L = cpu.layers, PH = gpuModel.gpu.layout.PH;
  cpu.initialize(Float64Array.from(surfaceT0), ice); gpu.initialize(Float64Array.from(surfaceT0), ice);
  const transmissionGap = async () => {
    const [[published], { h, eta }] = await Promise.all([physicsArrays(gpuModel, ['TRANS']), gpu.download()]);
    const column = new Float64Array(L);
    let engines = 0, own = 0, largest = 0;
    for (let i = 0; i < C; i++) {
      if (!gpu.cellOcean[i]) continue;
      for (let k = 0; k < L; k++) column[k] = h[k * C + i];
      const exact = penetrationShares(column, L, gpu.D[i] + eta[i]);
      largest = Math.max(largest, cpu.transmission[i]);
      if (cpu.transmission[i] > 0) engines = Math.max(engines, Math.abs(published[i] - cpu.transmission[i]) / cpu.transmission[i]); else assert.equal(published[i], 0);
      if (exact > 0) own = Math.max(own, Math.abs(published[i] - exact) / exact);
    }
    return { engines, own, largest };
  };
  const started = await transmissionGap();
  const flux = new Float64Array(C), surfaceT = Float64Array.from(surfaceT0);
  const light = Float32Array.from(cpu.transmission, (t) => 250 * dt * t);
  const pass = () => { cpu.penetrated.set(light); gpuModel.gpu.device.queue.writeBuffer(gpuModel.gpu.buffers.PH, 4 * PH.PEN, light); };
  const heat = async () => { const [q] = await readRanges(gpuModel.gpu.device, gpu.buffers.S, [{ offset: gpu.layout.OS.OQ, length: L * C }]); return Float64Array.from(q); };
  const cpuBefore = Float64Array.from(cpu.Q), gpuBefore = await heat();
  pass();
  cpu.readSurface(Float64Array.from(surfaceT), ice);
  gpu.readSurface(Float64Array.from(surfaceT), ice);
  const gpuAfter = await heat(), [waiting] = await physicsArrays(gpuModel, ['PEN']);
  let classes = 0, coarse = 0, fine = 0, worstFine = 0, worstClass = 0, worstColumn = 0;
  for (let i = 0; i < C; i++) {
    if (!(light[i] > 0)) continue;
    let change = 0, allowance = 0;
    for (let k = 1, z = cpu.h[i]; k < L && z < PENETRATION_DEPTH + cpu.h[k * C + i]; z += cpu.h[k * C + i], k++) if (cpu.h[k * C + i] > THIN) allowance = Math.max(allowance, spacing(gpuBefore[k * C + i]));
    for (let k = 1; k < L; k++) {
      const n = k * C + i, dc = cpu.Q[n] - cpuBefore[n], dg = gpuAfter[n] - gpuBefore[n];
      change += dg;
      if (dc === 0 && dg === 0) continue;
      classes++;
      if (Math.abs(dg - dc) > 1e-5 * Math.abs(dc)) coarse++;
      worstClass = Math.max(worstClass, Math.abs(dg - dc) / (1e-5 * Math.abs(dc) + allowance));
    }
    worstColumn = Math.max(worstColumn, Math.abs(cpu.rhoCp * change + waiting[i] - light[i]) / light[i]);
  }
  console.log(`one deposit of 250 W/m²'s passing light on the same columns: ${classes} class changes, ${coarse} beyond 1e-5 of the CPU's, all within ${worstClass.toFixed(2)} of 1e-5 plus the coarsest spacing of the column's receivers; each GPU column's stores and waiting light within ${worstColumn.toExponential(2)} of its light`);
  assert.ok(classes > 0 && worstClass <= 1 && worstColumn < 1e-6);

  const strong = Float32Array.from(light, (x) => 2e4 * x), heatBefore = await heat(), cpuHeld = Float64Array.from(cpu.Q);
  cpu.penetrated.set(strong); gpuModel.gpu.device.queue.writeBuffer(gpuModel.gpu.buffers.PH, 4 * PH.PEN, strong);
  cpu.readSurface(Float64Array.from(surfaceT), ice);
  gpu.readSurface(Float64Array.from(surfaceT), ice);
  const heatAfter = await heat();
  for (let i = 0; i < C; i++) for (let k = 1; k < L; k++) {
    const n = k * C + i, dc = cpu.Q[n] - cpuHeld[n], dg = heatAfter[n] - heatBefore[n];
    if (Math.abs(dc) > 1e5 * spacing(heatBefore[n])) { fine++; worstFine = Math.max(worstFine, Math.abs(dg - dc) / Math.abs(dc)); }
  }
  console.log(`twenty thousand times that light, ${fine} class changes of more than 1e5 spacings: within ${worstFine.toExponential(2)} of the CPU's`);
  assert.ok(fine > 100 && worstFine < 1e-5);
  cpu.initialize(Float64Array.from(surfaceT0), ice); gpu.initialize(Float64Array.from(surfaceT0), ice);

  cpu.advance(surfaceT, ice, flux, stress, dt);
  await gpu.advance(Float64Array.from(surfaceT0), ice, stress, dt);
  const stepped = await transmissionGap();
  console.log(`the transmission, up to ${stepped.largest.toFixed(4)}: uploaded within ${started.engines.toExponential(2)} of the CPU's; after a step the WGSL's within ${stepped.own.toExponential(2)} of the double-precision value for the GPU's own columns and ${stepped.engines.toExponential(2)} of the CPU's`);
  assert.ok(started.largest > 0.005 && started.engines < 1e-6 && stepped.own < 1e-6 && stepped.engines < 1e-4);

  for (let n = 0; n < 21; n++) {
    pass();
    surfaceT.set(surfaceT0);
    cpu.advance(surfaceT, ice, flux, stress, dt);
    await gpu.advance(Float64Array.from(surfaceT0), ice, stress, dt);
    if (n === 0 || n === 20) {
      const c = cpu.serialize(), g = await gpu.serialize(), wet = Array.from(c.h, (v) => v > 1);
      const gaps = ['h', 'u', 'T', 'S', 'eta'].map((f) => [f, rms(c[f], g[f], f === 'T' || f === 'S' ? wet : null)]);
      console.log(`${n + 2} steps passing light from the second: rms relative ${gaps.map(([f, r]) => `${f} ${r.toExponential(2)}`).join(', ')}`);
      for (const [f, r] of gaps) assert.ok(r < 2e-3, `${f} rms relative diff ${r} after ${n + 2} steps`);
    }
  }
  gpuModel.destroy();
});

test('on the GPU open water keeps out of its temperature the passing share of the sunlight it absorbs and adds it to PH_PEN, under ice none passes, and the net flux is the same as without penetration', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const build = async (penetration) => {
    const model = await createGpuModel(new Grid(8), { ocean: { penetration } });
    initializeState(model).forEach((values, s) => model.state[s].set(values));
    model.time = 3 * 3600;
    model.load({ ocean: false });
    model.ocean.initialize(model.state[3], model.state[6]);
    return model;
  };
  const plain = await build(null), lit = await build('IB'), C = lit.mesh.nCells, dt = 900;
  const ice = Float64Array.from(lit.state[6]);
  const [transmission, capacity] = await physicsArrays(lit, ['TRANS', 'CAP']);
  await plain.step(dt); await lit.step(dt);
  await plain.sync(); await lit.sync();
  const [[net, absorbed, atmosphere, pending], [plainNet, plainPending]] = await Promise.all([physicsArrays(lit, ['SFLUX', 'ABS', 'ATMSW', 'PEN']), physicsArrays(plain, ['SFLUX', 'PEN'])]);
  assert.deepEqual(net, plainNet, 'the net surface flux');
  for (const a of [0, 1, 2, 4, 5]) assert.deepEqual(lit.state[a], plain.state[a], `state ${a}`);
  assert.ok(plainPending.every((x) => x === 0));
  let open = 0, iced = 0, worst = 0, kept = 0;
  for (let i = 0; i < C; i++) {
    if (ice[i] > 0 || lit.state[6][i] > 0) { iced++; assert.equal(pending[i], 0, `cell ${i} under ice`); continue; }
    const pen = dt * (absorbed[i] - atmosphere[i]) * transmission[i];
    if (!(pen > 0)) continue;
    open++;
    worst = Math.max(worst, Math.abs(pending[i] - pen) / pen);
    kept = Math.max(kept, Math.abs(capacity[i] * (plain.state[3][i] - lit.state[3][i]) - pending[i]) / (capacity[i] * spacing(lit.state[3][i])));
  }
  console.log(`${open} sunlit open-water cells passed their share of the surface's absorbed sunlight within ${worst.toExponential(2)}, kept out of their temperatures within ${kept.toFixed(2)} of the spacing of a float32 temperature; ${iced} iced cells passed none`);
  assert.ok(open > 0 && iced > 0 && worst < 1e-5 && kept <= 2);
  plain.destroy(); lit.destroy();
});

test('the forced ocean passes the recorded absorbed shortwave\'s share below the mixed layer of open water, not the restoring, and deposits it at its ocean step', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const model = await createGpuModel(new Grid(8), { ocean: { penetration: 'IB' } }), C = model.mesh.nCells, E = model.mesh.nEdges, dt = 900, ocean = model.oceanEngine;
  initializeState(model).forEach((values, s) => model.state[s].set(values));
  model.load({ ocean: false });
  model.ocean.initialize(model.state[3], model.state[6]);
  const ice = Float64Array.from(model.state[6]), sst = Float64Array.from(model.state[3], (t) => t + 2);
  const forced = createForcedOcean(model), cells = (x) => new Float32Array(C).fill(x);
  forced.setDay({ netFlux: cells(80), shortwave: cells(200), shortwaveDown: cells(230), sst, snowfall: cells(0), stress: new Float32Array(E), evaporation: cells(0), rain: cells(0), runoff: cells(0) }, ocean.everySteps * dt);
  const [transmission] = await physicsArrays(model, ['TRANS']);
  forced.step(dt, 30);
  const [pending] = await physicsArrays(model, ['PEN']);
  let open = 0, worst = 0;
  for (let i = 0; i < C; i++) {
    if (ice[i] > 0) { assert.equal(pending[i], 0, `cell ${i} under ice`); continue; }
    open++;
    worst = Math.max(worst, Math.abs(pending[i] - dt * 200 * transmission[i]) / (dt * 200 * transmission[i]));
  }
  let steps = 1;
  while (!forced.step(dt, 30)) steps++;
  steps++;
  const [[left], [h, q]] = await Promise.all([physicsArrays(model, ['PEN']), readRanges(model.gpu.device, ocean.buffers.S, [{ offset: ocean.layout.OS.OH, length: ocean.layers * C }, { offset: ocean.layout.OS.OQ, length: ocean.layers * C }])]);
  const rhoCp = ocean.options.density * ocean.options.specificHeat;
  let before = 0, after = 0, most = 0, beyond = 0;
  for (let i = 0; i < C; i++) {
    if (ice[i] > 0) continue;
    const given = dt * 200 * transmission[i] * steps;
    before += model.mesh.areaCell[i] * given; after += model.mesh.areaCell[i] * Math.abs(left[i]); most = Math.max(most, Math.abs(left[i]) / given);
    let coarsest = 0;
    for (let k = 1, z = h[i]; k < ocean.layers && z < PENETRATION_DEPTH + h[k * C + i]; z += h[k * C + i], k++) if (h[k * C + i] > THIN) coarsest = Math.max(coarsest, spacing(q[k * C + i]));
    if (Math.abs(left[i]) > rhoCp * coarsest) beyond++;
  }
  console.log(`${open} open-water cells passed the recorded shortwave's share within ${worst.toExponential(2)}; the ocean step after ${steps} steps deposited it, leaving ${(after / before).toExponential(2)} of it waiting, ${(100 * most).toFixed(2)} % in the cell that kept the most, and ${beyond} cells more than a spacing of their coarsest receiver`);
  assert.ok(open > 0 && worst < 1e-6 && after < 0.1 * before && beyond === 0);
  model.destroy();
});

test('a GPU ocean restarted from a state without the pending light or the published transmission loads as restartable with none waiting and its transmission recomputed, and one with them takes them back word for word', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const model = await createGpuModel(new Grid(8), { ocean: { penetration: 'IB' } }), C = model.mesh.nCells;
  initializeState(model).forEach((values, s) => model.state[s].set(values));
  model.load({ ocean: false });
  model.ocean.initialize(model.state[3], model.state[6]);
  for (let n = 0; n < 6; n++) await model.step(900);
  await model.sync();
  const saved = await model.ocean.serialize({ restart: true });
  assert.ok(saved.penetrated.some((x) => x > 0) && saved.transmission.some((x) => x > 0), 'light is waiting and the share is published');
  const fresh = async (file) => {
    const other = await createGpuModel(new Grid(8), { ocean: { penetration: 'IB' } });
    other.state.forEach((a, s) => a.set(model.state[s]));
    other.load({ ocean: false });
    other.ocean.load(file, other.state[3], other.state[6]);
    const back = await other.ocean.serialize({ restart: true });
    other.destroy();
    return back;
  };
  const again = await fresh(saved);
  for (const name of ['penetrated', 'transmission', 'capacity', 'flux', 'previousT0']) assert.deepEqual(again[name], saved[name], name);
  const { penetrated, transmission, ...legacy } = saved;
  const old = await fresh(legacy);
  assert.ok(old.penetrated.every((x) => x === 0), 'nothing waits after a legacy load');
  for (const name of ['capacity', 'flux', 'previousT0', 'buoyancy']) assert.deepEqual(old[name], saved[name], `${name}: the legacy state is still restartable`);
  let worst = 0;
  for (let i = 0; i < C; i++) if (transmission[i] > 0) worst = Math.max(worst, Math.abs(old.transmission[i] - transmission[i]) / transmission[i]);
  console.log(`a legacy restart recomputes the transmission within ${worst.toExponential(2)} of the published one`);
  assert.ok(worst < 1e-5);
  model.destroy();
});
