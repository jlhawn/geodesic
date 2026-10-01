import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { topographyFromInt16, createGeography } from '../js/geography.module.js';
import { createOcean, atlasColumns, EPS, THIN, RESTORE_TOLERANCE, LAYER_DENSITIES, LAYER_SALINITIES } from '../js/ocean/layered.module.js';
import { decodeClimatology, encodeClimatology, loadClimatology, profileAt, CLIMATOLOGY_FILE } from '../js/ocean/climatology.module.js';
import { seawaterDensity, labelTemperature } from '../js/ocean/seawater.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const DEG = Math.PI / 180, KELVIN = 273.15;
const FILE = new URL(`../${CLIMATOLOGY_FILE}`, import.meta.url);
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);

/*
 * An 8 × 4 grid at 45° spacing over three depths, with temperature and
 * salinity linear in depth, latitude and longitude: no third level at
 * row 2, columns 3, 6 and 7, an all-land column at row 1, column 5, and
 * row 3 holding only the surface level.
 */
const GRID = { nLon: 8, nLat: 4, lon0: -157.5, dLon: 45, lat0: -67.5, dLat: 45, depths: [0, 100, 1000] };
const fixtureT = (z, lat, lon) => 20 - 0.01 * z + 0.1 * lat + 0.02 * lon;
const fixtureS = (z, lat) => 35 + 0.001 * z - 0.01 * lat;
function fixture() {
  const { nLon, nLat, lon0, dLon, lat0, dLat, depths } = GRID, T = [], S = [];
  depths.forEach((z, j) => {
    for (let r = 0; r < nLat; r++) for (let c = 0; c < nLon; c++) {
      const lat = lat0 + r * dLat, lon = lon0 + c * dLon;
      const missing = (r === 1 && c === 5) || (r === 2 && [3, 6, 7].includes(c) && j === 2) || (r === 3 && j > 0);
      T.push(missing ? NaN : fixtureT(z, lat, lon)); S.push(missing ? NaN : fixtureS(z, lat));
    }
  });
  return { ...GRID, T, S, source: 'fixture' };
}

test('a climatology packed from arrays decodes to the same grid, values within half a step and missing points kept', () => {
  const grid = fixture();
  const atlas = decodeClimatology(encodeClimatology(grid));
  for (const key of ['nLon', 'nLat', 'lon0', 'dLon', 'lat0', 'dLat', 'source']) assert.equal(atlas[key], grid[key], key);
  assert.deepEqual(Array.from(atlas.depths), grid.depths);
  for (let n = 0; n < grid.T.length; n++) {
    if (Number.isNaN(grid.T[n])) { assert.equal(atlas.T[n], atlas.missing); assert.equal(atlas.S[n], atlas.missing); continue; }
    assert.ok(Math.abs(atlas.tOffset + atlas.tScale * atlas.T[n] - grid.T[n]) <= 0.5 * atlas.tScale + 1e-9, `T ${n}`);
    assert.ok(Math.abs(atlas.sOffset + atlas.sScale * atlas.S[n] - grid.S[n]) <= 0.5 * atlas.sScale + 1e-9, `S ${n}`);
  }
  const node = atlas.columnAt(-22.5 * DEG, 22.5 * DEG);
  assert.equal(node.depths.length, 3);
  for (let j = 0; j < 3; j++) assert.ok(Math.abs(node.T[j] - KELVIN - fixtureT(GRID.depths[j], -22.5, 22.5)) < 1e-3 && Math.abs(node.S[j] - fixtureS(GRID.depths[j], -22.5)) < 1e-3, `level ${j} at a grid point`);
});

test('the repository file re-encodes byte for byte, so scripts/packWoa.py and encodeClimatology share one layout', () => {
  const bytes = new Uint8Array(readFileSync(FILE));
  const atlas = decodeClimatology(bytes);
  const unpack = (values, offset, scale) => Array.from(values, (v) => (v === atlas.missing ? NaN : offset + scale * v));
  const again = encodeClimatology({ ...atlas, depths: Array.from(atlas.depths), T: unpack(atlas.T, atlas.tOffset, atlas.tScale), S: unpack(atlas.S, atlas.sOffset, atlas.sScale) });
  assert.equal(again.length, bytes.length);
  assert.ok(again.every((b, n) => b === bytes[n]), 'identical bytes');
  assert.equal(atlas.nLon, 360); assert.equal(atlas.nLat, 180); assert.equal(atlas.depths.length, 25);
  assert.ok(bytes.length <= 8 * 1024 * 1024, `${bytes.length} bytes`);
});

test('columnAt interpolates bilinearly among four wet points, takes the nearest wet one where a corner is dry, ends where none is wet and reaches for the nearest column over land', () => {
  const atlas = decodeClimatology(encodeClimatology(fixture()));
  const inside = atlas.columnAt(-40 * DEG, -100 * DEG);
  assert.equal(inside.depths.length, 3);
  for (let j = 0; j < 3; j++) assert.ok(Math.abs(inside.T[j] - KELVIN - fixtureT(GRID.depths[j], -40, -100)) < 2e-3, `bilinear T at level ${j}: ${inside.T[j] - KELVIN}`);
  const nearDry = atlas.columnAt(10 * DEG, -30 * DEG);
  assert.equal(nearDry.depths.length, 3, 'three of the four corners reach the third level');
  assert.ok(Math.abs(nearDry.T[1] - KELVIN - fixtureT(100, 10, -30)) < 2e-3, 'bilinear where all four are wet');
  assert.ok(Math.abs(nearDry.T[2] - KELVIN - fixtureT(1000, -22.5, -22.5)) < 2e-3, `the nearest wet corner at the third level: ${nearDry.T[2] - KELVIN}`);
  const shallow = atlas.columnAt(40 * DEG, 130 * DEG);
  assert.equal(shallow.depths.length, 2, 'the column ends at the first level none of the corners reaches');
  assert.ok(Math.abs(shallow.T[1] - KELVIN - fixtureT(100, 22.5, 112.5)) < 2e-3, 'the second level from the nearest corner that has it');
  const wrapped = atlas.columnAt(-22.5 * DEG, 180 * DEG);
  const east = atlas.columnAt(-22.5 * DEG, 157.5 * DEG), west = atlas.columnAt(-22.5 * DEG, -157.5 * DEG);
  assert.ok(Math.abs(wrapped.T[1] - 0.5 * (east.T[1] + west.T[1])) < 1e-9, 'across the dateline the columns on either side are averaged');
  const island = fixture(), dryCorners = [[1, 5], [1, 6], [2, 5], [2, 6]];
  for (let n = 0; n < island.T.length; n++) if (dryCorners.some(([r, c]) => n % 32 === 8 * r + c)) island.T[n] = NaN;
  const packed = encodeClimatology(island), point = [-10, 90];
  let nearest = null, distance = Infinity;
  for (let r = 0; r < 4; r++) for (let c = 0; c < 8; c++) {
    if (!Number.isFinite(island.T[8 * r + c])) continue;
    const la = (GRID.lat0 + 45 * r) * DEG, lo = (GRID.lon0 + 45 * c) * DEG, pla = point[0] * DEG, plo = point[1] * DEG;
    const d = Math.acos(Math.sin(la) * Math.sin(pla) + Math.cos(la) * Math.cos(pla) * Math.cos(lo - plo)) / DEG;
    if (d < distance) { distance = d; nearest = [r, c]; }
  }
  const reached = decodeClimatology(packed, { reach: distance + 1 }).columnAt(point[0] * DEG, point[1] * DEG);
  assert.ok(reached, 'a point over land takes a column within reach');
  assert.ok(Math.abs(reached.T[0] - KELVIN - island.T[8 * nearest[0] + nearest[1]]) < 1e-3, 'the nearest wet column, whole');
  assert.equal(reached.depths.length, nearest[0] === 3 ? 1 : 3);
  assert.equal(decodeClimatology(packed, { reach: distance - 1 }).columnAt(point[0] * DEG, point[1] * DEG), null, 'and none beyond it');
  const [t50, s50] = profileAt(inside, 50), [tDeep] = profileAt(inside, 3000), [tTop] = profileAt(inside, -5);
  assert.ok(Math.abs(t50 - 0.5 * (inside.T[0] + inside.T[1])) < 1e-9 && Math.abs(s50 - 0.5 * (inside.S[0] + inside.S[1])) < 1e-9, 'linear between levels');
  assert.equal(tDeep, inside.T[2]); assert.equal(tTop, inside.T[0]);
});

test('a sea cell under an atlas column that ends far above its bottom takes the levels below from the nearest column that reaches them, and keeps the plug when none is within reach', () => {
  const depths = [0, 50, 75, 100, 200, 500, 1000, 1500, 2000], nLon = 6, nLat = 3, T = [], S = [];
  const below = (z) => Math.exp(-Math.max(0, z - 50) / 150), profileT = (z) => 5 + 24.5 * below(z), profileS = (z) => 34.6 - 0.6 * below(z);
  depths.forEach((z, j) => {
    for (let r = 0; r < nLat; r++) for (let c = 0; c < nLon; c++) { const dry = c < 3 && j > 1; T.push(dry ? NaN : profileT(z)); S.push(dry ? NaN : profileS(z)); }
  });
  const grid = { nLon, nLat, lon0: 0, dLon: 1, lat0: -1, dLat: 1, depths, T, S, source: 'fixture' };
  const rho = [1025, ...LAYER_DENSITIES], labelS = [35, ...LAYER_SALINITIES], labelT = rho.map((r, k) => Math.max(FREEZING_POINT, labelTemperature(r, labelS[k]))), L = rho.length;
  const column = (atlas) => {
    const h = new Float64Array(L), Q = new Float64Array(L), W = new Float64Array(L), T0 = new Float64Array(1);
    atlasColumns({ nCells: 1, latCell: [0], lonCell: [1 * DEG] }, atlas, { D: [1200], cellOcean: [1], ice: [0], rho, labelT, labelS, h, Q, W, T0 });
    const layers = [];
    let z = 0;
    for (let k = 0; k < L; k++) { if (k === 0 || h[k] > THIN) layers.push({ rho: k ? rho[k] : 'ML', top: z, bottom: z + h[k] }); z += h[k]; }
    return layers;
  };
  const reached = column(decodeClimatology(encodeClimatology(grid))), plugged = column(decodeClimatology(encodeClimatology(grid), { deepReach: 1 }));
  const show = (layers) => layers.map((l) => `${l.rho}: ${l.top.toFixed(0)}–${l.bottom.toFixed(0)}`).join(', ');
  console.log(`a 1200 m cell under a 50 m atlas column 2° from a 2000 m one: ${show(reached)}; with no deep column in reach ${show(plugged)}`);
  assert.ok(reached.filter((l) => l.rho !== 'ML' && l.rho > 1025 && l.bottom > 100).length >= 3, 'classes denser than 1025 below 100 m');
  assert.ok(!reached.some((l) => l.rho === 1020.5 && l.bottom > 75), 'no 1020.5 water below 75 m');
  assert.ok(plugged.some((l) => l.rho === 1020.5 && l.bottom > 1000), 'without a deep column within reach the surface water fills the column');
});

let built = null;
async function atlasStart() {
  if (built) return built;
  const mesh = buildMesh(new Grid(32)), C = mesh.nCells;
  const geography = createGeography(mesh, topography);
  const atlas = await loadClimatology(FILE.pathname);
  const ocean = createOcean(mesh, { geography, climatology: atlas });
  const ice = Float64Array.from(mesh.latCell, (lat, i) => (geography.land[i] ? 0 : lat > 72 * DEG ? 1.5 : lat < -68 * DEG ? 0.7 : 0));
  const before = Float64Array.from(mesh.latCell, (lat) => 272 + 30 * Math.cos(lat) ** 2), surfaceT = Float64Array.from(before);
  const started = ocean.initialize(surfaceT, ice);
  built = { mesh, C, geography, atlas, ocean, ice, before, surfaceT, started };
  return built;
}
const nearestSea = ({ mesh, ocean }, lat, lon) => {
  const x = [Math.cos(lat * DEG) * Math.cos(lon * DEG), Math.cos(lat * DEG) * Math.sin(lon * DEG), Math.sin(lat * DEG)];
  let best = -1, dot = -2;
  for (let i = 0; i < mesh.nCells; i++) {
    if (!ocean.cellOcean[i]) continue;
    const r = Math.hypot(mesh.xCell[3 * i], mesh.xCell[3 * i + 1], mesh.xCell[3 * i + 2]);
    const d = (x[0] * mesh.xCell[3 * i] + x[1] * mesh.xCell[3 * i + 1] + x[2] * mesh.xCell[3 * i + 2]) / r;
    if (d > dot) { dot = d; best = i; }
  }
  return best;
};
const classTop = ({ ocean, C }, label, i) => { const k = ocean.densities.indexOf(label, 1); let z = 0; for (let j = 0; j < k; j++) z += ocean.h[j * C + i]; return z; };

test('an N=32 start from the World Ocean Atlas: columns sum to the bathymetry, every class on its label or a token, stable, finite, and the sea surface taken from the mixed layer', async () => {
  const { mesh, C, geography, ocean, ice, before, surfaceT, started } = await atlasStart();
  const L = ocean.layers, rho = ocean.densities;
  const densest = rho.map((r, k) => r + 0.5 * (k < L - 1 ? rho[k + 1] - r : r - rho[k - 1]));
  let sea = 0, unstable = 0;
  for (let i = 0; i < C; i++) {
    if (!ocean.cellOcean[i]) { assert.equal(surfaceT[i], before[i], 'land keeps its surface'); continue; }
    sea++;
    let sum = 0;
    for (let k = 0; k < L; k++) {
      const n = k * C + i;
      assert.ok(Number.isFinite(ocean.h[n]) && Number.isFinite(ocean.Q[n]) && Number.isFinite(ocean.W[n]), `finite layer ${k} at cell ${i}`);
      sum += ocean.h[n];
      if (k === 0) continue;
      assert.ok(ocean.h[n] >= EPS - 1e-12, `layer ${k} at cell ${i} holds at least the token (${ocean.h[n]})`);
      if (ocean.h[n] === EPS) {
        assert.ok(Math.abs(seawaterDensity(ocean.Q[n] / EPS, ocean.W[n] / EPS) - rho[k]) < 1e-6, `a token of class ${k} at cell ${i} sits at its label`);
      } else if (ocean.h[n] > THIN) {
        assert.ok(Math.abs(seawaterDensity(ocean.Q[n] / ocean.h[n], ocean.W[n] / ocean.h[n]) - rho[k]) <= RESTORE_TOLERANCE + 1e-9, `class ${k} at cell ${i} within the restoring tolerance of its label`);
      }
    }
    assert.ok(Math.abs(sum - ocean.D[i] - ocean.eta[i]) < 1e-6, `column ${i} sums to its depth plus its sea level (${sum} vs ${ocean.D[i]} + ${ocean.eta[i]})`);
    assert.ok(Math.abs(ocean.eta[i]) < 3, `steric sea level at cell ${i} (${ocean.eta[i]} m)`);
    assert.ok(ocean.h[i] >= Math.min(50, ocean.D[i]) - 3 && ocean.h[i] <= 600, `mixed layer at cell ${i} (${ocean.h[i]} m)`);
    const t0 = ocean.Q[i] / ocean.h[i], rm = seawaterDensity(t0, ocean.W[i] / ocean.h[i]);
    if (ice[i] > 0) { assert.ok(Math.abs(t0 - FREEZING_POINT) < 1e-9, 'at the freezing point under ice'); assert.equal(surfaceT[i], before[i], 'the ice keeps its surface'); }
    else { assert.ok(t0 >= FREEZING_POINT - 1e-9 && t0 < 305, `mixed layer at cell ${i} is ${t0} K`); assert.ok(Math.abs(surfaceT[i] - t0) < 1e-9, 'the open sea surface is the mixed layer'); }
    let k = 1;
    while (k < L && ocean.h[k * C + i] <= THIN) k++;
    if (k < L && rm > densest[k] + 1e-9) {
      unstable++;
      assert.ok(ice[i] > 0 || k === L - 1, `cell ${i} at ${(mesh.latCell[i] / DEG).toFixed(0)}° has a mixed layer (${rm.toFixed(3)}) denser than class ${k} can hold, with neither ice nor only the densest class beneath`);
    }
  }
  assert.deepEqual(started, { atlas: sea, analytic: 0 });
  assert.ok(unstable < 0.01 * sea, `${unstable} of ${sea} mixed layers denser than the class beneath`);
  assert.ok(geography.landFraction.length === C);
});

test('the start places the 1024 class top near the 20 °C isotherm: 150–200 m in the west Pacific at 160E, 40–80 m in the east at 100W, outcropped poleward of 35°', async () => {
  const start = await atlasStart();
  const { mesh, C, ocean } = start;
  const west = classTop(start, 1024, nearestSea(start, 0, 160)), east = classTop(start, 1024, nearestSea(start, 0, -100));
  assert.ok(west >= 150 && west <= 200, `0N 160E: ${west.toFixed(0)} m`);
  assert.ok(east >= 40 && east <= 80, `0N 100W: ${east.toFixed(0)} m`);
  let n = 0, out = 0;
  for (let i = 0; i < C; i++) {
    const lat = Math.abs(mesh.latCell[i] / DEG);
    if (!ocean.cellOcean[i] || lat < 35 || lat > 70) continue;
    n++;
    if (classTop(start, 1024, i) <= ocean.h[i] + 1) out++;
  }
  assert.ok(out > 0.95 * n, `outcropped at ${out} of ${n} cells between 35° and 70°`);
});

test('in the subtropical gyres the 1025 class top lies deeper than 200 m (30S 90W, 30N 150W) where the 1024 class sits at the atlas\'s own 20 °C depth, near the mixed-layer base', async () => {
  const start = await atlasStart();
  const { atlas } = start;
  for (const [lat, lon] of [[-30, -90], [30, -150]]) {
    const i = nearestSea(start, lat, lon);
    const top1025 = classTop(start, 1025, i), top1024 = classTop(start, 1024, i);
    assert.ok(top1025 > 200 && top1025 < 260, `${lat} ${lon}: 1025 class top ${top1025.toFixed(0)} m`);
    const column = atlas.columnAt(start.mesh.latCell[i], start.mesh.lonCell[i]);
    let crossing = 0;
    for (let z = 0, running = -Infinity; z < 400; z += 1) { const [t, s] = profileAt(column, z); running = Math.max(running, seawaterDensity(t, s)); if (running >= 1023.875) { crossing = z; break; } }
    assert.ok(Math.abs(top1024 - Math.max(start.ocean.h[i], crossing)) < 3, `${lat} ${lon}: 1024 class top ${top1024.toFixed(0)} m against the atlas's ${crossing} m below a ${start.ocean.h[i].toFixed(0)} m mixed layer`);
  }
});

test('the Southern Ocean starts with Circumpolar Deep Water: 60–70S at 300–800 m between +0.5 and +2 °C, and every depth band within 0.5 K of the atlas at the same cells', async () => {
  const { mesh, C, atlas, ocean } = await atlasStart();
  const BANDS = [[0, 60], [60, 200], [200, 500], [500, 1000], [1000, 3000], [300, 800]];
  const model = BANDS.map(() => [0, 0]), observed = BANDS.map(() => [0, 0]);
  for (let i = 0; i < C; i++) {
    const lat = mesh.latCell[i] / DEG;
    if (!ocean.cellOcean[i] || lat < -70 || lat >= -60) continue;
    const a = mesh.areaCell[i], column = atlas.columnAt(mesh.latCell[i], mesh.lonCell[i]);
    let z = 0;
    for (let k = 0; k < ocean.layers; k++) {
      const hk = ocean.h[k * C + i], t = ocean.Q[k * C + i] / hk - KELVIN;
      BANDS.forEach(([z0, z1], b) => { const o = Math.max(0, Math.min(z + hk, z1) - Math.max(z, z0)); model[b][0] += o * a * t; model[b][1] += o * a; });
      z += hk;
    }
    BANDS.forEach(([z0, z1], b) => { for (let d = z0 + 5; d < Math.min(z1, ocean.D[i]); d += 10) { observed[b][0] += 10 * a * (profileAt(column, d)[0] - KELVIN); observed[b][1] += 10 * a; } });
  }
  const mean = ([q, w]) => q / w;
  const cdw = mean(model[5]);
  assert.ok(cdw > 0.5 && cdw < 2, `60–70S 300–800 m: ${cdw.toFixed(2)} °C`);
  BANDS.forEach(([z0, z1], b) => assert.ok(Math.abs(mean(model[b]) - mean(observed[b])) < 0.5, `60–70S ${z0}–${z1} m: ${mean(model[b]).toFixed(2)} °C against the atlas's ${mean(observed[b]).toFixed(2)}`));
});

test('climatology: null gives the analytic start byte for byte and leaves the surface alone; a climatology that is not decoded is refused', async () => {
  const { mesh, geography, atlas } = await atlasStart();
  const C = mesh.nCells, ice = new Float64Array(C);
  const plain = createOcean(mesh, { geography }), withAtlas = createOcean(mesh, { geography, climatology: atlas });
  const a = Float64Array.from(mesh.latCell, (lat) => 275 + 25 * Math.cos(lat) ** 2), b = Float64Array.from(a);
  assert.equal(plain.initialize(a, ice), null);
  assert.equal(withAtlas.initialize(b, ice, { climatology: null }), null);
  assert.deepEqual(Array.from(b), Array.from(a));
  for (const name of ['h', 'Q', 'W', 'eta', 'T0', 'S0', 'capacity']) assert.deepEqual(Array.from(withAtlas[name]), Array.from(plain[name]), name);
  assert.throws(() => createOcean(mesh, { geography, climatology: CLIMATOLOGY_FILE }), /decoded/);
  assert.equal(await loadClimatology(null), null);
  assert.equal(await loadClimatology(atlas), atlas);
});

test('the GPU ocean starts from the atlas as the CPU ocean does and sends the device the sea surface it sets', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const atlas = await loadClimatology(FILE.pathname);
  const model = await createGpuModel(new Grid(8), { topography, ocean: { climatology: FILE.pathname } });
  const { mesh, state } = model, C = mesh.nCells;
  for (let i = 0; i < C; i++) { state[3][i] = 272 + 30 * Math.cos(mesh.latCell[i]) ** 2; state[6][i] = !model.geography.land[i] && mesh.latCell[i] < -68 * DEG ? 0.7 : 0; }
  const surfaceT = Float64Array.from(state[3]);
  model.load();
  const gpuStart = model.ocean.initialize(state[3], state[6]);
  const cpu = createOcean(mesh, { geography: model.geography, climatology: atlas });
  const cpuStart = cpu.initialize(surfaceT, state[6]);
  assert.deepEqual(gpuStart, cpuStart);
  const gpu = await model.ocean.serialize(), expected = cpu.serialize();
  const worst = (name, scale) => { let m = 0; for (let n = 0; n < expected[name].length; n++) m = Math.max(m, Math.abs(expected[name][n] - gpu[name][n]) / scale(n)); return m; };
  assert.ok(worst('h', (n) => Math.max(1, Math.abs(expected.h[n]))) < 1e-6, 'thicknesses to single precision');
  assert.ok(worst('T', () => 300) < 1e-6 && worst('S', () => 35) < 1e-6, 'temperatures and salinities to single precision');
  assert.ok(worst('eta', () => 1) < 1e-5, 'sea level');
  assert.deepEqual(Array.from(state[3]), Array.from(surfaceT), 'the same sea surface');
  const [, , , deviceSurface] = await model.gpu.download();
  for (let i = 0; i < C; i++) assert.ok(Math.abs(deviceSurface[i] - surfaceT[i]) < 1e-4, `the device holds the sea surface at cell ${i}`);
  model.destroy();
});
