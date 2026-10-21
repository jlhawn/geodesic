import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { topographyFromInt16, createGeography } from '../js/geography.module.js';
import { createOcean, EPS, THIN, RESTORE_TOLERANCE, LAYER_DENSITIES, LAYER_SALINITIES } from '../js/ocean/layered.module.js';
import { decodeClimatology, encodeClimatology, loadClimatology, profileAt, CLIMATOLOGY_FILE, decodeMonthlySst, encodeMonthlySst, monthWeights, climatologyAt, MONTHLY_SST_FILE } from '../js/ocean/climatology.module.js';
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

/*
 * A 12 × 6 monthly grid at 30° spacing over four months, linear in month,
 * latitude and longitude, with land from pole to pole in columns 2–8
 * (105W–75E).
 */
const MONTHLY = { nLon: 12, nLat: 6, lon0: -165, dLon: 30, lat0: -75, dLat: 30, centers: [45, 136, 228, 319] };
const monthlyT = (m, lat, lon) => 10 + m + 0.1 * lat + 0.02 * lon;
function monthlyFixture() {
  const { nLon, nLat, lon0, dLon, lat0, dLat, centers } = MONTHLY, T = [];
  centers.forEach((_, m) => { for (let r = 0; r < nLat; r++) for (let c = 0; c < nLon; c++) T.push(c >= 2 && c <= 8 ? NaN : monthlyT(m, lat0 + r * dLat, lon0 + c * dLon)); });
  return { ...MONTHLY, T, source: 'monthly fixture' };
}

test('a monthly sea surface climatology packed from arrays decodes to the same grid, values within half a step and missing points kept, and the repository file re-encodes byte for byte, so scripts/packWoaMonthlySst.py and encodeMonthlySst share one layout', () => {
  const grid = monthlyFixture(), sst = decodeMonthlySst(encodeMonthlySst(grid));
  for (const key of ['nLon', 'nLat', 'lon0', 'dLon', 'lat0', 'dLat', 'source']) assert.equal(sst[key], grid[key], key);
  assert.deepEqual(Array.from(sst.centers), grid.centers);
  for (let n = 0; n < grid.T.length; n++) {
    if (Number.isNaN(grid.T[n])) { assert.ok(Number.isNaN(sst.T[n]), `missing ${n}`); continue; }
    assert.ok(Math.abs(sst.T[n] - grid.T[n]) <= 0.5 * sst.tScale + 1e-9, `T ${n}`);
  }
  const bytes = new Uint8Array(readFileSync(new URL(`../${MONTHLY_SST_FILE}`, import.meta.url))), woa = decodeMonthlySst(bytes), again = encodeMonthlySst(woa);
  assert.equal(again.length, bytes.length);
  assert.ok(again.every((b, n) => b === bytes[n]), 'identical bytes');
  assert.equal(woa.nLon, 360); assert.equal(woa.nLat, 180); assert.equal(woa.nMonth, 12);
  assert.deepEqual(Array.from(woa.centers), [15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5]);
  let wet = 0, low = Infinity, high = -Infinity;
  for (const t of woa.T) if (Number.isFinite(t)) { wet++; low = Math.min(low, t); high = Math.max(high, t); }
  assert.ok(wet > 12 * 40000 && low > -2.5 && high < 36, `${wet} values from ${low} to ${high} °C`);
});

test('monthWeights takes the two months whose centers bracket a day, linearly, across the turn of the year; climatologyAt blends them, bilinearly among wet points, renormalized over dry ones, else from the nearest wet point within two steps, else NaN', () => {
  const woa = decodeMonthlySst(readFileSync(new URL(`../${MONTHLY_SST_FILE}`, import.meta.url)));
  const close = (got, want, label) => { assert.equal(got.length, 2, label); got.forEach(([m, w], q) => { assert.equal(m, want[q][0], label); assert.ok(Math.abs(w - want[q][1]) < 1e-12, `${label}: ${w}`); }); };
  close(monthWeights(woa, 108), [[3, 27.5 / 30.5], [4, 3 / 30.5]], 'April 19');
  close(monthWeights(woa, 15.5), [[0, 1], [1, 0]], 'January 16, noon');
  close(monthWeights(woa, 0), [[11, 0.5], [0, 0.5]], 'January 1');
  close(monthWeights(woa, 360), [[11, 20.5 / 31], [0, 10.5 / 31]], 'December 27');
  for (let day = 0; day < 365; day += 0.25) { const [[, a], [, b]] = monthWeights(woa, day); assert.ok(a >= 0 && b >= 0 && Math.abs(a + b - 1) < 1e-12, `day ${day}`); }
  const sst = decodeMonthlySst(encodeMonthlySst(monthlyFixture())), months = [[0, 0.25], [1, 0.75]], blend = (lat, lon) => 0.25 * monthlyT(0, lat, lon) + 0.75 * monthlyT(1, lat, lon);
  close(monthWeights(sst, 113.25), [[0, 0.25], [1, 0.75]], 'the fixture');
  assert.ok(Math.abs(climatologyAt(sst, months, 0, 120 * DEG) - KELVIN - blend(0, 120)) < 1e-3, 'between four wet points');
  assert.ok(Math.abs(climatologyAt(sst, months, 0, 80 * DEG) - KELVIN - blend(0, 105)) < 1e-3, 'beside the coast, from the wet corners');
  assert.ok(Math.abs(climatologyAt(sst, months, 15 * DEG, -105 * DEG) - KELVIN - blend(15, -135)) < 1e-3, 'on land, the nearest wet point');
  assert.ok(Number.isNaN(climatologyAt(sst, months, 0, 15 * DEG)), 'more than two steps inland');
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
  const model = await createGpuModel(new Grid(16), { topography, ocean: { climatology: FILE.pathname } });
  const { mesh, state } = model, C = mesh.nCells;
  for (let i = 0; i < C; i++) { state[3][i] = 272 + 30 * Math.cos(mesh.latCell[i]) ** 2; state[6][i] = !model.geography.land[i] && mesh.latCell[i] < -68 * DEG ? 0.7 : 0; }
  const surfaceT = Float64Array.from(state[3]);
  model.load();
  const gpuStart = model.ocean.initialize(state[3], state[6]);
  const cpu = createOcean(mesh, { geography: model.geography, climatology: atlas });
  const cpuStart = cpu.initialize(surfaceT, state[6]);
  assert.deepEqual(gpuStart, cpuStart);
  const gpu = await model.ocean.serialize(), expected = cpu.serialize();
  const K = expected.h.length / C, profile = (h, T, i, z) => { let zz = 0; for (let k = 0; k < K; k++) { zz += h[k * C + i]; if (z < zz) return T[k * C + i]; } return T[(K - 1) * C + i]; };
  let interfaces = 0, close = 0, heat = 0, salt = 0, rmsSum = 0, rmsMax = 0, sea = 0;
  for (let i = 0; i < C; i++) {
    if (model.geography.land[i]) continue;
    sea++;
    let za = 0, zb = 0, heatA = 0, heatB = 0, saltA = 0, saltB = 0;
    for (let k = 0; k < K; k++) {
      const n = k * C + i, ha = expected.h[n], hb = gpu.h[n];
      heatA += ha * expected.T[n]; heatB += hb * gpu.T[n]; saltA += ha * expected.S[n]; saltB += hb * gpu.S[n];
      za += ha; zb += hb;
      if (k + 1 < K && expected.h[n + C] > 1 && Math.abs(expected.T[n] - expected.T[n + C]) > 0.05) { interfaces++; if (Math.abs(za - zb) < 2) close++; }
    }
    heat = Math.max(heat, Math.abs(heatA - heatB) / za); salt = Math.max(salt, Math.abs(saltA - saltB) / za);
    let ss = 0, m = 0;
    for (let z = 2.5; z < za; z += 5) { ss += (profile(expected.h, expected.T, i, z) - profile(gpu.h, gpu.T, i, z)) ** 2; m++; }
    const rms = Math.sqrt(ss / m); rmsSum += rms; rmsMax = Math.max(rmsMax, rms);
  }
  assert.ok(close > 0.95 * interfaces, `the stratified interfaces within 2 m of the CPU's but for the device reconstruction's knife edges in weakly stratified water (${close} of ${interfaces})`);
  assert.ok(rmsSum / sea < 0.05 && rmsMax < 0.5, `the temperature profiles alike (rms ${(rmsSum / sea).toFixed(3)} K on average, ${rmsMax.toFixed(3)} K at worst)`);
  assert.ok(heat < 0.02 && salt < 0.002, `each column's heat and salt within the device remap's rounding (${heat.toExponential(2)} K, ${salt.toExponential(2)} g/kg)`);
  let eta = 0; for (let i = 0; i < C; i++) eta = Math.max(eta, Math.abs(expected.eta[i] - gpu.eta[i]));
  assert.ok(eta < 1e-5, 'sea level');
  assert.deepEqual(Array.from(state[3]), Array.from(surfaceT), 'the same sea surface');
  const [, , , deviceSurface] = await model.gpu.download();
  for (let i = 0; i < C; i++) assert.ok(Math.abs(deviceSurface[i] - surfaceT[i]) < 1e-4, `the device holds the sea surface at cell ${i}`);
  model.destroy();
});
