import { test } from 'node:test';
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { existsSync, mkdtempSync, readFileSync, writeFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { AEROSOL_FILE, AEROSOL_FIELDS, AEROSOL_RANGES, WAVELENGTHS, decodeAerosol, aerosolFields, aerosolDay, aerosolDayOfYear, createAerosolSampler, loadAerosol, streamOptics, aerosolTable } from '../js/aerosol.module.js';
import { AEROSOL_SPECTRA } from '../js/physics/aerosolSpectra.module.js';
import { LAYERED_TABLE } from '../js/physics/layeredTable.module.js';
import { monthWeights } from '../js/ocean/climatology.module.js';
import { netcdfClassic } from '../scripts/netcdfClassic.mjs';
import { syntheticAerosol } from './helpers/aerosol.mjs';
import { writeClassic } from './helpers/netcdf.mjs';

const here = dirname(fileURLToPath(import.meta.url)), root = join(here, '..'), mesh = buildMesh(new Grid(16)), C = mesh.nCells, deg = 180 / Math.PI;
const shipped = join(root, AEROSOL_FILE);
const step = (field) => { const [low, high] = AEROSOL_RANGES[field]; return (high - low) / 255; };

test('a climatology file decodes to the fields it was written from, each within the byte\'s step of its range', () => {
  const { data, bytes } = syntheticAerosol({ nLon: 12, nLat: 6, nMonth: 4 });
  assert.equal(String.fromCharCode(...bytes.subarray(0, 4)), 'AERM');
  assert.deepEqual([data.nLon, data.nLat, data.nMonth, data.nField, data.source], [12, 6, 4, AEROSOL_FIELDS, 'synthetic']);
  assert.deepEqual(Array.from(data.centers), [45.625, 136.875, 228.125, 319.375]);
  assert.equal(data.values.length, AEROSOL_FIELDS * 4 * 6 * 12);
  const wanted = syntheticAerosol({ nLon: 12, nLat: 6, nMonth: 4, sample: (f, lat, lon, m) => AEROSOL_RANGES[f][0] + (AEROSOL_RANGES[f][1] - AEROSOL_RANGES[f][0]) * ((Math.sin(lat + 2 * lon + m + f) + 1) / 2) ** (AEROSOL_RANGES[f][2] ? 2 : 1) }).data;
  for (let f = 0; f < AEROSOL_FIELDS; f++) {
    let worst = 0;
    for (let n = 0; n < 4 * 6 * 12; n++) {
      const lat = -90 + 30 * (Math.floor(n / 12) % 6) + 15, lon = -180 + 30 * (n % 12) + 15, m = Math.floor(n / 72);
      const [low, high, square] = AEROSOL_RANGES[f], s = (Math.sin(lat + 2 * lon + m + f) + 1) / 2, value = low + (high - low) * (square ? s * s : s);
      worst = Math.max(worst, Math.abs(wanted.lookup[f * 256 + wanted.values[f * 288 + n]] - value));
    }
    assert.ok(worst <= (AEROSOL_RANGES[f][2] ? 2 * (AEROSOL_RANGES[f][1] - AEROSOL_RANGES[f][0]) / 255 : 0.5 * step(f) + 1e-6), `field ${f}: worst ${worst}`);
  }
  assert.throws(() => decodeAerosol(bytes.subarray(0, bytes.length - 1)), /bytes, not/);
  assert.throws(() => decodeAerosol(Uint8Array.from(bytes, (b, n) => (n === 0 ? 0x57 : b))), /not an aerosol climatology/);
});

test('two months blend linearly between the month centers, in the stored variables, and the year wraps from December to January', () => {
  const nMonth = 4, { data } = syntheticAerosol({ nLon: 12, nLat: 6, nMonth, sample: (f, lat, lon, m) => (f < 2 ? 0.1 + 0.2 * m : f < 4 ? 0.02 + 0.02 * m : f < 6 ? 0.5 + 0.05 * m : 1 + m) });
  const lookup = (f, month) => data.lookup[f * 256 + data.values[f * nMonth * 72 + month * 72]];
  const out = new Float32Array(AEROSOL_FIELDS * C);
  for (const day of [0, 10, 45.625, 90, 136.875, 200, 228.125, 280, 319.375, 340, 364.9]) {
    aerosolFields(mesh, data, day, out);
    const weights = monthWeights(data, day);
    for (const f of [0, 3, 5, 7]) {
      const expected = weights.reduce((s, [m, w]) => s + w * lookup(f, m), 0);
      for (const i of [0, 100, 1000, C - 1]) assert.ok(Math.abs(out[f * C + i] - expected) < 1e-6 * Math.max(1, expected), `day ${day} field ${f}: ${out[f * C + i]} against ${expected}`);
    }
  }
  const [[m0, w0], [m1, w1]] = monthWeights(data, 340);
  assert.deepEqual([m0, m1], [3, 0]);
  assert.ok(Math.abs(w0 + w1 - 1) < 1e-12 && w0 > 0 && w1 > 0, 'the last center and the first of the next year bracket December');
  const [[a], [b]] = monthWeights(data, 10);
  assert.deepEqual([a, b], [3, 0]);
});

test('the sampler\'s bilinear taps are the grid\'s, at a hundred cells, with the rows past the first and last centers held and the columns wrapping', () => {
  const { data } = syntheticAerosol({ nLon: 36, nLat: 18, nMonth: 12 }), nLon = 36, nLat = 18, plane = 12 * nLat * nLon;
  const day = 100.5, out = aerosolFields(mesh, data, day), weights = monthWeights(data, day);
  const at = (f, month, r, c) => data.lookup[f * 256 + data.values[f * plane + (month * nLat + r) * nLon + ((c % nLon) + nLon) % nLon]];
  let worst = 0, edge = 0;
  for (let i = 0; i < C; i += Math.floor(C / 100)) {
    const lat = mesh.latCell[i] * deg, lon = mesh.lonCell[i] * deg, y = Math.min(nLat - 1, Math.max(0, (lat - data.lat0) / data.dLat)), x = (lon - data.lon0) / data.dLon;
    const r0 = Math.min(nLat - 2, Math.floor(y)), c0 = Math.floor(x), fy = y - r0, fx = x - c0;
    for (let f = 0; f < AEROSOL_FIELDS; f++) {
      let expected = 0;
      for (const [m, w] of weights) expected += w * ((1 - fy) * ((1 - fx) * at(f, m, r0, c0) + fx * at(f, m, r0, c0 + 1)) + fy * ((1 - fx) * at(f, m, r0 + 1, c0) + fx * at(f, m, r0 + 1, c0 + 1)));
      worst = Math.max(worst, Math.abs(out[f * C + i] - expected));
    }
    if (Math.abs(lat) > 85 || Math.abs(Math.abs(lon) - 180) < 6) edge++;
  }
  assert.ok(worst < 2e-6, `the sampler differs from the direct evaluation by ${worst}`);
  assert.ok(edge >= 3, `${edge} sampled cells near a pole or the date line`);
  for (let i = 0; i < C; i++) for (let f = 0; f < AEROSOL_FIELDS; f++) assert.ok(Number.isFinite(out[f * C + i]) && out[f * C + i] >= 0, `cell ${i} field ${f}: ${out[f * C + i]}`);
});

test('the model day maps to the year as the repository\'s scripts do, and the sampler is refreshed once a model day', () => {
  assert.equal(aerosolDay(0), 0);
  assert.equal(aerosolDay(86399), 0);
  assert.equal(aerosolDay(86400), 1);
  assert.equal(aerosolDayOfYear(0), 78.5);
  assert.equal(aerosolDayOfYear(286), 364.5);
  assert.equal(aerosolDayOfYear(287), 0.5);
  assert.equal(aerosolDayOfYear(287 + 365), 0.5);
  assert.equal(aerosolDayOfYear(-79), 364.5);
  const { data } = syntheticAerosol(), sampler = createAerosolSampler(mesh, data);
  assert.equal(sampler.update(0), true);
  const before = Float32Array.from(sampler.fields);
  assert.equal(sampler.update(86399), false);
  assert.deepEqual(sampler.fields, before);
  assert.equal(sampler.update(86400), true);
  assert.equal(sampler.day, 1);
  assert.equal(sampler.update(40 * 86400), true);
  assert.notDeepEqual(sampler.fields, before, 'forty days later the fields differ');
  assert.equal(sampler.update(365 * 86400 + 86400), true);
  assert.equal(sampler.day, 366);
  const year = createAerosolSampler(mesh, data);
  year.update(0);
  assert.deepEqual(year.fields, before, 'a model year later is the same day of the year');
});

test('the data are given decoded, as bytes or, under node, from the repository\'s file, which a missing file names', () => {
  const { data, bytes } = syntheticAerosol();
  assert.equal(loadAerosol(data), data);
  assert.deepEqual(loadAerosol(bytes).values, data.values);
  assert.deepEqual(loadAerosol(bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength)).values, data.values);
  if (!existsSync(shipped)) assert.throws(() => loadAerosol(), /scripts\/packAerosol\.py/);
  else assert.equal(loadAerosol().nField, AEROSOL_FIELDS);
});

test('the spectral constants cover MACv2\'s 12 wavelengths, take the 550 nm triple to itself, and weigh the two streams and the visible slices to one', () => {
  const s = AEROSOL_SPECTRA;
  assert.equal(s.wavelengths.length, WAVELENGTHS);
  assert.deepEqual(s.wavelengths, [400, 550, 700, 1000, 1270, 1460, 1780, 2050, 2320, 2790, 3470, 8000]);
  const near = LAYERED_TABLE.bands.map((b) => b.band).sort((a, b) => a - b);
  assert.deepEqual(near, [16, 17, 18, 19, 20, 21, 22, 23, 24, 29]);
  for (const band of near) assert.ok(s.bands.includes(band), `RRTMG_SW band ${band} has a wavelength`);
  assert.deepEqual(s.bands.slice(0, 2), [26, 25]);
  for (const mode of [s.fine, s.coarse]) {
    for (const key of ['r', 'a', 'c', 'd', 'e']) assert.equal(mode[key].length, WAVELENGTHS, key);
    assert.deepEqual([mode.r[1], mode.a[1], mode.c[1], mode.d[1], mode.e[1]], [1, 0, 1, 0, 1], 'the 550 nm row is the identity');
    for (let b = 9; b < WAVELENGTHS; b++) assert.deepEqual([mode.c[b], mode.e[b]], [0, 0], `beyond 2.5 um the band means stand in (${s.wavelengths[b]} nm)`);
  }
  for (const weights of [s.visible, s.nearInfrared, s.slices]) assert.ok(Math.abs(weights.reduce((a, b) => a + b, 0) - 1) < 1e-5);
  assert.deepEqual(s.visible.map((w) => w > 0), s.wavelengths.map((nm) => nm <= 700));
  assert.deepEqual(s.nearInfrared.map((w) => w > 0), s.wavelengths.map((nm) => nm >= 700));
  assert.deepEqual(s.slices.map((w) => w > 0), s.wavelengths.map((nm) => nm <= 550));
  const table = streamOptics([0.1, 0.2], [0.05, 0.02], [0.65, 0.76], aerosolTable());
  assert.ok(Math.abs(table.tau[1] - 0.1) < 1e-12 && Math.abs(table.tau[WAVELENGTHS + 1] - 0.2) < 1e-12 && Math.abs(table.x[1] - 0.05) < 1e-12 && Math.abs(table.g[WAVELENGTHS + 1] - 0.76) < 1e-12);
});

test('netcdfClassic reads the variables and records of a classic file, whole or by record', () => {
  const lat = Float64Array.from([60, 20, -20, -60]), lon = Float64Array.from({ length: 8 }, (_, c) => -157.5 + 45 * c);
  const values = Float64Array.from({ length: 3 * 4 * 8 }, (_, n) => 0.5 + n / 100);
  const bytes = writeClassic({ dims: [['time', 0], ['lat', 4], ['lon', 8]], numrecs: 3, vars: [{ name: 'lat', dims: ['lat'], values: lat }, { name: 'lon', dims: ['lon'], values: lon }, { name: 'time', dims: ['time'], values: Float64Array.of(14, 45, 74) }, { name: 'aod', dims: ['time', 'lat', 'lon'], values }, { name: 'ssa', dims: ['time', 'lat', 'lon'], values: values.map((v) => v / 2) }] });
  const dir = mkdtempSync(join(tmpdir(), 'netcdf-')), path = join(dir, 'x.nc');
  writeFileSync(path, bytes);
  try {
    const file = netcdfClassic(path);
    assert.deepEqual(file.dims, { time: 0, lat: 4, lon: 8 });
    assert.equal(file.numrecs, 3);
    assert.deepEqual(Array.from(file.read('lat')), Array.from(lat));
    assert.deepEqual(Array.from(file.read('time', 2)), [74]);
    for (let r = 0; r < 3; r++) {
      assert.deepEqual(Array.from(file.read('aod', r)), Array.from(values.subarray(r * 32, (r + 1) * 32), Math.fround));
      assert.deepEqual(Array.from(file.read('ssa', r)), Array.from(values.subarray(r * 32, (r + 1) * 32), (v) => Math.fround(v / 2)));
    }
    file.close();
  } finally { rmSync(dir, { recursive: true, force: true }); }
});

const python = spawnSync('python3', ['-I', '-c', 'import sys; print(sys.version_info[:2] >= (3, 8))'], { encoding: 'utf8' }).stdout?.trim() === 'True';

/*
 * MACv2's four files on a 45 degree grid of 8 x 4 cells with three months, and the file the builder makes
 * of them at 2 x 2 cells to a cell: the 550 nm optics of each mode, and the mode's layer fractions with the
 * dataset's surface altitude or as a one-month record.
 */
function rawFiles(dir, { full }) {
  const nLon = 8, nLat = 4, months = 3, plane = nLat * nLon, lat = Float64Array.from([67.5, 22.5, -22.5, -67.5]), lon = Float64Array.from({ length: nLon }, (_, c) => -157.5 + 45 * c);
  const field = (mode, month, row, col) => {
    const m = mode === 0 ? 1 : 2;
    return { aod: 0.05 * m * (1 + 0.5 * Math.sin(1 + row + 2 * col + month)) + 0.01 * row, ssa: 0.97 - 0.03 * ((row + col + month + mode) % 4), asy: 0.55 + 0.06 * mode + 0.01 * ((2 * row + col) % 5) };
  };
  const optics = (mode) => {
    const out = { aod: new Float64Array(months * plane), ssa: new Float64Array(months * plane), asy: new Float64Array(months * plane) };
    for (let m = 0; m < months; m++) for (let r = 0; r < nLat; r++) for (let c = 0; c < nLon; c++) for (const key of ['aod', 'ssa', 'asy']) out[key][m * plane + r * nLon + c] = Math.fround(field(mode, m, r, c)[key]);
    return out;
  };
  const centers = Float64Array.of(14, 100, 200), layers = 20;
  const top = Float64Array.from([500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 15000, 20000]), bottom = Float64Array.from([0, ...top.subarray(0, layers - 1)]);
  const asl = Float64Array.from({ length: plane }, (_, n) => (n % nLon === 3 ? 1200 : 0));
  const fractions = (mode, month) => {
    const out = new Float64Array(layers * plane);
    for (let n = 0; n < plane; n++) {
      const lowest = asl[n] > 0 ? 2 : 0;
      let total = 0;
      for (let l = lowest; l < layers; l++) { const f = Math.exp(-(l - lowest) / (3 + mode + month)) * (1 + 0.2 * Math.sin(n + l)); out[l * plane + n] = Math.fround(f); total += Math.fround(f); }
      for (let l = lowest; l < layers; l++) out[l * plane + n] = Math.fround(out[l * plane + n] / total);
    }
    return out;
  };
  const expected = [];
  for (const mode of [0, 1]) {
    const name = mode === 0 ? 'f' : 'c', o = optics(mode);
    writeFileSync(join(dir, `gt_${name}_00550nm.nc`), writeClassic({ dims: [['time', 0], ['lat', nLat], ['lon', nLon]], numrecs: months, vars: [{ name: 'lon', dims: ['lon'], values: lon }, { name: 'lat', dims: ['lat'], values: lat }, { name: 'time', dims: ['time'], values: centers }, { name: 'aod', dims: ['time', 'lat', 'lon'], values: o.aod }, { name: 'ssa', dims: ['time', 'lat', 'lon'], values: o.ssa }, { name: 'asy', dims: ['time', 'lat', 'lon'], values: o.asy }] }));
    const file = mode === 0 ? 'aeropt_fin_alt_km20' : 'aeropt_coa_alt_km20', variable = mode === 0 ? 'zkm_aer_fine_mo' : 'zkm_aer_coarse_mo';
    const perMonth = Array.from({ length: months }, (_, m) => fractions(mode, full ? m : 0));
    if (full) {
      writeFileSync(join(dir, `${file}.nc`), writeClassic({ dims: [['time', 0], ['lay', layers], ['lat', nLat], ['lon', nLon]], numrecs: months, vars: [{ name: 'time', dims: ['time'], values: centers }, { name: 'ztop', dims: ['lay'], values: top }, { name: 'zbot', dims: ['lay'], values: bottom }, { name: 'asl', dims: ['lat', 'lon'], values: asl }, { name: variable, dims: ['time', 'lay', 'lat', 'lon'], values: Float64Array.from({ length: months * layers * plane }, (_, n) => perMonth[Math.floor(n / (layers * plane))][n % (layers * plane)]) }] }));
    } else {
      const record = Buffer.alloc(4 + 4 * layers * plane);
      record.writeFloatBE(20060832, 0);
      perMonth[0].forEach((v, n) => record.writeFloatBE(v, 4 + 4 * n));
      writeFileSync(join(dir, `${file}.record`), record);
    }
    expected.push({ o, perMonth });
  }
  return { nLon, nLat, months, lat, expected, asl, top, bottom, layers, plane };
}

function packed(dir, factor = 2) {
  const out = join(dir, 'out.bin'), run = spawnSync('python3', ['-I', join(root, 'scripts', 'packAerosol.py'), dir, out, String(factor)], { encoding: 'utf8' });
  assert.equal(run.status, 0, run.stderr || run.stdout);
  return { data: decodeAerosol(readFileSync(out)), log: run.stdout };
}

for (const full of [true, false]) test(`scripts/packAerosol.py makes the file of ${full ? 'complete files with the dataset\'s surface altitude' : 'one-month altitude records, whose ground is the lowest layer that holds aerosol'}: AOD by area, 1 - SSA by AOD, asymmetry by AOD times SSA, height by AOD, as bytes of their ranges`, { skip: !python && 'python3 not found' }, () => {
  const dir = mkdtempSync(join(tmpdir(), 'aerosol-raw-'));
  try {
    const source = rawFiles(dir, { full }), { data, log } = packed(dir);
    assert.deepEqual([data.nLon, data.nLat, data.nMonth, data.nField], [4, 2, 3, AEROSOL_FIELDS]);
    assert.deepEqual([data.lon0, data.dLon, data.lat0, data.dLat], [-135, 90, -45, 90]);
    assert.deepEqual(Array.from(data.centers), [14, 100, 200]);
    assert.match(data.source, /MACv2.*Kinne 2019.*no license stated/);
    assert.match(log, /clipped, per field: \[0, 0, 0, 0, 0, 0, 0, 0\]/);
    const middle = source.expected.map(() => Float64Array.from(source.top, (t, l) => 0.5 * (t + source.bottom[l]) / 1000));
    let worst = new Float64Array(AEROSOL_FIELDS);
    const plane = 4 * 2;
    for (let m = 0; m < 3; m++) for (let j = 0; j < 2; j++) for (let i = 0; i < 4; i++) {
      for (const mode of [0, 1]) {
        const { o, perMonth } = source.expected[mode];
        let area = 0, depth = 0, scatter = 0, tilt = 0, lift = 0;
        for (let dj = 0; dj < 2; dj++) for (let di = 0; di < 2; di++) {
          const r = 3 - (2 * j + dj), c = 2 * i + di, n = r * 8 + c, w = Math.cos(source.lat[r] * Math.PI / 180);
          const a = o.aod[m * 32 + n], s = o.ssa[m * 32 + n], g = o.asy[m * 32 + n], fractions = perMonth[m];
          let total = 0, weighted = 0, lowest = -1;
          for (let l = 0; l < source.layers; l++) { const f = fractions[l * 32 + n]; if (f > 0) { if (lowest < 0) lowest = l; total += f; weighted += f * middle[mode][l]; } }
          const ground = full ? source.asl[n] / 1000 : source.bottom[lowest] / 1000, height = Math.max(0, weighted / total - ground);
          area += w; depth += w * a; scatter += w * a * s; tilt += w * a * s * g; lift += w * a * height;
        }
        const wanted = [depth / area, 1 - scatter / depth, tilt / scatter, lift / depth];
        wanted.forEach((value, kind) => {
          const f = 2 * kind + mode, [low, high] = AEROSOL_RANGES[f], got = data.lookup[f * 256 + data.values[f * 3 * plane + m * plane + j * 4 + i]];
          worst[f] = Math.max(worst[f], Math.abs(got - value));
        });
      }
    }
    for (let f = 0; f < AEROSOL_FIELDS; f++) {
      const [low, high, isSquare] = AEROSOL_RANGES[f], limit = isSquare ? 2 * (high - low) * Math.sqrt(0.9) / 255 + 1e-5 : 0.5 * (high - low) / 255 + 1e-5;
      assert.ok(worst[f] <= limit, `field ${f}: the file differs from the coarsened data by ${worst[f]}, more than ${limit}`);
    }
  } finally { rmSync(dir, { recursive: true, force: true }); }
});

test('scripts/packAerosol.py refuses a factor that does not divide the grid and a missing altitude file', { skip: !python && 'python3 not found' }, () => {
  const dir = mkdtempSync(join(tmpdir(), 'aerosol-raw-'));
  try {
    rawFiles(dir, { full: false });
    const run = spawnSync('python3', ['-I', join(root, 'scripts', 'packAerosol.py'), dir, join(dir, 'x.bin'), '3'], { encoding: 'utf8' });
    assert.notEqual(run.status, 0);
    assert.match(run.stderr, /does not divide/);
    rmSync(join(dir, 'aeropt_coa_alt_km20.record'));
    const missing = spawnSync('python3', ['-I', join(root, 'scripts', 'packAerosol.py'), dir, join(dir, 'x.bin')], { encoding: 'utf8' });
    assert.notEqual(missing.status, 0);
    assert.match(missing.stderr, /aeropt_coa_alt_km20/);
  } finally { rmSync(dir, { recursive: true, force: true }); }
});

/*
 * The 24 numbers of scripts/aerosolFit.mjs against MACv2's band data (August): the Planck-weighted depth and
 * absorption depth of the visible stream and of the rest of the beam, in the area means of six boxes.
 */
const BANDS = {
  global: [0.1522, 0.01099, 0.08681, 0.005235], 'tropics 10S-10N': [0.1845, 0.01581, 0.08707, 0.00662], 'midlatitudes 40-60N': [0.1604, 0.009624, 0.07104, 0.004362],
  India: [0.3491, 0.02297, 0.1846, 0.01062], 'Sahara and Arabia': [0.5265, 0.04924, 0.3696, 0.02556], 'East Asia': [0.4471, 0.02727, 0.1579, 0.01165],
};
const BOXES = { global: [-90, 90, -180, 180], 'tropics 10S-10N': [-10, 10, -180, 180], 'midlatitudes 40-60N': [40, 60, -180, 180], India: [8, 28, 70, 90], 'Sahara and Arabia': [15, 32, -17, 60], 'East Asia': [22, 45, 100, 125] };

test('the shipped file reproduces MACv2\'s global means, fills its byte ranges without crowding their ends, and rebuilds the band data\'s stream depths in its boxes to a few per cent', { skip: !existsSync(shipped) && `${AEROSOL_FILE} is not in the repository (MPI-M states no license for it: scripts/packAerosol.py builds it, docs/aerosol-climatology.md)` }, () => {
  const data = decodeAerosol(readFileSync(shipped)), { nLon, nLat, nMonth, lat0, dLat, lon0, dLon } = data, plane = nMonth * nLat * nLon;
  assert.deepEqual([nLon, nLat, nMonth], [180, 90, 12]);
  const value = (f, m, r, c) => data.lookup[f * 256 + data.values[f * plane + (m * nLat + r) * nLon + c]];
  const mean = (f, months) => {
    let sum = 0, area = 0;
    for (const m of months) for (let r = 0; r < nLat; r++) { const w = Math.cos((lat0 + r * dLat) / deg); for (let c = 0; c < nLon; c++) { sum += w * value(f, m, r, c); area += w; } }
    return sum / area;
  };
  const year = Array.from({ length: 12 }, (_, m) => m);
  assert.ok(Math.abs(mean(0, year) - 0.063) < 0.004 && Math.abs(mean(1, year) - 0.058) < 0.004, `annual AOD ${mean(0, year)} fine, ${mean(1, year)} coarse against 0.063 and 0.058 (Kinne 2019, Table 2)`);
  assert.ok(Math.abs(mean(0, [7]) - 0.083) < 0.005 && Math.abs(mean(1, [7]) - 0.059) < 0.005, 'August');
  assert.ok(mean(2, year) > 0.05 && mean(2, year) < 0.09 && mean(3, year) > 0.015 && mean(3, year) < 0.04, '1 - SSA');
  assert.ok(mean(4, year) > 0.6 && mean(4, year) < 0.7 && mean(5, year) > 0.73 && mean(5, year) < 0.79, 'asymmetry');
  assert.ok(mean(6, year) > 1 && mean(6, year) < 3.5 && mean(7, year) > 1 && mean(7, year) < 2.5, 'mean heights');
  for (let f = 0; f < AEROSOL_FIELDS; f++) {
    let top = 0, count = 0;
    for (let n = 0; n < plane; n++) { const q = data.values[f * plane + n]; if (q === 255) top++; count++; }
    assert.ok(top < 0.001 * count, `field ${f}: ${top} of ${count} values at the top of the range`);
  }
  const table = aerosolTable(), tau = [0, 0], x = [0, 0], g = [0, 0];
  for (const [name, [south, north, west, east]] of Object.entries(BOXES)) {
    const sums = [0, 0, 0, 0];
    let area = 0;
    for (let r = 0; r < nLat; r++) {
      const lat = lat0 + r * dLat;
      if (lat < south || lat > north) continue;
      for (let c = 0; c < nLon; c++) {
        const lon = lon0 + c * dLon;
        if (lon < west || lon > east) continue;
        for (let m = 0; m < 2; m++) { tau[m] = value(m, 7, r, c); x[m] = value(2 + m, 7, r, c); g[m] = value(4 + m, 7, r, c); }
        streamOptics(tau, x, g, table);
        const w = Math.cos(lat / deg);
        area += w;
        let aV = 0, aN = 0;
        for (let m = 0; m < 2; m++) for (let b = 0; b < WAVELENGTHS; b++) { aV += AEROSOL_SPECTRA.visible[b] * table.tau[m * WAVELENGTHS + b] * table.x[m * WAVELENGTHS + b]; aN += AEROSOL_SPECTRA.nearInfrared[b] * table.tau[m * WAVELENGTHS + b] * table.x[m * WAVELENGTHS + b]; }
        sums[0] += w * table.tauV; sums[1] += w * aV; sums[2] += w * table.tauN; sums[3] += w * aN;
      }
    }
    const want = BANDS[name], got = sums.map((s) => s / area), errors = got.map((v, n) => v / want[n] - 1);
    console.log(`${name}: stream depths from the file ${got.map((v) => v.toFixed(4)).join(' ')} against the band data ${want.join(' ')} (${errors.map((e) => `${(100 * e).toFixed(1)}%`).join(' ')})`);
    errors.forEach((e, n) => assert.ok(Math.abs(e) < 0.1, `${name} quantity ${n}: ${(100 * e).toFixed(1)} % off`));
  }
});
