import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, existsSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { topographyFromInt16, syntheticTopography, decodeSubgrid, encodeSubgrid, meshSubgrid, subgridUrl, SUBGRID_FIELDS } from '../js/geography.module.js';
import { createModel, orographyFields } from '../js/model.module.js';
import { meshFields } from '../scripts/subgridTerrain.mjs';

const close = (actual, expected, tolerance, what) => assert.ok(Math.abs(actual - expected) <= tolerance * Math.abs(expected), `${what}: ${actual} against ${expected}`);
const deg = Math.PI / 180, R = 6371220;
const python = spawnSync('python3', ['-c', 'import numpy, PIL'], { encoding: 'utf8' });
const havePython = python.status === 0;

test('the 30″ filters of analytic orography: the 5 km smoothing and the 3–22 km band against the kernel’s response, and the spectral slope fit', { skip: !havePython && 'python3 with numpy and Pillow not found' }, () => {
  const run = spawnSync('python3', [new URL('../scripts/subgridTerrain.py', import.meta.url).pathname, 'selftest', '-'], { encoding: 'utf8', maxBuffer: 1 << 24 });
  assert.equal(run.status, 0, run.stderr);
  const { filters, spectrumSlope } = JSON.parse(run.stdout);
  for (const [name, f] of Object.entries(filters)) {
    close(f.mean, 1000, 1e-5, `${name} block mean`);
    close(f.h5Mean, 1000, 1e-5, `${name} smoothed mean`);
    close(f.h5Std, f.h5StdAnalytic, 0.005, `${name} amplitude after the 5 km smoothing`);
    close(f.sigmaFlt, f.sigmaFltAnalytic, 0.008, `${name} σ_flt`);
    assert.equal(f.land, 1);
  }
  close(filters.ridge8.sigmaFltAnalytic / (300 / Math.SQRT2), 0.816, 0.01, 'an 8 km ridge keeps 0.82 of its amplitude in the band');
  close(filters.ridge40.h5StdAnalytic / (300 / Math.SQRT2), 0.978, 0.002, 'a 40 km ridge keeps 0.98 of its amplitude at 5 km');
  close(spectrumSlope, -1.9, 0.03, 'the fit recovers a −1.9 spectrum');
});

function analyticFields(N, elevationAt, band, rows = 2160, cols = 4320) {
  const h5 = new Float32Array(rows * cols), flt2 = new Float32Array(rows * cols).fill(band * band);
  for (let r = 0; r < rows; r++) {
    const lat = Math.PI / 2 - (r + 0.5) * Math.PI / rows;
    for (let c = 0; c < cols; c++) h5[r * cols + c] = elevationAt(lat, -Math.PI + (c + 0.5) * 2 * Math.PI / cols);
  }
  const mesh = buildMesh(new Grid(N));
  const model = { mesh, geography: { land: new Uint8Array(mesh.nCells).fill(1) }, surfaceGeopotential: null, core: { diagnostics: { g: 9.80616 } } };
  return { mesh, fields: meshFields(model, h5, flt2, rows, cols), dLon: 2 * Math.PI / cols, dLat: Math.PI / rows };
}

test('an analytic ridge and an analytic isotropic field through the per-mesh script: μ, γ, θ, σ and σ_flt', () => {
  const amplitude = 300, waves = 400;
  const ridge = analyticFields(16, (lat, lon) => 1000 + amplitude * Math.cos(waves * lon), 20);
  const crate = analyticFields(16, (lat, lon) => 1000 + amplitude * Math.cos(waves * lon) * Math.cos(waves * lat), 35);
  let n = 0;
  for (let i = 0; i < ridge.mesh.nCells; i++) {
    const lat = ridge.mesh.latCell[i];
    if (Math.abs(lat) > 12 * deg) continue;
    n++;
    const kx = waves / (R * Math.cos(lat)), factor = Math.sin(waves * ridge.dLon) / (waves * ridge.dLon);
    close(ridge.fields.deviation[i], amplitude / Math.SQRT2, 0.02, `ridge μ at ${(lat / deg).toFixed(1)}°`);
    close(ridge.fields.slope[i], amplitude * kx * factor / Math.SQRT2, 0.03, `ridge σ at ${(lat / deg).toFixed(1)}°`);
    assert.ok(ridge.fields.anisotropy[i] < 0.03, `ridge γ ${ridge.fields.anisotropy[i]}`);
    assert.ok(Math.abs(Math.sin(ridge.fields.orientation[i])) < 0.02, `ridge θ ${ridge.fields.orientation[i]}: the steepest slope is east–west`);
    close(ridge.fields.filtered[i], 20, 1e-6, 'σ_flt of a uniform band variance');
    const k = waves / R, f = Math.sin(waves * crate.dLat) / (waves * crate.dLat);
    close(crate.fields.deviation[i], amplitude / 2, 0.02, `isotropic μ at ${(lat / deg).toFixed(1)}°`);
    close(crate.fields.slope[i], amplitude * k * f / 2 * Math.sqrt(1 / Math.cos(lat) ** 2), 0.04, `isotropic σ at ${(lat / deg).toFixed(1)}°`);
    assert.ok(crate.fields.anisotropy[i] > 0.95 * Math.cos(lat), `isotropic γ ${crate.fields.anisotropy[i]}`);
    close(crate.fields.filtered[i], 35, 1e-6, 'σ_flt');
  }
  assert.ok(n > 30, `${n} cells`);
});

test('the gradients reach 5 km each way: a north–south ridge of 20 km wavelength on 2′30″ rows', () => {
  const amplitude = 300, waves = 2000, rows = 4320, cols = 1440;
  const { mesh, fields } = analyticFields(16, (lat) => 1000 + amplitude * Math.cos(waves * lat), 20, rows, cols);
  const k = waves / R, step = R * Math.PI / rows, reach = 5000 / step, whole = Math.floor(reach), t = reach - whole;
  const factor = ((1 - t) * Math.sin(k * whole * step) + t * Math.sin(k * (whole + 1) * step)) / (k * reach * step);
  assert.ok(Math.sin(k * step) / (k * step) - factor > 0.08, 'the neighbouring rows would give a larger slope');
  let n = 0;
  for (let i = 0; i < mesh.nCells; i++) {
    if (Math.abs(mesh.latCell[i]) > 60 * deg) continue;
    n++;
    close(fields.slope[i], amplitude * k * factor / Math.SQRT2, 0.01, `σ at ${(mesh.latCell[i] / deg).toFixed(1)}°`);
    close(fields.deviation[i], amplitude / Math.SQRT2, 0.02, 'μ');
  }
  assert.ok(n > 1000, `${n} cells`);
});

/*
 * scripts/subgridTerrainHand.py's values for three N=64 cells, each
 * recomputed from GMTED2010's 30″ grid by direct sums (the cell's own 30″
 * points for σ_flt): the Great Plains, the Himalayan front, the Andes.
 */
const HAND_N64 = [
  { cell: 14486, where: 'the Great Plains, 39.4N 98.7W', deviation: 30.376, anisotropy: 0.53879, orientation: 1.38416, slope: 0.0032324, filtered: 15.014 },
  { cell: 4854, where: 'the Himalayan front, 28.5N 84.4E', deviation: 1421.24, anisotropy: 0.79910, orientation: 1.13378, slope: 0.087253, filtered: 461.36 },
  { cell: 34912, where: 'the Andes, 32.7S 70.2W', deviation: 940.09, anisotropy: 0.92716, orientation: 0.06030, slope: 0.059268, filtered: 359.34 },
];

test('the N=64 file holds three real cells as recomputed by hand from the 30″ raster', () => {
  const fields = decodeSubgrid(readFileSync(subgridUrl(64)));
  for (const hand of HAND_N64) {
    close(fields.deviation[hand.cell], hand.deviation, 0.005, `${hand.where} μ`);
    close(fields.slope[hand.cell], hand.slope, 0.005, `${hand.where} σ`);
    assert.ok(Math.abs(fields.anisotropy[hand.cell] - hand.anisotropy) < 0.005, `${hand.where} γ ${fields.anisotropy[hand.cell]}`);
    assert.ok(Math.abs(fields.orientation[hand.cell] - hand.orientation) < 0.005, `${hand.where} θ ${fields.orientation[hand.cell]}`);
    close(fields.filtered[hand.cell], hand.filtered, 0.01, `${hand.where} σ_flt`);
  }
});

const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);

test('every mesh the page offers has its file; sea cells hold nothing and the codec round-trips', () => {
  for (const N of [16, 32, 64, 128]) {
    assert.ok(existsSync(subgridUrl(N)), `data/subgrid_N${N}.bin`);
    const mesh = buildMesh(new Grid(N)), fields = meshSubgrid(mesh);
    assert.equal(fields.deviation.length, mesh.nCells);
    assert.equal(readFileSync(subgridUrl(N)).byteLength, 16 + 2 * SUBGRID_FIELDS.length * mesh.nCells);
    if (N > 32) continue;
    const model = createModel(new Grid(N), { physics: false, topography }), land = model.geography.land;
    let landWith = 0;
    for (let i = 0; i < mesh.nCells; i++) {
      if (!land[i]) for (const [name] of SUBGRID_FIELDS) if (name !== 'orientation') assert.equal(fields[name][i], 0, `N=${N} sea cell ${i} ${name}`);
      if (land[i] && fields.deviation[i] > 0 && fields.filtered[i] > 0) landWith++;
    }
    assert.ok(landWith > 0.95 * land.reduce((a, b) => a + b, 0), `N=${N}: ${landWith} land cells hold fields`);
    const again = decodeSubgrid(encodeSubgrid(fields));
    for (const [name, scale] of SUBGRID_FIELDS) for (let i = 0; i < mesh.nCells; i++) assert.ok(Math.abs(again[name][i] - fields[name][i]) <= 0.5 * scale + 1e-12, `${name} round trip`);
  }
  assert.equal(meshSubgrid(buildMesh(new Grid(8))), null, 'no file for N=8');
});

test('the files serve the bundled land mask with its terrain only: given fields lose their sea cells, and a run without terrain, on another mask or with the files off takes the raster’s', () => {
  const N = 16, mesh = buildMesh(new Grid(N)), file = meshSubgrid(mesh), g = 9.80616;
  const model = createModel(new Grid(N), { physics: false, topography }), { geography, surfaceGeopotential: phis } = model;
  const spread = Object.fromEntries(SUBGRID_FIELDS.map(([name]) => [name, Float64Array.from(file[name], (v, i) => (geography.land[i] ? v : 1))]));
  const fitted = orographyFields(mesh, topography, geography, phis, spread, g);
  assert.ok(!fitted.raster);
  for (let i = 0; i < mesh.nCells; i++) {
    if (geography.land[i]) assert.equal(fitted.deviation[i], file.deviation[i]);
    else for (const [name] of SUBGRID_FIELDS) assert.equal(fitted[name][i], 0, `sea cell ${i} ${name}`);
  }
  assert.ok(orographyFields(mesh, topography, geography, null, undefined, g).raster, 'no terrain');
  const said = [], log = console.log;
  console.log = (line) => said.push(line);
  const off = orographyFields(mesh, topography, geography, phis, false, g);
  console.log = log;
  assert.ok(off.raster && said.length === 1 && said[0].includes('turned off'), `files turned off: ${said.join(' | ')}`);
  const other = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
  const elsewhere = createModel(new Grid(N), { physics: false, topography: other });
  assert.ok(orographyFields(mesh, other, elsewhere.geography, elsewhere.surfaceGeopotential, undefined, g).raster, 'another land mask');
  assert.equal(createModel(new Grid(N), { topography: other }).boundaryLayer.formDrag, null);
  assert.equal(createModel(new Grid(N), { topography, terrain: false }).boundaryLayer.formDrag, null);
});
