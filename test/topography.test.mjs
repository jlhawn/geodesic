import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const COLS = 1440;
const ROWS = 720;
const PATH = join(dirname(fileURLToPath(import.meta.url)), '..', 'data', 'topography_0p25.bin');

const buffer = readFileSync(PATH);
const view = new DataView(buffer.buffer, buffer.byteOffset, buffer.byteLength);

function latOf(row) { return 89.875 - 0.25 * row; }
function lonOf(col) { return -179.875 + 0.25 * col; }
function cellAt(lat, lon) {
  const row = Math.round((89.875 - lat) / 0.25);
  const col = Math.round((lon + 179.875) / 0.25);
  return view.getInt16((row * COLS + col) * 2, true);
}

test('file is 1440x720 int16 cells', () => {
  assert.equal(buffer.length, COLS * ROWS * 2);
});

test('area-weighted land fraction is close to Earth\'s ~29%', () => {
  let landArea = 0, totalArea = 0;
  for (let row = 0; row < ROWS; row++) {
    const weight = Math.cos((latOf(row) * Math.PI) / 180);
    for (let col = 0; col < COLS; col++) {
      totalArea += weight;
      if (view.getInt16((row * COLS + col) * 2, true) > 0) landArea += weight;
    }
  }
  const fraction = landArea / totalArea;
  console.log(`land fraction ${fraction.toFixed(4)}`);
  assert.ok(fraction > 0.27 && fraction < 0.31);
});

test('extreme elevations are physically plausible', () => {
  let min = Infinity, max = -Infinity, minIdx = -1, maxIdx = -1;
  for (let i = 0; i < ROWS * COLS; i++) {
    const v = view.getInt16(i * 2, true);
    if (v < min) { min = v; minIdx = i; }
    if (v > max) { max = v; maxIdx = i; }
  }
  const maxLat = latOf(Math.floor(maxIdx / COLS));
  const maxLon = lonOf(maxIdx % COLS);
  console.log(`min ${min} m, max ${max} m at ${maxLat.toFixed(3)},${maxLon.toFixed(3)}`);
  assert.ok(max > 5000 && max < 8900);
  // The 0.25 degree area mean is highest over the Kunlun/Karakoram plateau,
  // not the Everest massif: valleys pull the block mean down next to sharp
  // peaks, so the check is against High Mountain Asia broadly rather than a
  // single summit.
  assert.ok(maxLat > 25 && maxLat < 40 && maxLon > 70 && maxLon < 100);
  assert.ok(min < -8000);
});

test('known land and ocean cells have the right sign and magnitude', () => {
  assert.ok(cellAt(40, -100) > 0); // central USA
  assert.ok(cellAt(0, -160) < -3000); // central Pacific abyssal plain
  assert.ok(cellAt(-75, 90) > 2000); // East Antarctica: ice surface, not bedrock
});

test('the Panama isthmus is closed at N=64, where the land fraction alone would open a seaway', async () => {
  const { Grid } = await import('../js/grid.module.js');
  const { buildMesh } = await import('../js/mesh.module.js');
  const { createGeography, topographyFromInt16 } = await import('../js/geography.module.js');
  const mesh = buildMesh(new Grid(64)), topography = topographyFromInt16(buffer.buffer.slice(buffer.byteOffset, buffer.byteOffset + buffer.byteLength));
  const deg = 180 / Math.PI, C = mesh.nCells, lat = (i) => mesh.latCell[i] * deg, lon = (i) => mesh.lonCell[i] * deg;
  const inBox = (i) => lat(i) >= 0 && lat(i) <= 15 && lon(i) >= -92 && lon(i) <= -70;
  const nearest = (la, lo) => { let best = -1, d = 1e9; for (let i = 0; i < C; i++) { const dd = (lat(i) - la) ** 2 + (lon(i) - lo) ** 2; if (dd < d) { d = dd; best = i; } } return best; };
  const connected = (geo) => {
    const seen = new Uint8Array(C), queue = [nearest(7.5, -80.5)], target = nearest(10.5, -79.5);
    while (queue.length) { const i = queue.pop(); if (i === target) return true; for (let m = 0; m < mesh.nEdgesOnCell[i]; m++) { const j = mesh.cellsOnCell[mesh.maxEdges * i + m]; if (j < 0 || seen[j] || geo.land[j] || !inBox(j)) continue; seen[j] = 1; queue.push(j); } }
    return false;
  };
  assert.equal(connected(createGeography(mesh, topography, { landBridges: {} })), true, 'without the bridge the coarse mask opens Panama');
  assert.equal(connected(createGeography(mesh, topography)), false, 'the default bridge closes it');
});

test('Hormuz and Bab-el-Mandeb stay open at N=64 through the default sea straits', async () => {
  const { Grid } = await import('../js/grid.module.js');
  const { buildMesh } = await import('../js/mesh.module.js');
  const { createGeography, topographyFromInt16 } = await import('../js/geography.module.js');
  const mesh = buildMesh(new Grid(64)), topography = topographyFromInt16(buffer.buffer.slice(buffer.byteOffset, buffer.byteOffset + buffer.byteLength));
  const deg = 180 / Math.PI, C = mesh.nCells, lat = (i) => mesh.latCell[i] * deg, lon = (i) => mesh.lonCell[i] * deg;
  const connected = (geo, a, b, box) => {
    const inBox = (i) => lat(i) >= box[0] && lat(i) <= box[1] && lon(i) >= box[2] && lon(i) <= box[3];
    const nearestSea = (la, lo) => { let best = -1, d = 1e9; for (let i = 0; i < C; i++) { if (geo.land[i]) continue; const dd = (lat(i) - la) ** 2 + (lon(i) - lo) ** 2; if (dd < d) { d = dd; best = i; } } return best; };
    const seen = new Uint8Array(C), queue = [nearestSea(...a)], target = nearestSea(...b);
    while (queue.length) { const i = queue.pop(); if (i === target) return true; for (let m = 0; m < mesh.nEdgesOnCell[i]; m++) { const j = mesh.cellsOnCell[mesh.maxEdges * i + m]; if (j < 0 || seen[j] || geo.land[j] || !inBox(j)) continue; seen[j] = 1; queue.push(j); } }
    return false;
  };
  const straits = [[[26.5, 52.0], [25.0, 57.5], [22, 31, 46, 62]], [[15.0, 41.5], [12.0, 45.5], [10, 30, 32, 52]]];
  const bare = createGeography(mesh, topography, { seaStraits: {} }), forced = createGeography(mesh, topography);
  for (const [a, b, box] of straits) {
    assert.equal(connected(bare, a, b, box), false, 'the coarse mask alone closes the strait');
    assert.equal(connected(forced, a, b, box), true, 'the default strait keeps it open');
  }
  for (let i = 0; i < C; i++) if (!forced.land[i] && bare.land[i]) assert.ok(forced.elevation[i] <= -50, `a forced strait cell is at least 50 m deep, got ${forced.elevation[i]}`);
});
