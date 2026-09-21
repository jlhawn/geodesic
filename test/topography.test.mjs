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
