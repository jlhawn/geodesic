// Where the stratocumulus deck's top sits against the grid over saved states
// (ocean columns; deck columns have mlmGate > 0.5). Per box: the deck's share,
// tops within 5 m of the ceiling, 100 m bins of the top, and the top's place in
// its layer by fifths. Per site: the 4 K and the largest θv jumps below 3 km,
// the deck's top, and the 4 K jump's height in deck columns by 50 m bins:
//   node scripts/deckTops.mjs <state.bin>...
import { readFileSync } from 'node:fs';
import { decodeState } from '../js/stateFile.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { stateGeography } from '../js/geography.module.js';
import { R_DRY, CP_DRY, P0, GRAVITY, VIRTUAL_FACTOR } from '../js/dynamics/sigmaCore.module.js';

const FILES = process.argv.slice(2), KAPPA = R_DRY / CP_DRY, JUMP = 4, TOP = 3000;
const box = (south, north, west, east) => (lat, lon) => lat >= south && lat <= north && lon >= west && lon <= east;
const BOXES = {
  'SE Pacific 25-10S 90-72W': box(-25, -10, -90, -72),
  'Namibia 25-10S 0-14E': box(-25, -10, 0, 14),
  'California 20-35N 135-118W': box(20, 35, -135, -118),
  'Canary 15-30N 30-15W': box(15, 30, -30, -15),
  'W Australia 35-20S 95-112E': box(-35, -20, 95, 112),
  'cold tongue 3S-3N 140-90W': box(-3, 3, -140, -90),
};
const SITES = {
  'DYCOMS-II 30-33N 124-119W': box(30, 33, -124, -119),
  '20S coast 22-18S 73-70W': box(-22, -18, -73, -70),
  '20S 75W 22-18S 77-73W': box(-22, -18, -77, -73),
  '20S 85W 22-18S 87-83W': box(-22, -18, -87, -83),
};
const boxes = {}, sites = {};
let mesh = null;
for (const file of FILES) {
  const saved = await decodeState(new Uint8Array(readFileSync(file)));
  mesh ??= stateGrid(saved, stateGeography(saved).topography).mesh;
  const C = mesh.nCells, K = saved.K, levels = saved.levels, L = saved.ocean.h.length / C;
  const z = new Float64Array(K + 1), tv = new Float64Array(K);
  for (let i = 0; i < C; i++) {
    const lat = mesh.latCell[i] * 180 / Math.PI, east = mesh.lonCell[i] * 180 / Math.PI, lon = east > 180 ? east - 360 : east;
    const inBoxes = Object.keys(BOXES).filter((name) => BOXES[name](lat, lon)), inSites = Object.keys(SITES).filter((name) => SITES[name](lat, lon));
    if (!inBoxes.length && !inSites.length) continue;
    let depth = 0;
    for (let k = 0; k < L; k++) depth += saved.ocean.h[k * C + i];
    if (!(depth > 1)) continue;
    const exner = (sigma) => sigma > 0 ? (saved.pi[i] * sigma / P0) ** KAPPA : 0;
    z[K] = 0;
    for (let k = K - 1; k >= 0; k--) {
      const n = k * C + i;
      tv[k] = saved.theta[n] * (1 + VIRTUAL_FACTOR * saved.q[n] - saved.qc[n]);
      z[k] = z[k + 1] + CP_DRY * tv[k] * (exner(levels[k + 1]) - exner(levels[k])) / GRAVITY;
    }
    const middle = (k) => 0.5 * (z[k] + z[k + 1]);
    let cap = -1, largest = -1, largestJump = 0;
    for (let k = K - 2; k >= 1 && z[k + 1] < TOP; k--) {
      const jump = tv[k] - tv[k + 1];
      if (cap < 0 && middle(k + 1) > saved.boundaryDepth[i] && jump >= JUMP) cap = k + 1;
      if (z[k + 1] >= 300 && jump > largestJump) { largestJump = jump; largest = k + 1; }
    }
    const h = saved.mlmHeight[i], deck = saved.mlmGate[i] > 0.5 && h > 0, ceiling = cap >= 0 ? middle(cap - 1) - 1 : NaN;
    let layer = -1;
    for (let k = K - 1; k >= 0; k--) if (h >= z[k + 1] && h < z[k]) { layer = k; break; }
    for (const name of inBoxes) {
      const b = boxes[name] ??= { columns: 0, decks: 0, atCeiling: 0, bins: {}, layers: {} };
      b.columns++;
      if (!deck) continue;
      b.decks++;
      if (Math.abs(h - ceiling) < 5) b.atCeiling++;
      const bin = Math.floor(h / 100) * 100;
      b.bins[bin] = (b.bins[bin] ?? 0) + 1;
      if (layer < 0) continue;
      const l = b.layers[layer] ??= { n: 0, bottom: 0, top: 0, fifths: [0, 0, 0, 0, 0] };
      l.n++; l.bottom += z[layer + 1]; l.top += z[layer];
      l.fifths[Math.min(4, Math.floor(5 * (h - z[layer + 1]) / (z[layer] - z[layer + 1])))]++;
    }
    for (const name of inSites) {
      const t = sites[name] ??= { columns: 0, largest: 0, largestJump: 0, capped: 0, cap: 0, decks: 0, top: 0, caps: {} };
      t.columns++;
      if (largest >= 0) { t.largest += z[largest]; t.largestJump += largestJump; }
      if (cap >= 0) { t.capped++; t.cap += z[cap]; }
      if (!deck) continue;
      t.decks++; t.top += h;
      const bin = cap >= 0 ? Math.round(z[cap] / 50) * 50 : 'none';
      t.caps[bin] = (t.caps[bin] ?? 0) + 1;
    }
  }
}
const pct = (a, b) => (100 * a / Math.max(1, b)).toFixed(0);
console.log(`${FILES.length} states; deck columns have mlmGate > 0.5`);
for (const [name, b] of Object.entries(boxes)) {
  console.log(`\n${name}: ${b.columns} columns, deck in ${pct(b.decks, b.columns)} %, ${pct(b.atCeiling, b.decks)} % of deck tops within 5 m of the ceiling`);
  console.log('  top, 100 m bins (% of deck tops): ' + Object.keys(b.bins).map(Number).sort((x, y) => x - y).filter((bin) => b.bins[bin] / b.decks > 0.01).map((bin) => `${bin}: ${pct(b.bins[bin], b.decks)}`).join(', '));
  console.log('  layer holding the top (bottom-top m, % of deck tops; fifths bottom to top, %):');
  for (const [k, l] of Object.entries(b.layers).sort((x, y) => y[0] - x[0])) if (l.n / b.decks > 0.03) console.log(`    k ${k}: ${(l.bottom / l.n).toFixed(0)}-${(l.top / l.n).toFixed(0)}, ${pct(l.n, b.decks)} %; ${l.fifths.map((v) => pct(v, l.n)).join(' ')}`);
}
console.log('');
for (const [name, t] of Object.entries(sites)) {
  console.log(`${name}: ${t.columns} columns, largest θv jump at ${(t.largest / t.columns).toFixed(0)} m (${(t.largestJump / t.columns).toFixed(1)} K), 4 K jump at ${(t.cap / Math.max(1, t.capped)).toFixed(0)} m in ${pct(t.capped, t.columns)} %, deck in ${pct(t.decks, t.columns)} % with its top at ${(t.top / Math.max(1, t.decks)).toFixed(0)} m`);
  console.log('  4 K jump in deck columns, 50 m bins (%): ' + Object.entries(t.caps).sort((x, y) => Number(x[0]) - Number(y[0])).filter(([, n]) => n / t.decks > 0.02).map(([bin, n]) => `${bin}: ${pct(n, t.decks)}`).join(', '));
}
