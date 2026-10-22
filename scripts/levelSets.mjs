// Level sets side by side in the US Standard Atmosphere 1976 (1013.25 hPa at
// the ground): the thickness of the layer that holds each height, layers per
// band by middle height, the ratio of each layer's thickness to the one below
// it under 1 hPa, and the deck's mixing index z(1 − z/h)² for a jump at z of
// 0.3–3.6 km with the deck's top h half a layer above it (per κu*, in m):
//   node scripts/levelSets.mjs [name...]       (default: every bl grid)
//   TABLE=1 node scripts/levelSets.mjs bl60    (every layer, top to ground)
import { SIGMA_GRIDS, sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';

const P_GROUND = 101325, G0 = 9.80665, R_AIR = 287.053;
const US76 = [[0, 288.15, -0.0065, 101325], [11000, 216.65, 0, 22632.06], [20000, 216.65, 0.001, 5474.889], [32000, 228.65, 0.0028, 868.0187],
  [47000, 270.65, 0, 110.9063], [51000, 270.65, -0.0028, 66.93887], [71000, 214.65, -0.002, 3.956420], [84852, 186.946, 0, 0.3734]];
function height(p) {
  if (!(p > 0)) return Infinity;
  let base = US76[0];
  for (const layer of US76) if (p <= layer[3]) base = layer;
  const [z, t, lapse, pBase] = base;
  return lapse === 0 ? z + R_AIR * t / G0 * Math.log(pBase / p) : z + t / lapse * ((pBase / p) ** (R_AIR * lapse / G0) - 1);
}

const NAMES = process.argv.length > 2 ? process.argv.slice(2) : SIGMA_GRIDS.filter((name) => name.startsWith('bl'));
const HEIGHTS = [500, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 8000, 12000, 18000, 25000, 35000, 45000];
const STRATOPAUSE = height(100), TROPOPAUSE = 16600;
const BANDS = [[0, 1000], [1000, 2500], [2500, 5000], [5000, TROPOPAUSE], [TROPOPAUSE, STRATOPAUSE], [STRATOPAUSE, Infinity]];
const meters = (x) => Number.isFinite(x) ? Math.round(x).toLocaleString('en-US') : 'top';
for (const name of NAMES) {
  const sigma = Array.from(sigmaInterfaces(name)), K = sigma.length - 1;
  const z = sigma.map((s) => height(s * P_GROUND)), middle = (k) => height(0.5 * (sigma[k] + sigma[k + 1]) * P_GROUND);
  const holding = (h) => { for (let k = 0; k < K; k++) if (z[k] > h && z[k + 1] <= h) return z[k] - z[k + 1]; };
  const bands = BANDS.map(([a, b]) => { let n = 0; for (let k = 0; k < K; k++) if (middle(k) >= a && middle(k) < b) n++; return n; });
  const ratios = [];
  for (let k = 1; k < K - 1; k++) if (z[k - 1] <= STRATOPAUSE + 1) ratios.push((z[k - 1] - z[k]) / (z[k] - z[k + 1]));
  const index = [];
  for (let k = K - 1; k >= 1; k--) {
    const jump = z[k], h = jump + (z[k - 1] - jump) / 2;
    if (jump >= 300 && jump <= 3600) index.push(`${Math.round(jump)}: ${(jump * (1 - jump / h) ** 2).toFixed(1)}`);
  }
  const under = sigma.findIndex((s) => s * P_GROUND >= 100);
  console.log(`\n${name}: ${K} layers`);
  console.log(`  thickness (m) at ${HEIGHTS.map((h) => h / 1000).join(', ')} km: ${HEIGHTS.map((h) => meters(holding(h))).join(', ')}`);
  console.log(`  layers 0-1 / 1-2.5 / 2.5-5 / 5-16.6 km / 16.6 km-1 hPa / above: ${bands.join(' / ')}`);
  console.log(`  thickness over the one below, under 1 hPa: ${Math.min(...ratios).toFixed(2)} to ${Math.max(...ratios).toFixed(2)}; layer under 1 hPa ${meters(z[under] - z[under + 1])} m, the one over it ${meters(z[under - 1] - z[under])} m`);
  console.log(`  deck mixing index by jump height (m): ${index.join(', ')}`);
  if (!process.env.TABLE) continue;
  console.log('  k | σ top - bottom | p top - bottom (hPa) | thickness (m) | middle (m)');
  for (let k = 0; k < K; k++) console.log(`  ${k} | ${sigma[k].toPrecision(4)} - ${sigma[k + 1].toPrecision(4)} | ${(sigma[k] * P_GROUND / 100).toPrecision(4)} - ${(sigma[k + 1] * P_GROUND / 100).toPrecision(4)} | ${meters(z[k] - z[k + 1])} | ${meters(middle(k))}`);
}
