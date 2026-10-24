#!/usr/bin/env node
/*
 * The controlling depth of each passage of scripts/straits.mjs straight
 * from ETOPO 2022's 60″ surface (CACHE, as for scripts/landIce.mjs): the
 * widest path between the deepest points within 0.3° of the two basin
 * points, through the points of the passage's box, stepping to the four or
 * the eight neighbors of a point.
 *
 *   node scripts/straitSills.mjs [CACHE]
 */
import { openTiff } from './landIce.mjs';
import { STRAITS } from './straits.mjs';

const cache = process.argv[2] ?? '/Users/jlhawn/git_repos/jlhawn/geodesic-terrain-cache/etopo2022';
const tiff = openTiff(`${cache}/ETOPO_2022_v1_60s_N90W180_surface.tif`), tiles = new Map();
const sample = (r, c) => {
  const tx = Math.floor(c / tiff.tileWidth), ty = Math.floor(r / tiff.tileLength), key = ty * tiff.across + tx;
  if (!tiles.has(key)) tiles.set(key, tiff.readTile(tx, ty));
  return tiles.get(key)[(r - ty * tiff.tileLength) * tiff.tileWidth + (c - tx * tiff.tileWidth)];
};
const row = (lat) => Math.floor((90 - lat) * 60), col = (lon) => Math.floor((lon + 180) * 60);

export function controllingDepth({ a, b, box: [south, north, west, east] }, neighbors) {
  const r0 = row(north), r1 = row(south), c0 = col(west), c1 = col(east), W = c1 - c0 + 1, H = r1 - r0 + 1;
  const depth = new Float64Array(W * H);
  for (let r = 0; r < H; r++) for (let c = 0; c < W; c++) depth[r * W + c] = Math.max(0, -sample(r0 + r, c0 + c));
  const deepest = ([lat, lon]) => {
    let best = -1;
    for (let r = Math.max(0, row(lat + 0.3) - r0); r <= Math.min(H - 1, row(lat - 0.3) - r0); r++)
      for (let c = Math.max(0, col(lon - 0.3) - c0); c <= Math.min(W - 1, col(lon + 0.3) - c0); c++) if (best < 0 || depth[r * W + c] > depth[best]) best = r * W + c;
    return best;
  };
  const s = deepest(a), t = deepest(b), width = new Float64Array(W * H), done = new Uint8Array(W * H), buckets = [];
  width[s] = Math.round(depth[s]);
  buckets[width[s]] = [s];
  for (let d = width[s]; d > 0; d--) {
    const list = buckets[d];
    if (!list) continue;
    while (list.length) {
      const p = list.pop();
      if (done[p]) continue;
      done[p] = 1;
      if (p === t) return { depth: d, at: [+(90 - (r0 + Math.floor(p / W) + 0.5) / 60).toFixed(3)] };
      const r = Math.floor(p / W), c = p - r * W;
      for (const [dr, dc] of neighbors) {
        const rr = r + dr, cc = c + dc;
        if (rr < 0 || rr >= H || cc < 0 || cc >= W) continue;
        const q = rr * W + cc, v = Math.min(d, Math.round(depth[q]));
        if (v > width[q]) { width[q] = v; (buckets[v] ??= []).push(q); }
      }
    }
  }
  return { depth: 0 };
}

const FOUR = [[-1, 0], [1, 0], [0, -1], [0, 1]], EIGHT = [...FOUR, [-1, -1], [-1, 1], [1, -1], [1, 1]];
for (const [name, strait] of Object.entries(STRAITS)) console.log(`${name}: eight-connected ${controllingDepth(strait, EIGHT).depth} m, four-connected ${controllingDepth(strait, FOUR).depth} m`);
