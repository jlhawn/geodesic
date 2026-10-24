import { Grid } from './grid.module.js';

/*
 * The geography of a fine relaxed mesh (N = 1024 by default, about 10.5
 * million cells of ~7 km), built once from ETOPO 2022's 60″ rasters by
 * scripts/fineGeography.mjs, and the geography of a model mesh made from
 * it. Each fine cell carries area-weighted means of its 60″ points: the
 * land fraction, the mean surface, the mean surface of its land and of its
 * sea, the bed, the deepest point, and the ice of scripts/landIce.mjs's
 * rules. Each fine edge carries its sill: the greatest depth to which the
 * two cells' deepest points stay joined through their own points
 * (eight-connected), 0 where land parts them. A model cell takes the fine
 * cells whose centers lie nearest it, area-weighted.
 */

/*
 * Grid(N, { relax }) on typed arrays: the same seeds, topology and Lloyd
 * passes, with the operations in Grid's order so that the centers come out
 * bit for bit the same. Grid itself holds an object per cell and vertex,
 * ~25 GB at N = 1024.
 */
export function relaxedMesh(N, relax = 30) {
  let grid = new Grid(N, { relax: 0 });
  const C = grid.size, V = grid.vertices.length;
  const center = new Float64Array(3 * C), vertex = new Float64Array(3 * V);
  const nEdgesOnCell = new Uint8Array(C), cellsOnCell = new Int32Array(6 * C).fill(-1), verticesOnCell = new Int32Array(6 * C).fill(-1);
  for (const cell of grid) {
    const i = cell.index, c = cell.centerVertex;
    center[3 * i] = c.x; center[3 * i + 1] = c.y; center[3 * i + 2] = c.z;
    nEdgesOnCell[i] = cell.neighbors.length;
    cell.neighbors.forEach((n, k) => { cellsOnCell[6 * i + k] = n.index; verticesOnCell[6 * i + k] = cell.vertices[k].index; });
  }
  grid.vertices.forEach((v, j) => { vertex[3 * j] = v.x; vertex[3 * j + 1] = v.y; vertex[3 * j + 2] = v.z; });
  grid = null;
  const excess = (ax, ay, az, bx, by, bz, cx, cy, cz) => {
    const tx = by * cz - bz * cy, ty = bz * cx - bx * cz, tz = bx * cy - by * cx;
    return 2 * Math.atan2(tx * ax + ty * ay + tz * az, 1 + (ax * bx + ay * by + az * bz) + (bx * cx + by * cy + bz * cz) + (cx * ax + cy * ay + cz * az));
  };
  const centroid = new Float64Array(3 * C);
  for (let pass = 0; pass < relax; pass++) {
    for (let i = 0; i < C; i++) {
      const n = nEdgesOnCell[i], cx = center[3 * i], cy = center[3 * i + 1], cz = center[3 * i + 2];
      let x = 0, y = 0, z = 0;
      for (let k = 0; k < n; k++) {
        const a = 3 * verticesOnCell[6 * i + k], b = 3 * verticesOnCell[6 * i + (k + 1) % n];
        let sx = cx + vertex[a] + vertex[b], sy = cy + vertex[a + 1] + vertex[b + 1], sz = cz + vertex[a + 2] + vertex[b + 2];
        const s = 1 / (Math.sqrt(sx * sx + sy * sy + sz * sz) || 1);
        sx *= s; sy *= s; sz *= s;
        const e = excess(cx, cy, cz, vertex[a], vertex[a + 1], vertex[a + 2], vertex[b], vertex[b + 1], vertex[b + 2]);
        x += sx * e; y += sy * e; z += sz * e;
      }
      const s = 1 / (Math.sqrt(x * x + y * y + z * z) || 1);
      centroid[3 * i] = x * s; centroid[3 * i + 1] = y * s; centroid[3 * i + 2] = z * s;
    }
    center.set(centroid);
    for (let i = 0; i < C; i++) {
      const n = nEdgesOnCell[i];
      for (let k = 0; k < n; k++) {
        const p = cellsOnCell[6 * i + k], q = cellsOnCell[6 * i + (k + 1) % n];
        if (!(p > i && q > i)) continue;
        const ax = center[3 * i], ay = center[3 * i + 1], az = center[3 * i + 2];
        const ux = center[3 * p] - ax, uy = center[3 * p + 1] - ay, uz = center[3 * p + 2] - az;
        const wx = center[3 * q] - ax, wy = center[3 * q + 1] - ay, wz = center[3 * q + 2] - az;
        let tx = uy * wz - uz * wy, ty = uz * wx - ux * wz, tz = ux * wy - uy * wx;
        const o = 3 * verticesOnCell[6 * i + k];
        if (tx * tx + ty * ty + tz * tz === 0) {
          tx = (ax + center[3 * p] + center[3 * q]) * (1 / 3); ty = (ay + center[3 * p + 1] + center[3 * q + 1]) * (1 / 3); tz = (az + center[3 * p + 2] + center[3 * q + 2]) * (1 / 3);
          const s = 1 / (Math.sqrt(tx * tx + ty * ty + tz * tz) || 1);
          vertex[o] = tx * s; vertex[o + 1] = ty * s; vertex[o + 2] = tz * s;
          continue;
        }
        const s = 1 / (Math.sqrt(tx * tx + ty * ty + tz * tz) || 1);
        tx *= s; ty *= s; tz *= s;
        if (tx * ax + ty * ay + tz * az < 0) { tx = -tx; ty = -ty; tz = -tz; }
        vertex[o] = tx; vertex[o + 1] = ty; vertex[o + 2] = tz;
      }
    }
  }
  const area = new Float64Array(C);
  for (let i = 0; i < C; i++) {
    const n = nEdgesOnCell[i], cx = center[3 * i], cy = center[3 * i + 1], cz = center[3 * i + 2];
    for (let k = 0; k < n; k++) {
      const a = 3 * verticesOnCell[6 * i + k], b = 3 * verticesOnCell[6 * i + (k + 1) % n];
      area[i] += excess(cx, cy, cz, vertex[a], vertex[a + 1], vertex[a + 2], vertex[b], vertex[b + 1], vertex[b + 2]);
    }
  }
  return { N, relax, nCells: C, xCell: center, nEdgesOnCell, cellsOnCell, area };
}

/*
 * The walk to the nearest center over a cell graph from `start`: x, y, z
 * a unit vector, xCell the centers, neighbors `stride` to a cell.
 */
export function nearestCell(x, y, z, start, xCell, nEdgesOnCell, cellsOnCell, stride) {
  let best = start, bestDot = x * xCell[3 * start] + y * xCell[3 * start + 1] + z * xCell[3 * start + 2];
  for (;;) {
    let next = best;
    for (let m = 0; m < nEdgesOnCell[best]; m++) {
      const j = cellsOnCell[stride * best + m];
      const dot = x * xCell[3 * j] + y * xCell[3 * j + 1] + z * xCell[3 * j + 2];
      if (dot > bestDot) { bestDot = dot; next = j; }
    }
    if (next === best) return best;
    best = next;
  }
}

/*
 * data/geography1024.bin: a 32-byte header ('FGE1', N, relax, cells,
 * edges, then three spare words), then per cell its center (float32 ×
 * 3), area on the unit sphere (float32), edge count (uint8) and neighbors
 * (int32 × 6, −1 past a pentagon’s five), the fields of FIELDS, and per
 * edge, in the order of fineEdges, its sill (uint16, meters).
 */
export const FINE_GEOGRAPHY_URL = new URL('../data/geography1024.bin', import.meta.url);
export const FIELDS = {
  points: Uint32Array, landFraction: Float32Array, elevation: Float32Array, landElevation: Float32Array, seaElevation: Float32Array, bed: Float32Array, maxDepth: Float32Array,
  iceFraction: Float32Array, iceThickness: Float32Array, floatFraction: Float32Array, floatThickness: Float32Array, floatSurface: Float32Array, floatBed: Float32Array,
};

export function fineEdges({ nCells: C, nEdgesOnCell, cellsOnCell }) {
  let E = 0;
  for (let i = 0; i < C; i++) for (let k = 0; k < nEdgesOnCell[i]; k++) if (cellsOnCell[6 * i + k] > i) E++;
  const cellsOnEdge = new Int32Array(2 * E), edgeOfSlot = new Int32Array(6 * C).fill(-1);
  for (let i = 0, e = 0; i < C; i++) for (let k = 0; k < nEdgesOnCell[i]; k++) if (cellsOnCell[6 * i + k] > i) { cellsOnEdge[2 * e] = i; cellsOnEdge[2 * e + 1] = cellsOnCell[6 * i + k]; edgeOfSlot[6 * i + k] = e++; }
  for (let i = 0; i < C; i++) {
    for (let k = 0; k < nEdgesOnCell[i]; k++) {
      const j = cellsOnCell[6 * i + k];
      if (j > i) continue;
      for (let m = 0; m < nEdgesOnCell[j]; m++) if (cellsOnCell[6 * j + m] === i) edgeOfSlot[6 * i + k] = edgeOfSlot[6 * j + m];
    }
  }
  return { nEdges: E, cellsOnEdge, edgeOfSlot };
}

export function encodeFineGeography(mesh, fields, edgeSill) {
  const C = mesh.nCells, E = edgeSill.length;
  const parts = [Float32Array.from(mesh.xCell), Float32Array.from(mesh.area), mesh.nEdgesOnCell, mesh.cellsOnCell, ...Object.entries(FIELDS).map(([name, Type]) => Type.from(fields[name])), Uint16Array.from(edgeSill, (s) => Math.min(65535, Math.round(s)))];
  const head = new DataView(new ArrayBuffer(32));
  [70, 71, 69, 49].forEach((b, j) => head.setUint8(j, b));
  head.setUint32(4, mesh.N, true); head.setUint32(8, mesh.relax, true); head.setUint32(12, C, true); head.setUint32(16, E, true);
  return [new Uint8Array(head.buffer), ...parts.map((p) => new Uint8Array(p.buffer, p.byteOffset, p.byteLength))];
}

export function decodeFineGeography(bytes) {
  const buffer = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength), head = new DataView(buffer);
  if (String.fromCharCode(...new Uint8Array(buffer, 0, 4)) !== 'FGE1') throw new Error('not a fine geography file');
  const N = head.getUint32(4, true), relax = head.getUint32(8, true), C = head.getUint32(12, true), E = head.getUint32(16, true);
  let at = 32;
  const take = (Type, n) => { const out = new Type(buffer.slice(at, at + n * Type.BYTES_PER_ELEMENT)); at += n * Type.BYTES_PER_ELEMENT; return out; };
  const xCell = take(Float32Array, 3 * C), area = take(Float32Array, C), nEdgesOnCell = take(Uint8Array, C);
  const cellsOnCell = take(Int32Array, 6 * C), fields = {};
  for (const [name, Type] of Object.entries(FIELDS)) fields[name] = take(Type, C);
  const edgeSill = take(Uint16Array, E);
  if (at !== buffer.byteLength) throw new Error('fine geography file of the wrong size');
  const mesh = { N, relax, nCells: C, xCell, area, nEdgesOnCell, cellsOnCell };
  return { ...mesh, ...fineEdges(mesh), ...fields, edgeSill };
}

export function readFineGeography(url = FINE_GEOGRAPHY_URL) {
  const fs = globalThis.process?.getBuiltinModule?.('node:fs');
  if (!fs) throw new Error('the fine geography is read under node only');
  return decodeFineGeography(fs.readFileSync(url));
}

/*
 * A model mesh's geography from the fine one. Every fine cell goes to the
 * model cell nearest its center; a cell's fractions and means are its fine
 * cells' weighted by area (the land's elevation by land area, the sea's by
 * sea area). edgeSill on every model edge is the greatest depth to which
 * the deepest fine cells of its two cells stay joined through fine edges
 * between fine cells of those two cells, 0 where no water crosses. The
 * class is land where the land fraction exceeds landThreshold.
 */
export function coarseGeography(mesh, fine, { landThreshold = 0.5 } = {}) {
  const { nCells: C, nEdges: E, xCell, nEdgesOnCell, cellsOnCell, cellsOnEdge, maxEdges } = mesh;
  const F = fine.nCells, owner = new Int32Array(F);
  let cell = 0;
  for (let f = 0; f < F; f++) owner[f] = cell = nearestCell(fine.xCell[3 * f], fine.xCell[3 * f + 1], fine.xCell[3 * f + 2], cell, xCell, nEdgesOnCell, cellsOnCell, maxEdges);
  const area = new Float64Array(C), landArea = new Float64Array(C), elevation = new Float64Array(C), landElevation = new Float64Array(C), seaElevation = new Float64Array(C), maxDepth = new Float64Array(C), deepest = new Int32Array(C).fill(-1);
  const iceArea = new Float64Array(C), iceVolume = new Float64Array(C), floatArea = new Float64Array(C), floatVolume = new Float64Array(C);
  for (let f = 0; f < F; f++) {
    const i = owner[f], a = fine.area[f], l = fine.landFraction[f];
    area[i] += a; landArea[i] += a * l; elevation[i] += a * fine.elevation[f];
    landElevation[i] += a * l * fine.landElevation[f]; seaElevation[i] += a * (1 - l) * fine.seaElevation[f];
    iceArea[i] += a * fine.iceFraction[f]; iceVolume[i] += a * fine.iceFraction[f] * fine.iceThickness[f];
    floatArea[i] += a * fine.floatFraction[f]; floatVolume[i] += a * fine.floatFraction[f] * fine.floatThickness[f];
    if (fine.maxDepth[f] > maxDepth[i] || deepest[i] < 0) { maxDepth[i] = fine.maxDepth[f]; deepest[i] = f; }
  }
  const land = new Uint8Array(C), landFraction = new Float64Array(C), iceFraction = new Float64Array(C), iceThickness = new Float64Array(C), floatFraction = new Float64Array(C), floatThickness = new Float64Array(C);
  for (let i = 0; i < C; i++) {
    landFraction[i] = landArea[i] / area[i];
    elevation[i] /= area[i];
    landElevation[i] = landArea[i] > 0 ? landElevation[i] / landArea[i] : 0;
    seaElevation[i] = area[i] > landArea[i] ? seaElevation[i] / (area[i] - landArea[i]) : 0;
    iceFraction[i] = iceArea[i] / area[i]; iceThickness[i] = iceArea[i] > 0 ? iceVolume[i] / iceArea[i] : 0;
    floatFraction[i] = floatArea[i] / area[i]; floatThickness[i] = floatArea[i] > 0 ? floatVolume[i] / floatArea[i] : 0;
    land[i] = landFraction[i] > landThreshold ? 1 : 0;
  }
  const start = new Int32Array(C + 1), list = new Int32Array(F);
  for (let f = 0; f < F; f++) start[owner[f] + 1]++;
  for (let i = 0; i < C; i++) start[i + 1] += start[i];
  { const fill = Int32Array.from(start.subarray(0, C)); for (let f = 0; f < F; f++) list[fill[owner[f]]++] = f; }
  const slotEdge = fine.edgeOfSlot, edgeSill = new Float64Array(E), parent = new Map();
  const find = (p) => { let r = p; while (parent.get(r) !== r) r = parent.get(r); while (parent.get(p) !== r) { const q = parent.get(p); parent.set(p, r); p = q; } return r; };
  for (let e = 0; e < E; e++) {
    const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
    if (maxDepth[a] <= 0 || maxDepth[b] <= 0) continue;
    const links = [];
    for (const i of [a, b]) {
      for (let p = start[i]; p < start[i + 1]; p++) {
        const f = list[p];
        for (let k = 0; k < fine.nEdgesOnCell[f]; k++) {
          const g = fine.cellsOnCell[6 * f + k];
          if (g < f || (owner[g] !== a && owner[g] !== b)) continue;
          const s = fine.edgeSill[slotEdge[6 * f + k]];
          if (s > 0) links.push(s, f, g);
        }
      }
    }
    const order = Array.from({ length: links.length / 3 }, (_, n) => n).sort((m, n) => links[3 * n] - links[3 * m]);
    parent.clear();
    const touch = (f) => { if (!parent.has(f)) parent.set(f, f); };
    touch(deepest[a]); touch(deepest[b]);
    for (const n of order) {
      const f = links[3 * n + 1], g = links[3 * n + 2];
      touch(f); touch(g);
      const x = find(f), y = find(g);
      if (x !== y) parent.set(x, y);
      if (find(deepest[a]) === find(deepest[b])) { edgeSill[e] = links[3 * n]; break; }
    }
  }
  return { owner, area, land, landFraction, elevation, landElevation, seaElevation, maxDepth, iceFraction, iceThickness, floatFraction, floatThickness, edgeSill };
}

/*
 * The connectivity rule on a coarse geography. Water is resolved at fine
 * scale where it crosses fine edges (sill > 0) between fine cells that are
 * mostly sea; land where fine cells that are mostly land touch. A model
 * ocean edge is open where water crosses it at fine scale (edgeSill > 0).
 * Straits: where the resolved water joins model ocean components that the
 * model's open edges keep apart, the model land cells on the fine paths
 * between them turn to sea, at least as deep as the path's sill through
 * them, and no edge of theirs deeper than that sill. The paths are the
 * widest, any sill of openSill or more counting the same, then those
 * through the model cells with the least land. Bridges: where resolved
 * land joins model land components kept apart by model sea cells, those
 * cells on the fine paths through the model cells with the most land turn
 * to land, unless a strait needs them.
 */
export function connectedGeography(mesh, coarse, fine, { landThreshold = 0.5, openSill = 150 } = {}) {
  const { nCells: C, maxEdges, nEdgesOnCell, cellsOnCell, edgesOnCell } = mesh, F = fine.nCells, E = fine.nEdges;
  const { owner } = coarse, land = Uint8Array.from(coarse.land), elevation = Float64Array.from(coarse.elevation);
  const edgeSill = Float64Array.from(coarse.edgeSill), flipped = new Int8Array(C), crossing = new Float64Array(C).fill(Infinity);
  const components = (isLand) => {
    const id = new Int32Array(C).fill(-1);
    let n = 0;
    for (let s = 0; s < C; s++) {
      if (id[s] >= 0 || land[s] !== isLand) continue;
      const stack = [s]; id[s] = n;
      while (stack.length) {
        const i = stack.pop();
        for (let m = 0; m < nEdgesOnCell[i]; m++) {
          const j = cellsOnCell[maxEdges * i + m];
          if (id[j] >= 0 || land[j] !== isLand || (!isLand && !(edgeSill[edgesOnCell[maxEdges * i + m]] > 0))) continue;
          id[j] = n; stack.push(j);
        }
      }
      n++;
    }
    return { id, n };
  };
  const join = (isLand, allowed, weightOf) => {
    const { id, n } = components(isLand), parent = new Int32Array(F + n);
    for (let p = 0; p < F + n; p++) parent[p] = p;
    const find = (p) => { while (parent[p] !== p) { parent[p] = parent[parent[p]]; p = parent[p]; } return p; };
    const node = (f) => (land[owner[f]] === isLand ? F + id[owner[f]] : f);
    const buckets = new Map();
    for (let e = 0; e < E; e++) {
      const f = fine.cellsOnEdge[2 * e], g = fine.cellsOnEdge[2 * e + 1];
      if (!allowed(f, g, e)) continue;
      const w = weightOf(f, g, e);
      if (!buckets.has(w)) buckets.set(w, []);
      buckets.get(w).push(e);
    }
    const tree = [];
    for (const w of [...buckets.keys()].sort((a, b) => b - a)) {
      for (const e of buckets.get(w)) {
        const a = node(fine.cellsOnEdge[2 * e]), b = node(fine.cellsOnEdge[2 * e + 1]);
        if (a === b) continue;
        const x = find(a), y = find(b);
        if (x !== y) { parent[x] = y; tree.push(e, a, b, w); }
      }
    }
    const degree = new Int32Array(F + n), links = new Map();
    for (let t = 0; t < tree.length; t += 4) {
      for (const [p, q] of [[tree[t + 1], tree[t + 2]], [tree[t + 2], tree[t + 1]]]) { degree[p]++; if (!links.has(p)) links.set(p, []); links.get(p).push(q, t); }
    }
    const removed = new Uint8Array(tree.length / 4), leaves = [];
    for (const [p] of links) if (p < F && degree[p] === 1) leaves.push(p);
    while (leaves.length) {
      const p = leaves.pop();
      for (let k = 0; k < links.get(p).length; k += 2) {
        const q = links.get(p)[k], t = links.get(p)[k + 1];
        if (removed[t / 4]) continue;
        removed[t / 4] = 1; degree[p]--; degree[q]--;
        if (q < F && degree[q] === 1) leaves.push(q);
      }
    }
    const kept = [];
    for (let t = 0; t < tree.length; t += 4) if (!removed[t / 4]) kept.push(tree[t], tree[t + 1], tree[t + 2], tree[t + 3]);
    return { kept, before: n };
  };
  const seaFine = (f) => fine.landFraction[f] <= landThreshold;
  const most = (f, g) => Math.max(coarse.landFraction[owner[f]], coarse.landFraction[owner[g]]);
  const strait = join(0, (f, g, e) => seaFine(f) && seaFine(g) && fine.edgeSill[e] > 0, (f, g, e) => 1000 * Math.min(fine.edgeSill[e], openSill) + Math.round(999 * (1 - most(f, g))));
  const straitCells = new Set();
  for (let t = 0; t < strait.kept.length; t += 4) {
    for (const p of [strait.kept[t + 1], strait.kept[t + 2]]) {
      if (p >= F) continue;
      const i = owner[p];
      straitCells.add(i); crossing[i] = Math.min(crossing[i], fine.edgeSill[strait.kept[t]]);
    }
  }
  for (let t = 0; t < strait.kept.length; t += 4) {
    const a = owner[fine.cellsOnEdge[2 * strait.kept[t]]], b = owner[fine.cellsOnEdge[2 * strait.kept[t] + 1]];
    if (a === b) continue;
    for (let m = 0; m < nEdgesOnCell[a]; m++) if (cellsOnCell[maxEdges * a + m] === b) { const e = edgesOnCell[maxEdges * a + m]; edgeSill[e] = Math.max(edgeSill[e], fine.edgeSill[strait.kept[t]]); }
  }
  for (const i of straitCells) {
    land[i] = 0; flipped[i] = -1; elevation[i] = Math.min(coarse.seaElevation[i], -crossing[i]);
    for (let m = 0; m < nEdgesOnCell[i]; m++) { const e = edgesOnCell[maxEdges * i + m]; edgeSill[e] = Math.min(edgeSill[e], crossing[i]); }
  }
  const bridge = join(1, (f, g) => !seaFine(f) && !seaFine(g), (f, g) => Math.round(1000 * Math.min(coarse.landFraction[owner[f]], coarse.landFraction[owner[g]])));
  const conflicts = [];
  for (let t = 0; t < bridge.kept.length; t += 4) {
    for (const p of [bridge.kept[t + 1], bridge.kept[t + 2]]) {
      if (p >= F) continue;
      const i = owner[p];
      if (straitCells.has(i)) { if (!conflicts.includes(i)) conflicts.push(i); continue; }
      if (land[i]) continue;
      land[i] = 1; flipped[i] = 1; elevation[i] = Math.max(coarse.landElevation[i], 1);
    }
  }
  const edgeOcean = new Uint8Array(mesh.nEdges);
  for (let e = 0; e < mesh.nEdges; e++) edgeOcean[e] = !land[mesh.cellsOnEdge[2 * e]] && !land[mesh.cellsOnEdge[2 * e + 1]] && edgeSill[e] > 0 ? 1 : 0;
  return { land, elevation, edgeOcean, edgeSill, flipped, conflicts, seaComponents: strait.before, landComponents: bridge.before };
}
