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
 * The depth to which each fine cell's water stays joined to the deepest
 * fine cell's, the world ocean's: the widest path over fine edges, its
 * value the shallowest sill on the way; 0 for water the ocean never
 * reaches.
 */
export function oceanDepth(fine) {
  const F = fine.nCells, depth = new Float64Array(F), done = new Uint8Array(F), buckets = [];
  let seed = 0;
  for (let f = 1; f < F; f++) if (fine.maxDepth[f] > fine.maxDepth[seed]) seed = f;
  depth[seed] = Math.round(fine.maxDepth[seed]);
  buckets[depth[seed]] = [seed];
  for (let d = depth[seed]; d > 0; d--) {
    const list = buckets[d];
    if (!list) continue;
    while (list.length) {
      const f = list.pop();
      if (done[f]) continue;
      done[f] = 1;
      for (let k = 0; k < fine.nEdgesOnCell[f]; k++) {
        const g = fine.cellsOnCell[6 * f + k], v = Math.min(d, fine.edgeSill[fine.edgeOfSlot[6 * f + k]]);
        if (v > depth[g]) { depth[g] = v; (buckets[v] ??= []).push(g); }
      }
    }
  }
  return depth;
}

/*
 * A model mesh's geography from the fine one. Every fine cell goes to the
 * model cell nearest its center; a cell's fractions and means are its fine
 * cells' weighted by area (the land's elevation by land area, the sea's by
 * sea area). A cell's water is represented by its fine cell joined deepest
 * to the ocean (oceanDepth), the deepest of them for water the ocean never
 * reaches. edgeSill on a model edge is the depth to which the
 * representatives of its two cells stay joined through fine edges among
 * the fine cells of its two cells and of the two cells at its ends, so that
 * a channel leaving through a third cell counts and land between them gives
 * 0. Where that is shallower than both representatives' depth to the
 * ocean, the fine cells of every neighbor of the two cells join the search,
 * which can only deepen the edge. By the ultrametric inequality no model
 * path is then deeper than the fine water between its ends. The class is land where the land fraction
 * exceeds landThreshold.
 */
export function coarseGeography(mesh, fine, { landThreshold = 0.5, ocean = oceanDepth(fine) } = {}) {
  const { nCells: C, nEdges: E, xCell, nEdgesOnCell, cellsOnCell, cellsOnEdge, verticesOnEdge, cellsOnVertex, maxEdges } = mesh;
  const F = fine.nCells, owner = new Int32Array(F);
  let cell = 0;
  for (let f = 0; f < F; f++) owner[f] = cell = nearestCell(fine.xCell[3 * f], fine.xCell[3 * f + 1], fine.xCell[3 * f + 2], cell, xCell, nEdgesOnCell, cellsOnCell, maxEdges);
  const area = new Float64Array(C), landArea = new Float64Array(C), elevation = new Float64Array(C), landElevation = new Float64Array(C), seaElevation = new Float64Array(C), maxDepth = new Float64Array(C), rep = new Int32Array(C).fill(-1);
  const iceArea = new Float64Array(C), iceVolume = new Float64Array(C), floatArea = new Float64Array(C), floatVolume = new Float64Array(C);
  for (let f = 0; f < F; f++) {
    const i = owner[f], a = fine.area[f], l = fine.landFraction[f];
    area[i] += a; landArea[i] += a * l; elevation[i] += a * fine.elevation[f];
    landElevation[i] += a * l * fine.landElevation[f]; seaElevation[i] += a * (1 - l) * fine.seaElevation[f];
    iceArea[i] += a * fine.iceFraction[f]; iceVolume[i] += a * fine.iceFraction[f] * fine.iceThickness[f];
    floatArea[i] += a * fine.floatFraction[f]; floatVolume[i] += a * fine.floatFraction[f] * fine.floatThickness[f];
    maxDepth[i] = Math.max(maxDepth[i], fine.maxDepth[f]);
    const r = rep[i];
    if (fine.maxDepth[f] > 0 && (r < 0 || ocean[f] > ocean[r] || (ocean[f] === ocean[r] && fine.maxDepth[f] > fine.maxDepth[r]))) rep[i] = f;
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
  const edgeSill = new Float64Array(E), local = new Int32Array(F).fill(-1), SHIFT = 2 ** 26;
  let parent = new Int32Array(1024), keys = new Float64Array(1024), ends = new Int32Array(2048);
  const find = (p) => { while (parent[p] !== p) { parent[p] = parent[parent[p]]; p = parent[p]; } return p; };
  const sillWithin = (a, b, region) => {
    let nodes = 0, links = 0, sill = 0;
    for (const i of region) {
      for (let p = start[i]; p < start[i + 1]; p++) {
        if (nodes >= parent.length) { const grown = new Int32Array(2 * parent.length); grown.set(parent); parent = grown; }
        local[list[p]] = nodes; parent[nodes] = nodes; nodes++;
      }
    }
    for (const i of region) {
      for (let p = start[i]; p < start[i + 1]; p++) {
        const f = list[p];
        for (let k = 0; k < fine.nEdgesOnCell[f]; k++) {
          const g = fine.cellsOnCell[6 * f + k];
          if (g < f || local[g] < 0) continue;
          const s = fine.edgeSill[fine.edgeOfSlot[6 * f + k]];
          if (!(s > 0)) continue;
          if (links >= keys.length) { const k2 = new Float64Array(2 * keys.length); k2.set(keys); keys = k2; const e2 = new Int32Array(2 * keys.length); e2.set(ends); ends = e2; }
          keys[links] = s * SHIFT + links; ends[2 * links] = local[f]; ends[2 * links + 1] = local[g]; links++;
        }
      }
    }
    const order = keys.subarray(0, links).sort(), ra = local[rep[a]], rb = local[rep[b]];
    for (let q = links - 1; q >= 0; q--) {
      const n = order[q] % SHIFT, x = find(ends[2 * n]), y = find(ends[2 * n + 1]);
      if (x !== y) parent[x] = y;
      if (find(ra) === find(rb)) { sill = Math.floor(order[q] / SHIFT); break; }
    }
    for (const i of region) for (let p = start[i]; p < start[i + 1]; p++) local[list[p]] = -1;
    return sill;
  };
  for (let e = 0; e < E; e++) {
    const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
    if (rep[a] < 0 || rep[b] < 0) continue;
    const region = [a, b];
    for (const v of [verticesOnEdge[2 * e], verticesOnEdge[2 * e + 1]]) for (let m = 0; m < 3; m++) { const c = cellsOnVertex[3 * v + m]; if (c >= 0 && !region.includes(c)) region.push(c); }
    edgeSill[e] = sillWithin(a, b, region);
    if (!(edgeSill[e] > 0 && edgeSill[e] < Math.min(ocean[rep[a]], ocean[rep[b]]))) continue;
    for (const i of [a, b]) for (let m = 0; m < nEdgesOnCell[i]; m++) { const c = cellsOnCell[maxEdges * i + m]; if (!region.includes(c)) region.push(c); }
    edgeSill[e] = Math.max(edgeSill[e], sillWithin(a, b, region));
  }
  return { owner, rep, area, land, landFraction, elevation, landElevation, seaElevation, maxDepth, iceFraction, iceThickness, floatFraction, floatThickness, edgeSill };
}

/*
 * The connectivity rule on a coarse geography. Water is resolved at fine
 * scale where it crosses fine edges (sill > 0) between fine cells that are
 * mostly sea; land where fine cells that are mostly land touch. A model
 * ocean edge is open where its edgeSill is above 0. Model ocean components
 * are the sea cells joined by open edges. One of fewer than minimumSea
 * cells that resolved water joins to one of minimumSea or more is a
 * feature the mesh does not resolve and turns to land; one the water joins
 * to no other component (a lake, or dry land below sea level) is left as
 * it is. Straits: where resolved water joins the remaining components,
 * the model land cells on the fine paths between them turn to sea, as deep
 * as their deepest edge at least. The paths are the widest, any sill of
 * openSill or more counting the same, then those through the model cells
 * with the least land. Bridges: where resolved land joins model land
 * components kept apart by model sea cells, those cells on the fine paths
 * through the model cells with the most land turn to land, unless a strait
 * needs them. Bridges keep no minimum size.
 */
export function connectedGeography(mesh, coarse, fine, { landThreshold = 0.5, openSill = 150, minimumSea = 3 } = {}) {
  const { nCells: C, maxEdges, nEdgesOnCell, cellsOnCell, edgesOnCell } = mesh, F = fine.nCells, E = fine.nEdges;
  const { owner, edgeSill } = coarse, land = Uint8Array.from(coarse.land), elevation = Float64Array.from(coarse.elevation);
  const flipped = new Int8Array(C);
  const components = (isLand) => {
    const id = new Int32Array(C).fill(-1), size = [];
    let n = 0;
    for (let s = 0; s < C; s++) {
      if (id[s] >= 0 || land[s] !== isLand) continue;
      const stack = [s]; id[s] = n; size.push(1);
      while (stack.length) {
        const i = stack.pop();
        for (let m = 0; m < nEdgesOnCell[i]; m++) {
          const j = cellsOnCell[maxEdges * i + m];
          if (id[j] >= 0 || land[j] !== isLand || (!isLand && !(edgeSill[edgesOnCell[maxEdges * i + m]] > 0))) continue;
          id[j] = n; size[n]++; stack.push(j);
        }
      }
      n++;
    }
    return { id, size, n };
  };
  const join = (isLand, allowed, weightOf, smallest = 0) => {
    const { id, size, n } = components(isLand), parent = new Int32Array(F + n);
    for (let p = 0; p < F + n; p++) parent[p] = p;
    const find = (p) => { while (parent[p] !== p) { parent[p] = parent[parent[p]]; p = parent[p]; } return p; };
    const node = (f) => (land[owner[f]] === isLand && size[id[owner[f]]] >= smallest ? F + id[owner[f]] : f);
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
        if (x !== y) { parent[x] = y; tree.push(e, a, b); }
      }
    }
    const degree = new Int32Array(F + n), links = new Map();
    for (let t = 0; t < tree.length; t += 3) {
      for (const [p, q] of [[tree[t + 1], tree[t + 2]], [tree[t + 2], tree[t + 1]]]) { degree[p]++; if (!links.has(p)) links.set(p, []); links.get(p).push(q, t); }
    }
    const removed = new Uint8Array(tree.length / 3), leaves = [];
    for (const [p] of links) if (p < F && degree[p] === 1) leaves.push(p);
    while (leaves.length) {
      const p = leaves.pop();
      for (let k = 0; k < links.get(p).length; k += 2) {
        const q = links.get(p)[k], t = links.get(p)[k + 1];
        if (removed[t / 3]) continue;
        removed[t / 3] = 1; degree[p]--; degree[q]--;
        if (q < F && degree[q] === 1) leaves.push(q);
      }
    }
    const kept = [];
    for (let t = 0; t < tree.length; t += 3) if (!removed[t / 3]) kept.push(tree[t + 1], tree[t + 2]);
    return { kept, id, size, n, root: (c) => find(F + c) };
  };
  const seaFine = (f) => fine.landFraction[f] <= landThreshold;
  const water = (f, g, e) => seaFine(f) && seaFine(g) && fine.edgeSill[e] > 0;
  const most = (f, g) => Math.max(coarse.landFraction[owner[f]], coarse.landFraction[owner[g]]);
  const widest = (f, g, e) => 1000 * Math.min(fine.edgeSill[e], openSill) + Math.round(999 * (1 - most(f, g)));
  const cellsOf = (id, c) => { const out = []; for (let i = 0; i < C; i++) if (id[i] === c) out.push(i); return out; };
  const first = join(0, water, widest), withLarge = new Set(), seen = new Map(), unresolved = [], isolated = [];
  for (let c = 0; c < first.n; c++) {
    const r = first.root(c);
    seen.set(r, (seen.get(r) ?? 0) + 1);
    if (first.size[c] >= minimumSea) withLarge.add(r);
  }
  for (let c = 0; c < first.n; c++) {
    if (first.size[c] >= minimumSea) { if (seen.get(first.root(c)) === 1) isolated.push(cellsOf(first.id, c)); continue; }
    const cells = cellsOf(first.id, c);
    if (!withLarge.has(first.root(c))) { isolated.push(cells); continue; }
    unresolved.push(cells);
    for (const i of cells) { land[i] = 1; flipped[i] = 1; elevation[i] = Math.max(coarse.landElevation[i], 1); }
  }
  const strait = join(0, water, widest, minimumSea), straitCells = new Set();
  for (const p of strait.kept) if (p < F && land[owner[p]]) straitCells.add(owner[p]);
  for (const i of straitCells) {
    let deepest = 0;
    for (let m = 0; m < nEdgesOnCell[i]; m++) deepest = Math.max(deepest, edgeSill[edgesOnCell[maxEdges * i + m]]);
    land[i] = 0; flipped[i] = flipped[i] === 1 ? 0 : -1; elevation[i] = Math.min(coarse.seaElevation[i], -deepest);
  }
  const bridge = join(1, (f, g) => !seaFine(f) && !seaFine(g), (f, g) => Math.round(1000 * Math.min(coarse.landFraction[owner[f]], coarse.landFraction[owner[g]])));
  const conflicts = [];
  for (const p of bridge.kept) {
    if (p >= F) continue;
    const i = owner[p];
    if (straitCells.has(i)) { if (!conflicts.includes(i)) conflicts.push(i); continue; }
    if (land[i]) continue;
    land[i] = 1; flipped[i] = 1; elevation[i] = Math.max(coarse.landElevation[i], 1);
  }
  const edgeOcean = new Uint8Array(mesh.nEdges);
  for (let e = 0; e < mesh.nEdges; e++) edgeOcean[e] = !land[mesh.cellsOnEdge[2 * e]] && !land[mesh.cellsOnEdge[2 * e + 1]] && edgeSill[e] > 0 ? 1 : 0;
  return { land, elevation, edgeOcean, edgeSill, flipped, conflicts, unresolved, isolated, seaComponents: first.n, landComponents: bridge.n };
}
