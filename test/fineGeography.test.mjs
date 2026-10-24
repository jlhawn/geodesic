import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { gridMesh } from '../js/stateGrid.module.js';
import { relaxedMesh, nearestCell, fineEdges, encodeFineGeography, decodeFineGeography, coarseGeography, connectedGeography, maskGeography, FIELDS } from '../js/fineGeography.module.js';
import { encodeMask, decodeMask, createGeography } from '../js/geography.module.js';
import { lakeFeatures, rasterizeLakes, lakeLevels } from '../scripts/lakes.mjs';
import { landRaster } from '../scripts/landWater.mjs';

test('relaxedMesh reproduces Grid bit for bit at N=16', () => {
  const grid = new Grid(16, { relax: 30 }), mesh = relaxedMesh(16, 30);
  for (const cell of grid) {
    const i = cell.index, c = cell.centerVertex;
    assert.deepEqual([mesh.xCell[3 * i], mesh.xCell[3 * i + 1], mesh.xCell[3 * i + 2], mesh.area[i]], [c.x, c.y, c.z, cell.area]);
    assert.deepEqual(Array.from(mesh.cellsOnCell.subarray(6 * i, 6 * i + cell.neighbors.length)), cell.neighbors.map((n) => n.index));
  }
});

test('a fine geography survives its file, and a coarse mesh takes its land by area and its sills from the fine edges', () => {
  const fine = relaxedMesh(64, 30), C = fine.nCells, { nEdges, cellsOnEdge } = fineEdges(fine), z = (f) => fine.xCell[3 * f + 2];
  const fields = Object.fromEntries(Object.keys(FIELDS).map((name) => [name, new Float64Array(C)]));
  for (let f = 0; f < C; f++) {
    const land = z(f) > 0.3 ? 1 : 0;
    fields.points[f] = 14; fields.landFraction[f] = land; fields.elevation[f] = land ? 500 : -4000;
    fields.landElevation[f] = land ? 500 : 0; fields.seaElevation[f] = land ? 0 : -4000; fields.bed[f] = fields.elevation[f]; fields.maxDepth[f] = land ? 0 : 4000;
  }
  const edgeSill = Float64Array.from({ length: nEdges }, (_, e) => (fields.maxDepth[cellsOnEdge[2 * e]] > 0 && fields.maxDepth[cellsOnEdge[2 * e + 1]] > 0 ? 3000 : 0));
  const parts = encodeFineGeography(fine, fields, edgeSill), bytes = new Uint8Array(parts.reduce((n, p) => n + p.byteLength, 0));
  parts.reduce((at, p) => { bytes.set(p, at); return at + p.byteLength; }, 0);
  const read = decodeFineGeography(bytes);
  assert.equal(read.nCells, C);
  assert.equal(read.nEdges, nEdges);
  assert.deepEqual(Array.from(read.landFraction), Array.from(fields.landFraction));
  const mesh = gridMesh(new Grid(16)), coarse = coarseGeography(mesh, read);
  let fineLand = 0, coarseLand = 0, total = 0;
  for (let f = 0; f < C; f++) { fineLand += read.area[f] * read.landFraction[f]; total += read.area[f]; }
  for (let i = 0; i < mesh.nCells; i++) coarseLand += coarse.area[i] * coarse.landFraction[i];
  assert.ok(Math.abs(coarseLand - fineLand) < 1e-9 * total);
  for (let i = 0; i < mesh.nCells; i++) {
    const zi = mesh.xCell[3 * i + 2];
    if (zi > 0.45) assert.equal(coarse.land[i], 1);
    if (zi < 0.15) assert.equal(coarse.land[i], 0);
  }
  for (let e = 0; e < mesh.nEdges; e++) {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
    if (coarse.landFraction[a] === 0 && coarse.landFraction[b] === 0) assert.equal(coarse.edgeSill[e], 3000);
    if (coarse.landFraction[a] === 1 || coarse.landFraction[b] === 1) assert.equal(coarse.edgeSill[e], 0);
  }
  const rule = connectedGeography(mesh, coarse, read);
  assert.equal(rule.flipped.filter(Boolean).length, 0);
});

test('lakes fill the raster points whose centers lie inside an odd number of their rings, intermittent ones are left out, and a lake\'s level is its flat surface or else its shore\'s low quantile', () => {
  const square = (w, s, e, n) => [[w, s], [e, s], [e, n], [w, n], [w, s]];
  const geojson = { features: [
    { properties: { featurecla: 'Lake', name: 'ring' }, geometry: { type: 'Polygon', coordinates: [square(10, 10, 20, 20), square(14, 14, 16, 16)] } },
    { properties: { featurecla: 'Reservoir', name: 'pair' }, geometry: { type: 'MultiPolygon', coordinates: [[square(-40, -10, -35, -5)], [square(-30, -10, -25, -5)]] } },
    { properties: { featurecla: 'Playa', name: 'pan' }, geometry: { type: 'Polygon', coordinates: [square(100, 0, 110, 10)] } },
  ] };
  const { lakes, left } = lakeFeatures(geojson);
  assert.deepEqual([lakes.map((l) => l.name), left.map((l) => l.name)], [['ring', 'pair'], ['pan']]);
  const rows = 180, cols = 360, id = rasterizeLakes(lakes, rows, cols), at = (lat, lon) => id[Math.floor(90 - lat) * cols + Math.floor(lon + 180)];
  let ring = 0, pair = 0;
  for (const k of id) { if (k === 1) ring++; if (k === 2) pair++; }
  assert.deepEqual([ring, pair], [100 - 4, 50]);
  assert.deepEqual([at(12.5, 12.5), at(15.5, 15.5), at(-7.5, -37.5), at(-7.5, -32.5), at(5.5, 105.5)], [1, 0, 2, 0, 0]);
  const surface = Float64Array.from({ length: rows * cols }, (_, n) => (id[n] === 1 ? 200 : id[n] === 2 ? (n % 7) * 10 - 100 : 250 + (n % 10)));
  const [flat, bed] = lakeLevels(id, surface, rows, cols, lakes.length);
  assert.deepEqual([flat.level, flat.flat], [200, true]);
  assert.equal(bed.flat, false);
  assert.ok(bed.level >= 250 && bed.level <= 252, `the shore's level ${bed.level}`);
});

const meshes = new Map();
function synthetic(N, depthOf, sillOf, lakes = null) {
  if (!meshes.has(N)) {
    const coarse = gridMesh(new Grid(N)), fine = relaxedMesh(4 * N, 30), edges = fineEdges(fine), owner = new Int32Array(fine.nCells);
    for (let f = 0, c = 0; f < fine.nCells; f++) owner[f] = c = nearestCell(fine.xCell[3 * f], fine.xCell[3 * f + 1], fine.xCell[3 * f + 2], c, coarse.xCell, coarse.nEdgesOnCell, coarse.cellsOnCell, coarse.maxEdges);
    meshes.set(N, { coarse, fine: { ...fine, ...edges }, owner });
  }
  const { coarse, fine, owner } = meshes.get(N), C = fine.nCells;
  const depth = Float64Array.from({ length: C }, (_, f) => depthOf(f, owner[f], fine));
  const fields = { points: new Float64Array(C).fill(14), landFraction: depth.map((d) => (d > 0 ? 0 : 1)), elevation: depth.map((d) => (d > 0 ? -d : 300)), landElevation: depth.map((d) => (d > 0 ? 0 : 300)), seaElevation: depth.map((d) => -d), bed: depth.map((d) => (d > 0 ? -d : 300)), maxDepth: depth };
  fields.dryMaxDepth = Float64Array.from(depth);
  fields.wetDepth = Float64Array.from(depth);
  for (const name of Object.keys(FIELDS)) fields[name] ??= new Float64Array(C);
  if (lakes) {
    for (let f = 0; f < C; f++) {
      const lake = lakes(f, owner[f], fine, owner);
      if (!lake) continue;
      fields.lakeFraction[f] = 1; fields.lakeLevel[f] = lake.level; fields.lakeDepth[f] = lake.depth; fields.lakeId[f] = lake.id; fields.dryMaxDepth[f] = 0;
    }
  }
  const edgeSill = Float64Array.from({ length: fine.nEdges }, (_, e) => {
    const f = fine.cellsOnEdge[2 * e], g = fine.cellsOnEdge[2 * e + 1];
    return depth[f] > 0 && depth[g] > 0 ? sillOf(f, g, owner, Math.min(depth[f], depth[g])) : 0;
  });
  const geography = { ...fine, ...fields, edgeSill }, raw = coarseGeography(coarse, geography);
  return { coarse, owner, fine: geography, raw, rule: connectedGeography(coarse, raw, geography), mask: maskGeography(coarse, geography) };
}
const cellNear = (mesh, lat, lon) => {
  const r = Math.PI / 180;
  return nearestCell(Math.cos(lat * r) * Math.cos(lon * r), Math.cos(lat * r) * Math.sin(lon * r), Math.sin(lat * r), 0, mesh.xCell, mesh.nEdgesOnCell, mesh.cellsOnCell, mesh.maxEdges);
};
const neighbors = (mesh, i) => Array.from({ length: mesh.nEdgesOnCell[i] }, (_, m) => mesh.cellsOnCell[mesh.maxEdges * i + m]);
function joinedAt(mesh, geography, a, b, depth) {
  const seen = new Uint8Array(mesh.nCells), stack = [a];
  seen[a] = 1;
  while (stack.length) {
    const i = stack.pop();
    if (i === b) return true;
    for (let m = 0; m < mesh.nEdgesOnCell[i]; m++) {
      const j = mesh.cellsOnCell[mesh.maxEdges * i + m];
      if (!seen[j] && !geography.land[j] && geography.edgeSill[mesh.edgesOnCell[mesh.maxEdges * i + m]] >= depth) { seen[j] = 1; stack.push(j); }
    }
  }
  return false;
}
const controls = (mesh, geography, a, b, depth) => joinedAt(mesh, geography, a, b, depth) && !joinedAt(mesh, geography, a, b, depth + 1);
const along = (fine, [lat0, lon0], [lat1, lon1]) => {
  const r = Math.PI / 180, p = [Math.cos(lat0 * r) * Math.cos(lon0 * r), Math.cos(lat0 * r) * Math.sin(lon0 * r), Math.sin(lat0 * r)], q = [Math.cos(lat1 * r) * Math.cos(lon1 * r), Math.cos(lat1 * r) * Math.sin(lon1 * r), Math.sin(lat1 * r)];
  const cells = new Set(), steps = 4000;
  for (let n = 0, c = 0; n <= steps; n++) {
    const t = n / steps, v = p.map((x, k) => (1 - t) * x + t * q[k]), s = Math.hypot(...v);
    c = nearestCell(v[0] / s, v[1] / s, v[2] / s, c, fine.xCell, fine.nEdgesOnCell, fine.cellsOnCell, 6);
    cells.add(c);
  }
  return cells;
};

for (const N of [16, 64]) {
  test(`N=${N}: a two-cell estuary that a channel joins to the ocean turns to land, a three-cell sea is joined`, () => {
    const { coarse } = synthetic(N, () => 0, () => 0), i3 = cellNear(coarse, 40, 0), i2 = cellNear(coarse, 40, 120);
    const sea = [i3, ...neighbors(coarse, i3).slice(0, 2)], estuary = [i2, neighbors(coarse, i2)[0]], ocean = cellNear(coarse, -70, 0);
    let channels = null;
    const { raw, rule } = synthetic(N, (f, o, fine) => {
      channels ??= new Set([...along(fine, [40, 0], [-45, 0]), ...along(fine, [40, 120], [-45, 120])]);
      if (fine.xCell[3 * f + 2] < -0.5) return 4000;
      if (sea.includes(o) || estuary.includes(o)) return 1000;
      return channels.has(f) ? 100 : 0;
    }, (f, g, o, d) => d);
    for (const i of [...sea, ...estuary]) assert.equal(raw.land[i], 0);
    for (const i of estuary) assert.equal(rule.land[i], 1, 'the estuary is land');
    for (const i of sea) assert.equal(rule.land[i], 0, 'the sea stays sea');
    assert.ok(controls(coarse, rule, sea[0], ocean, 100), 'the sea joins the ocean through its 100 m channel');
    assert.ok(!joinedAt(coarse, rule, estuary[0], ocean, 1));
    for (let i = 0; i < coarse.nCells; i++) if (rule.flipped[i] < 0) assert.ok(Math.abs(coarse.lonCell[i] * 180 / Math.PI) < 20, 'only the sea\'s channel is cut');
  });

  test(`N=${N}: a sill inside one cell holds along the model path, a channel through a third cell is found, and a wall stays closed`, () => {
    const { coarse } = synthetic(N, () => 0, () => 0), a = cellNear(coarse, 30, 0), x = coarse.xCell;
    const s = neighbors(coarse, a)[0], dot = (i, j) => x[3 * i] * x[3 * j] + x[3 * i + 1] * x[3 * j + 1] + x[3 * i + 2] * x[3 * j + 2];
    const b = neighbors(coarse, s).filter((j) => j !== a && !neighbors(coarse, a).includes(j)).sort((i, j) => dot(i, a) - dot(j, a))[0];
    const axis = [0, 1, 2].map((k) => x[3 * b + k] - x[3 * a + k]), side = (fine, f) => Math.sign([0, 1, 2].reduce((sum, k) => sum + axis[k] * (fine.xCell[3 * f + k] - x[3 * s + k]), 0));
    const row = [a, s, b];
    let fineRef = null;
    const one = synthetic(N, (f, o, fine) => { fineRef = fine; return row.includes(o) ? 1000 : 0; }, (f, g, o, d) => (o[f] === s && o[g] === s && side(fineRef, f) !== side(fineRef, g) ? 50 : d));
    assert.ok(controls(coarse, one.rule, a, b, 50), 'the 50 m sill inside the middle cell controls');

    const c = neighbors(coarse, a).find((j) => neighbors(coarse, s).includes(j));
    const vertex = [0, 1, 2].map((k) => x[3 * a + k] + x[3 * s + k] + x[3 * c + k]), norm = Math.hypot(...vertex), reach = 1.6 * Math.sqrt(4 * Math.PI / (10 * (4 * N) ** 2));
    const near = (fine, f) => Math.acos(Math.min(1, [0, 1, 2].reduce((sum, k) => sum + fine.xCell[3 * f + k] * vertex[k] / norm, 0))) < reach;
    const wall = (f, g, o, d) => ((o[f] === a && o[g] === s) || (o[f] === s && o[g] === a) ? 0 : d);
    const through = synthetic(N, (f, o, fine) => (o === a || o === s ? 1000 : o === c && near(fine, f) ? 200 : 0), wall);
    assert.equal(through.rule.land[c], 1);
    assert.ok(controls(coarse, through.rule, a, s, 200), 'the 200 m channel through the third cell joins the two cells');
    const closed = synthetic(N, (f, o) => (o === a || o === s ? 1000 : 0), wall);
    assert.ok(!joinedAt(coarse, closed.rule, a, s, 1), 'the wall keeps them apart');
  });

  test(`N=${N}: a cell is a lake where lakes cover most of it, water at the lake's level and depth open only to its own lake, and land where they cover less`, () => {
    const { coarse } = synthetic(N, () => 0, () => 0), i = cellNear(coarse, 40, 60), around = neighbors(coarse, i);
    const lake = [i, around[0], around[1]], most = around[3], least = around[4];
    const x = coarse.xCell, dot = (fine, f, j) => [0, 1, 2].reduce((sum, k) => sum + fine.xCell[3 * f + k] * x[3 * j + k], 0);
    let share = null;
    const { mask, raw } = synthetic(N, (f, o, fine) => (fine.xCell[3 * f + 2] < -0.5 ? 4000 : 0), (f, g, o, d) => d, (f, o, fine, owner) => {
      if (!share) {
        share = new Set();
        for (const [j, part] of [[most, 0.6], [least, 0.4]]) {
          const mine = [];
          for (let g = 0; g < fine.nCells; g++) if (owner[g] === j) mine.push(g);
          mine.sort((a, b) => dot(fine, b, i) - dot(fine, a, i));
          const total = mine.reduce((sum, g) => sum + fine.area[g], 0);
          let taken = 0;
          for (const g of mine) { if (taken >= part * total) break; share.add(g); taken += fine.area[g]; }
        }
      }
      return lake.includes(o) || share.has(f) ? { level: 300, depth: 80, id: 1 } : null;
    });
    assert.ok(Math.abs(raw.lakeFraction[most] - 0.6) < 0.1 && Math.abs(raw.lakeFraction[least] - 0.4) < 0.1, `lake fractions ${raw.lakeFraction[most]}, ${raw.lakeFraction[least]}`);
    for (const j of [...lake, most]) assert.deepEqual([mask.lake[j], mask.land[j], Math.round(mask.elevation[j]), Math.round(mask.surface[j])], [1, 0, -80, 300], `cell ${j} is the lake`);
    assert.deepEqual([mask.lake[least], mask.land[least]], [0, 1], 'the cell the lake covers less than half of is land');
    for (let e = 0; e < coarse.nEdges; e++) {
      const a = coarse.cellsOnEdge[2 * e], b = coarse.cellsOnEdge[2 * e + 1];
      if (mask.lake[a] && mask.lake[b]) assert.ok(mask.edgeOcean[e] && Math.abs(mask.edgeSill[e] - 80) < 1e-9);
      else if (mask.lake[a] || mask.lake[b]) assert.equal(mask.edgeOcean[e], 0);
    }
    assert.equal(mask.depressions.length, 0);
  });

  test(`N=${N}: a closed basin no deeper than any dry land turns to land at its own elevation, a deeper one stays sea, and one the lakes mark stays water`, () => {
    const { coarse } = synthetic(N, () => 0, () => 0), basin = (lat, lon) => { const i = cellNear(coarse, lat, lon); return [i, ...neighbors(coarse, i).slice(0, 2)]; };
    const dry = basin(40, 0), deep = basin(40, 120), wet = basin(40, -120);
    const { mask } = synthetic(N, (f, o, fine) => (fine.xCell[3 * f + 2] < -0.5 ? 4000 : dry.includes(o) ? 100 : deep.includes(o) ? 1000 : wet.includes(o) ? 200 : 0), (f, g, o, d) => d,
      (f, o) => (wet.includes(o) ? { level: -30, depth: 170, id: 1 } : null));
    for (const i of dry) assert.deepEqual([mask.land[i], mask.depression[i], Math.round(mask.elevation[i]), Math.round(mask.surface[i])], [1, 1, -100, -100], `dry cell ${i}`);
    for (const i of deep) assert.deepEqual([mask.land[i], mask.depression[i]], [0, 0], `deep cell ${i}`);
    for (const i of wet) assert.deepEqual([mask.land[i], mask.lake[i], Math.round(mask.elevation[i]), Math.round(mask.surface[i])], [0, 1, -170, -30], `lake cell ${i}`);
    assert.deepEqual([mask.depressions.length, mask.closed.length], [1, 1]);
    const parts = encodeMask(coarse, mask, mask), bytes = new Uint8Array(parts.reduce((n, p) => n + p.byteLength, 0));
    parts.reduce((at, p) => { bytes.set(p, at); return at + p.byteLength; }, 0);
    const read = createGeography(coarse, null, { mask: 'n1024', maskData: decodeMask(bytes) });
    assert.deepEqual([Array.from(read.land), Array.from(read.lake), Array.from(read.depression), Array.from(read.edgeOcean)], [Array.from(mask.land), Array.from(mask.lake), Array.from(mask.depression), Array.from(mask.edgeOcean)]);
    for (let i = 0; i < coarse.nCells; i++) assert.ok(Math.abs(read.surface[i] - mask.surface[i]) <= 0.5 && Math.abs(read.elevation[i] - mask.elevation[i]) <= 0.5, `cell ${i} read back`);
  });
}

test('GSHHG land fills a raster point only where its whole cell is land: inside the land, not the lake on it, but the island in that lake, across the date line, and never by Antarctica\'s levels', () => {
  const square = (w, s, e, n, level) => ({ level, xy: Float64Array.of(w, s, e, s, e, n, w, n) });
  const polygons = [square(0, 0, 60, 60, 1), square(10, 10, 50, 50, 2), square(20, 20, 40, 40, 3), square(170, -30, 190, -10, 1), square(-60, -80, 60, -70, 5), square(-30, 2, -21, 8, 1)];
  const rows = 18, cols = 36, land = landRaster(polygons, rows, cols);
  const at = (lat, lon) => land[Math.floor((90 - lat) / 10) * cols + Math.floor((lon + 180) / 10)];
  assert.equal(land.reduce((n, v) => n + v, 0), 24 + 4);
  assert.deepEqual([at(5, 5), at(15, 15), at(25, 25), at(55, 55), at(-15, 175), at(-25, -175), at(-75, 0), at(5, -25)], [1, 0, 1, 1, 1, 1, 0, 0]);
});
