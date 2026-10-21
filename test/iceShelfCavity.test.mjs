import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { gridMesh } from '../js/stateGrid.module.js';
import { bathymetryFrom } from '../js/ocean/layered.module.js';
import { syntheticTopography, createGeography, encodeLandIce, decodeLandIce, readTopography, geographyOptions, FRESH_GEOGRAPHY, ICE_DENSITY, SEAWATER_DENSITY } from '../js/geography.module.js';

const deg = 180 / Math.PI;
const sparse = ({ rows, cols, points, count, ...fields }) => {
  const index = Uint32Array.from(Array.from(count.keys()).filter((p) => count[p] > 0)), take = (Kind, field) => Kind.from(index, (p) => (field ? field[p] : 0));
  return decodeLandIce(encodeLandIce({
    rows, cols, points, index, count: take(Uint8Array, count), floating: take(Uint8Array, fields.floating), thickness: take(Uint16Array, fields.thickness),
    floatThickness: take(Uint16Array, fields.floatThickness), floatSurface: take(Uint16Array, fields.floatSurface), floatBed: take(Int16Array, fields.floatBed),
  }));
};

test('a cavity is land to the atmosphere and sea beneath its shelf, from the draft of the floating ice down to the bed, and stays land where less than cavityDepth of water lies under the shelf', () => {
  const mesh = gridMesh(new Grid(8)), rows = 90, cols = 180, points = 225, freeboard = 50, H = 400, draft = H * ICE_DENSITY / SEAWATER_DENSITY, deep = -800, shallow = -Math.round(draft + 30);
  const zone = (lat) => (lat < -70 ? 'grounded' : lat < -54 ? 'shallow' : lat < -38 ? 'deep' : 'sea');
  const raster = syntheticTopography(rows, cols, (lat) => ({ grounded: 500, shallow: freeboard, deep: freeboard, sea: -4000 })[zone(lat * deg)]);
  const pixels = (f) => Array.from({ length: rows * cols }, (_, n) => f(zone(90 - (Math.floor(n / cols) + 0.5) * 180 / rows))), afloat = (z) => z === 'shallow' || z === 'deep';
  const ice = sparse({
    rows, cols, points, count: pixels((z) => (z === 'sea' ? 0 : points)), thickness: pixels((z) => (z === 'grounded' ? 1500 : afloat(z) ? H : 0)),
    floating: pixels((z) => (afloat(z) ? points : 0)), floatThickness: pixels((z) => (afloat(z) ? H : 0)), floatSurface: pixels((z) => (afloat(z) ? freeboard : 0)),
    floatBed: pixels((z) => (z === 'deep' ? deep : z === 'shallow' ? shallow : 0)),
  });
  const options = { landBridges: {}, seaStraits: {}, landIce: ice }, shelf = createGeography(mesh, raster, options), cavity = createGeography(mesh, raster, { ...options, iceShelves: 'cavity' });
  for (const name of ['land', 'landFraction', 'elevation', 'iceSheet', 'iceFraction', 'iceThickness', 'shelfFraction', 'shelfThickness', 'wetDepth', 'edgeOcean', 'edgeSill', 'coastEdges']) {
    assert.deepEqual(Array.from(cavity[name]), Array.from(shelf[name]), `the atmosphere's ${name}`);
  }
  const { ocean } = cavity, D = bathymetryFrom(mesh, ocean), pure = new Int8Array(mesh.nCells).fill(-1), seen = { grounded: 0, shallow: 0, deep: 0, sea: 0 };
  for (let i = 0; i < mesh.nCells; i++) {
    const lat = mesh.latCell[i] * deg, z = zone(lat);
    if (zone(lat - 5) !== z || zone(lat + 5) !== z) continue;
    seen[z]++; pure[i] = 1;
    const wet = { grounded: 0, shallow: 0, deep: 1, sea: 1 }[z];
    assert.deepEqual([cavity.cavity[i], 1 - ocean.land[i]], [z === 'deep' ? 1 : 0, wet], `${z} cell ${i} at ${lat.toFixed(1)}°`);
    if (z === 'deep') {
      assert.deepEqual([cavity.land[i], cavity.iceSheet[i]], [1, 1]);
      assert.ok(Math.abs(cavity.iceDraft[i] - draft) < 1e-9 && Math.abs(ocean.elevation[i] - deep) < 1e-9 && D[i] >= -deep, `cavity ${i}: draft ${cavity.iceDraft[i]} m, bed ${ocean.elevation[i]} m, depth ${D[i]} m`);
    } else assert.equal(cavity.iceDraft[i], 0);
    if (z === 'sea') assert.equal(ocean.elevation[i], -4000);
  }
  assert.ok(Object.values(seen).every((n) => n >= 5), JSON.stringify(seen));
  let inner = 0;
  for (let e = 0; e < mesh.nEdges; e++) {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
    assert.equal(ocean.edgeOcean[e], !ocean.land[a] && !ocean.land[b] ? 1 : 0);
    if (!ocean.edgeOcean[e]) assert.equal(ocean.edgeSill[e], 0);
    if (pure[a] > 0 && pure[b] > 0 && cavity.cavity[a] && cavity.cavity[b]) { inner++; assert.equal(ocean.edgeSill[e], -deep); }
  }
  assert.ok(inner >= 5, `${inner} edges inside the cavity`);
  const thin = createGeography(mesh, raster, { ...options, iceShelves: 'cavity', cavityDepth: 20 });
  for (let i = 0; i < mesh.nCells; i++) if (pure[i] > 0 && zone(mesh.latCell[i] * deg) === 'shallow') assert.equal(thin.cavity[i], 1, `shallow cell ${i} holds 30 m under its shelf`);
  assert.throws(() => createModel(new Grid(8), { topography: raster, geography: { ...options, iceShelves: 'cavity' } }), /does not take ice-shelf cavities yet/);
});

test('on ETOPO 2022 the cavities lie under the Antarctic shelves near their observed extent, Ross and Filchner–Ronne the largest, the Amundsen Sea shelves resolved at N=128', () => {
  const topography = readTopography('etopo2022'), { geography: base } = geographyOptions(FRESH_GEOGRAPHY);
  const regions = {
    Ross: (la, lo) => la < -77.3 && (lo > 158 || lo < -145),
    'Filchner–Ronne': (la, lo) => la < -74.4 && lo > -85 && lo < -28,
    Amundsen: (la, lo) => la > -76.5 && la < -72.5 && lo > -136 && lo < -98,
    Amery: (la, lo) => la > -73.8 && la < -68.5 && lo > 66 && lo < 75,
    Larsen: (la, lo) => la > -70 && la < -64 && lo > -66 && lo < -58,
  };
  const least = { 64: { all: [100, 1.3], Ross: [30, 0.38], 'Filchner–Ronne': [30, 0.38], Amundsen: [3, 0], Amery: [4, 0] }, 128: { all: [400, 1.2], Ross: [120, 0.38], 'Filchner–Ronne': [110, 0.36], Amundsen: [15, 0], Amery: [15, 0] } };
  for (const N of [64, 128]) {
    const mesh = gridMesh(new Grid(N)), shelf = createGeography(mesh, topography, base), g = createGeography(mesh, topography, { ...base, iceShelves: 'cavity' }), D = bathymetryFrom(mesh, g.ocean);
    for (const name of ['land', 'elevation', 'iceSheet']) assert.deepEqual(Array.from(g[name]), Array.from(shelf[name]), `N=${N}: the atmosphere's ${name}`);
    const out = {}, add = (k, i) => { const o = (out[k] ??= { n: 0, area: 0, draft: 0, water: 0 }); o.n++; o.area += mesh.areaCell[i] / 1e12; o.draft += g.iceDraft[i] * mesh.areaCell[i]; o.water += (D[i] - g.iceDraft[i]) * mesh.areaCell[i]; };
    for (let i = 0; i < mesh.nCells; i++) {
      if (!g.cavity[i]) continue;
      const la = mesh.latCell[i] * deg, lo = mesh.lonCell[i] * deg;
      assert.deepEqual([g.land[i], g.ocean.land[i], g.iceSheet[i]], [1, 0, 1], `N=${N}: cavity ${i} at ${la.toFixed(1)}°, ${lo.toFixed(1)}°`);
      assert.ok(la < -60 || la > 75, `N=${N}: cavity at ${la.toFixed(1)}°`);
      assert.ok(g.iceDraft[i] > 0 && D[i] - g.iceDraft[i] >= 50, `N=${N}: cavity ${i}: draft ${g.iceDraft[i].toFixed(0)} m over ${D[i].toFixed(0)} m`);
      add(Object.keys(regions).find((r) => regions[r](la, lo)) ?? (la < 0 ? 'other Antarctic' : 'Greenland'), i); add('all', i);
    }
    for (const [k, o] of Object.entries(out)) console.log(`N=${N} ${k}: ${o.n} cavities, ${o.area.toFixed(3)} million km², mean draft ${(o.draft / o.area / 1e12).toFixed(0)} m, mean water beneath ${(o.water / o.area / 1e12).toFixed(0)} m`);
    for (const [k, [n, area]] of Object.entries(least[N])) assert.ok(out[k] && out[k].n >= n && out[k].area >= area && out[k].area < 1.7, `N=${N} ${k}: ${JSON.stringify(out[k])}`);
  }
});
